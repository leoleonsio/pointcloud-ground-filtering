#include "GroundFilter.h"

// -- LAS reading and writing
#include <lasreader.hpp>
#include <laswriter.hpp>

// -- CGAL delaunay triangulation
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>

// -- CGAL kd-tree
//#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>

typedef CGAL::Projection_traits_xy_3<Kernel>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> DT;

bool x_comparator(const Point &p1, const Point &p2){
    return p1.x() < p2.x();
}

bool y_comparator(const Point &p1, const Point &p2){
    return p1.y() < p2.y();
}

bool z_comparator(const Point &p1, const Point &p2){
    return p1.z() < p2.z();
}

// square used for the square grid
struct Square {
    double min_x, min_y, max_x, max_y, length;

    Square(){;}

    Square(double min_x, double min_y, double length){
        this->min_x = min_x;
        this->min_y = min_y;
        this->max_x = min_x + length;
        this->max_y = min_y + length;
        this->length = length;
    }
};

std::ostream& operator<<(std::ostream& os, const Square& s) {
    os << "Square {min_x=" << s.min_x << ", min_y=" << s.min_y << ", length=" << s.length << "}";
    return os;
}

// adding cells to the grid row - variant type
template <class T>
void add_cell(std::vector<T>& grid_row, double x, double y, double cellsize, double height);

// add square
template<>
void add_cell<Square>(std::vector<Square>& grid_row, double x, double y, double cellsize, double height){
    grid_row.emplace_back(Square(x, y, cellsize));
}

// grid can be made of squares (tin) or particles (csf)
template<typename cell_type>
struct Grid{
    unsigned int nrows = 0, ncols = 0;
    double min_x, max_x, min_y, max_y;
    std::vector<std::vector<cell_type>> values;

    Grid() {;}

    // initialize a 2d grid based on pointcloud extent
    Grid(const std::vector<Point>& pointcloud, double cellsize, double height=0){
        auto bb_x = std::minmax_element(pointcloud.begin(), pointcloud.end(), x_comparator);
//        std::cout << "XMIN, XMAX " << bb_x.first->x() << ", " << bb_x.second->x() << std::endl;
        auto bb_y = std::minmax_element(pointcloud.begin(), pointcloud.end(), y_comparator);
//        std::cout << "YMIN, YMAX " << bb_y.first->y() << ", " << bb_y.second->y() << std::endl;
        min_x = bb_x.first->x();
        max_x = bb_x.second->x();
        min_y = bb_y.first->y();
        max_y = bb_y.second->y();

        // fill the grid
        auto row = min_y;
        unsigned int nvalues = 0;
        while (row <= max_y){
            auto col = min_x;
            nrows++;
            std::vector<cell_type> grid_row;
            while (col <= max_x){
                nvalues++;
                add_cell<cell_type>(grid_row, col, row, cellsize, height);
                col += cellsize;
            }
            values.emplace_back(grid_row);
            row += cellsize;
        }
        ncols = nvalues / nrows;
    }

    cell_type &operator()(unsigned int row, unsigned int col) {
        assert(col >= 0 && col < ncols);
        assert(row >= 0 && row < nrows);
        return values[row][col];
    }

    cell_type operator()(unsigned int row, unsigned int col) const {
        assert(col >= 0 && col < ncols);
        assert(row >= 0 && row < nrows);
        return values[row][col];
    }

    // neighbors for a cell
    std::vector<cell_type> get_neighbors(unsigned int row, unsigned int col){
        std::vector<cell_type> neighbors;
        const std::vector<std::pair<int, int>> deltas {
            std::make_pair(0, 1),
            std::make_pair(1, 0),
            std::make_pair(-1, 0),
            std::make_pair(0, -1)
        };

        for (const auto &delta: deltas){
            int dx = delta.first;
            int dy = delta.second;

            unsigned int new_row = row + dy;
            unsigned int new_col = col + dx;

            // check if valid neighbor
            if (new_row >= 0 && new_row < nrows && new_col >= 0 && new_col < ncols){
                neighbors.push_back(values[new_row][new_col]);
            }
        }
        return neighbors;
    }
};

struct GroundFilterTin{
    Grid<Square> grid;
    std::vector<Point> pointcloud;
    std::vector<bool> ground;
    DT dt;

    GroundFilterTin(const std::vector<Point>& pointcloud, double cellsize){
        this->pointcloud = pointcloud;
        this->grid = Grid<Square>(pointcloud, cellsize);
    }

    void init_triangulation(){
        std::vector<Point_2d> points_2d;
        std::vector<int> indices;
        for (int i = 0; i < pointcloud.size(); i++){
            auto p = pointcloud[i];
            points_2d.emplace_back(Point_2d(p.x(), p.y()));
            indices.push_back(i);

            //mark all as non-ground
            ground.push_back(false);
        }

        // create a tree, query by fuzzy_iso_box or ball
        // tree should be 2 dimensional, because it will be searched with each square of the grid in 2d
        Tree tree(boost::make_zip_iterator(boost::make_tuple( points_2d.begin(),indices.begin())),
                  boost::make_zip_iterator(boost::make_tuple( points_2d.end(),indices.end())));

        // iterate grid of squares and search the pointcloud with each square and select the lowest point
        for (const auto &row: grid.values){
            for (const auto &square: row){
//                std::cout << square << std::endl;
                auto p1 = Point_2d(square.min_x, square.min_y);
                auto p2 = Point_2d(square.max_x, square.max_y);
                Fuzzy_iso_box search_box(p1, p2, 0.1);
                std::vector<boost::tuples::tuple<CGAL::Point_2<CGAL::Epick>, int>> points_in_square;
                tree.search(std::insert_iterator<std::vector<boost::tuples::tuple<CGAL::Point_2<CGAL::Epick>, int>>>(points_in_square, points_in_square.end()), search_box);

                int lowest_point_idx;
                double current_z = std::numeric_limits<double>::infinity();
                for (auto it: points_in_square){
                    int point_idx = boost::get<1>(it);
                    auto current_point = pointcloud[point_idx];
                    if (current_point.z() < current_z){
                        current_z = current_point.z();
                        lowest_point_idx = point_idx;
                    }
                }

                // insert into DT and mark as ground
                if (not ground[lowest_point_idx]) {
                    dt.insert(pointcloud[lowest_point_idx]);
                    ground[lowest_point_idx] = true;
                }
            }
        }

        //calculate convex hull
        std::vector<int> convex_hull;
        CGAL::convex_hull_2(indices.begin(), indices.end(), std::back_inserter(convex_hull),
                            Convex_hull_traits_2(CGAL::make_property_map(points_2d)));

        // construct a 2d tree of initial ground points
        std::vector<Point_2d> initial_ground;
        std::vector<int> ground_indices;
        for (int i = 0; i < ground.size(); i++){
            if (ground[i]){
                auto p = pointcloud[i];
                initial_ground.emplace_back(Point_2d(p.x(), p.y()));
                ground_indices.push_back(i);
            }
        }
        Tree2 tree_g(boost::make_zip_iterator(boost::make_tuple( initial_ground.begin(),ground_indices.begin())),
                   boost::make_zip_iterator(boost::make_tuple( initial_ground.end(),ground_indices.end())));

        // add convex hull points to ground with the height of the nearest ground point.
        // this is done to make sure no tested point will lie outside of any of the triangles
        for (auto idx: convex_hull){
            if (not ground[idx]) {
                // find the closest ground point and insert the convex hull point with its height
                auto p = pointcloud[idx];
                auto p2d = Point_2d(p.x(), p.y());
                Neighbor_search search(tree_g, p2d, 1);

                // height of the nearest point
                auto z = pointcloud[boost::get<1>(search.begin()->first)].z();
                dt.insert(Point(p.x(), p.y(), z));
            }
        }
    }

    std::vector<int> refine_tin(double distance, double angle){
        int counter = 1;
        // stop when no new points ground points were detected
        while (counter > 0 ) {
            counter = 0;
            for (unsigned int i = 0; i < ground.size(); i++) {
                if (not ground[i]) {
                    auto tested_point = pointcloud[i];

                    // Find triangle that intersects the given point
                    DT::Face_handle triangle = dt.locate(tested_point);
                    // get the 3 vertices of the triangle:
                    DT::Vertex_handle v0 = triangle->vertex(0);
                    DT::Vertex_handle v1 = triangle->vertex(1);
                    DT::Vertex_handle v2 = triangle->vertex(2);
                    Point vertices[] = {v0->point(), v1->point(), v2->point()};

                    // calculate d
                    Plane plane = Plane(vertices[0], vertices[1], vertices[2]);
                    double d = sqrt(CGAL::squared_distance(plane, tested_point));

                    // calculate alpha
                    double max_alpha = -1.0;
                    for (const auto &vertex_point: vertices) {
                        //dist from vertex to p
                        double d_vp = sqrt(CGAL::squared_distance(vertex_point, tested_point));
                        double alpha = asin(d / d_vp) * (180 / M_PI);
                        if (alpha > max_alpha) {
                            max_alpha = alpha;
                        }
                    }

                    // ground test
                    if (d < distance && max_alpha < angle) {
                        ground[i] = true;
                        dt.insert(tested_point);
                        counter++;
                    }
                }
            }
        }

        std::vector<int> result;
        for (auto is_ground : ground){
            result.push_back(is_ground ? 2 : 1);
        }

        return result;
    }
};

// Particle of the cloth
struct Particle{
    double x;
    double y;
    double z_min;
    double z_prev;
    double z_cur;
    bool moved = false;

    //initial displacement has to be bigger than epsilon_zmax
    constexpr static const double displacement = 0.1;

    Particle(){;}

    Particle(double x, double y, double z, double z_min){
        this->x = x;
        this->y = y;
        this->z_cur = z;
        z_prev = z + displacement;
        this->z_min = z_min;
    }

    Point point() const{
        return {x, y, z_cur};
    }

    Point_2d point2d() const{
        return {x, y};
    }

    bool movable() {
        return z_cur > z_min;
    }

    // update by gravity
    void update(int sign = -1){
        // default gravity vector downwards
        auto tmp = z_cur;
        z_cur += sign * fabs(z_cur - z_prev);
        z_prev = tmp;

        if (z_cur < z_min)
            z_cur = z_min;
    }

    // update by internal forces
    void update(int sign, double magnitude){
        z_cur += sign * magnitude;

        if (z_cur < z_min)
            z_cur = z_min;
    }
};

std::ostream& operator<<(std::ostream& os, const Particle& s) {
    os << "Particle {x=" << s.x << ", y=" << s.y << ", z_prev=" << s.z_prev << ", z_cur=" << s.z_cur << ", z_min=" << s.z_min << "}";
    return os;
}

// add particle to grid
template<>
void add_cell<Particle>(std::vector<Particle>& grid_row, double x, double y, double cellsize, double height){
    grid_row.emplace_back(Particle(x, y, height, 0));
}

struct GroundFilterCSF{
    Grid<Particle> grid;
    std::vector<Point> pointcloud;
    std::vector<bool> ground;

    GroundFilterCSF(const std::vector<Point>& pointcloud, double cellsize){
        std::vector<Point_2d> points_2d;
        std::vector<unsigned int> indices;

        for (unsigned int i = 0; i < pointcloud.size(); i++){
            auto p = pointcloud[i];

            // copy pointcloud with inverted z
            this->pointcloud.emplace_back(p.x(), p.y(), -p.z());
            points_2d.emplace_back(Point_2d(p.x(), p.y()));
            indices.push_back(i);

            //mark all as non-ground
            ground.push_back(false);
        }

        auto max_z_point = std::max_element(this->pointcloud.begin(), this->pointcloud.end(), z_comparator);
        auto z0 = max_z_point->z() + 1.0;

        // init the cloth grid slightly above the max height of the inverted pointcloud
        grid = Grid<Particle>(pointcloud, cellsize, z0);

        // 2d tree to find the closest point to particle in 2d
        Tree2 tree(boost::make_zip_iterator(boost::make_tuple( points_2d.begin(),indices.begin())),
                   boost::make_zip_iterator(boost::make_tuple( points_2d.end(),indices.end())));

        // for every particle of the grid find its closest point in 2d and set its height as the zmin of the particle
        for (auto &row: grid.values){
            for (auto &particle : row){
                auto p = particle.point2d();
                Neighbor_search search(tree, p, 1);

                auto z_min = this->pointcloud[boost::get<1>(search.begin()->first)].z();
                particle.z_min = z_min;
            }
        }
    }

    std::vector<int> run(double epsilon_zmax, double epsilon_ground, double rigidity, int n_iterations=5000) {
        // scale is 1 if rigidity is 2 or 3 because then the update is done 2 or 3 times so the vector is not actually scaled
        double scale = rigidity <= 1 ? rigidity : 1;
        unsigned int counter = 0;
        while (counter < n_iterations) {
            counter++;
//            std::cout << counter << std::endl;

            //external forces
            for (auto &row: grid.values) {
                for (Particle &p: row) {
                    if (p.movable()) {
                        p.update();
                    }
                }
            }

            //internal forces
            for (unsigned int r = 0; r < grid.nrows; r++) {
                for (unsigned int c = 0; c < grid.ncols; c++) {
                    Particle &p = grid(r, c);
                    if (!p.movable())
                        continue;

                    for (auto &neighbor: grid.get_neighbors(r, c)) {
                        int direction = p.z_cur > neighbor.z_cur ? -1 : 1;
                        double magnitude = scale * fabs(p.z_cur - neighbor.z_cur) / 2.0;

                        if (neighbor.movable()) {
                            neighbor.update(-direction, magnitude);
                            p.update(direction, magnitude);
                        }
                        else{
                            // move the particle N times
                            for (int j = 0; j < ceil(rigidity); j++){
                                p.update(direction, magnitude);
                                magnitude = fabs(p.z_cur - neighbor.z_cur) / 2.0;
                            }
                        }
                    }
                }
            }

            double delta_z = 0;
            Particle max_p;
            for (auto &row: grid.values) {
                for (Particle &p: row) {
                    if (fabs(p.z_cur - p.z_prev) > delta_z && p.movable()) {
                        delta_z = fabs(p.z_cur - p.z_prev);
                        max_p = p;
                    }
                }
            }

//            std::cout<< "delta: " << delta_z <<std::endl;
//            std::cout<<max_p<<std::endl;
            if (delta_z <= epsilon_zmax)
                break;
        }

        // make a kdtree of cloth points
        std::vector<Point> cloth_points;
        for (auto& row: grid.values)
            for (auto& p: row){
                cloth_points.emplace_back(p.point());
            }

        Tree3 tree(cloth_points.begin(), cloth_points.end());

        // find the distance of each pointcloud point to the nearest cloth point
        std::vector<int> result;
        for (auto &point: pointcloud){
            Neighbor_search3 search(tree, point, 1);
            auto dist = std::sqrt(search.begin()->second);
            result.push_back(dist < epsilon_ground ? 2 : 1);
        }

        return result;
    }
};

void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams) {
    /*
      Inputs:
        pointcloud: input point cloud (an Nx3 numpy array),
        jparams: a dictionary jparams with all the parameters that are to be used in this function:
          - resolution:    resolution (cellsize) for the initial grid that is computed as part of the ground filtering algorithm,
          - distance:      distance threshold used in the ground filtering algorithm,
          - angle:         angle threshold used in the ground filtering algorithm in degrees,
          - output_las:    path to output .las file that contains your ground classification,
    */
    double distance = jparams["distance"];
    double angle = jparams["angle"];
    double res = jparams["resolution"];

    auto gf = GroundFilterTin(pointcloud, res);
    gf.init_triangulation();
    std::vector<int> class_labels = gf.refine_tin(distance, angle);
    write_lasfile(jparams["output_las"], pointcloud, class_labels);

    std::vector<Point> grid;
    std::vector<int> grid_labels;
    for (auto &row : gf.grid.values)
        for (auto &s : row){
            grid.emplace_back(Point(s.min_x, s.min_y, 0));
            grid.emplace_back(Point(s.min_x, s.max_y, 0));
            grid.emplace_back(Point(s.max_x, s.min_y, 0));
            grid.emplace_back(Point(s.max_x, s.max_y, 0));
            grid_labels.push_back(3);
            grid_labels.push_back(3);
            grid_labels.push_back(3);
            grid_labels.push_back(3);
//            break;
        }
    write_lasfile("grid.laz", grid, grid_labels);
}

void groundfilter_csf(const std::vector<Point>& pointcloud, const json& jparams) {
  /*
  Inputs:
    pointcloud: input point cloud (an Nx3 numpy array),
    jparams: a dictionary with all the parameters that are to be used in this function:
      - resolution:     resolution of the cloth grid,
      - epsilon_zmax:   tolerance to stop the iterations,
      - epsilon_ground: threshold used to classify ground points,
      - output_las:     path to output .las file that contains your ground classification
  */

    std::cout << std::setprecision(15);
   double resolution = jparams["resolution"];
   double epsilon_zmax = jparams["epsilon_zmax"];
   double epsilon_ground = jparams["epsilon_ground"];
   std::string output_las = jparams["output_las"];

   auto gf = GroundFilterCSF(pointcloud, resolution);
   auto class_labels = gf.run(epsilon_zmax, epsilon_ground, 0.3, 5000);

   std::vector<Point> cloth;
   std::vector<int> cloth_labels;
   for (auto &row : gf.grid.values)
       for (auto &particle : row){
           cloth.emplace_back(Point(particle.x, particle.y, -particle.z_cur));
           cloth_labels.push_back(3);
       }

   write_lasfile(jparams["output_las"], pointcloud, class_labels);
   write_lasfile("cloth.laz", cloth, cloth_labels);
}



std::vector<Point> read_lasfile(const json& jparams) {
  /*
  Function to read points from a LAS file

  Inputs:
    jparams["filename"]:   the filename to read the LAS file to

  Returns:
    a std::vector<Point> with the points from the LAS file
  */
  std::string filename = jparams["filename"];
	LASreadOpener lasreadopener;
	lasreadopener.set_file_name(filename.c_str());
	LASreader* lasreader = lasreadopener.open();
	
	if (!lasreader){
		std::cerr << "cannot read las file: " << filename << "\n";
		exit(1);
	}

  //-- store each point in a CGAL Point_3 object
  //-- https://doc.cgal.org/latest/Kernel_23/classCGAL_1_1Point__3.html
	std::vector<Point> points;
	while (lasreader->read_point()) {
		points.push_back( 
			Point(
				lasreader->point.get_x(),
				lasreader->point.get_y(),
				lasreader->point.get_z()
			)
		);
	}
	lasreader->close();
	delete lasreader;

	return points;
}



void write_lasfile(const std::string filename, const std::vector<Point>& pointcloud, const std::vector<int>& class_labels) {
  /*
  Function to write a new LAS file with point labels (for the LAS classification field)

  Inputs:
    filename:   the filename to write the LAS file to
    pointcloud: input point cloud (a vector of Points),
    Labels:     Contains point labels. Should be a vector of ints of the same size as pointcloud (ie. one label for each point in the same order as pointcloud). Uses LAS classification codes, ie 2 = ground. 1 = unclassified.
  */
  LASwriteOpener laswriteopener;
  laswriteopener.set_file_name(filename.c_str());

  LASheader lasheader;
  lasheader.x_scale_factor = 0.01;
  lasheader.y_scale_factor = 0.01;
  lasheader.z_scale_factor = 0.01;
  lasheader.x_offset = 0.0;
  lasheader.y_offset = 0.0;
  lasheader.z_offset = 0.0;
  lasheader.point_data_format = 0;
  lasheader.point_data_record_length = 20;

  LASpoint laspoint;
  laspoint.init(&lasheader, lasheader.point_data_format, lasheader.point_data_record_length, 0);

  LASwriter* laswriter = laswriteopener.open(&lasheader);
  if (laswriter == 0)
  {
    std::cerr << "ERROR: could not open laswriter\n";
    exit(1);
  }

	if (pointcloud.size()!=class_labels.size()) {
		std::cerr << "ERROR: points has a different size than class_labels\n";
		exit(1);
	}

  for (size_t i=0; i<pointcloud.size(); ++i) {
		const Point& p = pointcloud[i];
		const int& label = class_labels[i];

    laspoint.set_x(p[0]);
    laspoint.set_y(p[1]);
    laspoint.set_z(p[2]);
		laspoint.set_classification(label);

    laswriter->write_point(&laspoint);
    laswriter->update_inventory(&laspoint);    
  } 

  laswriter->update_header(&lasheader, TRUE);
  laswriter->close();
  delete laswriter;
}