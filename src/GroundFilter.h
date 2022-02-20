#include <string>
#include <vector>

// -- json parsing
#include <json.hpp>

// -- CGAL kernel
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/property_map.h>

using json = nlohmann::json;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Point_2 Point_2d;

#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Plane_3.h>
typedef CGAL::Plane_3<Kernel> Plane;
typedef CGAL::Search_traits_2<Kernel>  Traits;
typedef boost::tuple<Point_2d,int> Point_and_int;
typedef CGAL::Search_traits_adapter<Point_and_int,
        CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
        Traits>   Traits_derived;
typedef CGAL::Kd_tree<Traits_derived> Tree;
typedef CGAL::Fuzzy_iso_box<Traits_derived> Fuzzy_iso_box;


typedef CGAL::Orthogonal_k_neighbor_search<Traits_derived> Neighbor_search;
typedef Neighbor_search::Tree Tree2;
typedef Neighbor_search::Distance Distance;


typedef CGAL::Search_traits_3<Kernel> Traits3;
typedef CGAL::Orthogonal_k_neighbor_search<Traits3> Neighbor_search3;
typedef Neighbor_search3::Tree Tree3;
typedef Neighbor_search3::Distance Distance3;

typedef CGAL::Convex_hull_traits_adapter_2<Kernel,
        CGAL::Pointer_property_map<Point_2d>::type > Convex_hull_traits_2;

std::vector<Point> read_lasfile(const json& jparams);
void write_lasfile(const std::string filename, const std::vector<Point>& pointcloud, const std::vector<int>& class_labels);

void groundfilter_tin(const std::vector<Point>& pointcloud, const json& jparams);
void groundfilter_csf(const std::vector<Point>& pointcloud, const json& jparams);