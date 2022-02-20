# pointcloud-ground-filtering

![cloth](https://user-images.githubusercontent.com/38667147/154866844-c375f032-5ae4-4656-b151-49a0f83fe217.png)

**Implementation of point cloud ground filtering algorithms**:
 - filtering by TIN refinement
 - filtering by Cloth Simulation Filtering (CSF https://www.cloudcompare.org/doc/wiki/index.php?title=CSF_(plugin))

This implementation uses the knowledge described in the 3D terrain book: https://github.com/tudelft3d/terrainbook

The code runs based on the parameters in params.json file, which specifies:
 - input file (.laz format)
- for TIN refinement:
  - output pointcloud file (.laz format)
  - resolution of the initial grid used for selecting the lowest point in each cell
  - distance [m] and angle[°] for the ground test

- for CSF
  - output pointcloud file (.laz format)
  - resolution of the cloth grid
  - epsilon_zmax - if the maximum displacement of all the cloth particles in an iteration becomes smaller than this value, the cloth stops falling. 
  - epsilon_ground - used for classifying the points as ground if they are within epsilon_ground distance of the cloth in its final position


For more details look at example _params.json_

