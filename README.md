# Globes
Creating globes on cube / dodecahedron ready to print using the Python arcpy library. Will not work without an active ArcGIS account.

## Project Description
This project aims to combine GIS, programming, cartography and art to generate globes on two different shapes: a dodecahedron and a cube. The faces of each of these shapes are saved to a pdf file ready to print. It is possible to choose the background map from the ESRI basemap list (described further in the documentation) and whether or not you want to have meridians and parallels visualized on the globe. After printed, the globe can be glued together and used as a Christmas tree decoration, for example.

## Installation Guide & How to Run
The code utilizes the arcpy library, it is necessary to have arcpy installed and running on your computer. The arcpy installation guide can be found here: https://pro.arcgis.com/en/pro-app/latest/arcpy/get-started/installing-arcpy.htm. 

The following libraries need to be installed:
- arcpy
- geopandas
- numpy
- os
- pathlib
- math

Further on, it is necessary to create a blank ArcGIS project called "project.aprx" and save it to your root folder where you will be running this code. Once you have arcpy and all the further mentioned libraries installed, and you created a project, you can proceed to running the code. 

The lines 11-21 in globes_main.py include the variables to be changed by the user. The variables that can be changed are:
- _globe\_type_: '12' or '6',
- _basemap_: defaults to "Charted Territory Map" but feel free to change it from the list of ESRI Basemaps found here: https://www.arcgis.com/home/group.html?id=702026e41f6641fb85da88efe79dc166#overview,
- _arcpy.env.workspace_: absolute path to your workspace,
- _project\_path_: relative path to your project (defaults to "project.aprx" but feel free to change if you have the project saved somewhere else/you named the project differently).

## Theoretical Background
- geographic vs projected coordinates
- gnomonic projection

## Documentation
The project consists of 3 Python files: **globes_main.py**, **globe_definitions.py** and **map_projections.py**.
- **globes_main.py**: the main script containing user defined variables and called functions that produce the final results,
- **globe_definitions.py**: script that contains hardcoded numbers that define the globe faces (vertices of each face, cartographic poles),
- **map_projections**: script that contains dictionaries of defined projections for each globe face.

The **map projections** are in the ArcGIS format. They are gnomonic projections with a different cartographic pole for each of the globe faces and are used for projecting the basemap for each of the globe faces.

The **globe definitions** contain values necessary for generating parallels and meridians for each of the faces. The values are stored as lists of integers that are fetched based on index (index of face 1 is 0 etc.).The values it contains are based on the supplementary material we obtained from doc. Bayer in the course Mathematical Cartography. The supplementary material can be downloaded here: https://web.natur.cuni.cz/~bayertom/images/courses/mmk/mmk_cv_2_navod.pdf. All inputs are in radians.

![image](https://github.com/eliskasieglova/globes/assets/57137417/3e0a4a40-2e7b-40e7-b0a1-21c0a73dba34)
![image](https://github.com/eliskasieglova/globes/assets/57137417/9aa8f0cd-4a7f-4b54-aead-8d6d5a25c028)

- ualpha:  
- u_min, u_max, v_min, v_max: x- and y-coordinates of a bounding box for each face (which area on the globe we should create meridians and parallels for, so that we do not make the whole globe for each face) 
- u_k, v_k: x- and y-coordinate of cartographic pole
- u_b: 
- v_b: 

- us1: (only for dodecahedron),
- uj1: (only for dodecahedron),
- us2: (only for dodecahedron),
- uj2: (only for dodecahedron),
- us: (only for cube),
- uj: -us, (only for cube).

The next section will be a simple walk-through of the code and explanation of its sections.

### 1. Set user defined variables. 
Described above.

### 2. Create directories 
Create directories to store temporary data ("/temp") and results ("/results").

### 3. Layout settings. 
DPI, layout width and layout height. Possible to change DPI if you want to print on larger paper, would not recommend changing width and height as the position of the elements on the paper is absolute.

### 4. Positioning of globe faces on layout. 
The faces are positioned on the layout so that they need as little cutting and glueing as possible in the end. The values of the globe face positions are absolute so that they connect on the correct edges and are turned around accordingly. Values and projections are set based on user-defined variable _globe\type_.

- _map_rotations_, _frame_rotations_: pre-defined rotations of the map frames and the maps inside the map frames for the faces to fit each other on the layout,
- _zooms_: the zoom needs to be adjusted for the pole faces on the cube, all other values = 1,
- _x_biases_, _y_biases_: movements of position of each face on layout,
- _anchor_point_x_, _anchor_point_y_: anchor points for each globe face (only for cube).

### 5. Functions 
This section contains the functions used to produce shapefiles of the face boundaries and the according parallels and meridians in the gnomonic projection specific for each globe face. Each function should have its documentation in the code, but we will go through them quickly here as well.

- _mapframe_vertices()_: returns vertices of the map frame for each globe face based on the center point (anchor point), edge length and rotation,
- _uv_to_sd()_: converts cartographic coordinates to projected coordinates,
- _gnom()_: definition of gnomonic projection,
- _graticule()_: returns lists of coordinates in gnomonic projection for meridians and parallels for given face based on variables from _globe_definitions.py_,
- _boundary()_: returns lists of coordinates in gnomonic projection for face boundary of given face based on variables from _globe_definitions.py_,
- _rotate()_: rotates list of input points by input angle,
- _splitMeridians()_: splits a list of points for all meridians to list of lists separated for each individual meridian,
- _splitParallels()_: splits a list of points for all parallels to list of lists separated for each individual parallel.

## Acknowledgements and Credits
The code for generating a dodecahedron was written together with Anna Brázdová and David Martínek based on code that was written in classes of the Mathematical Cartography course at the Charles University in Prague. This was used as examination for this course. The part for generating a cube was written as an extension to the already existing code for examination in the course Programming for GIS at the Charles University in Prague. 
  
