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


## Documentation
The project consists of 3 Python files: **globes_main.py**, **globe_definitions.py** and **map_projections.py**.
- **globes_main.py**: the main script containing user defined variables and called functions that produce the final results,
- **globe_definitions.py**: script that contains hardcoded numbers that define the globe faces (vertices of each face, cartographic poles)
- **map_projections**: script that contains dictionaries of defined projections for each globe face

The next section will be a simple walk-through of the code and explanation of its sections.

### 1. Set user defined variables. 
Described above.
### 2. Create directories 
Create directories to store temporary data ("/temp") and results ("/results").
### 3. Layout settings. 
DPI, layout width and layout height. Possible to change DPI if you want to print on larger paper, would not recommend changing width and height as the position of the elements on the paper is absolute.
### 4. Globe definitions. 
Values and projections will be fetched based on user-defined variable _globe\type_.

- _map_rotations_, _frame_rotations_:
- _zooms_:
- _x_biases_, _y_biases_: 
### 


## Acknowledgements and Credits
The code for generating a dodecahedron was written together with Anna Brázdová and David Martínek based on code that was written in classes of the Mathematical Cartography course at the Charles University in Prague. This was used as examination for this course. The part for generating a cube was written as an extension to the already existing code for examination in the course Programming for GIS at the Charles University in Prague. 
  
