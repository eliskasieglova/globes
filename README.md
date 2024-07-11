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
- _display_: True or False whether or not the parallels and meridians should be present in the output PDF file,
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

The **globe definitions** contain values necessary for generating parallels and meridians for each of the faces. The values are stored as lists of integers that are fetched based on index (index of face 1 is 0 etc.).The values it contains are based on the supplementary material we obtained from doc. Bayer in the course Mathematical Cartography. The supplementary material can be downloaded here: https://web.natur.cuni.cz/~bayertom/images/courses/mmk/mmk_cv_2_navod.pdf.

![image](https://github.com/eliskasieglova/globes/assets/57137417/3e0a4a40-2e7b-40e7-b0a1-21c0a73dba34 | width = 100)

![image](https://github.com/eliskasieglova/globes/assets/57137417/9aa8f0cd-4a7f-4b54-aead-8d6d5a25c028 | width = 100)

Images from doc. Bayer's materials.

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
Create directories to store temporary data ("/weloveshapefiles") and results ("/results").

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

### Creating Globes
The actual script!!
The script loops through the globe faces (number of iterations based on variable _globe_type_) and fetches the variables from the script _globe_definitions.py_. All the variables are converted to radians. 3 shapefiles are created for each face and saved to the temporary folder _/weloveshapefiles_: _boundary.shp_ (face boundary), _meridians.shp_ and _parallels.shp_. All the shapefiles are named based on the globe face they represent (e.g. _boundary_face1.shp_). This process is cached, therefore if the wanted output file (for example _meridians_face1.shp_) already exists, this computationally expensive part of the code is skipped. This is also where arcpy comes into the game. We initially have all the coordinates saved in lists by x- and y-coordinates (for example _XM_ and _YM_), which we zip into array of arcpy Point objects and then convert those to a polyline such.
```
# Create array of points from lists of coordinates
boundary_pts = arcpy.Array([arcpy.Point(x, y) for x, y in zip(XM, YM)])

# Convert to polyline
output = arcpy.Polyline(boundary_pts)

# Save as shapefile
arcpy.CopyFeatures_management(output, f'temp/meridians_{i + 1}{globe_type}.shp')
```

After this, the parallels and meridians are all in one file and it is necessary to split them to avoid the ends and beginnings from connecting. This was done using the _arcpy.management.SplitLine_ function which splits the line by points, and then using _arcpy.da.UpdateCursor()_ and _cursor.deleteRow()_ the longest segments (= the ones connecting the end of one meridian/parallel and the beginning of another) are removed.

```
arcpy.management.SplitLine(f'temp/meridians_face{i + 1}_notsplit{globe_type}.shp', f'temp/meridians_face{i + 1}_split{globe_type}.shp')

with arcpy.da.UpdateCursor(f"temp/meridians_face{i + 1}_split{globe_type}.shp", ["FID", "SHAPE@LENGTH"]) as cursor:
    for row in cursor:
        if row[1] > 1000000:
            cursor.deleteRow()
```

After we have the shapefiles for the boundary, meridians and parallels of given face, a new map in the ArcGIS Pro project is created, the basemapp is added and the the target projection is set based on the defined projections in the _map_projections.py_ script. Why do we have to create a new map for each face, one might ask? It is because each face has to have its own set projection, and also each face is rotated in a different direction so that the edges connect.

```
# Create Map
new_map = aprx.createMap(f"Face{i + 1}")

# Add basemap territory
new_map.addBasemap(basemap)

# set target spatial reference
spatial_reference = projs[i]
```

Then each of the shapefiles is projected to the defined target spatial reference. There were tries to do this using the _arcpy.Project_management()_ function, however, that returned a fatal error. Because debugging, looking for help online or asking ChatGPT was not much help, the reprojection was done using the Python library _geopandas_. All the above created and edited shapefiles are stored in the temporary folder, the results of the projection operation are saved into the results folder and are plotted on the resulting globe.

```
# reproject meridians to given spatial_reference
if not Path(f'results/meridians_face{i + 1}{globe_type}.shp').is_file():
    # open and project shapefile with geopandas - arcpy gave us fatalerror on arcpy.Project_management()
    m_gdf = gpd.read_file(f"temp/meridians_face{i + 1}_split{globe_type}.shp").set_crs(spatial_reference)
    m_gdf.to_file(f'results/meridians_face{i + 1}{globe_type}.shp')
```

Feature layers are created from these shapefiles and they are added to the map.

```
    # make feature layers from projected shapefiles
    m_path = f"results/meridians_face{i + 1}{globe_type}.shp"
    p_path = f"results/parallels_face{i + 1}{globe_type}.shp"
    b_path = f"results/boundary_face{i + 1}{globe_type}.shp"

    meridians_shp = arcpy.management.MakeFeatureLayer(m_path, "meridians")[0]
    parallels_shp = arcpy.management.MakeFeatureLayer(p_path, "parallels")[0]
    boundary_shp = arcpy.management.MakeFeatureLayer(b_path, "boundary")[0]

    # add meridians and parallels to map
    new_map.addLayer(meridians_shp)
    new_map.addLayer(parallels_shp)
    new_map.addLayer(boundary_shp)
```
Then, the face is placed on the layout based on the vertices defined at the beginning of the script.
```
# Set boundary points of map frame
map_boundary_points = arcpy.Array([
    arcpy.Point(p1[0], p1[1]),
    arcpy.Point(p2[0], p2[1]),
    arcpy.Point(p3[0], p3[1]),
    arcpy.Point(p4[0], p4[1]),
    arcpy.Point(p1[0], p1[1])
])

# Create the polygon of map frame boundary
map_boundary = arcpy.Polygon(map_boundary_points)

# Create the map frame
map_frame = layout.createMapFrame(map_boundary, new_map)
map_frame.map = new_map
map_frame.setAnchor = "CENTER_POINT"
```

Then comes the finishing touches of the given map frame (we are still working face by face). If ```display == True```, the parallels and meridians will be visualized on the globe. For that it is necessary to set the symbology of these layers. The color is set to slightly transparent black with the width 1. If ```display == False``` alpha is turned to 0 making all the layers invisible.

```
layer1 = new_map.listLayers('meridians')[0]
layer2 = new_map.listLayers('parallels')[0]

symbology1 = layer1.symbology
symbology2 = layer2.symbology

# Set the symbol color to black
symbology1.renderer.symbol.color = {'RGB': [0, 0, 0, 20]}
symbology1.renderer.symbol.width = 1
symbology2.renderer.symbol.color = {'RGB': [0, 0, 0, 20]}
symbology2.renderer.symbol.width = 1

# Update the layer symbology
layer1.symbology = symbology1
layer2.symbology = symbology2
```
The projection of the background map is set.
```
# set projection of background map
spatial_ref = arcpy.SpatialReference(text=projs[i])
map_frame.map.spatialReference = spatial_ref
```
The map frame is set to fit the extent of the boundary layer (the extent is extracted using arcpy.Describe()) and then the map frame is rotated based on the predefined rotations. Note that by default, the map frame does not rotate together with the elements inside it, therefore the map element has to be rotated separately. A special handling is needed for the pole faces on the cube - setting the extent made the map element rotated by 45° compared to the rest of the faces, the map element had to be rotated by 45° and zoomed in to fit the map frame and connect to the neighboring faces.

```
# set extent to fit the boundary face
desc = arcpy.Describe(f'results/boundary_face{i + 1}{globe_type}.shp')
extent = desc.extent

# zoom in on extent (for poles on cube)
map_frame.camera.setExtent(arcpy.Extent(extent.XMin * zoom,
                           extent.YMin * zoom,
                           extent.XMax * zoom,
                           extent.YMax * zoom))

# rotate map frame
map_frame.camera.heading = rot_map
map_frame.elementRotation = rot_frame
```
Once this part of the script was run for each of the faces, the layout is exported to PDF, named {globe_name}.pdf which is defined by the selected globe, at a resolution of 300 dpi by default.

## Acknowledgements and Credits
The code for generating a dodecahedron was written together with Anna Brázdová and David Martínek based on code that was written in classes of the Mathematical Cartography course at the Charles University in Prague. This was used as examination for this course. The part for generating a cube was written as an extension to the already existing code for examination in the course Programming for GIS at the Charles University in Prague. 
  
