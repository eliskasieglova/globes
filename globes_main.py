import arcpy
import numpy as np
import geopandas as gpd
import os
from pathlib import Path
from map_projections import projs_12, projs_6
from globe_definitions import globe_12, globe_6
import math
import sys

# input when running from terminal
globe_type = str(sys.argv[1])  # '12' or '6'

# USER DEFINED VARIABLES
display = False
dpi = 300
basemap = "Charted Territory Map"  # choose map to be displayed on globe

# Set workspace - SET YOUR OWN WORKSPACE
arcpy.env.workspace = 'C:/path/to/your/workspace'
arcpy.env.overwriteOutput = 1

# Create Project - MAKE PROJECT IN ARCGIS PRO CALLED "project.aprx" AND SAVE IN CURRENT DIRECTORY
project_path = "project.aprx"
aprx = arcpy.mp.ArcGISProject(project_path)

# CREATE DIRECTORIES
if not Path('results').exists():
    os.mkdir('results')

if not Path('weloveshapefiles').exists():
    os.mkdir('weloveshapefiles')

# LAYOUT SETTING
# set page width and page height (mm)
if globe_type == '12':
    pw = 430
    ph = 450

elif globe_type == '6':
    pw = 430
    ph = 300
edge_length = 75

# Create Layout
layout = aprx.createLayout(pw, ph, 'MILLIMETER')
layout.name = "Layout"

# Globe definitions
if globe_type == '12':
    globe = globe_12
    projs = projs_12
elif globe_type == '6':
    globe = globe_6
    projs = projs_6

if globe_type == '12':
    # antarctica, arctic, africa, australia, australia2, oceania, south america, africa, asia, japan, pacific, america
    map_rotations = [-180, -216, 144, -288, 0, -72, -144, -216, -144, -72, 0, 72]
    frame_rotations = [-180, -216, 144, -288, 0, -72, -144, -216, -144, -72, 0, 72]
    zooms = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    x_biases = [271, 202, 242.8, 185, 190.1, 251, 283.65, 247.5, 116, 30.4, 108.65, 242.7]
    y_biases = [370, 159, 283.8, 307.9, 370.2, 384.5, 331.3, 229.5, 277, 167, 51.5, 90.2]

elif globe_type == '6':
    anchor_point_x = pw * 1/4
    anchor_point_y = ph * 1/2
    # antarctica, arctic, africa, asia/australia, oceania, america
    map_rotations = [225, 0, 0, 0, 0, 135]
    frame_rotations = [0, 0, 0, 0, 0, 0]
    zooms = [0.7, 1, 1, 1, 1, 0.7]
    x_biases = [anchor_point_x + edge_length, anchor_point_x, anchor_point_x + edge_length, anchor_point_x + edge_length * 2, anchor_point_x + edge_length * 3, anchor_point_x + edge_length]
    y_biases = [anchor_point_y + edge_length, anchor_point_y, anchor_point_y, anchor_point_y, anchor_point_y, anchor_point_y - edge_length]


# FUNCTIONS
def mapframe_vertices(x_center, y_center, edge_length, theta_0):
    """
    Create vertices of map frame.

    :param x_center: x-coord of center point of pentagon
    :param y_center: y-coord of center point of pentagon
    :param edge_length: length of one edge
    :param theta_0: rotation

    :return: list of tuples (x,y)-coordinates of each vertex
    """

    # Calculate the circumradius R
    R = edge_length / (2 * math.sin(math.pi / 5))

    # Calculate the vertices
    vertices = []
    for i in range(5):
        theta = theta_0 + i * 2 * math.pi / 5
        x = x_center + R * math.cos(theta)
        y = y_center + R * math.sin(theta)
        vertices.append((x, y))

    return vertices


def uv_to_sd(u, v, uk, vk):
    """
    Transforms geographic coordinates (u, v) to projected coordinates (s, d).

    :param u: longitude of input point
    :param v: latitude of input point
    :param uk: longitude of cartographic pole
    :param vk: latitude of cartographic pole

    :return: Input points u, v in projected coordinates as s, d.
    """
    dv = vk - v

    # Transformed latitude
    s = np.arcsin(np.sin(u) * np.sin(uk) + np.cos(u) * np.cos(uk) * np.cos(dv))

    # Transformed longitude
    d = -1 * np.arctan2(np.cos(u) * np.sin(dv), np.cos(u) * np.sin(uk) * np.cos(dv) - np.sin(u) * np.cos(uk))

    return s, d


def gnom(R, s, d):
    """
    Projects point using the gnomonic projection.

    :param R: radius of the globe,
    :param s: longitude of point,
    :param d: latitude of point

    :return: projected x-, y-c oordinates to gnomonic projection
    """

    # defined projection functions
    x = R * np.tan(np.pi / 2 - s) * np.cos(d)
    y = R * np.tan(np.pi / 2 - s) * np.sin(d)

    return x, y


def graticule(u_min, u_max, v_min, v_max, D_u, D_v, d_u, d_v, R, uk, vk, u0, proj):
    """
    Creates a graticule of parallels and meridians based on input information in radians.

    :param u_min: the
    :param u_max:
    :param v_min:
    :param v_max:
    :param D_u: step between X-coordinates of meridians/parallels
    :param D_v: step between Y-coordinates of meridians/parallels
    :param d_u: step between points to save (level of detail) for each meridian/parallel
                if too large --> lines can be pointy
                if too small --> files will be large
    :param d_v: step between points to save (level of detail) for each meridian/parallel
                if large --> lines can be pointy
                if small --> files will be large
    :param R:
    :param uk: X-coordinate of cartographic pole
    :param vk: Y-coordinate of cartographic pole
    :param u0: not used
    :param proj: function of projection (gnom)

    Returns lists of X and Y coordinates of meridians and parallels in gnomonic projection
    in the following order:
    XM (X-coordinates of meridians)
    YM (Y-coordinates of meridians)
    XP (X-coordinates of parallels)
    YP (Y-coordinates of parallels)
    """

    # Create meridians
    XM = []
    YM = []

    for v in np.arange(v_min, v_max + D_v, D_v):
        # Meridian
        um = np.arange(u_min, u_max + d_u, d_u)
        m = len(um)
        vm = np.repeat(v, m)

        # Convert to oblique aspect
        sm, dm = uv_to_sd(um, vm, uk, vk)

        # Project a meridian
        xm, ym = proj(R, sm, dm)

        # Add meridian
        XM.extend(xm)
        YM.extend(ym)

    # Create parallels
    XP = []
    YP = []
    for u in np.arange(u_min, u_max + D_u, D_u):
        # Longitude of parallel points
        vp = np.arange(v_min, v_max + d_v, d_v)
        n = len(vp)
        up = np.repeat(u, n)

        # Convert to oblique aspect
        sp, dp = uv_to_sd(up, vp, uk, vk)

        # Project a parallel
        xp, yp = proj(R, sp, dp)

        # Add parallel
        XP.extend(xp)
        YP.extend(yp)

    return XM, YM, XP, YP


def boundary(u, v, R, uk, vk, u0, proj):
    """
    Creates a boundary of globe face in gnomonic projection coordinates
    based on input information in radians.

    :param u: longitude
    :param v: latitude
    :param R: radius
    :param uk: longitude of cartographic pole
    :param vk: latitude of cartographic pole
    :param u0: not used
    :param proj: projection function (gnom)

    :return: Lists of x-coordinates and y-coordinates in gnomonic projection.
    """

    # initiate empty lists
    XB, YB = [], []

    # loop through points
    for u, v in zip(u, v):
        # convert them to projected coordinates
        s, d = uv_to_sd(u, v, uk, vk)
        dv = vk - v
        # project to gnomonic projection
        XB_, YB_ = proj(R, s, d)

        # append converted points to lists
        XB.append(XB_)
        YB.append(YB_)

    return XB, YB


def rotate(X, Y, angle):
    """
    Rotate input points by given angle.

    :param X: list of X-coordinates
    :param Y: list of Y-coordinates
    :param angle: angle to rotate by

    :return: new value for data in lists X, Y rotated by given angle
    """

    # Convert angle to radians
    angle_rad = np.radians(angle)

    # Define the rotation matrix
    rotation_matrix = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                                [np.sin(angle_rad), np.cos(angle_rad)]])

    # Convert boundary points to numpy array for easier manipulation
    boundary_points = np.array([X, Y])

    # Apply the rotation matrix to the boundary points
    rotated_boundary_points = np.dot(rotation_matrix, boundary_points)

    # Extract rotated X and Y coordinates
    rotated_x = rotated_boundary_points[0]
    rotated_y = rotated_boundary_points[1]

    return rotated_x, rotated_y


def splitMeridians(XM, YM):
    """
    Splits the whole meridians shapefile to individual meridians (for plotting so that the ends
    of meridians don't connect with the beginnings of the following meridians).

    :param XM: list of X-coords of meridians
    :param YM: list of Y-coords of meridians

    :return: list containing lists of points for each meridians
    """

    meridians = []
    temp_meridian = []
    xm_current = XM[0]
    for i in range(0, len(XM)):
        if (xm_current != XM[i]) & (i != 0):
            meridians.append(temp_meridian)
            temp_meridian = []
            temp_meridian.append((XM[i], YM[i]))
        else:
            temp_meridian.append((XM[i], YM[i]))
            xm_current = XM[i]

    return meridians


def splitParallels(XP, YP):
    """
    Splits the whole parallels shapefile to individual parallels (for plotting so that the ends
    of parallels don't connect with the beginnings of the following parallels).

    :param XP: list of X-coords of parallels
    :param YP: list of Y-coords of parallels

    :return: list containing lists of points for each parallels
    """

    parallels = []

    temp_parallel = []
    yp_current = YP[0]
    for i in range(0, len(YP)):
        if (yp_current != YP[i]) & (i != 0):
            parallels.append(temp_parallel)
            temp_parallel = []
            temp_parallel.append((XP[i], YP[i]))
        else:
            temp_parallel.append((XP[i], YP[i]))
            yp_current = YP[i]

    return parallels


# CREATE GLOBES
for i in range(int(globe_type)):
    print(i)  # print progress
    # set the variables
    umin = globe['u_min'][i] * np.pi / 180
    umax = globe['u_max'][i] * np.pi / 180
    vmin = globe['v_min'][i] * np.pi / 180
    vmax = globe['v_max'][i] * np.pi / 180
    Du = 10 * np.pi / 180
    Dv = Du
    du = np.pi / 180
    dv = du
    R = 6378000
    uk = globe['u_k'][i] * np.pi / 180
    vk = globe['v_k'][i] * np.pi / 180
    proj = gnom
    u0 = 0
    ub = [u * np.pi / 180 for u in globe['u_b'][i]]
    vb = [v * np.pi / 180 for v in globe['v_b'][i]]

    # create boundary face shapefile if it does not exist
    if not Path(f'temp/rotated_boundary_face{i + 1}{globe}.shp').is_file():
        XB, YB = boundary(ub, vb, R, uk, vk, u0, proj)
        XM, YM, XP, YP = graticule(umin, umax, vmin, vmax, Du, Dv, du, dv, R, uk, vk, u0, proj)

        # Rotate boundary points by -90 degrees
        rotated_XB, rotated_YB = rotate(XB, YB, -90)

        # Create polyline from rotated boundary points
        boundary_pts = arcpy.Array([arcpy.Point(x, y) for x, y in zip(rotated_XB, rotated_YB)])

        # Save the rotated boundary as a shapefile
        output = arcpy.Polyline(boundary_pts)
        arcpy.CopyFeatures_management(output, f'temp/rotated_boundary_face{i + 1}{globe_type}.shp')

        rotated_XM, rotated_YM = rotate(XM, YM, -90)
        rotated_XP, rotated_YP = rotate(XP, YP, -90)

    # create shapefile of meridians if it does not exist
    if not Path(f'temp/meridians_face{i + 1}_split{globe_type}.shp').is_file():
        # generate meridians
        meridians = arcpy.Array()
        for x in range(len(rotated_XM)):
            meridians.add(arcpy.Point(rotated_XM[x], rotated_YM[x]))

        output = arcpy.Polyline(meridians)
        arcpy.CopyFeatures_management(output, f'temp/meridians_face{i + 1}_notsplit{globe_type}.shp')

        arcpy.management.SplitLine(f'temp/meridians_face{i + 1}_notsplit{globe_type}.shp', f'temp/meridians_face{i + 1}_split{globe_type}.shp')

        with arcpy.da.UpdateCursor(f"temp/meridians_face{i + 1}_split{globe_type}.shp", ["FID", "SHAPE@LENGTH"]) as cursor:
            for row in cursor:
                if row[1] > 1000000:
                    cursor.deleteRow()

    # create shapefile of parallels if it does not exist
    if not Path(f'temp/parallels_face{i + 1}_split{globe_type}.shp').is_file():

        # generate parallels
        parallels = arcpy.Array()
        for y in range(len(rotated_XP)):
            parallels.add(arcpy.Point(rotated_XP[y], rotated_YP[y]))

        output = arcpy.Polyline(parallels)
        arcpy.CopyFeatures_management(output, f'temp/parallels_face{i + 1}_notsplit{globe_type}.shp')

        arcpy.management.SplitLine(f'temp/parallels_face{i + 1}_notsplit{globe_type}.shp', f'temp/parallels_face{i + 1}_split{globe_type}.shp')

        with arcpy.da.UpdateCursor(f"temp/parallels_face{i + 1}_split{globe_type}.shp", ["FID", "SHAPE@LENGTH"]) as cursor:
            for row in cursor:
                if row[1] > 1000000:
                    cursor.deleteRow()

    # Create Map
    new_map = aprx.createMap(f"Face{i + 1}")

    # Add basemap territory
    new_map.addBasemap(basemap)

    # set target spatial reference
    spatial_reference = projs[i]

    # reproject meridians to given spatial_reference
    if not Path(f'results/meridians_face{i + 1}{globe_type}.shp').is_file():
        # open and project shapefile with geopandas - arcpy gave us fatalerror on arcpy.Project_management()
        m_gdf = gpd.read_file(f"temp/meridians_face{i + 1}_split{globe_type}.shp").set_crs(spatial_reference)
        m_gdf.to_file(f'results/meridians_face{i + 1}{globe_type}.shp')

    # reproject parallels to given spatial_reference
    if not Path(f'results/parallels_face{i + 1}{globe_type}.shp').is_file():
        p_gdf = gpd.read_file(f"temp/parallels_face{i + 1}_split{globe_type}.shp").set_crs(spatial_reference)
        p_gdf.to_file(f'results/parallels_face{i + 1}{globe_type}.shp')

    # reproject boundary to given spatial_reference
    if not Path(f'results/boundary_face{i + 1}{globe_type}.shp').is_file():
        p_gdf = gpd.read_file(f"temp/rotated_boundary_face{i + 1}{globe_type}.shp").set_crs(spatial_reference)
        p_gdf.to_file(f'results/boundary_face{i + 1}{globe_type}.shp')

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

    # input parameters for counting vertices of map frame
    dx = x_biases[i]  # x-coordinate of center point (set at the beginning of code)
    dy = y_biases[i]  # y-coordinate of center point (set at the beginning of code)
    rot_map = map_rotations[i]
    rot_frame = frame_rotations[i]
    zoom = zooms[i]

    if globe_type == '12':
        # count vertices
        edge_length = 50  # edge length
        center = (dx, dy)
        theta_0 = math.pi / 10  # rotation of map frames for northern hemisphere
        if i > 6:
            theta_0 = 19 * math.pi / 10  # rotation of map frames for southern hemisphere

        # count vertices of map frame
        vertices = mapframe_vertices(center[0], center[1], edge_length, theta_0)

        p1 = (vertices[0])
        p2 = (vertices[1])
        p3 = (vertices[2])
        p4 = (vertices[3])
        p5 = (vertices[4])

        # Set boundary points of map frame
        map_boundary_points = arcpy.Array([
            arcpy.Point(p1[0], p1[1]),
            arcpy.Point(p2[0], p2[1]),
            arcpy.Point(p3[0], p3[1]),
            arcpy.Point(p4[0], p4[1]),
            arcpy.Point(p5[0], p5[1]),
            arcpy.Point(p1[0], p1[1])
        ])

        # Create the polygon of map frame boundary
        map_boundary = arcpy.Polygon(map_boundary_points)

    elif globe_type == '6':
        edge_length = 75  # edge length
        center = (dx, dy)

        # count vertices of map frame
        p1 = (dx - edge_length / 2, dy - edge_length / 2)
        p2 = (dx - edge_length / 2, dy + edge_length / 2)
        p3 = (dx + edge_length / 2, dy + edge_length / 2)
        p4 = (dx + edge_length / 2, dy - edge_length / 2)

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

    # Edit symbology
    if display:
        layer1 = new_map.listLayers('meridians')[0]
        layer2 = new_map.listLayers('parallels')[0]
        layer3 = new_map.listLayers('boundary')[0]

        symbology1 = layer1.symbology
        symbology2 = layer2.symbology
        symbology3 = layer3.symbology

        # Set the symbol color to black
        symbology1.renderer.symbol.color = {'RGB': [0, 0, 0, 20]}
        symbology1.renderer.symbol.width = 1
        symbology2.renderer.symbol.color = {'RGB': [0, 0, 0, 20]}
        symbology2.renderer.symbol.width = 1
        symbology3.renderer.symbol.color = {'RGB': [0, 0, 0, 0]}
        symbology3.renderer.symbol.width = 0

        # Update the layer symbology
        layer1.symbology = symbology1
        layer2.symbology = symbology2
        layer3.symbology = symbology3
    else:
        layer1 = new_map.listLayers('meridians')[0]
        layer2 = new_map.listLayers('parallels')[0]
        layer3 = new_map.listLayers('boundary')[0]

        symbology1 = layer1.symbology
        symbology2 = layer2.symbology
        symbology3 = layer3.symbology

        # Set the symbol color to black
        symbology1.renderer.symbol.color = {'RGB': [0, 0, 0, 0]}
        symbology1.renderer.symbol.width = 1
        symbology2.renderer.symbol.color = {'RGB': [0, 0, 0, 0]}
        symbology2.renderer.symbol.width = 1
        symbology3.renderer.symbol.color = {'RGB': [0, 0, 0, 0]}
        symbology3.renderer.symbol.width = 0

        # Update the layer symbology
        layer1.symbology = symbology1
        layer2.symbology = symbology2
        layer3.symbology = symbology3

    # set projection of background map
    spatial_ref = arcpy.SpatialReference(text=projs[i])
    map_frame.map.spatialReference = spatial_ref

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


# Save the changes to the project
# aprx.save()

# Export the layout to a PDF
print('exporting')
if globe_type == '12':
    globe_name = 'dodecahedron'
elif globe_type == '6':
    globe_name = 'cube'

layout.exportToPDF(f"{globe_name}.pdf", resolution=dpi)

