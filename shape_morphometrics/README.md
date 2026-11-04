# Cell shape and orientation analysis
This code is associated with the paper "Claudin-dependent clonal expansion couples epicardial activation to regenerative fate acquisition".

The script shape_morphometrics.py measures the area, perimeter, orientation, aspect ratio and other shape descriptors of every cell in a segmented label image.

# Requirements
Python 3.11.15 (expected to work on >= 3.9)

Packages:
    numpy
    scipy
    scikit-image
    matplotlib
    pandas
    openpyxl

Install them with:
    pip install numpy scipy scikit-image matplotlib pandas openpyxl

# Input
A segmented label image (TIFF) exported from ilastik, where each cell has its own value. See the example image veh_R_no1.tif

# How to use
To analyse one image, give the file name and the pixel size in micrometres:

    python shape_morphometrics.py image.tif 0.375

To analyse all .tif and .tiff files in the same folder as the script:

    python shape_morphometrics.py 0.375

If no pixel size is given, the value of PIXEL_SIZE_UM at the top of the script is used (default 0.375).

# Output
For each image, a folder named output_ followed by the image name is created next to the image, for example output_image for image.tif. It contains:

    cell_metrics.xlsx            measurements, one row per cell
    01_cell_area_map.png         cells coloured by area
    02_cell_area_labels.png      area and cell ID written on each cell
    03_small_cells_orange.png    cells smaller than SMALL_AREA_UM2 (default 250 µm²) in orange
    04_small_cells_magma.png     cells smaller than SMALL_AREA_UM2 coloured by area
    05_orientation_field.png     long axis of each cell drawn as a line
    06_orientation_rose.png      rose histogram of cell orientations (10 degree bins)

# Measurements in cell_metrics.xlsx

    Cell_ID           cell number
    Area_um2          cell area in µm²
    Perimeter_um      perimeter in µm, Crofton estimate (skimage.measure.perimeter_crofton)
    Feret_Max_um      maximum Feret diameter in µm
    Circularity       4 x pi x area / perimeter squared
    Shape_Index_q     perimeter / square root of area
    Eccentricity      eccentricity of the ellipse with the same second moments as the cell
    Aspect_Ratio      major axis length / minor axis length of that ellipse (empty if the minor axis is close to 0)
    Solidity          cell area / convex hull area
    Orientation_deg   angle of the cell long axis, from -90 to 90 degrees
    Centroid_X_um     cell centre, distance from the left edge of the image
    Centroid_Y_um     cell centre, distance from the top edge of the image
    N_Neighbours      number of cells touching this cell
    Free_Edge         True if part of the cell border touches background instead of another cell

Orientation is calculated from the second moments of the cell outline about the cell centre. It is given in degrees counter-clockwise from the x axis, with the y axis pointing up. 0 means the cell is horizontal, and -90 or 90 means it is vertical.

Other settings can be changed at the top of the script. Each setting is explained in a comment in the code.


# Orientation statistics
The script orientation_stats.py tests whether cells in each field of view (FOV) are aligned, and compares the proportion of aligned FOVs between conditions.

# Input
An Excel file named Orientation.xlsx in the folder you run the script from. Each column is one FOV, and contains the cell orientation angles of that FOV in degrees.

The column names are used to group the FOVs:

    Condition   the name must contain naive, veh (vehicle) or ccpe (cCpE), upper or lower case
    Animal      FOVs from the left and right side of the same animal, or repeated fields of the same animal, are grouped as one animal

# How to use

    python orientation_stats.py


# Contact
Tuba Majid, tuba.majid@ki.se
