# Junctional protein quantification
This code is associated with the paper "Claudin-dependent clonal expansion couples epicardial activation to regenerative fate acquisition".

It contains two Python scripts:
    make_skeletons.py        makes 1 px junction skeletons from ilastik Multicut label images
    junctional_protein.py    measures junctional protein intensity using the skeleton and the matching fluorescence image

# Requirements
Python 3.11.15 (expected to work on >= 3.9)

Packages:
    numpy
    scipy
    scikit-image
    matplotlib
    pandas
    tifffile

Install them with:
    pip install numpy scipy scikit-image matplotlib pandas tifffile

# How to use
Step 1. Segment the cells in ilastik using Multicut (see paper for details) and export the segmented image as a TIFF file.

Step 2. Make the skeleton. Put the label TIFF files in a folder together with make_skeletons.py, open in that folder and run:
    python make_skeletons.py
This makes a skeleton for every .tif and .tiff file in the folder, so the folder should contain only the label images. To run it on one file only use:
    python make_skeletons.py image.tif
For each file, a new file is saved with _skeleton added to the name, for example image.tif gives image_skeleton.tif. If you run the script again, remove the _skeleton.tif files from the folder first, otherwise they are processed too.

Step 3. Measure the protein intensity. Open junctional_protein.py and change these lines at the top of the script:
    SKELETON_PATH          the skeleton file from step 2
    PROTEIN_PATH           the fluorescence image (TIFF) of the same field
    BACKGROUND_INTENSITY   background value subtracted from each junction

Then run:
    python junctional_protein.py

Note: Cells touching the image border are excluded. See the image veh_R_no1.tif to see what my input image this code is build for looks like. 

# Output
A folder named after the skeleton file is created. It contains:
    junction_intensities.csv     results for each junction
    cell_intensities.csv         results for each cell
    junction_heatmap_bgsub.png   junctions coloured by background-subtracted intensity
    cell_heatmap_bgsub.png       cells coloured by background-subtracted mean intensity
    cv.png                       cells coloured by CV
    QC.png                       fluorescence image with junction and cell IDs, for checking

Other settings can be changed at the top of each script. Each setting is explained in a comment in the code.

# Contact
Tuba Majid, tuba.majid@ki.se
