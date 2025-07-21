# filament-reconstuction-tool
A tool for reconstructing filamentous structures considering a Z-stack from a monocular light microscope as input.

A repository with the code and data for reproducing the results in: 3D graph reconstruction of filamentous structures from z-stack images.
Article authors: Oscar Sten, Emanuela del Dottore, Marilena Ronzan, 
Nicola Pugno, and Barbara Mazzolai.

[Open access article](https://doi.org/10.1016/j.ecoinf.2025.103317)

Code author: Oscar Sten, unless otherwise indicated.


Src:
Contains:

``mycelium_3D_reconstruction.m`` for reproducing the results in Figure 7. 

``overlap_detection_in_2D.m`` for reproducing the results in Table 3.

``printed_filament_3D_reconstruction.m`` and ``printed_filaments_3D_reconstruction_more_complex.m`` for reproducing the results in Figure 6 and Table 2. The first script is for Prototypes #1-2, and the second for Prototypes #3-5. The distinction is due to differences in image aquisition setup.

``qualitative_results.m`` for reproducing the results in Figures 8 and 9.

``time_lapse_plotting.m`` for reproducing the results in Figure 10

``fcn/`` contains all necessary functions.

``fcn/external/`` functions originally authored by external authors (see each specific function).

 

### Result_tables:
The excel spreadsheet contains the data and computations in Table 3.


### Data:
Contains all data including images and annotations.

### GT_overlap_tests
Contains the annotation files for reproducing the results in Table 3.

Copyright (2025) Oscar Sten.

### Image_Z_stacks_fungi
Contains Z-stack images of samples of the fungus Rhizophagus Irregolaris
used to reproduce the results in Figures 7, 8, and 9.

Copyright (2025) Oscar Sten.

### Images_2D_fungi
Contains 2D images of samples of the fungus Rhizophagus Irregolaris.

Original publication: Cardini, A., Pellegrino, E., Del Dottore, E. et al. 
HyLength: a semi-automated digital image analysis tool for measuring the 
length of roots and fungal hyphae of dense mycelia. Mycorrhiza 30, 
229–242 (2020). https://doi.org/10.1007/s00572-020-00956-w 

Copyright (2020) Cardini, A., Pellegrino, E., Del Dottore, E. et al.

### Prototype(1-5)
Contains the images of Prototypes used to reproduce the results in  Figure 6 and Table 2.

Copyright (2025) Oscar Sten.

### Imgs4TimeLapse

Contains images for reproducing the graphs in Figure 10.

### CAD-files

Contains the stl-files used for fabricating the 3D-printed prototypes.

 

### Remarks

The branch/crossing detection and branch/crossing matching tests are not fully automated, since a couple of cases were too complex to annotate. In these cases the person performing the test must determine visually if the algorithm's assessment is correct or not.

In image "o.jpg", a filament in the top right corner traverses multiple structures and must be assessed maunally.  Running overlap_detection_in_2D.m. with the parameters reported in the manuscript makes this crossing not detected, so one
 "false negative" must be added.


The reader who is currious about the functionaliy of the graph merging algorithm many run the function mergeGraphs.m in debug mode. Just open the function and set debug_mode = 1 instead of debug_mode = 0.
