# 2D-NMR-Stacking-Bruker

This GUI is made to assist 2D characterisation, by aligning multiple 2D spectra via a common F1 dimension.

Choose file folder from the data browser e.g. "usr/NMR_data/Sample_A/1" for the F1 1D spectrum, select nucleus and tweak axis limits as desired, then click "Load F1".
Do the same process for the F2 dimension.
Choose the 2D data using the file browser in the 2D panel, and choose number of contour levels, and the threshold for the contours (as % of maximum signal).

Once the F1, F2, and 2D data are all input, click plot first and the first 2D spectrum will be plotted, if you change the F2 and 2D inputs and click "Plot Next", the next 2D spectrum will be plotted, aligned via the F1 axis.

The 1D data can also have a magnified section plotted above the spectrum using the "Plot Magnified Section" button, this uses 3 mouse-click inputs: one to select the line, one to select the start point, and the third to select the end point.

The save section can be used to save a pdf image of the spectra to the specified directory (default is home directory). To change this, input the desired folder path in line 250.
The default paths for the data browsers can also be set manually in lines 48, 89, and 108.
