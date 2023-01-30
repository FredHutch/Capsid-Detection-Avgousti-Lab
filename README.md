# Capsid-Detection-Avgousti-Lab
detect and classify capsids in Electron Microscopy images

This repository contains MATLAB scripts used by the Avgousti lab to count and classify HSV capsids in nuclei of infected cells from Electron Microscopy images. See related manuscript for details about image acquisition and experimental set up.

The pipeline follows three steps: 1- nuclear boundaries identification, 2- capsid detection, 3- capsid classification. Nuclear boundaries are identified with user input, by outlining a freehand contour of the nucleus of interest, from which a binary mask is extracted. Capsid detection is performed on the median-filtered complement of the original image, using a Circular Hough Transform based algorithm, with phase coding for radii estimation, and search radius ranging from 10.6 to 31.7 nm. Capsids residing within the nuclear mask are then counted and classified. Detected capsids are classified in three categories (Empty, Intermediate and Full), depending on the distribution of the pixel grayscale intensities within the capsids relative to a normal distribution. Moreover, the distance of the capsids from the nearest nuclear membrane pixel is measured using a distance transform, and capsids within a 200nm range from the nuclear membrane are counted.

REQUIREMENTS:

-MATLAB R2020a or later

-Image Processing ToolBox

  
