# Capsid-Detection-Avgousti-Lab
detect and classify capsids in Electron Microscopy images

This repository contains MATLAB scripts used by the Avgousti lab to count HSV capsids in nuclei of infected cells from Electron Microscopy images. It also quantifies the thickness of dense heterochromatin abutting the nuclear enveloppe. See related manuscript for details about image acquisition and experimental set up: https://doi.org/10.1101/2022.05.31.494218

The pipeline counting the capsids follows two steps: 1- nuclear boundaries identification, 2- capsid detection. Nuclear boundaries are identified with user input, by outlining a freehand contour of the nucleus of interest, from which a binary mask is extracted. Capsid detection is performed on the median-filtered complement of the original image, using a Circular Hough Transform based algorithm, with phase coding for radii estimation, and search radius ranging from 10.6 to 31.7 nm. Capsids residing within the nuclear mask are then counted and classified. Moreover, the distance of the capsids from the nearest nuclear membrane pixel is measured using a distance transform, and capsids within a 200nm range from the nuclear membrane are counted.

Width of heterochromatin abutting the nuclear envelope is quantified by measuring the length of the binarized chromatin from 1D intensity profiles along the normal of the nuclear perimeter, sampled at every 10 perimeter pixels. Dense heterochromatin is first binarized using global Otsu’s thresholding method applied to the background-corrected complement of the contrast-adjusted original image. Noise from the binarized image is further reduced by applying a 2D order statistic filter using the minimum value of a varying domain interactively defined by the user, with default value of 8-by-8 pixels. The resulting heterochromatin density distribution is normalized to the total length of the nucleus’ perimeter. 

REQUIREMENTS:

-MATLAB R2020a or later

-Image Processing ToolBox

  
