---
title: "shinyForce3D: 3D Visualisation of Ground Reaction Forces"
tags:
  - R
  - Shiny
  - Biomechanics
  - Visualisation
  - Force vector
authors:
  - name: Daniel Kadlec
    orcid: 0000-0002-6986-4689
    affiliation: "1"
  - name: Shayne Vial
    orcid: 0000-0002-9235-8979
    affiliation: "1"


affiliations:
  - name: Edith Cowan University, School of Medical and Health Sciences
    index: 1

date: "`r format(Sys.Date(), '%d %B %Y')`"
bibliography: paper.bib
csl: joss.csl
---

# Summary
`shinyForce3D` is an open-source R Shiny application that enables dynamic, interactive visualisation of three-dimensional ground reaction force data using animated force vector trajectories. The app supports real-time rendering of animated 3D vectors, static 3D visualisations with planar projections, two-dimensional time-aligned Pedotti plots of selected vector components, and decomposition of traditional Cartesian components (Fx, Fy, Fz) into polar coordinate representations (magnitude *r*, inclination *θ*, and azimuth *φ*) to aid interpretation.

Users can upload CSV files with trial-level force data, overlay multiple data sets for comparative analysis (e.g., between-limb, between-day, or between-cohort visualisations), and export publication-quality figures and animations in common formats such as `.csv` and `.gif`.

In addition to its research applications, `shinyForce3D` is designed as a pedagogical tool for use in kinesiology, biomechanics and motor control education. By visualising how GRF vectors evolve over time and across conditions, students can develop a more intuitive understanding of ground reaction forces during athletic tasks than is typically afforded by traditional 2D plots.


# Statement of Need
Ground reaction force vector analysis is a foundational tool in biomechanics for evaluating the cause of human locomotion and to determine human performance characteristics in sport and rehabilitation contexts [@Winter2009]. Traditionally, ground reaction force data are presented as three separate 2D time-series (Fx, Fy, Fz), which obscures the vectorial nature of force and fails to capture inter-component dependencies [@Pataky2013]. The representation of the individual forec components (Fx, Fy, Fz) may hinder interpretation of how force vector orientation, magnitude, and timing evolve dynamically across the a given time domain. These limitations hinder not only exploratory visualisation and adequately formulating mechanically sound research questions but also clear communication of findings, particularly in interdisciplinary or educational settings.

In contrast, ground reaction force is ontologically a vector quantity, defined by magnitude and direction in three-dimensional space, and is best preserved when represented in its native vectorial form. Visualising ground reaction force in three dimensions enables a more accurate depiction of the inter‐component dependencies, revealing how force reorientation and amplitude co‐evolve across a given time domain. By integrating 3D vector trajectories, Pedotti plots, and polar coordinate visualisations, `shinyForce3D` provides users with a richer and more coherent understanding of the emergence of ground reaction force.

Unlike traditional two-dimensional static plotting scripts, `shinyForce3D` offers a responsive, user-friendly interface tailored for biomechanics researchers, sports scientists, and educators. It is particularly well-suited for teaching vector-based reasoning and enhancing student comprehension of force dynamics through immediate, visual feedback.


# Features
•	Dynamic 3D Plot
	This tab provides an interactive three-dimensional rendering of the full GRF vector trajectory across the stance phase. Users can rotate, zoom, and explore the vector evolution in space, with optional overlays such as butterfly lines (vector paths) and peak force indicators. This mode offers a structural depiction of force, crucial for understanding how direction and magnitude co-evolve during movement.
	
•	Static 3D Plot
	The static plot generates dual orthographic projections (e.g., frontal and sagittal views) with optional 3D shadowing on anatomical planes, which is commonly used for most atheltic tasks (i.e. sprinting, jumping, decelerating, etc). These planar projections facilitate side-by-side comparison and allow users to examine differences between trials. This view supports static 3D analysis and figure export while preserving spatial relationships between components. These figure are suitable for publications.
	
•	2D Plot (Cartesian Coordinate)
	This tab shows conventional time-series plots of Fx, Fy, Fz, and the resultant vector across stance (i.e. time). Users can toggle between trial-level and mean ± SD visualisations, highlight peak values, and customise which components are displayed.
	
•	2D Plot (Polar Coordinate)
	Force vectors are transformed into polar coordinates, decomposing each time point into magnitude (r), inclination (θ), and azimuth (φ), with optional angle display in degrees. This representation allows users to disentangle directional and magnitude changes, helping to reveal rotational strategies or orientation shifts that Cartesian time-series can obscure. It also facilitates intuitive reasoning about vector directionality.
	
•	Pedotti Plot
	This plot maps any two selected GRF components (e.g., Fy vs. Fz) as time-evolving vectors from the origin, creating a continuous 2D vector field over stance [@Pedotti1977]. Animated and static options are available, with time encoded via color progression. Originally developed to visualise neuromotor control strategies, Pedotti plots are valuable for identifying coordination patterns, symmetry, and variability in multidimensional force application.
	
	
# Data Upload
`shinyForce3D` allows users to upload CSV files containing ground reaction force data in Cartesian form (Fx, Fy, Fz) for one or multiple trials. Each trial must consist of three consecutive columns representing the force components - in order: Fx, Fy, Fz - without any header rows or additional metadata. Files should contain only numeric data, and all trials must be aligned in time (i.e., each trial must have the same number of rows, corresponding to normalised time or percent stance).

When multiple trials are included, `shinyForce3D` automatically segments the data by triplets of columns and computes the trial-wise mean and standart deviation for visualisation. No specific column names are required, as the app infers structure based on column order.

To ensure smooth parsing:

- Do not include subject IDs, time stamps, or non-numeric rows.
- Ensure the total number of columns is divisible by three.
- Confirm that all rows are complete (no missing data in the middle of triplets).

# References
