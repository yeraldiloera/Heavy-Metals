# Heavy-Metals
Heavy Metals in the Amazon

# Data and Code for: “Heavy metal contamination in birds from protected regions in the Amazon”
Wedge-billed woodcreeper (Glyphorynchus spirurus) feather samples were collected across control sites in the Andean foothills of Ecuador and anticipated contaminated sites from oil activities in Ecuador (Tiputini) and gold mining in French Guiana (Nouragues) between the years 2000-2009 were sent for heavy metal detection of 23 metals at Dartmouth’s Trace Element Analysis Core. 

# Heavy metal concentrations were reported in the Excel file labeled: FULL_Metal_Results_Data.xlsx. 
This Excel document contains three workbooks, with “Samples Metadata” providing information about the individual feather samples, “proc results” reporting the actual heavy metal concentrations for each feather, and “QC” containing the quality control data for these measurements.

# Heavy metal analyses (WBWC_ChemicalAnalysis_Code.R)
Heavy metal analyses were conducted R programming language: R Core Team (2021). R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/. Packages used for analysis included: readxl, dplyr, tidyverse, corrr, corrplot, ggcorrplot, ggplot2, ggpubr, ggthemes, RColorBrewer, AICcmodavg, broom, FactoMineR, factoextra, rstatix, sjstats, survey, sp, amt, rgdal, raster, rasterVis, rgeos. 

Data analyses included normality tests (density and qq plots, F- test, Shapiro-Wilkes test), AIC tests for interactions by feather weight, PCA plots, weighted mann-whitney tests of full contamination load across the two experimental groups and the control group, and weighted mann-whitney tests for each metal between the control group and a single experimental group (Tiputini and Nourages separately). Plots of these analyses were created with pairwise Wilcoxon rank-sums tests with Holm-Bonferroni p-adjustment for multiple comparisons. 

# Mapping
Mapping and supplementary kernel density estimates were performed using known oil sites in Ecuador according to the Global Energy Monitor (GEM) organization’s global gas and oil network database for Latina America (restricted access: Portal Energético para América Latina- Global Oil and Gas Extraction Tracker, August 2022 release), as well as of known mining sites filtered for gold mining activities from the French Republic Camino Digital Mining Cadastre (accessed November 2023). These Excel databases were subsequently filtered for extraction sites with known locations and activity beginning before feather sampling ended in 2009. The locations of these oil and gold mining activities were mapped in an ArcGIS base map overlaid with protected areas according to the UN’s Environmental World Conservation Monitoring Centre (September 2023 update: https://services5.arcgis.com/Mj0hjvkNtV7NRhA7/arcgis/rest/services/WDPA_v0/FeatureServer/1) to identify our reference and possibly contaminated localities (ArcGIS- Princeton University License [GIS software]. Version 10.0. Redlands, CA: Environmental Systems Research Institute, Inc., 2010).
