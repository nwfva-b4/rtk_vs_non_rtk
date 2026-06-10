#-------------------------------------------------------------------------------
# Name:         stand_heights_vs_ALS_heights.R
# Description:  Comparison of terrestrial-based stand heights with ALS-based
#               height percentiles in forest inventory plots. 
#               Influence of plot position (GNSS RTK-corrected vs. just GNSS) on
#               the mathcing of stand heights and height percentiles is analyzed.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01: data reading
#-------------------------------------------------------------------------------

# stand (reference) heights
stand_heights <- readRDS(file.path(processed_data_dir, 'forest_inventory', 'referenceHeights.RDS'))

# BI plots (RTK and non-RTK)
bi_plots_rtk <- sf::st_read(file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots_rtk.gpkg'))
bi_plots_non_rtk <- sf::st_read(file.path(processed_data_dir, 'forest_inventory', 'inv_attr_plots_non_rtk.gpkg'))
                            
