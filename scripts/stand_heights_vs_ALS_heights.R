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

# ALS point clouds (leaf-on = 2023, leaf-off = 2024)
pc_lon   <- list.files(
  file.path(raw_data_dir, 'pc_leafon_2023'),
  pattern = '2023', full.names = T
  )
pc_loff  <- list.files(
  file.path(raw_data_dir, 'pc_leafoff_2024'),
  pattern = '2024', full.names = T
  )



# 02: data preparation
#-------------------------------------------------------------------------------

# join stand heights (Hm, Hl, H100) to the BI plots via kspnr
bi_plots_rtk <- dplyr::left_join(
  bi_plots_rtk,
  dplyr::select(stand_heights, kspnr, Hm, Hl, H100),
  by = 'kspnr'
)

bi_plots_non_rtk <- dplyr::left_join(
  bi_plots_non_rtk,
  dplyr::select(stand_heights, kspnr, Hm, Hl, H100),
  by = 'kspnr'
)



# 03: ALS height percentile calculation
#-------------------------------------------------------------------------------

# metrics to compute: height percentiles 90, 95, 99 and maximum height
metrics_to_calc <- c('HAG_p90', 'HAG_p95', 'HAG_p99', 'HAG_max')

# helper: compute plot-level height metrics using lasR
calc_plot_metrics <- function(plots_sf, pc_files, metrics, radius = 13) {

  coords <- sf::st_coordinates(plots_sf)

  read <- lasR::reader(xc = coords[, 1], yc = coords[, 2], r = radius)
  summ <- lasR::summarise(metrics = metrics)

  pipeline <- read + summ
  ans <- lasR::exec(pipeline, on = pc_files)

  cbind(plots_sf, ans$metrics)
}

# calculate metrics for each combination of point cloud and plot positions
plot_metrics_lon_rtk      <- calc_plot_metrics(bi_plots_rtk,     pc_lon,  metrics_to_calc)
plot_metrics_loff_rtk     <- calc_plot_metrics(bi_plots_rtk,     pc_loff, metrics_to_calc)
plot_metrics_lon_non_rtk  <- calc_plot_metrics(bi_plots_non_rtk, pc_lon,  metrics_to_calc)
plot_metrics_loff_non_rtk <- calc_plot_metrics(bi_plots_non_rtk, pc_loff, metrics_to_calc)

# write results to disk
sf::st_write(
  plot_metrics_lon_rtk,
  file.path(processed_data_dir, 'metrics', 'height_percentiles_lon_rtk.gpkg'),
  delete_dsn = T
  )
sf::st_write(
  plot_metrics_loff_rtk, 
  file.path(processed_data_dir, 'metrics', 'height_percentiles_loff_rtk.gpkg'),
  delete_dsn = T
  )

sf::st_write(
  plot_metrics_lon_non_rtk,
  file.path(processed_data_dir, 'metrics', 'height_percentiles_lon_non_rtk.gpkg'),
  delete_dsn = T
  )
sf::st_write(
  plot_metrics_loff_non_rtk,
  file.path(processed_data_dir, 'metrics', 'height_percentiles_loff_non_rtk.gpkg'),
  delete_dsn = T
  )

