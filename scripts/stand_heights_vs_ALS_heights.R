#-------------------------------------------------------------------------------
# Name:         stand_heights_vs_ALS_heights.R
# Description:  Comparison of terrestrial-based stand heights with ALS-based
#               height metrics in forest inventory plots. 
#               Influence of plot position (GNSS RTK-corrected vs. just GNSS) on
#               the mathcing of stand heights and height metrics is analyzed.
# Author:       Florian Franz
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01: data reading
#-------------------------------------------------------------------------------

# stand (reference) heights
stand_heights <- readRDS(file.path(
  processed_data_dir, 'forest_inventory', 'referenceHeights.RDS')
  )

# BI plots (RTK and non-RTK)
bi_plots_rtk <- sf::st_read(file.path(
  processed_data_dir, 'forest_inventory', 'inv_attr_plots_rtk.gpkg')
  )
bi_plots_non_rtk <- sf::st_read(file.path(
  processed_data_dir, 'forest_inventory', 'inv_attr_plots_non_rtk.gpkg')
  )

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

# cache: per-plot ALS height metrics for each combination of point cloud season
# and plot positioning
metrics_files <- list(
  lon_rtk      = file.path(processed_data_dir, 'metrics', 'height_percentiles_lon_rtk.gpkg'),
  loff_rtk     = file.path(processed_data_dir, 'metrics', 'height_percentiles_loff_rtk.gpkg'),
  lon_non_rtk  = file.path(processed_data_dir, 'metrics', 'height_percentiles_lon_non_rtk.gpkg'),
  loff_non_rtk = file.path(processed_data_dir, 'metrics', 'height_percentiles_loff_non_rtk.gpkg')
)

if (all(file.exists(unlist(metrics_files)))) {

  # load cached results
  cat('Loading cached ALS height percentiles from disk...\n')
  plot_metrics_lon_rtk      <- sf::st_read(metrics_files$lon_rtk, quiet = T)
  plot_metrics_loff_rtk     <- sf::st_read(metrics_files$loff_rtk, quiet = T)
  plot_metrics_lon_non_rtk  <- sf::st_read(metrics_files$lon_non_rtk, quiet = T)
  plot_metrics_loff_non_rtk <- sf::st_read(metrics_files$loff_non_rtk, quiet = T)

} else {

  cat('No cached percentiles found, computing ALS height metrics...\n')

  # metrics to compute: height percentiles 90, 95, 99 and maximum height
  metrics_to_calc <- c('HAG_p90', 'HAG_p95', 'HAG_p99', 'HAG_max')

  # function to compute plot-level height metrics
  calc_plot_metrics <- function(plots_sf, pc_files, metrics, radius = 13) {

    coords <- sf::st_coordinates(plots_sf)

    read <- lasR::reader(xc = coords[, 1], yc = coords[, 2], r = radius)
    summ <- lasR::summarise(metrics = metrics)

    pipeline <- read + summ
    ans <- lasR::exec(pipeline, on = pc_files)

    cbind(plots_sf, ans$metrics)
  }

  # calculate metrics for each combination of point cloud and plot positions
  cat('  leaf-on / RTK...\n')
  plot_metrics_lon_rtk      <- calc_plot_metrics(bi_plots_rtk, pc_lon, metrics_to_calc)
  cat('  leaf-off / RTK...\n')
  plot_metrics_loff_rtk     <- calc_plot_metrics(bi_plots_rtk, pc_loff, metrics_to_calc)
  cat('  leaf-on / non-RTK...\n')
  plot_metrics_lon_non_rtk  <- calc_plot_metrics(bi_plots_non_rtk, pc_lon, metrics_to_calc)
  cat('  leaf-off / non-RTK...\n')
  plot_metrics_loff_non_rtk <- calc_plot_metrics(bi_plots_non_rtk, pc_loff, metrics_to_calc)

  # write results to disk for reuse
  cat('Writing ALS height metrics to disk...\n')
  sf::st_write(plot_metrics_lon_rtk, metrics_files$lon_rtk, delete_dsn = T)
  sf::st_write(plot_metrics_loff_rtk, metrics_files$loff_rtk, delete_dsn = T)
  sf::st_write(plot_metrics_lon_non_rtk, metrics_files$lon_non_rtk, delete_dsn = T)
  sf::st_write(plot_metrics_loff_non_rtk, metrics_files$loff_non_rtk, delete_dsn = T)
  
}



# 04: comparison analysis
#-------------------------------------------------------------------------------

# columns to compare
stand_height_cols <- c('Hm', 'Hl', 'H100')
als_metric_cols   <- c('HAG_p90', 'HAG_p95', 'HAG_p99', 'HAG_max')

# pairwise pearson correlations between stand heights and ALS metrics,
# returned in long format
calc_height_cors <- function(plot_metrics_sf, dataset_name) {

  df <- sf::st_drop_geometry(plot_metrics_sf)

  cor_mat <- cor(df[stand_height_cols], df[als_metric_cols],
                 use = 'pairwise.complete.obs')

  cor_df <- as.data.frame(as.table(cor_mat))
  names(cor_df) <- c('stand_height', 'als_metric', 'r')
  cor_df$dataset <- dataset_name

  cor_df[, c('dataset', 'stand_height', 'als_metric', 'r')]
}

# correlations for each combination of point cloud and plot positions
height_cors <- dplyr::bind_rows(
  calc_height_cors(plot_metrics_lon_rtk, 'lon_rtk'),
  calc_height_cors(plot_metrics_loff_rtk, 'loff_rtk'),
  calc_height_cors(plot_metrics_lon_non_rtk, 'lon_non_rtk'),
  calc_height_cors(plot_metrics_loff_non_rtk, 'loff_non_rtk')
)

# fix display order of stand heights and ALS metrics
height_cors$stand_height <- factor(
  height_cors$stand_height, levels = stand_height_cols
  )
height_cors$als_metric <- factor(
  height_cors$als_metric,   levels = als_metric_cols
  )

# correlation matrices in wide format for inspection
height_cors_wide <- tidyr::pivot_wider(
  height_cors,
  names_from = als_metric, values_from = r
)
print(height_cors_wide, n = Inf)

# best fitting ALS metric per stand height and dataset
best_fits <- height_cors |>
  dplyr::group_by(dataset, stand_height) |>
  dplyr::slice_max(r, n = 1) |>
  dplyr::ungroup()
print(best_fits, n = Inf)

# visual comparison of all correlations
ggplot(height_cors, aes(x = stand_height, y = r, fill = als_metric)) +
  geom_col(position = 'dodge') +
  facet_wrap(~ dataset) +
  coord_cartesian(ylim = c(min(height_cors$r) - 0.01, 1)) +
  labs(x = 'Stand height', y = "Pearson's r",
       fill = 'ALS height metric') +
  theme_bw()

# pair each stand height with each ALS metric in long format
prep_scatter_data <- function(plot_metrics_sf, dataset_name) {

  sf::st_drop_geometry(plot_metrics_sf) |>
    dplyr::select(kspnr, dominant_leaf_type,
                  dplyr::all_of(c(stand_height_cols, als_metric_cols))) |>
    tidyr::pivot_longer(dplyr::all_of(stand_height_cols),
                        names_to = 'stand_height', values_to = 'h_stand') |>
    tidyr::pivot_longer(dplyr::all_of(als_metric_cols),
                        names_to = 'als_metric', values_to = 'h_als') |>
    dplyr::mutate(dataset = dataset_name)
}

scatter_data <- dplyr::bind_rows(
  prep_scatter_data(plot_metrics_lon_rtk, 'lon_rtk'),
  prep_scatter_data(plot_metrics_loff_rtk, 'loff_rtk'),
  prep_scatter_data(plot_metrics_lon_non_rtk, 'lon_non_rtk'),
  prep_scatter_data(plot_metrics_loff_non_rtk, 'loff_non_rtk')
)

# fix panel order
scatter_data$stand_height <- factor(
  scatter_data$stand_height, levels = stand_height_cols
  )
scatter_data$als_metric   <- factor(
  scatter_data$als_metric,   levels = als_metric_cols
  )

# split dataset into point cloud season and positioning method
scatter_data <- scatter_data |>
  dplyr::mutate(
    season = factor(ifelse(grepl('^lon', dataset), 'leaf-on', 'leaf-off'),
                         levels = c('leaf-on', 'leaf-off')),
    positioning = ifelse(grepl('non_rtk', dataset), 'non-RTK', 'RTK')
  )

height_cors <- height_cors |>
  dplyr::mutate(
    season = factor(ifelse(grepl('^lon', dataset), 'leaf-on', 'leaf-off'),
                         levels = c('leaf-on', 'leaf-off')),
    positioning = ifelse(grepl('non_rtk', dataset), 'non-RTK', 'RTK')
  )

# scatterplots ALS metric vs. Hm with 1:1 line and pearson's r,
# RTK and non-RTK in one figure (leaf-on and leaf-off as facets),
# one figure per ALS metric
for (m in als_metric_cols) {

  scatter_data_m <- dplyr::filter(scatter_data,
                                  stand_height == 'Hm', als_metric == m)

  # identical limits on both axes
  axis_lims <- range(c(scatter_data_m$h_als, scatter_data_m$h_stand),
                     na.rm = T)

  p <- ggplot(
    scatter_data_m,
    aes(x = h_als, y = h_stand, colour = positioning)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = 1, intercept = 0,
                linetype = 'dashed', colour = 'black') +
    geom_text(data = dplyr::filter(height_cors,
                                   stand_height == 'Hm', als_metric == m),
              aes(label = paste0('r = ', round(r, 2)),
                  colour = positioning,
                  vjust = ifelse(positioning == 'RTK', 1.5, 3)),
              x = -Inf, y = Inf, hjust = -0.2,
              size = 3.5, show.legend = F, inherit.aes = F) +
    facet_wrap(~ season) +
    coord_equal(xlim = axis_lims, ylim = axis_lims) +
    labs(x = paste(m, '[m]'), y = 'Hm [m]', colour = 'Positioning') +
    theme_bw()

  print(p)
}

# summary table: best performing ALS metric per season and stand height,
# and whether RTK or non-RTK positioning performs better
rtk_vs_non_rtk <- height_cors |>
  tidyr::pivot_wider(id_cols = c(season, stand_height, als_metric),
                     names_from = positioning, values_from = r) |>
  dplyr::rename(r_rtk = RTK, r_non_rtk = `non-RTK`) |>
  dplyr::group_by(season, stand_height) |>
  dplyr::slice_max(pmax(r_rtk, r_non_rtk), n = 1) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    better = ifelse(r_rtk > r_non_rtk, 'RTK', 'non-RTK'),
    r_diff = r_rtk - r_non_rtk
  )
print(rtk_vs_non_rtk, n = Inf)



# 05: comparison analysis by dominant leaf type
#-------------------------------------------------------------------------------

# pearson's r per leaf type, season, positioning,
# stand height, and ALS metric
height_cors_lt <- scatter_data |>
  dplyr::group_by(dominant_leaf_type, season, positioning,
                  stand_height, als_metric) |>
  dplyr::summarise(r = cor(h_stand, h_als, use = 'pairwise.complete.obs'),
                   .groups = 'drop')

# correlation matrices in wide format for inspection
height_cors_lt_wide <- tidyr::pivot_wider(
  height_cors_lt,
  names_from = als_metric, values_from = r
)
print(height_cors_lt_wide, n = Inf)

# visual comparison of all correlations
ggplot(height_cors_lt, aes(x = stand_height, y = r, fill = als_metric)) +
  geom_col(position = 'dodge') +
  facet_grid(dominant_leaf_type ~ season + positioning) +
  coord_cartesian(ylim = c(min(height_cors_lt$r) - 0.01, 1)) +
  labs(x = 'Stand height', y = "Pearson's r",
       fill = 'ALS height metric') +
  theme_bw()

# scatterplots ALS metric vs. Hm with 1:1 line and pearson's r,
# RTK and non-RTK in one figure (leaf types as rows, seasons as columns),
# one figure per ALS metric
for (m in als_metric_cols) {

  scatter_data_m <- dplyr::filter(scatter_data,
                                  stand_height == 'Hm', als_metric == m)

  # identical limits on both axes
  axis_lims <- range(c(scatter_data_m$h_als, scatter_data_m$h_stand),
                     na.rm = T)

  p <- ggplot(
    scatter_data_m,
    aes(x = h_als, y = h_stand, colour = positioning)) +
    geom_point(alpha = 0.5, size = 1.5) +
    geom_abline(slope = 1, intercept = 0,
                linetype = 'dashed', colour = 'black') +
    geom_text(data = dplyr::filter(height_cors_lt,
                                   stand_height == 'Hm', als_metric == m),
              aes(label = paste0('r = ', round(r, 2)),
                  colour = positioning,
                  vjust = ifelse(positioning == 'RTK', 1.5, 3)),
              x = -Inf, y = Inf, hjust = -0.2,
              size = 3.5, show.legend = F, inherit.aes = F) +
    facet_grid(dominant_leaf_type ~ season) +
    coord_equal(xlim = axis_lims, ylim = axis_lims) +
    labs(x = paste(m, '[m]'), y = 'Hm [m]', colour = 'Positioning') +
    theme_bw()

  print(p)
}

# summary table: best performing ALS metric per leaf type, season,
# and stand height, and whether RTK or non-RTK positioning performs better
rtk_vs_non_rtk_lt <- height_cors_lt |>
  tidyr::pivot_wider(id_cols = c(dominant_leaf_type, season,
                                 stand_height, als_metric),
                     names_from = positioning, values_from = r) |>
  dplyr::rename(r_rtk = RTK, r_non_rtk = `non-RTK`) |>
  dplyr::group_by(dominant_leaf_type, season, stand_height) |>
  dplyr::slice_max(pmax(r_rtk, r_non_rtk), n = 1) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    better = ifelse(r_rtk > r_non_rtk, 'RTK', 'non-RTK'),
    r_diff = r_rtk - r_non_rtk
  )
print(rtk_vs_non_rtk_lt, n = Inf)



# 06: plots deviating most from the 1:1 line
#-------------------------------------------------------------------------------

# signed deviation from the 1:1 line (stand height - ALS metric); large absolute
# values mark plots where terrestrial and ALS heights disagree most
scatter_data <- scatter_data |>
  dplyr::mutate(dev_from_1to1 = h_stand - h_als)

# the n plots furthest from the 1:1 line for a given stand height and
# ALS metric, per season and positioning
plots_far_from_line <- function(stand_h = 'Hm', als_m = 'HAG_p99', n = 10) {
  scatter_data |>
    dplyr::filter(stand_height == stand_h, als_metric == als_m) |>
    dplyr::group_by(season, positioning) |>
    dplyr::slice_max(abs(dev_from_1to1), n = n) |>
    dplyr::ungroup() |>
    dplyr::arrange(season, positioning, dplyr::desc(abs(dev_from_1to1))) |>
    dplyr::select(kspnr, dominant_leaf_type, season, positioning,
                  h_stand, h_als, dev_from_1to1)
}

# example: plots furthest from the 1:1 line for Hm vs. HAG_p99
print(plots_far_from_line('Hm', 'HAG_p99', n = 10), n = Inf)

