#-------------------------------------------------------------------------------
# Name:         dat_prep.R
# Description:  Preparation of RTK-remeasured forest inventory plots (BI)
#               in the Solling region.
# Contact:      florian.franz@nw-fva.de
#-------------------------------------------------------------------------------



# source setup script
source('src/setup.R', local = TRUE)



# 01 - data reading
#-------------------------------------------------------------------------------

# read RTK-remeasured forest inventory plot centers
# from 2023/24 and 2025/26
bi_plots_1 <- sf::st_read(
  file.path(raw_data_dir, 'forest_inventory', 'bi_center_points_solling_2023_24.gpkg')
  )

# read RTK-remeasured forest inventory plot centers
bi_plots_2 <- sf::st_read(
  file.path(raw_data_dir, 'forest_inventory', 'bi_center_points_solling_2025_26.gpkg')
)

# read metadata of RTK-remeasured plots from 2023/24
measurement_details <- read.table(
  file.path(meta_data_dir, 'Neueinmessung Solling-etrs89-messdetails.txt'),
  skip = 10,
  header = F,
  fill = T,
  stringsAsFactors = F
)

colnames(measurement_details) <- c(
  'point_nr', 'antenna_height',
  'time', 'day', 'month', 'year',
  'epoch', 'hz_accuracy', 'v_accuracy', 'sv', 'rms'
)

lc_time_old <- Sys.getlocale('LC_TIME')
Sys.setlocale('LC_TIME', 'C')

measurement_details$datetime <- as.POSIXct(
  paste(
    measurement_details$year, measurement_details$month,
    measurement_details$day, measurement_details$time
  ),
  format = '%Y %b %d %H:%M:%S'
)

Sys.setlocale('LC_TIME', lc_time_old)

measurement_details <- measurement_details[
  , !(names(measurement_details) %in% c('time', 'day', 'month', 'year', 'rms'))
]

# read comments table to identify estimated plot centers 
# of RTK-remeasured plots from 2023/24
comments <- readxl::read_excel(
  file.path(meta_data_dir, '2024_01_15_Hauck_GNSS_BI-Punkte_Liste_Kommentare.xlsx')
)



# 02 - data cleaning
#-------------------------------------------------------------------------------

estimated_kspnr <- comments$Pktnr[grepl('ohne Nagel gemessen',
                                         comments$`Kommentar 1`,
                                         ignore.case = T)]

# clean up bi_plots_1:
bi_plots_1$kspnr <- bi_plots_1$KSPNR
bi_plots_1$center_point_estimated <- ifelse(
  bi_plots_1$KSPNR %in% estimated_kspnr, 'yes', 'no'
)
sf::st_geometry(bi_plots_1) <- 'geometry'

# clean up bi_plots_2:
bi_plots_2 <- bi_plots_2[!sf::st_is_empty(bi_plots_2), ]
bi_plots_2$center_point_estimated <- ifelse(grepl('e$', bi_plots_2$Name), 'yes', 'no')
bi_plots_2$kspnr <- as.integer(gsub('[a-zA-Z]$', '', bi_plots_2$Name))
bi_plots_2 <- sf::st_zm(bi_plots_2, drop = T, what = 'ZM')
bi_plots_2 <- bi_plots_2 %>%
  dplyr::arrange(kspnr, dplyr::desc(center_point_estimated)) %>%
  dplyr::filter(!duplicated(kspnr))

# harmonize geometry column naming for compatibility
bi_plots_2 <- dplyr::rename(bi_plots_2, geometry = geom)
sf::st_geometry(bi_plots_2) <- 'geometry'



# 03 - solution status
#-------------------------------------------------------------------------------

# bi_plots_1: derive solution status from hz_accuracy via measurement details
accuracy_per_point <- measurement_details %>%
  arrange(point_nr, hz_accuracy, desc(datetime)) %>%
  distinct(point_nr, .keep_all = T) %>%
  select(point_nr, hz_accuracy, v_accuracy, datetime)

bi_plots_1$point_nr <- as.numeric(substring(bi_plots_1$kspnr, 2))
bi_plots_1 <- dplyr::left_join(bi_plots_1, accuracy_per_point, by = 'point_nr')
bi_plots_1$point_nr <- NULL

bi_plots_1$solution_status <- dplyr::case_when(
  is.na(bi_plots_1$hz_accuracy)    ~ 'single',
  bi_plots_1$hz_accuracy <= 0.1    ~ 'fix',
  TRUE                             ~ 'float'
)

bi_plots_1$measurement_date <- as.Date(bi_plots_1$datetime)

# bi_plots_2: solution status is already in the data
bi_plots_2$solution_status <- tolower(bi_plots_2$Solution)
bi_plots_2$measurement_date <- as.Date(bi_plots_2$Averaging.start)


# exclude single-solution plots from campaign 1 for analysis
bi_plots_1 <- bi_plots_1 %>%
  dplyr::filter(solution_status != 'single')

# merge both RTK datasets
bi_plots <- rbind(
  bi_plots_1[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')],
  bi_plots_2[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')]
)



# 04 - spatial filtering
#-------------------------------------------------------------------------------

# 04.1. extent filter
#---------------------

# read leaf-off point cloud extent
ext_loff <- sf::st_read(
  file.path(raw_data_dir, 'pc_leafoff_2024', 'leafoff.vpc')
)

ext_loff <- sf::st_transform(ext_loff, sf::st_crs(bi_plots))

# keep only plots whose 13 m buffered footprint lies fully within the extent
ext_geom <- sf::st_union(sf::st_buffer(ext_loff, 0.5))

filter_within_extent <- function(plots, ext, radius = 13) {
  buf <- sf::st_buffer(plots, radius)
  plots[lengths(sf::st_within(buf, ext)) > 0, ]
}

bi_plots_1 <- filter_within_extent(bi_plots_1, ext_geom)
bi_plots_2 <- filter_within_extent(bi_plots_2, ext_geom)
bi_plots   <- filter_within_extent(bi_plots,   ext_geom)

# 04.2. angular sector coverage filter
# drops plots whose point cloud has a large empty bin (data gap on one side),
# even though the plot lies inside the nominal extent
#-------------------------------------------------------------------------------

laz_dir <- file.path(raw_data_dir, 'pc_leafoff_2024')

laz_files_2024 <- list.files(
  laz_dir, pattern = "_2024\\.laz$", full.names = T
)

ctg <- lidR::readLAScatalog(laz_files_2024)

angular_sector_coverage <- function(ctg, x, y, radius = 13, n_sectors = 8) {
  
  # clip circular plot from the catalog around (x, y) with the given radius
  las <- lidR::clip_circle(ctg, x, y, radius)
  if (lidR::is.empty(las)) {
    return(list(min_sector_frac = 0))
  }
  
  # total number of points inside the circle
  n_pts <- lidR::npoints(las)
  
  # angle of each point relative to the plot center, in radians (-pi..pi)
  ang <- atan2(las$Y - y, las$X - x)
  
  # assign each point to one of n_sectors equal-width angular bins
  sec <- cut(ang, breaks = seq(-pi, pi, length.out = n_sectors + 1),
             include.lowest = T)
  
  # count how many points fell into each sector
  counts <- as.integer(table(sec))
  
  # ratio of the emptiest sector to the expected share under uniform coverage
  # (1.0 = perfectly uniform, 0 = at least one bin is completely empty)
  min_sector_frac <- min(counts) / (n_pts / n_sectors)
  
  list(min_sector_frac = min_sector_frac)
}

filter_by_angular_coverage <- function(
    plots, ctg, radius = 13, n_sectors = 8, min_frac = 0
    ) 
  {
  
  # for each plot, measure how evenly points are distributed around the center
  # a value of 0 means at least one bin is empty -> data gap on one side
  plots$min_sector_frac <- vapply(seq_len(nrow(plots)), function(i) {
    co <- sf::st_coordinates(plots[i, ])
    angular_sector_coverage(ctg, co[1], co[2],
                            radius = radius,
                            n_sectors = n_sectors)$min_sector_frac
  }, numeric(1))
  
  plots[plots$min_sector_frac > min_frac, ]
}

bi_plots_1 <- filter_by_angular_coverage(bi_plots_1, ctg)
bi_plots_2 <- filter_by_angular_coverage(bi_plots_2, ctg)
bi_plots   <- filter_by_angular_coverage(bi_plots,   ctg)

# write to disk
sf::st_write(
  bi_plots_1[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')],
  file.path(processed_data_dir, 'forest_inventory', 'bi_center_points_2023_24.gpkg'),
  delete_dsn = TRUE
)

sf::st_write(
  bi_plots_2[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')],
  file.path(processed_data_dir, 'forest_inventory', 'bi_center_points_2025_26.gpkg'),
  delete_dsn = TRUE
)

sf::st_write(
  bi_plots,
  file.path(processed_data_dir, 'forest_inventory', 'bi_center_points_merged.gpkg'),
  delete_dsn = TRUE
)
