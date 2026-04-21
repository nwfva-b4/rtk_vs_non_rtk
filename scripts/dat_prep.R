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

# merge both RTK datasets
bi_plots <- rbind(
  bi_plots_1[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')],
  bi_plots_2[, c('kspnr', 'center_point_estimated', 'solution_status', 'measurement_date')]
)



# 04 - spatial cropping to point cloud extent
#-------------------------------------------------------------------------------

# read leaf-off point cloud extent
ext_loff <- sf::st_read(
  file.path(raw_data_dir, 'pc_leafoff_2024', 'leafoff.vpc')
)

ext_loff <- sf::st_transform(ext_loff, sf::st_crs(bi_plots))

bi_plots_1 <- sf::st_crop(bi_plots_1, sf::st_bbox(ext_loff))
bi_plots_2 <- sf::st_crop(bi_plots_2, sf::st_bbox(ext_loff))
bi_plots   <- sf::st_crop(bi_plots,   sf::st_bbox(ext_loff))

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
