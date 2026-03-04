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

# read metadata of RTK-remeasured plots from 2023/2024
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

# read RTK-remeasured forest inventory plot centers
bi_plots <- sf::st_read(
  file.path(raw_data_dir, 'forest_inventory', 'bi_mittelpunkt_solling_2023.gpkg')
  )



# 02 - data joining
#-------------------------------------------------------------------------------

# deduplicate: keep row with min hz_accuracy, ties broken by newest datetime
accuracy_per_point <- measurement_details %>%
  arrange(point_nr, hz_accuracy, desc(datetime)) %>%
  distinct(point_nr, .keep_all = TRUE) %>%
  select(point_nr, hz_accuracy, v_accuracy)

# create join key by removing the first digit from KSPNR
bi_plots$point_nr <- as.numeric(substring(bi_plots$KSPNR, 2))

bi_plots <- dplyr::left_join(bi_plots, accuracy_per_point, by = 'point_nr')

bi_plots$point_nr <- NULL

sf::st_write(
  bi_plots,
  file.path(processed_data_dir, 'forest_inventory', 'bi_center_points_pos_acc.gpkg')
  )
