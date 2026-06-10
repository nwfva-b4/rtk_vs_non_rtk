#----------------------------------------------------------------------------------------------------------------------------------------------
# Name:         inv_attr_plots.R
# Description:  Calculation of forest attributes in inventory plots (Betriebsionventur (BI) Lower Saxony).
#               Inventory data is first pre-processed and then volume (total and merchantable)
#               and above ground biomass (AGB) are calculated for individual trees. 
#               The tree volumes and AGB are aggregated per sample plot to obtain the growing stock volume (GSV) [m³/ha] and AGB [t/ha].
#               Other attributes which are calculated per sample plot include tree density [n/ha],
#               basal area [m³/ha], and quadratic mean diameter (QMD) [cm].
# Author:       Christoph Fischer, Georgia Reeves, Florian Franz
# Contact:      christoph.fischer@nw-fva.de
#               georgia.reeves@nw-fva.de
#               florian.franz@nw-fva.de
#----------------------------------------------------------------------------------------------------------------------------------------------


# source setup script
source('src/setup.R', local = TRUE)



# 01: data reading
#-------------------------------------------------------------------------------

# input paths
bi_path <- file.path(raw_data_dir, 'forest_inventory')
bi_rtk_path <- file.path(processed_data_dir, 'forest_inventory')

# read forest inventory (BI) data
bi_points <- read.table(
  file.path(bi_path, 'tblDatPh2_ZE.txt'),
  header = T, sep = ';'
)

bi_trees <- read.table(
  file.path(bi_path, 'tblDatPh2_Vorr_ZE.txt'),
  header = T, sep = ';'
)

# select desired forestry offices (Solling --> Neuhaus, Dassel)
bi_points <- bi_points[
  bi_points$DatOrga_Key == '268-2022-002' | 
    bi_points$DatOrga_Key == '254-2022-002',
]

bi_trees <- bi_trees[
  bi_trees$DatOrga_Key == '268-2022-002' |
    bi_trees$DatOrga_Key == '254-2022-002',
]

head(bi_points)
head(bi_trees)

# read remeasured plots (RTK-GNSS)
# already prepared in script dat_prep.R
bi_plots_rtk <- sf::st_read(
  file.path(bi_rtk_path, 'bi_center_points_merged_filtered.gpkg')
  )

head(bi_plots_rtk)



# 02: data preparation
#-------------------------------------------------------------------------------

# source and apply function for data formatting
source('src/format_data.R', local = TRUE)

bi_points <- format_data(bi_points)
bi_trees <- format_data(bi_trees)

head(bi_points)
str(bi_points)
head(bi_trees)
str(bi_trees)

## delete deadwood and used trees
bi_trees <- bi_trees[!bi_trees$ba %in% seq(100,800,100),]

bi_trees <- bi_trees[bi_trees$'1' < 3 & bi_trees$'2' < 3,]

bi_trees <- bi_trees[bi_trees$art != 1 & bi_trees$art != 2 & bi_trees$bhd > 0,]

# select needed columns
bi_trees <- bi_trees[,c(1:12)]

# assign tree species groups
bi_trees$bagr <- 
  ifelse(bi_trees$ba > 0 & bi_trees$ba < 200, "EI",
  ifelse(bi_trees$ba > 199 & bi_trees$ba < 300, "BU",
  ifelse(bi_trees$ba > 299 & bi_trees$ba < 400, "ALH",	
  ifelse(bi_trees$ba > 399 & bi_trees$ba < 500, "ALN",
  ifelse(bi_trees$ba > 499 & bi_trees$ba < 600, "FI",
  ifelse(bi_trees$ba > 599 & bi_trees$ba < 700, "DGL",
  ifelse(bi_trees$ba > 699 & bi_trees$ba < 800, "KI",	"LAE")
  ))))))

# DBH correction
# if not measured at 1.3 m (deviating measuring height),
# then correction to 1.3 m
bi_trees$ba1 <- bi_trees$ba

# red oak to oak, fir to spruce, hornbeam to beech
source('src/d_corr_func.R', local = TRUE)

bi_trees$ba <- input_d_korr(bi_trees$bagr)

# average DBH from 'Kreuzkluppung (bhdklup)',
# convert to cm
bi_trees$bhd <- ifelse(
  bi_trees$bhdklup > 0,
  (0.5 * (bi_trees$bhd + bi_trees$bhdklup)) / 10,
  bi_trees$bhd / 10
)

# separate data:
# trees with diameter at deviating measurement height
# and without deviating measurement height
bi_trees_2 <- bi_trees[bi_trees$bhddiff > 0,]
bi_trees <- bi_trees[bi_trees$bhddiff == 0,]

# convert diameter to DBH in case of deviating measuring height
d <- d_korr(du = bi_trees_2$bhd, abwmh = bi_trees_2$bhddiff, ba = bi_trees_2$ba)
bi_trees_2$bhd <- d

# bind the two data frames together again
# (correctly back to original number)
bi_trees <- rbind(bi_trees, bi_trees_2)
rm(bi_trees_2, d)

# add original tree species again
bi_trees$ba <- bi_trees$ba1
bi_trees$ba1 <- NULL

head(bi_trees)

# calculate number of stems per ha

# concentric sample circles:
#	r = 6 m all trees 
#	r = 13 m all trees with DBH >= 30 cm
# radius must be projected into the plane

# correct sample circle sizes with slope, calculate N_ha
# inclination from degrees in rad
bi_points$hang_rad <- (pi / 180) * bi_points$hang

bi_points_trees <- merge(
  bi_trees, bi_points[,c("key", "kspnr", "abt", "hang_rad", "rw", "hw")],
  by = c("key", "kspnr")
)

# r_plane = r_slope * cos(slope_rad)
bi_points_trees$nha <- ifelse(
  bi_points_trees$bhd < 30, 
  10000 / (pi * 6**2 * cos(bi_points_trees$hang_rad)),
  10000 / (pi * 13**2 * cos(bi_points_trees$hang_rad))
)

## add heights

# heights in m
bi_points_trees$hoehe <- bi_points_trees$hoehe / 10

# new ID consisting of key + sample point number
bi_points_trees$id2 <- paste(bi_points_trees$key, bi_points_trees$kspnr, sep = "_")

# data format for heights adding
source('src/ehk_func.R', local = TRUE)

dat <- input_ehk(
  id = bi_points_trees$id2, bnr = bi_points_trees$id,
  bs = bi_points_trees$bestschicht, bhd = bi_points_trees$bhd,
  hoe = bi_points_trees$hoehe, nha = bi_points_trees$nha,
  bagr = bi_points_trees$bagr
)

head(dat)

# assigning appropriate quantile percentages in order
# to properly remove outliers (unrealistic DBH height value pairs)
# in a statistically sound way by building a scam model with height and DBH
dat <- dat[dat$hoe > 0, ]
summary(dat)

m <- scam::scam(
  hoe ~ s(bhd, bs = 'mpi'),
  data = dat, 
  family = Gamma(link = 'log')
)

nd <- data.frame('bhd' = floor(min(dat$bhd)):ceiling(max(dat$bhd)))
nd$hoe <- predict(m, newdata = nd, type = 'response')

p <- ggplot(data = dat, aes(x = bhd, y  = hoe)) + 
  geom_point(color = rgb(.5, .5, .5, alpha = .2)) + 
  geom_line(dat = nd, color = 1, linewidth = 2)

tmp <- NULL
for (x in seq(10, 110, by = 10)) {
  nd2 <- data.frame('bhd' = x)
  nd2$hoe <- predict(m, newdata = nd2, type = 'response')
  v <- 1/m$sig2
  d <- stats::dgamma(1:60, shape = (nd2$hoe[1]^2)/v, scale = v/nd2$hoe[1])
  tmp <- rbind(
    tmp, 
    data.frame(
      'bhd' = x - (9 * d / max(d)), 
      'hoe' = 1:60, 
      'x' = x
    )
  )
}

p1 <- p + 
  geom_path(dat = tmp, aes(group = x, color = factor(x)), show.legend = F) +
  geom_vline(
    data = data.frame('bhd' = seq(10, 110, by = 10)), 
    aes(xintercept = bhd, color = factor(bhd)), show.legend = F, 
    linetype = 2
  )

v <- 1/m$sig2
dat$p <- stats::pgamma(q = dat$hoe, shape = (fitted(m)^2)/v, scale = v/fitted(m))
summary(dat$p)

dat$lab <- ''
ix_label <- sort(c(which(dat$p > .99), which(dat$p < .01)))
dat$lab[ix_label] <- paste0(round(dat$p[ix_label] * 100, 2), '%')

p2 <- ggplot(data = dat, aes(x = bhd,y  = hoe)) + 
  geom_point(color = rgb(.5, .5, .5, alpha = .2)) + 
  geom_line(dat = nd, color = 1, linewidth = 2) + 
  ggrepel::geom_text_repel(aes(label = lab))

dat$lab <- round(dat$p*100, 2)
dat_extremes <- dat[dat$lab >= 99.98 | dat$lab <= 0.02,]
dat_without_extremes <- dat[dat$lab < 99.98 & dat$lab > 0.02,]

cowplot::plot_grid(p1, p2, ncol = 2)

# re-assigning dat to its original value
dat <- input_ehk(
  id = bi_points_trees$id2, bnr = bi_points_trees$id,
  bs = bi_points_trees$bestschicht, bhd = bi_points_trees$bhd,
  hoe = bi_points_trees$hoehe, nha = bi_points_trees$nha,
  bagr = bi_points_trees$bagr
)

dat2 <- dat[dat$hoe == 0,]

dat3 <- subset(dat_without_extremes, select = -c(p, lab))
dat <- rbind(dat2, dat3)
rm(dat_without_extremes)
rm(dat2)
rm(dat3)
rm(dat_extremes)

# uniform height curve
dat2 <- ehk(dat)
plot(dat2$bhd, dat2$hoe_mod)

# merge modeled heights with original table
bi_points_trees <- merge(
  bi_points_trees,
  dat2[,c("id", "bnr", "hoe_mod")],
  by.x = c("id2","id"), 
  by.y = c("id", "bnr")
)

# remove unneeded data frames
rm(dat2, dat)



# 03: calculate individual tree volume and AGB
#-------------------------------------------------------------------------------

vol_agb <- bi_points_trees

# recode tree species before applying rBDAT
# mapping from BI tree species codes to rBDAT species codes
vol_agb$ba1 <- ifelse(vol_agb$ba == 110 | vol_agb$ba == 111 | vol_agb$ba == 112, 17,
               ifelse(vol_agb$ba == 113, 18,
               ifelse(vol_agb$ba == 211, 15,
               ifelse(vol_agb$ba == 221, 16,
               ifelse(vol_agb$ba == 311, 21,
               ifelse(vol_agb$ba == 320, 22,
               ifelse(vol_agb$ba == 321, 23,
               ifelse(vol_agb$ba == 322, 24,
               ifelse(vol_agb$ba == 323, 25,
               ifelse(vol_agb$ba == 330 | vol_agb$ba == 331 | vol_agb$ba == 332, 30,
               ifelse(vol_agb$ba == 341 | vol_agb$ba == 342, 27,
               ifelse(vol_agb$ba == 351, 31,
               ifelse(vol_agb$ba == 352 | vol_agb$ba == 442, 33,
               ifelse(vol_agb$ba == 353 | vol_agb$ba == 355, 35,
               ifelse(vol_agb$ba == 354 | vol_agb$ba == 452, 29,
               ifelse(vol_agb$ba == 357, 32,
               ifelse(vol_agb$ba == 410 | vol_agb$ba == 411 | vol_agb$ba == 412 | vol_agb$ba == 414, 26,
               ifelse(vol_agb$ba == 420 | vol_agb$ba == 421 | vol_agb$ba == 422, 28,
               ifelse(vol_agb$ba == 430 | vol_agb$ba == 431, 19,
               ifelse(vol_agb$ba == 441, 34,
               ifelse(vol_agb$ba == 451, 36,
               ifelse(vol_agb$ba == 511 | vol_agb$ba == 513 | vol_agb$ba == 551, 1,
               ifelse(vol_agb$ba == 512, 2,
               ifelse(vol_agb$ba == 521, 3,
               ifelse(vol_agb$ba == 523, 4,
               ifelse(vol_agb$ba == 541, 13,
               ifelse(vol_agb$ba == 542, 12,
               ifelse(vol_agb$ba == 611, 8,
               ifelse(vol_agb$ba == 711, 5,
               ifelse(vol_agb$ba == 712, 6,
               ifelse(vol_agb$ba == 731, 7,
               ifelse(vol_agb$ba == 810, 9,
               ifelse(vol_agb$ba == 811, 10,
               ifelse(vol_agb$ba == 812, 11, vol_agb$ba)
               )))))))))))))))))))))))))))))))))

# calculation of volume with rBDAT::getVolume
# https://gitlab.com/vochr/rbdat
print('Calculating total volume...')
vol_agb$total_vol <- rBDAT::getVolume(
  vol_agb,
  bark = T,
  mapping = c('ba1' = 'spp', 'bhd' = 'D1', 'hoe_mod' = 'H')
)

print('Calculating merchantable volume...')
vol_agb$merch_vol <- rBDAT::getVolume(
  vol_agb,
  bark = F,
  mapping = c('ba1' = 'spp', 'bhd' = 'D1', 'hoe_mod' = 'H')
)

# calculation of AGB with rBDAT::getBiomass
print('Calculating AGB...')
vol_agb$agb <- rBDAT::getBiomass(
  vol_agb,
  mapping = c('ba1' = 'spp', 'bhd' = 'D1', 'hoe_mod' = 'H')
)

vol_agb$ba1 <- NULL

# add results to bi_points_trees
bi_points_trees$total_vol <- vol_agb$total_vol
bi_points_trees$merch_vol <- vol_agb$merch_vol
bi_points_trees$agb <- vol_agb$agb

summary(bi_points_trees$total_vol)
summary(bi_points_trees$merch_vol)
summary(bi_points_trees$agb)

par(mfrow = c(1,3))
plot(bi_points_trees$bhd, bi_points_trees$total_vol)
plot(bi_points_trees$bhd, bi_points_trees$merch_vol)
plot(bi_points_trees$bhd, bi_points_trees$agb)
par(mfrow = c(1,1))

rm(vol_agb)

# preserve tree-level table before any aggregation
trees_base <- bi_points_trees



# 04: prepare plot geometries and include remeasured RTK-GNSS plots
#-------------------------------------------------------------------------------

# unique plot table used for RTK replacement and later filtering
plot_base <- unique(
  trees_base[, c('key', 'kspnr', 'abt', 'rw', 'hw')]
)

# conversion to sf object (DHDN / 3-degree Gauss-Kruger zone 3)
plot_base_gk <- sf::st_as_sf(
  plot_base, coords = c('rw', 'hw'), crs = 31467
  )

# transformation to ETRS89 / UTM zone 32N
plot_base_utm <- sf::st_transform(plot_base_gk, crs = 25832)

# merge remeasured plots into plot_base_utm
plot_base_utm$remeasured <- 'no'

# identify matching plots based on kspnr column
matching_plots <- plot_base_utm$kspnr %in% bi_plots_rtk$kspnr

# mark remeasured plots
plot_base_utm$remeasured[matching_plots] <- 'yes'

# add further RTK-related information
# NA for plots without RTK position
plot_base_utm <- dplyr::left_join(
  plot_base_utm,
  sf::st_drop_geometry(
    bi_plots_rtk[, c('kspnr', 'center_point_estimated',
                     'solution_status', 'measurement_date',
                     'height_loff', 'height_lon', 'height_diff')]
    ),
  by = 'kspnr'
)

# create sf object with non-RTK geometries (without RTK position)
remeasured_plots_non_rtk <- plot_base_utm[matching_plots, ]

# for plots that were remeasured,
# update their geometry with the more accurate RTK positions
if (any(matching_plots)) {
  
  # update geometry for remeasured plots
  for (i in which(matching_plots)) {
    kspnr_val <- plot_base_utm$kspnr[i]
    rtk_row <- which(bi_plots_rtk$kspnr == kspnr_val)
    if (length(rtk_row) > 0) {
      sf::st_geometry(plot_base_utm)[i] <- sf::st_geometry(bi_plots_rtk)[rtk_row[1]]
    }
  }
  
  # create sf object with RTK geometries (after RTK update)
  remeasured_plots_rtk <- plot_base_utm[matching_plots, ]
  
  cat('Updated', sum(matching_plots), 'plots with RTK-GNSS coordinates\n')
  cat('Created remeasured_plots_non_rtk:', nrow(remeasured_plots_non_rtk), 'plots with original geometries\n')
  cat('Created remeasured_plots_rtk:', nrow(remeasured_plots_rtk), 'plots with RTK geometries\n')
  
} else {
  
  cat('No matching plots found between inv_attr_plots_utm and bi_plots_rtk\n')
  remeasured_plots_non_rtk <- NULL
  remeasured_plots_rtk <- NULL
  
}



# 05: build tree-level table and aggregate per sample plot
#-------------------------------------------------------------------------------

# keep only trees from the RTK plot domain
plot_ids_rtk <- sf::st_drop_geometry(
  remeasured_plots_rtk[, c('key', 'kspnr')]
)
trees_rtk <- dplyr::semi_join(
  trees_base,
  plot_ids_rtk,
  by = c('key', 'kspnr')
)

# add final RTK coordinates and quality metadata to tree-level dataset
plot_coords_rtk <- sf::st_drop_geometry(
  dplyr::mutate(
    remeasured_plots_rtk,
    rtk_x = sf::st_coordinates(remeasured_plots_rtk)[, 1],
    rtk_y = sf::st_coordinates(remeasured_plots_rtk)[, 2]
  )[, c(
    'key', 'kspnr', 'remeasured', 'center_point_estimated',
    'solution_status', 'measurement_date', 'height_loff', 'height_lon',
    'height_diff', 'rtk_x', 'rtk_y'
  )]
)
inv_attr_trees_rtk <- dplyr::left_join(
  trees_rtk,
  plot_coords_rtk,
  by = c('key', 'kspnr')
)

# add original (non-RTK) coordinates and metadata to tree-level dataset
plot_ids_non_rtk <- sf::st_drop_geometry(
  remeasured_plots_non_rtk[, c('key', 'kspnr')]
)
trees_non_rtk <- dplyr::semi_join(
  trees_base,
  plot_ids_non_rtk,
  by = c('key', 'kspnr')
)
plot_coords_non_rtk <- sf::st_drop_geometry(
  dplyr::mutate(
    remeasured_plots_non_rtk,
    non_rtk_x = sf::st_coordinates(remeasured_plots_non_rtk)[, 1],
    non_rtk_y = sf::st_coordinates(remeasured_plots_non_rtk)[, 2]
  )[, c(
    'key', 'kspnr', 'remeasured', 'center_point_estimated',
    'solution_status', 'measurement_date', 'height_loff', 'height_lon',
    'height_diff', 'non_rtk_x', 'non_rtk_y'
  )]
)
inv_attr_trees_non_rtk <- dplyr::left_join(
  trees_non_rtk,
  plot_coords_non_rtk,
  by = c('key', 'kspnr')
)

# add column of leaf type
unique(trees_rtk$bagr)
trees_rtk$leaf_type <- ifelse(
  trees_rtk$bagr %in% c('EI', 'ALN', 'BU', 'ALH'),
  'deciduous',
  'coniferous'
)

# group sums of volume, AGB,
# and other forest inventory attributes like
# tree density, basal area, and QMD
trees_rtk <- trees_rtk %>%
  dplyr::group_by(key, kspnr) %>%
  dplyr::mutate(
    total_vol_ha = sum(total_vol * nha),
    merch_vol_ha = sum(merch_vol * nha),
    agb_ha = sum(agb * nha) / 1000,
    tree_density = mean(nha),
    basal_area_tree = (pi / 4) * (bhd / 100)^2,
    basal_area_ha = sum(basal_area_tree * nha, na.rm = T),
    dg = sqrt(sum(bhd^2 * nha, na.rm = T) / sum(nha, na.rm = T)),
    # assign dominant leaf type to each plot based on the basal area
    # share of deciduous vs. coniferous trees, considering only trees
    # from layer 1 (Hauptbestand) and 4 (Ueberhaelter)
    total_deciduous = sum(dplyr::if_else(
      leaf_type == 'deciduous' & bestschicht %in% c(1, 4),
      basal_area_tree * nha, 0, missing = 0), na.rm = T),
    total_coniferous = sum(dplyr::if_else(
      leaf_type == 'coniferous' & bestschicht %in% c(1, 4),
      basal_area_tree * nha, 0, missing = 0), na.rm = T),
    dominant_leaf_type = dplyr::case_when(
      total_deciduous > total_coniferous ~ 'deciduous',
      total_coniferous > total_deciduous ~ 'coniferous',
      TRUE                               ~ 'mixed'
    )) %>%
  dplyr::ungroup()

# extract unique forest inventory variables for all sample plots
inv_attr_plots <- unique(
  trees_rtk[, c(
    'key', 'kspnr', 'abt', 'total_vol_ha', 'merch_vol_ha',
    'agb_ha', 'tree_density', 'basal_area_ha', 'dg', 'dominant_leaf_type'
  )]
)
inv_attr_plots[is.na(inv_attr_plots)] <- 0

# attach aggregated attributes back to plot geometries
remeasured_plots_rtk <- dplyr::left_join(
  remeasured_plots_rtk,
  inv_attr_plots,
  by = c('key', 'kspnr', 'abt')
)
remeasured_plots_non_rtk <- dplyr::left_join(
  remeasured_plots_non_rtk,
  inv_attr_plots,
  by = c('key', 'kspnr', 'abt')
)

# summary statistics
summary(remeasured_plots_rtk)
summary_df <- as.data.frame(
  do.call(cbind, lapply(remeasured_plots_rtk, summary))
)
write.csv(
  summary_df,
  file.path(
    processed_data_dir,
    'forest_inventory',
    'summary_stats_inv_attr_plots_rtk.csv')
)
table(remeasured_plots_rtk$dominant_leaf_type)
table_df <- as.data.frame(table(remeasured_plots_rtk$dominant_leaf_type))
write.csv(
  table_df,
  file.path(
    processed_data_dir,
    'forest_inventory',
    'n_plots_dom_leaf_type.csv'),
  row.names = F
)

# boxplots
par(mfrow = c(2,3))
boxplot(remeasured_plots_rtk$total_vol_ha)
boxplot(remeasured_plots_rtk$merch_vol_ha)
boxplot(remeasured_plots_rtk$agb_ha)
boxplot(remeasured_plots_rtk$tree_density)
boxplot(remeasured_plots_rtk$basal_area_ha)
boxplot(remeasured_plots_rtk$dg)
par(mfrow = c(1,1))



# 06: save BI plots with the forest inventory attributes per sample plot
#-------------------------------------------------------------------------------

out_path <- file.path(processed_data_dir, 'forest_inventory')

# define file formats and corresponding save functions
file_formats <- list(
  rds = list(
    ext = '.RDS',
    save_func = saveRDS
  ),
  txt = list(
    ext = '.txt',
    save_func = function(data, file) 
      write.table(data, file, sep = '\t', row.names = F)
  ),
  gpkg = list(
    ext = '.gpkg',
    save_func = function(data, file) 
      sf::st_write(data, file, delete_dsn = T, quiet = T)
  )
)

# define datasets to save
datasets <- list(
  list(
    name = 'inv_attr_plots_rtk',
    data = remeasured_plots_rtk,
    drop_geom = T
  ),
  list(
    name = 'inv_attr_plots_non_rtk',
    data = remeasured_plots_non_rtk,
    drop_geom = T
  ),
  list(
    name = 'inv_attr_trees_rtk',
    data = inv_attr_trees_rtk,
    drop_geom = F
  ),
  list(
    name = 'inv_attr_trees_non_rtk',
    data = inv_attr_trees_non_rtk,
    drop_geom = F
  )
)

# loop through each file format
for (format_name in names(file_formats)) {
  
  format_info <- file_formats[[format_name]]
  
  # check if all files exist for this format
  all_files_exist <- all(sapply(datasets, function(ds) {
    file.exists(file.path(out_path, paste0(ds$name, format_info$ext)))
  }))
  
  if (!all_files_exist) {
    
    cat('Saving files in', format_name, 'format...\n')
    
    # save each dataset
    for (ds in datasets) {
      
      file_path <- file.path(out_path, paste0(ds$name, format_info$ext))
      
      # prepare data based on format
      if (format_name %in% c('rds', 'txt') && ds$drop_geom) {
        data_to_save <- sf::st_drop_geometry(ds$data)
      } else {
        data_to_save <- ds$data
      }
      
      # save using appropriate function
      format_info$save_func(data_to_save, file_path)
      cat('  Saved:', basename(file_path), '\n')
      
    }
    
  } else {
    
    cat('All', format_name, 'files already exist. Skipping.\n')
    
  }
  
}
