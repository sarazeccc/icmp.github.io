

corr_analysis <- function(df, title){
  subs_corr_df <- df |> 
    dplyr::select(VALUE_VALEUR, VARIABLE, SITE_NO) |> 
    spread(key = VARIABLE, value = VALUE_VALEUR) 
  corr_subs <- cor(subs_corr_df[,-1]|> drop_na())
  corrplot(corr_subs, method = 'color', tl.cex = 0.5,
           mar=c(0,0,1,0))
  return(subs_corr_df)
  
}

corr_analysis_final <- function(df, title){
  subs_corr_df <- df |> 
    dplyr::select(VALUE_VALEUR, VARIABLE, SITE_NO) |> 
    spread(key = VARIABLE, value = VALUE_VALEUR) 
  corr_subs <- cor(subs_corr_df[, c("O,P'-DDE (2,4'-DDE)",
                                    "P,P'-DDT (4,4'-DDT)",
                                    "MIREX", "ALDRIN",
                                    "ENDOSULFAN SULPHATE",
                                    "MAGNESIUM",
                                    "POTASSIUM", 
                                    "MOLYBDENUM",
                                    "ZINC")] |> drop_na())
  # Define categories for your variables
  categories <- c("PESTICIDE", "INSECTICIDE", "INSECTICIDE", 
                  "INSECTICIDE", "INSECTICIDE", "MAJOR ION", 
                  "MAJOR ION", "METAL", "METAL")
  

  corrplot(corr_subs, method = "circle", type = "lower",
           diag = FALSE,  tl.col = "black", tl.srt = 45)  
  # corrplot(corr_subs, method = 'color', tl.cex = 0.5,
  #         mar=c(0,0,1,0), order = 'original')
  return(subs_corr_df)
  
}

extract_data <- function(substance_name){
  # Load data
  water_quality <-
    read.csv('data/Great-Lakes-Sediment-Monitoring-Plan-St-Clair-River-2014 Data.csv')
  monitoring_sites <-
    read.csv('data/Great-Lakes-Sediment-Monitoring-Plan-Sites.csv')
  substance_types <-
    read.csv('data/Great-Lakes-Sediment-Monitoring-Plan-VMVs.csv')
  tentative_sites <- read.csv('data/tentative_1k_points.csv')
  
  substances_vector <- c(substance_name,  
                         'MAGNESIUM TOTAL RECOVERABLE',
                         'POTASSIUM TOTAL RECOVERABLE',
                         'MOLYBDENUM TOTAL RECOVERABLE',
                         'ZINC TOTAL RECOVERABLE',
                         "MIREX", "ALDRIN", "P,P'-DDT (4,4'-DDT)",
                         "ENDOSULFAN SULPHATE")
  
  
  
  # Separate data
  substances_types_unique <-
    unique(substance_types[, c('VARIABLE', 'VARIABLE_TYPE')])
  locs_quality <-
    merge(water_quality, monitoring_sites[, c('SITE_NO', 'LATITUDE', 'LONGITUDE')],
          by = 'SITE_NO')
  wq_type <-
    merge(
      locs_quality,
      substances_types_unique,
      by = 'VARIABLE',
      all.x = TRUE,
      all.y = FALSE
    )
  
  # Separate substance of interest
  substance_df <-
    wq_type |> dplyr::filter(VARIABLE %in% substances_vector)
  substance_df$VARIABLE <- gsub('\\s*TOTAL.*', "", substance_df$VARIABLE)
  
  substance_df_filtered <-
    substance_df |> dplyr::select(VARIABLE, SITE_NO, VALUE_VALEUR, LATITUDE, LONGITUDE)
  
  variable_name <- substance_df$VARIABLE
  
  substance_df_spread <-
    spread(substance_df_filtered, key = VARIABLE, value = VALUE_VALEUR) |>
    dplyr::select(SITE_NO, LONGITUDE, LATITUDE, all_of(variable_name))
  
  substances_filtered <- substance_df_spread |> drop_na()
  print(nrow(substances_filtered))
  
  
  return(substances_filtered)
  
}
add_border <- function(df){
  canadian_number_sites <- c('5060', '5062', '5063', '5061', '5075', '5081', 
                             '5077', '5070', '5073', '5066', '5080', '5072',
                             '5071', '5068', '5067', '5069', '5074', '5076',
                             '5079', '5078', '5097', '5092', '5094', '5095', 
                             '5093', '5098', '5100', '5101', '5091', '5096',
                             '5102')
  canadian_number_sites <- paste0('ON02GG', canadian_number_sites)
  
  df[df$SITE_NO %in% canadian_number_sites, 'Border'] <- 1
  df[df$SITE_NO %notin% canadian_number_sites, 'Border'] <- 0
  
  return(df)
  
}
munge_coordinates <- function(data_spread){
  
  coordinates(data_spread) <- ~ LONGITUDE + LATITUDE
  proj4string(data_spread) <- CRS("+proj=longlat")
  substance_utm <-
    spTransform(data_spread, CRS("+proj=utm +zone=17 +datum=WGS84"))
  substance_df <- as.data.frame(substance_utm)
  colnames(substance_df)[colnames(substance_df) == 'coords.x1'] <- 'X'
  colnames(substance_df)[colnames(substance_df) == 'coords.x2'] <- 'Y'
  
  # Change from meters to kilometers
  substance_df$X_km <- substance_df$X / 1000
  substance_df$Y_km <- substance_df$Y / 1000
  coords_scaled <- scale(substance_df[, c('X_km', 'Y_km')])
  
  dist_matrix <- dist(substance_df[, c('X_km', 'Y_km')]) |> as.matrix()
  
  # Scale covariates
  coords_scaled <-
    scale(substance_df[, c('X_km', 'Y_km')]) |> as.matrix()
  
  
  return(list('dist_matrix' = dist_matrix, 
              'coords_scaled' = coords_scaled,
              'substance_df' = substance_df))
  
  
}
data_preprocessing <- function(substance_of_interest, substances_filtered){
 
  output_coordinates <- munge_coordinates(substances_filtered)
  dist_matrix <- output_coordinates[['dist_matrix']]
  coords_scaled <- output_coordinates[['coords_scaled']]
  substances <- output_coordinates[['substance_df']]
  # Filter out the rows where the substance of interest is NA
  #  substances_interest <- substances[!is.na(substances[ ,substances_column]), ]
  
  d_max <- max(dist_matrix)
  
  phi_mean <-  (d_max / 6)
  phi_lower <- (5.7333 / 2) /3
  phi_upper <-  (25.5692 / 2) / 3
  
  
  ##  Get explanatory variables
  major_ions_df <- substances_filtered[, c('MAGNESIUM', 'POTASSIUM')]
  metals_df <- substances_filtered[,c('MOLYBDENUM', 'ZINC')]
  insecticides_df <- substances_filtered[,c("MIREX", "ALDRIN", 
                                            "P,P'-DDT (4,4'-DDT)",
                                            "ENDOSULFAN SULPHATE")]
  
  metals_scaled <- scale(metals_df) |> as.data.frame()
  insecticides_scaled <- scale(insecticides_df) |> as.data.frame()
  major_ions_scaled <- scale(major_ions_df) |> as.data.frame()
  n_var <- ncol(major_ions_scaled) + ncol(metals_scaled)  + ncol(insecticides_scaled) + 4
  n_obs <- nrow(substances_filtered)
  print(n_obs)
  
  
  ## Make a vector with the values of interest
  y_log <- log(substances_filtered[,substance_of_interest])
  
  # Initialize lists for nimble model
  
  const <- list(n_obs = n_obs, n_var = n_var)
  
  data <- list(
    dist_matrix = dist_matrix,
    y_log = y_log,
    phi_mean = phi_mean,
    coords = coords_scaled,
    border_dummy = substances_spread_border[,'Border'],
    major_ions = major_ions_scaled,
    metals = metals_scaled,
    insects = insecticides_scaled
  )
  
  init <- list(sigma2 = 0.5,
               tau2 = 0.1,
               phi = 0.1)
  
  return(list('const' = const,
              'data' = data,
              'init' = init))
}

prepare_prediction_sites <- function(prediction_sites, coords_scaled){
  coordinates(prediction_sites) <- ~ LONGITUDE + LATITUDE
  proj4string(prediction_sites) <- CRS("+proj=longlat")
  prediction_utm <- spTransform(prediction_sites,
                                CRS("+proj=utm +zone=17 +datum=WGS84"))
  prediction_df <- as.data.frame(prediction_utm)
  colnames(prediction_df)[colnames(prediction_df) == 'coords.x1'] <- 'X'
  colnames(prediction_df)[colnames(prediction_df) == 'coords.x2'] <- 'Y'
  
  # Change from meters to kilometers
  prediction_df$X_km <- prediction_df$X/1000
  prediction_df$Y_km <- prediction_df$Y/1000
  
  # Scale parameters with mean and scale from data
  means_to_scale <- attr(coords_scaled,"scaled:center")
  sd_to_scale <- attr(coords_scaled,"scaled:scale")
  prediction_coords_scaled <- scale(prediction_df[,c('X_km', 'Y_km')], 
                                    center = means_to_scale, 
                                    scale = sd_to_scale )
  
  
  return(list('predictions_coords_scaled' = prediction_coords_scaled,
              'pred_df' = prediction_df))
  
}

# Prepare distance matrix

prepare_dist_matrix <- function(prediction_df, substance_df){
  
  dist_11 <- dist(prediction_df[,c('X_km', 'Y_km')]) |> as.matrix()
  dist_12 <- cdist(prediction_df[,c('X_km', 'Y_km')], 
                   substance_df[,c('X_km', 'Y_km')]) |> as.matrix()
  dist_22 <- dist(substance_df[,c('X_km', 'Y_km')]) |> as.matrix()
  
  return(list('dist_11' = dist_11,
              'dist_12' = dist_12,
              'dist_22' = dist_22))
  
  
}



