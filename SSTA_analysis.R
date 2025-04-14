
# Data load and modification ----------------------------------------------

ssta_data = read.csv("/Users/megha/Documents/MIZZOU/SPATIAL/SSTA data analysis//datafile-2.csv")
#View(ssta_data)

# Cleanup the data #

nchar(ssta_data[,1])

# For longitude
simple_fun <- function(x){
  len <- nchar(x)
  last_char <- substr(x, start=len, stop = len)
  if(last_char!="E" & last_char!="W"){
    lon_tmp <- as.numeric(substr(x, start=1, stop = len))
    lon_tmp <- lon_tmp
  }else{
    lon_tmp <- as.numeric(substr(x, start=1, stop = len-1))
    if(last_char=="E") sign <- 1 else {sign <- -1; 
    lon_tmp <- sign*lon_tmp}
  }
  
  return(lon_tmp)
}

Simple_fun <- Vectorize(simple_fun)
Longitude <- Simple_fun(ssta_data[,1])


# For latitude
simple_fun_lat <- function(x){
  len <- nchar(x)
  lat_tmp <- as.numeric(substr(x, start=1, stop = len-1))
  if(substr(x, start=len, stop = len)=="N") sign <- 1 else  sign <- -1
  lat_tmp <- sign*lat_tmp
  return(lat_tmp)
}

Simple_fun_lat <- Vectorize(simple_fun_lat)
Latitude <- Simple_fun_lat(ssta_data[,2])


ssta_data$longitude <- Longitude
ssta_data$latitude <- Latitude

#View(ssta_data)
dim(ssta_data)
head(ssta_data)


# Plot of the data --------------------------------------------------------

library(transformr)
library(ggplot2)
library(gganimate)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(fields)

str(ssta_data)

# Load world map data
world_map <- ne_countries(scale = "medium", returnclass = "sf")

# Filter data to ensure it covers the desired time range
ssta_data$Time <- as.Date(paste0(ssta_data$Time, "-01"))  # Convert Time to Date format
ssta_data <- ssta_data[ssta_data$Time >= as.Date("1970-01-01") & ssta_data$Time <= as.Date("2003-03-01"), ]


# Load land polygons
land <- ne_countries(scale = "medium", returnclass = "sf")

library(sf)
# Convert your data into an sf object
ssta_sf <- st_as_sf(ssta_data, coords = c("longitude", "latitude"), crs = 4326)

# Remove points that intersect with land
ocean_points <- st_difference(ssta_sf, st_union(land))

# Convert back to a data frame and extract coordinates
ssta_ocean <- ocean_points %>%
  st_drop_geometry() %>%
  cbind(st_coordinates(ocean_points))

# Rename columns for longitude and latitude
colnames(ssta_ocean)[c(ncol(ssta_ocean) - 1, ncol(ssta_ocean))] <- c("longitude", "latitude")


# Create the plot
animation_plot_1 <- ggplot() +
  geom_sf(data = land, fill = "gray90", color = "black") + # Add land polygons
  geom_tile(data = ssta_ocean, aes(x = longitude, y = latitude, fill = Sea.Surface.Temperature.Anomaly), inherit.aes = FALSE) + # Ocean data
  scale_fill_viridis_c(option = "plasma") + # Vibrant color scale
  coord_sf(xlim = c(-180, 180), ylim = c(-90, 90), expand = FALSE) + # Map limits
  theme_minimal() +
  labs(
    title = "Sea Surface Temperature Anomaly: {frame_time}",
    x = "Longitude",
    y = "Latitude",
    fill = "SSTA (°C)"
  ) +
  transition_time(Time) + # Animate over time
  ease_aes('linear')      # Smooth transitions

# Render the animation
animate(animation_plot_1, nframes = 600, fps = 6, width = 1000, height = 800)



# Spatial EOF analysis ----------------------------------------------------

library(tidyr)
library(dplyr)

T = length(unique(ssta_data$Time)) 
m = nrow(unique(data.frame(ssta_data$latitude,ssta_data$longitude)))

# Comment : 399 time points and 2520 locations

# Create a location identifier (e.g., combining lat and lon into one string)
ssta_data <- ssta_data %>%
  mutate(location = paste(latitude, longitude, sep = "_"))


# Pivot data to wide format
ssta_data_spatial_wide <- ssta_data %>%
  pivot_wider(names_from = Time, values_from = Sea.Surface.Temperature.Anomaly)


#View(ssta_data)

ssta_data_spatial_wide = data.frame(ssta_data_spatial_wide)

ssta_data_spatial_wide[,6:ncol(ssta_data_spatial_wide)]

#View(ssta_data)

Z = t(as.data.frame(ssta_data_spatial_wide[,6:ncol(ssta_data_spatial_wide)]))
#View(Z_mat)
dim(Z) # Txm matrix

mu <- colMeans(Z, na.rm = TRUE)

one_mat <- matrix(rep(1,T),nrow = T, ncol = 1) # Tx1 vector of 1s
Z_tilde <- (Z - one_mat%*%t(mu))/(sqrt(T-1))
dim(Z_tilde)
C0_Z <- t(Z_tilde)%*%Z_tilde
C0_Z[1,1] #0.1976449

# Perform eigen decomposition
eigen_result <- eigen(C0_Z)

# Extract eigenvalues and eigenvectors
eigenvalues <- eigen_result$values
eigenvectors<- eigen_result$vectors
dim(eigenvectors)

# Sort eigenvalues and corresponding eigenvectors in decreasing order
sorted_indices <- order(eigenvalues, decreasing = TRUE)
sorted_eigenvalues <- eigenvalues[sorted_indices]
sorted_eigenvectors <- eigenvectors[, sorted_indices]

# Select the first 10 EOFs (eigenvectors)
num_eofs <- 10
EOFs <- sorted_eigenvectors[, 1:num_eofs]
dim(EOFs)

# Comment: Same EOF values as approach 1

# Now I will use these EOF to observe Transient Growth

## Project the data onto the first 10 EOFs

dim(EOFs)
dim(Z)

at <- Z %*% EOFs
dim(at)
#View(at)




# Var-like model fit and estimation of M with different lags --------------

# Step 1: Compute C0_a (covariance matrix at lag 0)

mu_hat <- colMeans(at)
length(mu_hat)
a_centered <- scale(at, center = mu_hat, scale = FALSE)
C0_a <- var(a_centered)
dim(C0_a)
#View(C0_a)

# Step 2: 

# Compute C1_a (lag-1 covariance matrix)
C1_a <- cov(a_centered[-1, ], a_centered[-T, ])
dim(C1_a)
#View(C1_a)

# Compute C3_a (lag-3 covariance matrix)
C3_a <- cov(a_centered[-c(1,2,3), ], a_centered[-c(T,T-1,T-2), ])
dim(C3_a)

# Compute C6_a (lag-6 covariance matrix)
C6_a <- cov(a_centered[-c(1:6), ], a_centered[-c(T,T-1,T-2,T-3,T-4,T-5), ])
dim(C6_a)

# Step 3: 

# Calculate the method-of-moments estimator M(lag-1)
M_freq_lag1 <- C1_a %*% solve(C0_a)
M_freq_lag1 

library(xtable)
latex_code_M_freq_lag1 <- print(xtable(M_freq_lag1), include.rownames=TRUE, include.colnames=TRUE, floating=FALSE, hline.after=NULL, print.results = FALSE)

# Calculate the method-of-moments estimator M(lag-3)
M_freq_lag3 <- C3_a %*% solve(C0_a)
M_freq_lag3

# Calculate the method-of-moments estimator M(lag-6)
M_freq_lag6 <- C6_a %*% solve(C0_a)
M_freq_lag6


# Some metrices of Non-normality ------------------------------------------
library(matrixcalc)

## By considering Frobenius metric
mu_2 <- function(A){
  A_star <- t(A)
  AA_star <- A %*% A_star
  A_star_A <- A_star %*% A
  diff_matrix <- AA_star - A_star_A
  frobenius_norm <- norm(diff_matrix, type = "F")
  sqrt_frobenius_norm <- sqrt(frobenius_norm)
  return(sqrt_frobenius_norm)
}

## By considering Henrici metric
mu_3 <- function(A) {
  frobenius_norm <- norm(A, type = "F")
  eigen_vals <- eigen(A)$values
  h_A <- sqrt(frobenius_norm^2 - sum(abs(eigen_vals)^2))
  return(h_A)
}

## By considering Ruhe's metric
mu_4 <- function(A) {
  singular_values <- svd(A)$d
  eigenvalues <- abs(eigen(A)$values)
  sigma_i <- singular_values
  lambda_i <- eigenvalues
  difference <- max(abs(sigma_i - lambda_i))
  return(difference)
}

# lag 1

mu_2(M_freq_lag1)
mu_3(M_freq_lag1)
mu_4(M_freq_lag1)

# lag 3

mu_2(M_freq_lag3)
mu_3(M_freq_lag3)
mu_4(M_freq_lag3)

# lag 6

mu_2(M_freq_lag6)
mu_3(M_freq_lag6)
mu_4(M_freq_lag6)

# Comment : Non-normal = evidence of transient growth


# Moving average estimation of M ------------------------------------------

library(stats)
at_ts <- ts(at, frequency=12) 

### 1-year MA

at_ma_1 <- stats::filter(at_ts, rep(1/12, 12), sides=1)
# Remove NA 
at_ma_1 <- at_ma_1[12:nrow(at_ma_1), ]

mu_hat_ma_1 <- colMeans(at_ma_1)
a_centered_ma_1 <- scale(at_ma_1, center = mu_hat_ma_1, scale = FALSE)
T_ma_1 <- nrow(a_centered_ma_1)
C0_a_ma_1 <- var(a_centered_ma_1)  

## Lag-1
C1_a_ma_1 <- cov(a_centered_ma_1[-1, ], a_centered_ma_1[-T_ma_1, ])
M1_hat_freq_ma_1 <- C1_a_ma_1 %*% solve(C0_a_ma_1)
print(M1_hat_freq_ma_1)

## Lag-3
C3_a_ma_1 <- cov(a_centered_ma_1[-c(1,2,3), ], a_centered_ma_1[-c(T_ma_1,T_ma_1-1,T_ma_1-2), ])
M3_hat_freq_ma_1 <- C3_a_ma_1 %*% solve(C0_a_ma_1)
print(M3_hat_freq_ma_1)

## Lag-6
C6_a_ma_1 <- cov(a_centered_ma_1[-c(1:6), ], a_centered_ma_1[-c(T_ma_1,T_ma_1-1,T_ma_1-2,T_ma_1-3,T_ma_1-4,T_ma_1-5), ])
M6_hat_freq_ma_1 <- C6_a_ma_1 %*% solve(C0_a_ma_1)
print(M6_hat_freq_ma_1)


### 4-year MA

at_ma_4 <- stats::filter(at_ts, rep(1/48, 48), sides=1)
# Remove NA 
at_ma_4 <- at_ma_4[48:nrow(at_ma_4), ]

mu_hat_ma_4 <- colMeans(at_ma_4)
a_centered_ma_4 <- scale(at_ma_4, center = mu_hat_ma_4, scale = FALSE)
T_ma_4 <- nrow(a_centered_ma_4)
C0_a_ma_4 <- var(a_centered_ma_4)  

## Lag-1
C1_a_ma_4 <- cov(a_centered_ma_4[-1, ], a_centered_ma_4[-T_ma_4, ])
M1_hat_freq_ma_4 <- C1_a_ma_4 %*% solve(C0_a_ma_4)
print(M1_hat_freq_ma_4)

## Lag-3
C3_a_ma_4 <- cov(a_centered_ma_4[-c(1,2,3), ], a_centered_ma_4[-c(T_ma_4,T_ma_4-1,T_ma_4-2), ])
M3_hat_freq_ma_4 <- C3_a_ma_4 %*% solve(C0_a_ma_4)
print(M3_hat_freq_ma_4)

## Lag-6
C6_a_ma_4 <- cov(a_centered_ma_4[-c(1:6), ], a_centered_ma_4[-c(T_ma_4,T_ma_4-1,T_ma_4-2,T_ma_4-3,T_ma_4-4,T_ma_4-5), ])
M6_hat_freq_ma_4 <- C6_a_ma_4 %*% solve(C0_a_ma_4)
print(M6_hat_freq_ma_4)

### 5-year MA

at_ma_5 <- stats::filter(at_ts, rep(1/60, 60), sides=1)
# Remove NA 
at_ma_5 <- at_ma_5[60:nrow(at_ma_5), ]

mu_hat_ma_5 <- colMeans(at_ma_5)
a_centered_ma_5 <- scale(at_ma_5, center = mu_hat_ma_5, scale = FALSE)
T_ma_5 <- nrow(a_centered_ma_5)
C0_a_ma_5 <- var(a_centered_ma_5)  

## Lag-1
C1_a_ma_5 <- cov(a_centered_ma_5[-1, ], a_centered_ma_5[-T_ma_5, ])
M1_hat_freq_ma_5 <- C1_a_ma_5 %*% solve(C0_a_ma_5)
print(M1_hat_freq_ma_5)

## Lag-3
C3_a_ma_5 <- cov(a_centered_ma_5[-c(1,2,3), ], a_centered_ma_5[-c(T_ma_5,T_ma_5-1,T_ma_5-2), ])
M3_hat_freq_ma_5 <- C3_a_ma_5 %*% solve(C0_a_ma_5)
print(M3_hat_freq_ma_5)

## Lag-6
C6_a_ma_5 <- cov(a_centered_ma_5[-c(1:6), ], a_centered_ma_5[-c(T_ma_5,T_ma_5-1,T_ma_5-2,T_ma_5-3,T_ma_5-4,T_ma_5-5), ])
M6_hat_freq_ma_5 <- C6_a_ma_5 %*% solve(C0_a_ma_5)
print(M6_hat_freq_ma_5)

### 7-year MA
at_ma_7 <- stats::filter(at_ts, rep(1/84,84), sides=1)
# Remove NA 
at_ma_7 <- at_ma_7[84:nrow(at_ma_7), ]

mu_hat_ma_7 <- colMeans(at_ma_7)
a_centered_ma_7 <- scale(at_ma_7, center = mu_hat_ma_7, scale = FALSE)
T_ma_7 <- nrow(a_centered_ma_7)
C0_a_ma_7 <- var(a_centered_ma_7)  

## Lag-1
C1_a_ma_7 <- cov(a_centered_ma_7[-1, ], a_centered_ma_7[-T_ma_7, ])
M1_hat_freq_ma_7 <- C1_a_ma_7 %*% solve(C0_a_ma_7)
print(M1_hat_freq_ma_7)

## Lag-3
C3_a_ma_7 <- cov(a_centered_ma_7[-c(1,2,3), ], a_centered_ma_7[-c(T_ma_7,T_ma_7-1,T_ma_7-2), ])
M3_hat_freq_ma_7 <- C3_a_ma_7 %*% solve(C0_a_ma_7)
print(M3_hat_freq_ma_7)

## Lag-6
C6_a_ma_7 <- cov(a_centered_ma_7[-c(1:6), ], a_centered_ma_7[-c(T_ma_7,T_ma_7-1,T_ma_7-2,T_ma_7-3,T_ma_7-4,T_ma_7-5), ])
M6_hat_freq_ma_7 <- C6_a_ma_7 %*% solve(C0_a_ma_7)
print(M6_hat_freq_ma_7)


# Replication of GAMM (Nicol-Harper) as GAM model -------------------------

library(matrixcalc)
library(ggplot2)
library(lubridate)

# Define time range: Jan 1970 to Mar 2003
start_date <- as.Date("1970-01-01")
end_date <- as.Date("2003-03-01")
time_points <- seq(start_date, end_date, by = "1 month")

# Define slice length (19 months per segment)
segment_length <- 19
num_segments <- ceiling(nrow(at) / segment_length)

# Create index groups for slicing
segment_indices <- split(1:nrow(at), ceiling(seq_along(1:nrow(at)) / segment_length))

# List to store M matrices
M_list_freq <- list()

# Loop through each 19-month slice and compute M
for (i in seq_along(segment_indices)) {
  at_i <- at[segment_indices[[i]], ]  # Extract sub-matrix
  T_i <- nrow(at_i)  # Adjust time length for this segment
  
  # Step 1: Compute C0_a (covariance matrix at lag 0)
  mu_hat <- colMeans(at_i)
  a_centered <- scale(at_i, center = mu_hat, scale = FALSE)
  C0_a <- var(a_centered)
  
  # Step 2: Compute C1_a (lag-1 covariance matrix)
  if (T_i > 1) {
    C1_a <- cov(a_centered[-1, ], a_centered[-T_i, ])
  } else {
    C1_a <- matrix(0, ncol = ncol(at), nrow = ncol(at))  # Handle single-row case
  }
  
  # Step 3: Compute the method-of-moments estimator M1_hat
  if (det(C0_a) != 0) {
    M_freq_lag1 <- C1_a %*% solve(C0_a)
  } else {
    M_freq_lag1 <- matrix(NA, ncol = ncol(at), nrow = ncol(at))  # Handle singular matrix
  }
  
  # Store result
  M_list_freq[[paste0("M", i)]] <- M_freq_lag1
}

# Compute mu_2 for each M matrix
mu_2_values <- sapply(M_list_freq, function(M) {
    return(mu_2(M))
})


# Compute mu_3 for each M matrix
mu_3_values <- sapply(M_list_freq, function(M) {
    return(mu_3(M))
  })

# Compute mu_4 for each M matrix
mu_4_values <- sapply(M_list_freq, function(M) {
    return(mu_4(M))
})


# Find the midpoint of each time slice (10th month)
mid_indices <- sapply(segment_indices, function(idx) {
    return(idx[10])  #  10th month is the midpoint
})

# Convert midpoint indices to actual date values
mid_dates <- time_points[mid_indices]

# Convert results to a data frame
mu_2_df <- data.frame(Date = mid_dates, mu_2_Value = mu_2_values)
mu_3_df <- data.frame(Date = mid_dates, mu_3_Value = mu_3_values)
mu_4_df <- data.frame(Date = mid_dates, mu_4_Value = mu_4_values)

library(mgcv)   # For GAM modeling
library(ggplot2)

## using mu_2
# Check Date is numeric 
mu_2_df$Year <- as.numeric(format(mu_2_df$Date, "%Y")) + (as.numeric(format(mu_2_df$Date, "%m")) / 12)
# Remove any NA values in mu_2_df
mu_2_df <- na.omit(mu_2_df)
# Fit a GAM model (without random effect for Location)
gam_model_mu_2 <- gam(mu_2_values ~ s(Year),  # Smooth function for time
                        data = mu_2_df)
# Add predictions to the dataframe
mu_2_df$Fitted_mu_2 <- predict(gam_model_mu_2, newdata = mu_2_df)


## using mu_3
# Check Date is numeric 
mu_3_df$Year <- as.numeric(format(mu_3_df$Date, "%Y")) + (as.numeric(format(mu_3_df$Date, "%m")) / 12)
# Remove any NA values in mu_2_df
mu_3_df <- na.omit(mu_3_df)
# Fit a GAM model (without random effect for Location)
gam_model_mu_3 <- gam(mu_3_values ~ s(Year),  # Smooth function for time
                        data = mu_3_df)
# Add predictions to the dataframe
mu_3_df$Fitted_mu_3 <- predict(gam_model_mu_3, newdata = mu_3_df)


## using mu_4
# Check Date is numeric 
mu_4_df$Year <- as.numeric(format(mu_4_df$Date, "%Y")) + (as.numeric(format(mu_4_df$Date, "%m")) / 12)
# Remove any NA values in mu_2_df
mu_4_df <- na.omit(mu_4_df)
# Fit a GAM model (without random effect for Location)
gam_model_mu_4 <- gam(mu_4_values ~ s(Year),  # Smooth function for time
                        data = mu_4_df)
# Add predictions to the dataframe
mu_4_df$Fitted_mu_4 <- predict(gam_model_mu_4, newdata = mu_4_df)


# Combine all data frames into one
library(dplyr)
library(ggplot2)

mu_2_df$Metric <- "Frobenius"
mu_3_df$Metric <- "Henrici"
mu_4_df$Metric <- "Ruhe"

# Rename columns to have consistent names
mu_2_df <- mu_2_df %>% rename(Value = mu_2_Value, Fitted = Fitted_mu_2)
mu_3_df <- mu_3_df %>% rename(Value = mu_3_Value, Fitted = Fitted_mu_3)
mu_4_df <- mu_4_df %>% rename(Value = mu_4_Value, Fitted = Fitted_mu_4)

# Combine all data frames
combined_df <- bind_rows(mu_2_df, mu_3_df, mu_4_df)

# Plot 
ggplot(combined_df, aes(x = Date, y = Value)) +
  geom_line(color = "blue", alpha = 0.5) +  # Original data points
  geom_line(aes(y = Fitted), color = "darkred", size = 1.2) +  # Fitted GAM trend
  labs(title = "Non-normality Metrics Over Time",
       x = "Year",
       y = "Non-normality Value",
       color = "Metric") +
  theme_minimal() +
  facet_wrap(~ Metric, scales = "free_y")


# Load NINO 3.4 and ANOM 3.4 index data and merge with ssta data (at)------------------------------------------------------

library(tibble)
library(zoo)
nino_url <- "https://www.cpc.ncep.noaa.gov/data/indices/ersst5.nino.mth.91-20.ascii"
nino_data <- readLines(nino_url, warn = FALSE)
# Convert into a structured dataframe
nino_df <- read.table(text = paste(nino_data, collapse = "\n"), header = FALSE, fill = TRUE)
dim(nino_df)

# Assign correct column names based on dataset headers
names(nino_df) <- c("Year", "Month", "NINO1+2", "ANOM1+2", "NINO3", "ANOM3", 
                    "NINO4", "ANOM4", "NINO3.4", "ANOM3.4")

names(nino_df)


nino_df <- nino_df %>%
  mutate(Time = as.yearmon(paste(Year, Month), format = "%Y %m")) %>%
  dplyr::select(Time, NINO3.4, ANOM3.4)

head(nino_df)


# Convert row names to a column for the first data frame
at <- data.frame(at) %>% rownames_to_column(var = "Time")
head(at)

# Standardize the format of the "Time" column
at$Time <- gsub("^X", "", at$Time)  # Remove the 'X' prefix
at$Time <- gsub("\\.", "-", at$Time)  # Replace dots with hyphens

# Convert the Time to Date class 
at$Time <- as.Date(at$Time, format = "%Y-%m-%d")
nino_df$Time <- as.Date(paste0(nino_df$Time, "-01"), format="%b %Y-%d")

# Ensure "Time" is in the same format in both data frames
at$Time <- as.character(at$Time)  
nino_df$Time <- as.character(nino_df$Time)

# Merge the data frames based on "Time"
merged_data <- left_join(at, nino_df[242:640,], by = "Time")

head(merged_data)
dim(merged_data)

# Plot of spatial-EOFs  ------------------------------------------------------------
head(ssta_data)
ssta_data <- ssta_data %>%
  mutate(location = paste(latitude, longitude, sep = "_"))

library(dplyr)

# Extract unique locations
unique_locations <- ssta_data %>%
  dplyr::select(longitude, latitude, location) %>%  
  distinct()  

dim(unique_locations)  
head(unique_locations)

spatial_EOF_data <- unique_locations
colnames(EOFs) <-  paste0("EOF", 1:num_eofs)
spatial_EOF_data <- cbind(spatial_EOF_data,EOFs)
head(spatial_EOF_data)

# Plot of EOF
library(ggplot2)
ggplot(spatial_EOF_data %>% mutate(longitude = ifelse(longitude < 0, longitude + 360, longitude))) + 
  geom_tile(aes(x = longitude, y = latitude, fill = EOF1)) +
  scale_fill_viridis_c(name = "degC") +  
  theme_bw() +
  labs(title = "Spatial EOF plot",
       x = "Longitude (deg)",
       y = "Latitude (deg)")+
  coord_cartesian(expand = FALSE)  


# Correlation between non-normality and ANOM 3.4 index  -------------------

mid_data_ANOM <- merged_data$ANOM3.4[mid_indices]

mu_2_df <- cbind(mu_2_df, ANOM3.4 = mid_data_ANOM)
mu_3_df <- cbind(mu_3_df, ANOM3.4 = mid_data_ANOM )
mu_4_df <- cbind(mu_4_df, ANOM3.4 = mid_data_ANOM )

library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)

mu_2_df <- mu_2_df[, !colnames(mu_2_df) %in% "ANOM3.4"] # Remove if exists
mu_3_df <- mu_3_df[, !colnames(mu_3_df) %in% "ANOM3.4"] # Remove if exists
mu_4_df <- mu_4_df[, !colnames(mu_4_df) %in% "ANOM3.4"] # Remove if exists

mu_2_df$ANOM3.4 <- mid_data_ANOM
mu_3_df$ANOM3.4 <- mid_data_ANOM
mu_4_df$ANOM3.4 <- mid_data_ANOM

# Normalize ANOM3.4
normalize_anom <- function(df) {
  df$ANOM3.4 <- as.numeric(df$ANOM3.4)
  df$ANOM3.4_normalized <- (df$ANOM3.4 - min(df$ANOM3.4, na.rm = TRUE)) / 
    (max(df$ANOM3.4, na.rm = TRUE) - min(df$ANOM3.4, na.rm = TRUE)) * 
    (max(df$Value, na.rm = TRUE) - min(df$Value, na.rm = TRUE)) + 
    min(df$Value, na.rm = TRUE)
  return(df)
}

mu_2_df <- normalize_anom(mu_2_df)
mu_3_df <- normalize_anom(mu_3_df)
mu_4_df <- normalize_anom(mu_4_df)

# Correlation
calculate_correlation <- function(df) {
  cor_value <- cor(df$Value, df$ANOM3.4, use = "complete.obs")
  return(paste("Correlation: ", round(cor_value, 4)))
}

cor_text_mu2 <- calculate_correlation(mu_2_df)
cor_text_mu3 <- calculate_correlation(mu_3_df)
cor_text_mu4 <- calculate_correlation(mu_4_df)

# Plotting function to create plots with dual axes
create_plot <- function(df, cor_text, title_suffix) {
  ggplot(df, aes(x = Date)) +
    geom_line(aes(y = Value, color = "Value"), size = 1) +
    geom_point(aes(y = Value, color = "Value"), size = 2) +
    geom_line(aes(y = ANOM3.4_normalized, color = "ANOM3.4"), size = 1) +
    geom_point(aes(y = ANOM3.4_normalized, color = "ANOM3.4"), size = 2) +
    labs(title = paste("mu", title_suffix, "(Frobenius metric) and ANOM Values Over Time"),
         x = "Year", y = "Value / Normalized ANOM3.4", color = "Metric") +
    scale_x_date(date_labels = "%Y", date_breaks = "5 years") +
    annotate("text", x = min(df$Date), y = max(df$Value, na.rm = TRUE),
             label = cor_text, hjust = 0, vjust = 1.5, size = 5, color = "black") +
    theme_minimal() +
    scale_color_manual(values = c("Value" = "lightblue4", "ANOM3.4" = "pink3"))
}

# Create plots for each dataset
p1 <- create_plot(mu_2_df, cor_text_mu2, 2)
p2 <- create_plot(mu_3_df, cor_text_mu3, 3)
p3 <- create_plot(mu_4_df, cor_text_mu4, 4)

# Arrange and display all plots together
grid.arrange(p1, p2, p3, ncol = 1)

# Time EOF analysis -------------------------------------------------------

head(ssta_data)
dim(ssta_data)
ssta_data_time_wide  <- ssta_data[,-c(1,2)]  %>%
  pivot_wider(names_from = c(latitude,longitude,location) , values_from = Sea.Surface.Temperature.Anomaly)

dim(ssta_data_time_wide)
ssta_data_time_wide = data.frame(ssta_data_time_wide)

final_ssta_data <- ssta_data_time_wide[,2:ncol(ssta_data_time_wide)]
dim(final_ssta_data)

Z_time_EOF = t(as.data.frame(final_ssta_data))
#View(Z_time_EOF)
dim(Z_time_EOF) # mxT matrix

mu_time_EOF <- colMeans(Z_time_EOF, na.rm = TRUE)
mu_time_EOF # Tx1 vector

one_mat_time_EOF <- matrix(rep(1,m),nrow = m, ncol = 1) # mx1 vector of 1s
Z_tilde_time_EOF <- (Z_time_EOF - one_mat_time_EOF%*%t(mu_time_EOF))/(sqrt(m-1))
dim(Z_tilde_time_EOF)

# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(FactoMineR)

# Function to perform EOF analysis using SVD
compute_time_EOF <- function(data) {
  svd_result_time_EOF <- svd(data)
  Psi_time_EOF <- t(svd_result_time_EOF$u)
  eigen_values_time_EOF <- svd_result_time_EOF$d^2
  total_variance_time_EOF <- sum(eigen_values_time_EOF)
  expl_var_time_EOF <- eigen_values_time_EOF / total_variance_time_EOF
  return(list(Psi_time_EOF = Psi_time_EOF, expl_var_time_EOF = expl_var_time_EOF))
}

# Compute EOFs and number of EOFs
Psi_time_EOF <- compute_time_EOF(Z_tilde_time_EOF)$Psi_time_EOF
actual_expl_var_time_EOF <- compute_time_EOF(Z_tilde_time_EOF)$expl_var_time_EOF

# Set up permutation test
num_permutations_time_EOF <- 100
permuted_var_time_EOF <- matrix(nrow = num_permutations_time_EOF, ncol = length(actual_expl_var_time_EOF))

# Generate permuted datasets and compute their EOFs
set.seed(123)  
for (i in 1:num_permutations_time_EOF) {
  # Time-wise permutation
  permuted_data_time_EOF <- apply(Z_tilde_time_EOF, 1, sample)
  permuted_var_time_EOF[i, ] <- compute_time_EOF(permuted_data_time_EOF)$expl_var_time_EOF
}

# Mean variance explained by permuted data
mean_permuted_var_time_EOF <- apply(permuted_var_time_EOF, 2, mean)

# Plotting actual vs permuted explained variance
var_df_time_EOF <- data.frame(Index = seq_len(ncol(Z_tilde_time_EOF)),
                          Actual = actual_expl_var_time_EOF,
                          Permuted = mean_permuted_var_time_EOF)

ggplot(var_df_time_EOF, aes(x = Index)) +
  geom_point(aes(y = Actual, color = "Actual Data")) +
  geom_point(aes(y = Permuted, color = "Permuted Data")) +
  scale_color_manual(values = c("Actual Data" = "blue", "Permuted Data" = "red")) +
  labs(title = "Time EOF Analysis: Actual vs. Permuted Data",
       x = "Index of EOF", y = "Proportion of Variance Explained") +
  theme_minimal()

# Determine the number of significant EOFs
num_time_EOFs <- which(actual_expl_var_time_EOF < mean_permuted_var_time_EOF)[1]
num_time_EOFs

time_EOFs <- Psi_time_EOF[,1:num_time_EOFs]
dim(time_EOFs)

# Clustering after Time-EOF -----------------------------------------------

dim(Z_time_EOF)
dim(time_EOFs)

am <- Z_time_EOF %*% time_EOFs
dim(am)
am_df <- data.frame(am) 

# Convert row names to a column and rename it as "Location"
am_df <- cbind(am_df , locations = unique_locations$location)
# Check the dimensions
dim(am_df)  #  2520 x 26


library(NbClust)   # Optimal number of clusters
library(ClustGeo)  # Spatially constrained clustering

# Step 1: Standardize EOF Features 
# Convert to tibble to ensure dplyr functions work properly
am_df <- as_tibble(am_df)

# Extract EOF features
time_EOF_features <- am_df %>%
  dplyr::select(X1:X25) 

# Normalize EOF features
time_EOF_scaled <- scale(time_EOF_features)

# Step 2: Compute Feature-Based Dissimilarity Matrix 
D0 <- dist(time_EOF_scaled)  # Feature-based dissimilarity

# Step 3: Compute Geographic Dissimilarity Matrix 
# Extract Latitude & Longitude
am_df <- am_df %>%
  mutate(latitude = as.numeric(sub("_.*", "", locations)),  
         longitude = as.numeric(sub(".*_", "", locations)))

# Convert negative longitudes to [0,360] for proper plotting
am_df <- am_df %>%
  mutate(longitude = ifelse(longitude < 0, longitude + 360, longitude))
print(colnames(am_df))
# Compute geographic dissimilarity matrix
D1 <- dist(am_df %>% dplyr::select(latitude, longitude))


#  Step 4: Use `NbClust` 
nb_results <- NbClust(data = time_EOF_scaled, distance = "minkowski", 
                      min.nc = 2, max.nc = 10, method = "centroid", index = "all")

# Get the most frequently suggested number of clusters
optimal_k <- as.integer(names(sort(table(nb_results$Best.nc[1,]), decreasing = TRUE)[1]))
print(paste("Optimal number of clusters suggested:", optimal_k))

par(mfrow=c(1,1))
# Plot the frequency of the best cluster choices
barplot(table(nb_results$Best.nc[1,]), 
        main = "Optimal Number of Clusters ", 
        xlab = "Number of Clusters", ylab = "Frequency")

#  Step 5: Perform Spatially-Constrained Clustering 
alpha <- 0.6  # Adjust this value for spatial influence in `hclustgeo`
hclust_result <- hclustgeo(D0, D1, alpha = alpha)

# Assign clusters based on improved `NbClust`
clusters <- cutree(hclust_result, optimal_k )
am_df$Cluster <- as.factor(clusters)

# Step 6: Visualize the Clustered Locations 
colors <- colorRampPalette(c("red", "green", "blue", "yellow", "hotpink1", "purple", "cyan"))(optimal_k)

ggplot(am_df, aes(x = longitude, y = latitude, color = Cluster)) +
  geom_point(size = 3) +
  theme_minimal() +
  scale_color_manual(values = colors) +
  labs(title = paste("Clustered Locations (Alpha =", alpha, ", k =", optimal_k, ")"), 
       x = "Longitude", y = "Latitude") +
  theme(legend.position = "right")

#  Step 7: Print Summary of Clustering 
cat("Clustering completed with", optimal_k, "clusters using spatially improved NbClust\n")
print(table(am_df$Cluster)) 

# View final dataset
head(am_df)
dim(am_df)
am_df <- data.frame(am_df)
rownames(am_df) <- 1:nrow(am_df)

am_df <- am_df %>%
  mutate(latitude = as.numeric(sub("_.*", "", locations)),  
         longitude = as.numeric(sub(".*_", "", locations)))

head(am_df)
dim(am_df)
View(am_df)


# Spatial EOFs for each clusters ------------------------------------------


cluster_list <- split(am_df, am_df$Cluster)

# Assigning separate data frames for each cluster
cluster_1_df <- cluster_list[["1"]]
cluster_2_df <- cluster_list[["2"]]
cluster_3_df <- cluster_list[["3"]]
cluster_4_df <- cluster_list[["4"]]
cluster_5_df <- cluster_list[["5"]]
cluster_6_df <- cluster_list[["6"]]
cluster_7_df <- cluster_list[["7"]]
cluster_8_df <- cluster_list[["8"]]
cluster_9_df <- cluster_list[["9"]]
cluster_10_df <- cluster_list[["10"]]


dim(cluster_1_df)
dim(cluster_2_df)
dim(cluster_3_df)
dim(cluster_4_df)
dim(cluster_5_df)
dim(cluster_6_df)
dim(cluster_7_df)
dim(cluster_8_df)
dim(cluster_9_df)
dim(cluster_10_df)



### Spatial EOF for 5 clusters ###

# Create a list to store SST anomaly data for each cluster
cluster_ssta_list <- list()

# Loop through clusters 1 to optimal_k
for (i in 1:optimal_k) {
  
  # Extract cluster-specific latitude and longitude
  cluster_locations <- unique(get(paste0("cluster_", i, "_df"))[, c("latitude", "longitude")])
  
  # Merge with 'data' to retain only matching lat/lon
  cluster_ssta_list[[i]] <- merge(ssta_data, cluster_locations, by = c("latitude", "longitude"))
}


# Now, each cluster's SST data is stored in cluster_ssta_list[[1]], cluster_ssta_list[[2]], etc.

library(ggplot2)
library(gridExtra)

# Create empty lists to store results
EOF_list <- list()
num_eofs_list <- list()
Z_list <- list()
plot_list <- list()
Psi_list <- list()
actual_expl_var_list <- list()


for (i in 1:optimal_k) {
 
  # Retrieve SST anomaly data for each cluster
  cluster_data_i <- cluster_ssta_list[[i]]
  cluster_data_i <- cluster_data_i[,-c(1,2,3,4)]
  
  # Handle different data formats
  if ("Time" %in% names(cluster_data_i)) {
    # Long-format → Pivot wider
    cluster_data_wide_i <- cluster_data_i %>%
      pivot_wider(names_from = Time, values_from = Sea.Surface.Temperature.Anomaly)
    
    # Get numeric data from time columns
    cluster_final_data_i <- cluster_data_wide_i %>% select(where(is.numeric))
  } else {
    # Already wide-format → Directly select "X"-type columns
    cluster_final_data_i <- cluster_data_i %>% select(starts_with("X"))
  }
  dim(cluster_final_data_i)
  # Remove fully NA rows (e.g., spatial points with missing data)
  cluster_final_data_i <- cluster_final_data_i[rowSums(is.na(cluster_final_data_i)) < ncol(cluster_final_data_i), ]
  
  # Create T x m matrix
  Z_i <- t(as.matrix(cluster_final_data_i))  # T x m
  Z_list[[i]] <- Z_i
  dim(Z_i)
  # Compute mean vector
  mu_vec_i <- colMeans(Z_i, na.rm = TRUE)
  length(mu_vec_i)
  # Mean-centered anomaly matrix
  T_i <- nrow(Z_i)
  m_i <- ncol(Z_i)
  one_mat_i <- matrix(1, nrow = T_i, ncol = 1)
  Z_tilde_i <- (Z_i - one_mat_i %*% t(mu_vec_i)) / sqrt(T_i - 1)
  dim(Z_tilde_i)
  # EOF function
  compute_EOF <- function(data) {
    svd_result <- svd(data, nu = T_i , nv = m_i )
    Psi <- t(svd_result$v)
    eigen_values <- svd_result$d^2
    total_variance <- sum(eigen_values)
    expl_var <- eigen_values / total_variance
    return(list(Psi = Psi, expl_var = expl_var))
  }

  
  # Compute EOFs
  eof_result_i <- compute_EOF(Z_tilde_i)
  Psi_i <- eof_result_i$Psi
  expl_var_i <- eof_result_i$expl_var
  Psi_list[[i]] <- Psi_i
  actual_expl_var_list[[i]] <- expl_var_i
  dim(Psi_list[[i]])

  # Permutation test
  num_permutations <- 100
  permuted_var_i <- matrix(NA,nrow = num_permutations, ncol = length(expl_var_i))
  set.seed(123)
  
  for (j in 1:num_permutations) {
    # Permute over time (i.e., shuffle each row independently)
    permuted_data_i <- apply(Z_tilde_i, 2 , sample)  # still T x m
    permuted_result_i <- compute_EOF(permuted_data_i)
    permuted_var_i[j, ] <- permuted_result_i$expl_var
  }
  
  # Mean permuted variance
  mean_permuted_var_i <- apply(permuted_var_i, 2, mean)
  
  # Plot actual vs permuted explained variance
  var_df_i <- data.frame(Index = seq_len(length(expl_var_i)),
                         Actual = expl_var_i,
                         Permuted = mean_permuted_var_i)
  
  p <- ggplot(var_df_i, aes(x = Index)) +
    geom_point(aes(y = Actual, color = "Actual Data")) +
    geom_point(aes(y = Permuted, color = "Permuted Data")) +
    scale_color_manual(values = c("Actual Data" = "black", "Permuted Data" = "red")) +
    labs(title = paste("Time EOF Analysis: Cluster", i),
         x = "Index of EOF", y = "Proportion of Variance Explained") +
    theme_minimal()
  
  # Determine number of significant EOFs
  num_eofs_i <- which(expl_var_i < mean_permuted_var_i)[1]
  if (is.na(num_eofs_i)) num_eofs_i <- length(expl_var_i)
  num_eofs_list[[i]] <- num_eofs_i
  EOF_list[[i]] <- Psi_i[, 1:num_eofs_i, drop = FALSE]
  plot_list[[i]] <- p
}



# Arrange all plots in a grid
grid.arrange(grobs = plot_list, nrow = 5, ncol = 2)


dim(EOF_list[[1]])
dim(EOF_list[[2]])
dim(EOF_list[[3]])
dim(EOF_list[[4]])
dim(EOF_list[[5]])
dim(EOF_list[[6]])
dim(EOF_list[[7]])
dim(EOF_list[[8]])
dim(EOF_list[[9]])
dim(EOF_list[[10]])


# calculate at for 10 clusters and merge them with nino data and estimate M for each cluster --------------

# Load necessary libraries
library(dplyr)
library(tibble)
library(zoo)

# Create lists to store transformed data for all clusters
at_list <- list()
merged_data_list <- list()
C0_a_list <- list()
C1_a_list <- list()
M_hat_list <- list()

# Loop through all 10 clusters
for (i in 1:10) {

  # Compute the transformed matrix at_c_i (Projection of Z_list[[i]] onto EOFs)
  at_c_i <- Z_list[[i]] %*% EOF_list[[i]]

    # Convert to a data frame and add row names as "Time"
  at_c_i <- data.frame(at_c_i) %>% rownames_to_column(var = "Time")
  
  # Standardize the "Time" format
  at_c_i$Time <- gsub("\\.", " ", at_c_i$Time)  # Replace dots with spaces
  at_c_i$Time <- as.character(at_c_i$Time)  # Ensure it's a character
  
  # Ensure "Time" is in the same format in both data frames
  nino_df$Time <- as.character(nino_df$Time)
  
  # Merge with the corresponding section of nino_df
  merged_data_i <- left_join(at_c_i, nino_df[242:640,], by = "Time")
  
  # Convert "Time" to a sortable date format
  merged_data_i$Time <- as.yearmon(merged_data_i$Time, "%b %Y")
  
  # Sort the dataframe by "Time"
  merged_data_i <- merged_data_i[order(merged_data_i$Time), ]
  
  # Convert back to "Month Year" format if needed
  merged_data_i$Time <- format(merged_data_i$Time, "%b %Y")
  
  # Store merged data
  at_list[[i]] <- at_c_i[, -1]  # Exclude Time column before covariance calculations
  merged_data_list[[i]] <- merged_data_i
  
  # =========================== Covariance & M_hat Calculation ============================
  at_i <- as.matrix(at_list[[i]])  # Convert to matrix for covariance calculations
  
  # Step 1: Compute C0_a (covariance matrix at lag 0)
  mu_hat_i <- colMeans(at_i, na.rm = TRUE)  # Compute mean vector
  a_centered_i <- scale(at_i, center = mu_hat_i, scale = FALSE)  # Center the data

  C0_a_i <- var(a_centered_i, use = "pairwise.complete.obs")  # Covariance matrix
  C0_a_list[[i]] <- C0_a_i  # Store in list
  
  # Step 2: Compute C1_a (lag-1 covariance matrix)
  T_i <- nrow(a_centered_i)  # Get the number of time points
  C1_a_i <- cov(a_centered_i[-1, ], a_centered_i[-T_i, ])
  C1_a_list[[i]] <- C1_a_i  # Store in list

  # Step 3: Compute the method-of-moments estimator M(lag-1)
  M_hat_i <- C1_a_i %*% solve(C0_a_i)
  M_hat_list[[i]] <- M_hat_i  # Store in list
  
}


# Dimensions of ats
dim(at_list[[1]])
dim(at_list[[2]])
dim(at_list[[3]])
dim(at_list[[4]])
dim(at_list[[5]])
dim(at_list[[6]])
dim(at_list[[7]])
dim(at_list[[8]])
dim(at_list[[9]])
dim(at_list[[10]])

dim(M_hat_list[[1]])
dim(M_hat_list[[2]])
dim(M_hat_list[[3]])
dim(M_hat_list[[4]])
dim(M_hat_list[[5]])
dim(M_hat_list[[6]])
dim(M_hat_list[[7]])
dim(M_hat_list[[8]])
dim(M_hat_list[[9]])
dim(M_hat_list[[10]])


# Different metrices of non-normality for each cluster --------------------

## Frobenius norm of Ms
mu_2(M_hat_list[[1]])
mu_2(M_hat_list[[2]])
mu_2(M_hat_list[[3]])
mu_2(M_hat_list[[4]])
mu_2(M_hat_list[[5]])
mu_2(M_hat_list[[6]])
mu_2(M_hat_list[[7]])
mu_2(M_hat_list[[8]])
mu_2(M_hat_list[[9]])
mu_2(M_hat_list[[10]])

## Henrici norm of Ms
mu_3(M_hat_list[[1]])
mu_3(M_hat_list[[2]])
mu_3(M_hat_list[[3]])
mu_3(M_hat_list[[4]])
mu_3(M_hat_list[[5]])
mu_3(M_hat_list[[6]])
mu_3(M_hat_list[[7]])
mu_3(M_hat_list[[8]])
mu_3(M_hat_list[[9]])
mu_3(M_hat_list[[10]])


## Ruhe norm of Ms
mu_4(M_hat_list[[1]])
mu_4(M_hat_list[[2]])
mu_4(M_hat_list[[3]])
mu_4(M_hat_list[[4]])
mu_4(M_hat_list[[5]])
mu_4(M_hat_list[[6]])
mu_4(M_hat_list[[7]])
mu_4(M_hat_list[[8]])
mu_4(M_hat_list[[9]])
mu_4(M_hat_list[[10]])

frobenius_norms <- c()
henrici_norms <- c()
ruhe_norms <- c()

# Loop through all 10 clusters to compute the norms
for (i in 1:10) {
  frobenius_norms[i] <- mu_2(M_hat_list[[i]])  # Frobenius norm
  henrici_norms[i] <- mu_3(M_hat_list[[i]])    # Henrici norm
  ruhe_norms[i] <- mu_4(M_hat_list[[i]])       # Ruhe norm
}

# Combine the results into a data frame for plotting
norm_data <- data.frame(
  Cluster = factor(1:10),  # Cluster numbers as factor for proper x-axis
  Frobenius = frobenius_norms,
  Henrici = henrici_norms,
  Ruhe = ruhe_norms
)

# Load ggplot2 for visualization
library(ggplot2)

# Convert data to long format for ggplot
library(tidyr)
norm_data_long <- pivot_longer(norm_data, cols = c("Frobenius", "Henrici", "Ruhe"), 
                               names_to = "Norm_Type", values_to = "Value")

# Plot the three norms for each cluster
ggplot(norm_data_long, aes(x = Cluster, y = Value, color = Norm_Type, group = Norm_Type)) +
  geom_line(size = 1) +  # Line plot
  geom_point(size = 3) +  # Add points
  labs(title = "Comparison of Non-Normality Metrics Across Clusters",
       x = "Cluster",
       y = "Norm Value",
       color = "Norm Type") +
  theme_minimal() +
  theme(legend.position = "top")

# BHM (gibbs sampler) ------------------------------------------------------

library(MASS)
library(Matrix)
library(mvtnorm)
library(MCMCpack)  
library(progress)


gibbs_sampler_simulated <- function(n_iter) {
  set.seed(1000)
  
  # Data (EOF)
  Z <- t(at) 
  
  n <- nrow(Z)
  T <- ncol(Z)
  
  # Initialize hyperparameters 
  mu_0 <- rep(0, n)
  Sigma_0 <- diag(1, n) 
  
  nu_R <- 10 + n
  nu_q <- 10 + n
  
  C_R <- diag(0.1, n) 
  C_q <- diag(0.1, n)
  
  mu_m <- rep(0, n^2)
  Sigma_m <- diag(1, n^2)  
  
  H_list <- lapply(1:T, function(x) matrix(runif(n * n), n, n))
  
  # Initial values
  Y_1 <- matrix(rnorm(n * (T + 1), 0, 1), nrow = n, ncol = (T + 1))
  M_1 <- matrix(rnorm(n * n, 0, 4), nrow = n)
  R_1 <- diag(runif(n, 0.5, 1.5), n)
  Q_1 <- diag(runif(n, 0.5, 1.5), n)
  
  # Storage
  Y_0_updated_1 <- matrix(NA, n, n_iter)
  Y_t <- matrix(NA, n, T - 1)
  Y_t_updated_1 <- list()
  Y_T_updated_1 <- matrix(NA, n, n_iter)
  Y_updated_1 <- list()
  R_updated_1 <- list()
  Q_updated_1 <- list()
  M_updated_1 <- list()
  v_t <- vector("list", T - 1)
  a_t <- matrix(NA, n, T - 1)
  diff_R <- matrix(NA, n, T)
  diff_Q <- matrix(NA, n, T)
  
  time_records_1 <- numeric(n_iter)
  
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  start_time_total <- proc.time()
  
  for (i in 1:n_iter) {
    iter_start <- proc.time()
    
    ## ===== Update Y_0 =====
    v_0 <- solve(t(M_1) %*% solve(Q_1 + diag(1e-6, n), M_1) + solve(Sigma_0 + diag(1e-6, n)))
    a_0 <- t(M_1) %*% solve(Q_1, Y_1[, 2]) + solve(Sigma_0, mu_0)
    Y_0 <- mvrnorm(1, v_0 %*% a_0, v_0)
    Y_0_updated_1[, i] <- Y_0
    
    ## ===== Update Y_t (1 to T-1) =====
    for (t in 1:(T - 1)) {
      v_t[[t]] <- solve(t(H_list[[t]]) %*% solve(R_1, H_list[[t]]) +
                          solve(Q_1) + t(M_1) %*% solve(Q_1, M_1))
      a_t[, t] <- t(H_list[[t]]) %*% solve(R_1, Z[, t]) +
        solve(Q_1, M_1 %*% Y_1[, t]) +
        t(M_1) %*% solve(Q_1, Y_1[, t + 1])
      Y_t[, t] <- mvrnorm(1, v_t[[t]] %*% a_t[, t], v_t[[t]])
    }
    Y_t_updated_1[[i]] <- Y_t
    
    ## ===== Update Y_T =====
    H_T <- H_list[[T]]
    v_T <- solve(t(H_T) %*% solve(R_1, H_T) + solve(Q_1))
    a_T <- t(H_T) %*% solve(R_1, Z[, T]) + solve(Q_1, M_1 %*% Y_1[, T])
    Y_T <- mvrnorm(1, v_T %*% a_T, v_T)
    Y_T_updated_1[, i] <- Y_T
    
    # Full latent trajectory
    Y_1 <- cbind(Y_0, Y_t, Y_T)
    Y_updated_1[[i]] <- Y_1
    
    ## ===== Update R =====
    scale_R <- matrix(0, n, n)
    for (t in 1:T) {
      diff_R[, t] <- Z[, t] - H_list[[t]] %*% Y_1[, t + 1]
      scale_R <- scale_R + diff_R[, t] %*% t(diff_R[, t])
    }
    scale_R <- solve(scale_R + nu_R * C_R)
    R_1_inv <- rwish(nu_R + T, scale_R)
    R_1 <- solve(R_1_inv)
    R_updated_1[[i]] <- R_1
    
    ## ===== Update Q =====
    scale_Q <- matrix(0, n, n)
    for (t in 2:(T + 1)) {
      diff_Q[, t - 1] <- Y_1[, t] - M_1 %*% Y_1[, t - 1]
      scale_Q <- scale_Q + diff_Q[, t - 1] %*% t(diff_Q[, t - 1])
    }
    scale_Q <- solve(scale_Q + nu_q * C_q)
    Q_1_inv <- rwish(nu_q + T, scale_Q)
    Q_1 <- solve(Q_1_inv)
    Q_updated_1[[i]] <- Q_1
    
    ## ===== Update M =====
    I_T <- Diagonal(T, x = 1)
    Y_kron <- kronecker(t(Y_1[, 1:T]), Diagonal(n, x = 1))
    Q_inv <- solve(kronecker(I_T, Q_1))
    v_m <- solve(t(Y_kron) %*% Q_inv %*% Y_kron + solve(Sigma_m))
    a_m <- t(Y_kron) %*% Q_inv %*% as.vector(Y_1[, 2:(T + 1)]) + solve(Sigma_m, mu_m)
    m <- mvrnorm(1, v_m %*% a_m, v_m)
    M_1 <- matrix(m, nrow = n)
    M_updated_1[[i]] <- M_1
    
    # Save iteration time
    time_records_1[i] <- (proc.time() - iter_start)[["elapsed"]]
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  end_time_total <- proc.time()
  
  return(list(
    Y_0_updated_1 = Y_0_updated_1,
    Y_t_updated_1 = Y_t_updated_1,
    Y_T_updated_1 = Y_T_updated_1,
    Y_updated_1 = Y_updated_1,
    R_updated_1 = R_updated_1,
    Q_updated_1 = Q_updated_1,
    M_updated_1 = M_updated_1,
    execution_time_total_1 = end_time_total - start_time_total,
    time_records_1 = time_records_1
  ))
}


result_fixed <- gibbs_sampler_simulated(10000)


# Print total execution time
print(result_fixed$execution_time_total_1)

# Extract M_updated_1 from the result list
M_sim_est_1 <- result_fixed$M_updated_1  

# Convert the list of matrices into a 3D array
M_sim_est_array_1 <- simplify2array(M_sim_est_1)
M_sim_est_array_1

# Compute the element-wise mean across all 100 matrices
M_sim_est <- apply(M_sim_est_array_1 , c(1, 2), mean)
M_sim_est

# stochastic optimals and eofs -----------------------------------------------------

solve_discrete_lyapunov <- function(M) {
  I <- diag(nrow(M))
  n <- nrow(M)
  A <- kronecker(t(M), t(M)) - diag(n^2)
  b <- as.vector(-I)
  x <- solve(A, b)
  B <- matrix(x, nrow = n, byrow = TRUE)
  return(B)
}

solve_discrete_lyapunov_eof <- function(M) {
  n <- nrow(M)
  I_big <- diag(n^2)
  
  # Identity for Q_eta
  Q_eta <- diag(n)
  
  # Construct the matrix A = I - (M ⊗ M)
  A <- I_big - kronecker(M, M)
  
  # Vectorized Q_eta
  b <- as.vector(Q_eta)
  
  # Solve for vec(Sigma_Y)
  x <- solve(A, b)
  
  # Reshape vector into matrix form
  Sigma_Y <- matrix(x, nrow = n, byrow = FALSE)
  
  return(Sigma_Y)
}

## Checking

M_s <- matrix(c(0.4,0.01,4.0,0.1),nrow = 2,byrow = T)
M_s

B_result <- solve_discrete_lyapunov(M_s)
B_result

Sigma_Y_result <- solve_discrete_lyapunov_eof(M_s)
Sigma_Y_result

# Eigen decomposition
eigen(B_result)
eigen(Sigma_Y_result)

## With Freq M

M_freq_lag1

B_result_freq <- solve_discrete_lyapunov(M_freq_lag1)
B_result_freq

Sigma_Y_result_freq <- solve_discrete_lyapunov_eof(M_freq_lag1)
Sigma_Y_result_freq

# Eigen decomposition
eigen(B_result_freq)
eigen(Sigma_Y_result_freq)




# Balanced truncation -----------------------------------------------------

# Sigma_Y = P and B = Q by FARRELL & IOANNOU

P_freq <- Sigma_Y_result_freq
Q_freq <- B_result_freq

eig_P_freq <- eigen(P_freq)
P_sqrt_freq <- eig_P_freq$vectors %*% diag(sqrt(eig_P_freq$values)) %*% t(eig_P_freq$vectors)

R_freq <- P_sqrt_freq %*% Q_freq %*% P_sqrt_freq
R_freq

eig_R_freq <- eigen(R_freq)
U_freq <- eig_R_freq$vectors  # Unitary matrix U
Sigma_squared_freq <- diag(eig_R_freq$values)  # Eigenvalues = Hankel singular values squared
Sigma_squared_freq


all.equal(t(U_freq) %*% R_freq %*% U_freq, Sigma_squared_freq) #TRUE


P_inv_sqrt_freq <- eig_P_freq$vectors %*% diag(1 / sqrt(eig_P_freq$values)) %*% t(eig_P_freq$vectors)

Sigma_vals_freq <- sqrt(diag(Sigma_squared_freq))  
Sigma_freq <- diag(Sigma_vals_freq)
Sigma_half_freq <- diag(sqrt(Sigma_vals_freq))  # S^{1/2}

T_freq <- Sigma_half_freq %*% t(U_freq) %*% P_inv_sqrt_freq
T_freq

T_inv_freq <- solve(T_freq)
A_tilde_freq <- T_freq %*% M_freq_lag1 %*% T_inv_freq
A_tilde_freq

sigma_freq <- sort(Sigma_vals_freq, decreasing = TRUE) # Hankel singular values
sigma_freq

sigma_df_freq <- data.frame(SingularValue = sigma_freq, Index = 1:length(sigma_freq))

# Plot the scree plot
library(ggplot2)
plot <- ggplot(sigma_df_freq, aes(x = Index, y = SingularValue)) +
  geom_line() +  
  geom_point() +  
  scale_x_continuous(breaks = 1:length(sigma_freq)) +  
  labs(title = "Scree Plot of Hankel Singular Values (Freq)",
       x = "Index",
       y = "Singular Value") +
  theme_minimal()

# Detecting the elbow point using the 'elbow' method
elbow_point_freq <- which(diff(diff(sigma_freq)) == min(diff(diff(sigma_freq)))) + 1 
plot + geom_vline(xintercept = elbow_point_freq, linetype = "dashed", color = "blue") +
  annotate("text", x = elbow_point_freq, y = max(sigma_freq), label = paste("Elbow at index", elbow_point_freq), vjust = -1)



k_freq <- elbow_point_freq  # Choose the number of modes to retain

A_tilde_11_freq <- A_tilde_freq[1:k_freq, 1:k_freq] # k=5 x k=5 matrix
A_tilde_11_freq
latex_code_A_tilde_11_freq <- print(xtable(A_tilde_11_freq), include.rownames=FALSE, include.colnames=FALSE, floating=FALSE, hline.after=NULL, print.results = FALSE)



T1_inv_freq <- T_inv_freq[, 1:k_freq]  # n=10 x k=5 matrix
dim(T1_inv_freq)



# Function that evolves the reduced system 
reduced_step <- function(z_t) {
  return(A_tilde_11_freq %*% z_t)
}

# Function to recover full state from reduced coordinates
recover_full_state <- function(z_t) {
  return(T1_inv_freq %*% z_t)
}


zt_reduce <- reduced_step(at[1:k_freq,1:10])
dim(zt_reduce)

at_recover <- recover_full_state(zt_reduce)
dim(at_recover)
