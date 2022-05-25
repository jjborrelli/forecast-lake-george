
# CHOOSE FORECAST DATE (default is for TODAY)
initial_fc_date <- Sys.Date()

# Check DODS database for info
dodsmod <- rNOMADS::GetDODSModelRuns(paste0("https://nomads.ncep.noaa.gov:443/dods/gefs/gefs",
                                            lubridate::year(initial_fc_date), 
                                            stringr::str_pad(lubridate::month(initial_fc_date), 2, "left", "0"),
                                            stringr::str_pad(lubridate::day(initial_fc_date), 2, "left", "0")))

# If there is nothing there then stop the forecast
if(length(dodsmod$model.run.info) == 0){stop("Today has no forecast available")}



# DOWNLOAD GEFS DATA
## using function from Rnoaa4cast
noaa_gefs_point_download_downscale(
  read_from_path = FALSE,
  lat_list = 43.56711729817296,
  lon_list = -73.64878182822284,
  site_list = "DFWI",
  forecast_time = "00",
  forecast_date = initial_fc_date,
  downscale = TRUE,
  overwrite = TRUE,
  model_name = "gefs",
  model_name_ds = "gefs_ds",
  output_directory = here::here("drivers")
)


# LOAD DATA BACK INTO R
folder <- paste0(here::here("drivers"), "/DFWI/", initial_fc_date, "/00/")
lf <- list.files(folder)


if(length(lf) == 0){stop("The forecast failed to download properly")}

df <- list()
for(i in 1:length(lf)){
  ncfile <- ncdf4::nc_open(paste0(folder, lf[i]))
  times <- hours(ncvar_get(ncfile, varid = "time")) + ymd_h(paste(initial_fc_date, "00"))
  airt <- ncvar_get(ncfile, varid = "air_temperature")
  airp <- ncvar_get(ncfile, varid = "air_pressure")
  relhum <- ncvar_get(ncfile, varid = "relative_humidity")
  long <- ncvar_get(ncfile, varid = "surface_downwelling_longwave_flux_in_air")
  short <- ncvar_get(ncfile, varid = "surface_downwelling_shortwave_flux_in_air")
  precip <- ncvar_get(ncfile, varid = "precipitation_flux")
  spechum <- ncvar_get(ncfile, varid = "specific_humidity")
  wnd <- ncvar_get(ncfile, varid = "wind_speed")
  
  
  df[[i]] <- data.frame(times, airt, airp, relhum, long, short, precip, spechum, wnd,ens = i)
}

# COMPILE DATA FOR LakeEnsemblR
ens_meteo <- do.call(rbind, df) %>%
  mutate(airt = airt - 273.15) %>%
  mutate(precip = (precip * 3600)) %>%
  rename(datetime = times, Surface_Level_Barometric_Pressure_pascal = airp,
         Relative_Humidity_percent = relhum,
         Shortwave_Radiation_Downwelling_wattPerMeterSquared = short,
         Longwave_Radiation_Downwelling_wattPerMeterSquared = long,
         Air_Temperature_celsius = airt, Ten_Meter_Elevation_Wind_Speed_meterPerSecond = wnd,
         Rainfall_millimeterPerHour = precip) %>%
  mutate(Snowfall_millimeterPerHour = 0) %>%
  filter(ymd_hms(datetime) <= (ymd_hms(paste(initial_fc_date, "00:00:00"))+days(15)))

# ARCHIVE MET PREDICTION DATA
metpredpath <- paste0("C:/Users/borre/Documents/lg-forecast/drivers/met_archive/meteo_fc_", initial_fc_date, ".csv")
data.table::fwrite(ens_meteo, metpredpath)

# SPLIT ENSEMBLES
ens_sp <- split(ens_meteo, ens_meteo$ens)

# SAVE EACH ENSEMBLE FOR LakeEnsemblR
for(i in 1:length(ens_sp)){
  write.csv(dplyr::select(ens_sp[[i]], -ens, -spechum),
            paste0("C:/Users/borre/Documents/lg-forecast/drivers/met_files/LakeEnsemblr_meteo_gefs_ens", i, ".csv"),
            row.names = FALSE)
}

