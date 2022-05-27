
fc_initial <- read.csv("drivers/initial_4fc.csv")

glm_nml <- read_nml("GLMAED/glm3.nml")
glm_nml <- set_nml(glm_nml, arg_name = "start", arg_val = paste(initial_fc_date, "00:00:00"))
glm_nml <- set_nml(glm_nml, arg_name = "stop", arg_val = paste(initial_fc_date+14, "23:00:00"))

glm_nml <- set_nml(glm_nml, arg_name = "the_depths", arg_val = fc_initial$depth)
glm_nml <- set_nml(glm_nml, arg_name = "the_temps", arg_val = round(fc_initial$temp_c,3))

write_nml(glm_nml, "GLMAED/glm3.nml")

t1 <- Sys.time()
for(i in 1:31){
  if(i < 10){ensrun <- paste0("0",i)}else{ensrun <- paste0(i)}
  metfile <- read.csv(paste0("drivers/met_files/LakeEnsemblr_meteo_gefs_ens",i,".csv"))
  metfile2 <- dplyr::select(metfile, datetime, Air_Temperature_celsius, Relative_Humidity_percent,
                            Longwave_Radiation_Downwelling_wattPerMeterSquared, 
                            Shortwave_Radiation_Downwelling_wattPerMeterSquared,
                            Rainfall_millimeterPerHour,Ten_Meter_Elevation_Wind_Speed_meterPerSecond,
                            Snowfall_millimeterPerHour)
  colnames(metfile2) <- c("Date", "AirTemp", "RelHum", "LongWave", "ShortWave", "Rain", "WindSpeed", "Snow")
  metfile2 <- dplyr::select(metfile2, Date, ShortWave, LongWave, AirTemp, RelHum, WindSpeed, Rain, Snow)
  write.csv(metfile2, "GLMAED/meteo_file_2.csv", row.names = FALSE, quote = FALSE)
  
  
  glm_nml <- set_nml(glm_nml, arg_name = "out_fn", arg_val = paste0("output_", ensrun))
  write_nml(glm_nml, "GLMAED/glm3.nml")
  run_glm(verbose = T, sim_folder = "GLMAED")
}
t2 <- Sys.time()
t2 - t1



get_var_df <- function(var_name = "temp", out_fn){
  var_out <- get_var(file = out_fn, 
                     var_name = var_name,
                     reference = 'surface', 
                     z_out = 0:58)
  var_df_tidy <- var_out %>% 
    tidyr::gather(key = var, value = value, starts_with(var_name)) %>% 
    tidyr::separate(var, into = c("var", "var2" ,"depth"), sep = "_") %>% 
    mutate(depth = round(as.numeric(depth)), DateTime = ymd_hms(DateTime))
  return(var_df_tidy)
}




chllist <- list()
phslist <- list()
oxylist <- list()
nitlist <- list()
for(i in 1:31){
  if(i < 10){ensrun <- paste0("0",i)}else{ensrun <- paste0(i)}
  out_file <- paste0("GLMAED/output/output_", ensrun, ".nc")
  
  tryCatch({
    chllist[[i]] <- get_var_df("PHY_TCHLA", out_file) %>% 
      mutate(ens = ensrun)
    
    phslist[[i]] <- get_var_df("PHS_frp", out_file) %>% 
      mutate(ens = ensrun)
    
    oxylist[[i]] <- get_var_df("OXY_oxy", out_file) %>% 
      mutate(ens = ensrun)
    
    nitlist[[i]] <- get_var_df("NIT_nit", out_file) %>% 
      mutate(ens = ensrun)
  }, error = function(e){cat("error");chllist[[i]] <- NA;
  phslist[[i]] <- NA;phslist[[i]] <- NA;nitlist[[i]] <- NA})
  
  
  
  print(i)
}


if(!dir.exists(paste0(here::here("GLMAED"), "/archive_fc"))){
  dir.create(paste0(here::here("GLMAED"), "/archive_fc"))
}


data.table::rbindlist(chllist) %>% 
  write.csv("GLMAED/archive_fc/chla.csv", row.names = FALSE)
# group_by(depth, DateTime) %>% summarize(value = median(value)) %>% 
# ggplot(aes(x = DateTime, y = depth, fill = value)) + geom_raster(interpolate = TRUE) + 
# scale_fill_viridis_c() + scale_y_reverse() +
# labs(title = "chla")

data.table::rbindlist(phslist) %>% 
  write.csv("GLMAED/archive_fc/phs.csv", row.names = FALSE)
# group_by(depth, DateTime) %>% summarize(value = median(value)) %>% 
# ggplot(aes(x = DateTime, y = depth, fill = value)) + geom_raster(interpolate = TRUE) + 
# scale_fill_viridis_c() + scale_y_reverse() +
# labs(title = "phosphorus")

data.table::rbindlist(oxylist) %>% 
  write.csv("GLMAED/archive_fc/oxy.csv", row.names = FALSE)
# group_by(depth, DateTime) %>% summarize(value = median(value)) %>% 
# ggplot(aes(x = DateTime, y = depth, fill = value)) + geom_raster(interpolate = TRUE) + 
# scale_fill_viridis_c() + scale_y_reverse() +
# labs(title = "oxygen")

data.table::rbindlist(nitlist) %>% 
  write.csv("GLMAED/archive_fc/nit.csv", row.names = FALSE)
# group_by(depth, DateTime) %>% summarize(value = median(value)) %>% 
# ggplot(aes(x = DateTime, y = depth, fill = value)) + geom_raster(interpolate = TRUE) + 
# scale_fill_viridis_c() + scale_y_reverse() +
# labs(title = "nitrogen (nitrate?)")

# p <- list()
# for(i in 1:31){
#   if(i < 10){ensrun <- paste0("0",i)}else{ensrun <- paste0(i)}
#   out_file <- paste0("GLMAED/output/output_", ensrun, ".nc")
#   # tryCatch(plot_var(nc_file = out_file, var_name = 'temp', reference = 'surface'),
#   #          finally = next)
#   p[[i]] <- plot_var(nc_file = out_file, var_name = 'temp', reference = 'surface', legend.title = "°C")
#   
#   print(i)
# }
# 
# p[[1]] + p[[2]] + p[[3]] + p[[4]] + p[[5]] + p[[6]] + 
#   p[[7]] + p[[8]] + p[[9]] + p[[10]] + p[[11]] + p[[12]] + 
#   p[[11]] + p[[12]] + p[[13]] + p[[14]] + plot_layout(ncol = 4, nrow = 4)
# 
# p[[15]] + p[[16]] + p[[17]] + p[[18]] + p[[19]] + p[[20]] +
#  p[[21]] + p[[22]] + p[[23]] + p[[24]] + p[[25]] + p[[26]] + 
#   p[[27]] + p[[28]] + p[[29]] + p[[30]] + p[[31]] + plot_layout(ncol = 4, nrow = 5)

# aed_nml <- read_nml('aed2/aed2.nml')
# 
# # heat maps of phosphate, nitrate and silica
# plot_var(nc_file = out_file, var_name = 'PHS_frp',reference = 'bottom')
# plot_var(nc_file = out_file, var_name = 'NIT_nit',reference = 'bottom')
# plot_var(nc_file = out_file, var_name = 'SIL_rsi',reference = 'bottom')
# 
# 
# aed_nml <- read_nml('GLMAED/aed2/aed2.nml')
# phyto_list <- get_nml_value(aed_nml,arg_name = 'aed2_phytoplankton::dbase')
# 
# path_phyto <- paste0("GLMAED/", phyto_list)
# phyto_nml <- read_nml(path_phyto)
# phyto_nam <- get_nml_value(phyto_nml,arg_name = 'pd%p_name')
# names <- unlist(strsplit(phyto_nam, ","))
# 
# lim_attributes <- c('fI', 'fNit', 'fPho', 'fSil', 'fT', 'fSal')
# plist <- list()
# pindex <- 1
# for (ii in seq_len(length(names))){
#   for (jj in seq_len(length(lim_attributes))){
#     
#     p1 <- plot_var(nc_file = out_file, var_name = paste0('PHY_',names[ii],'_',
#                                                          lim_attributes[jj]),
#                    legend.title = paste(names[ii], lim_attributes[jj]))
#     
#     plist[[pindex]] <- p1
#     pindex <- pindex + 1
#   }
# }
# 
# # limitation functions for cyanobacteria and diatoms
# p_cyano <- plist[[1]] / plist[[2]] / plist[[3]] / plist[[4]] / plist[[5]] / plist[[6]] 
# p_diatom <- plist[[7]] / plist[[8]] / plist[[9]] / plist[[10]] / plist[[11]] / plist[[12]] 
# 
# p_cyano
# p_diatom
# 
