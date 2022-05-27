# Required Libraries
## simulations
library(Rnoaa4cast)
library(LakeEnsemblR)
library(gotmtools)
library(rLakeAnalyzer)
library(glmtools)
library(GLM3r)

library(deSolve)
library(mgcv)

## data wrangling
library(readxl)
library(ncdf4)
library(lubridate)
library(dplyr)
library(tidyr)

## plotting
library(ggplot2)
library(patchwork)

## reporting
library(emayili)
library(rmarkdown)



brk <- "##################################################################"

cat(paste0(brk, "\n", brk, "\n", brk, "\n\n", "STARTING FORECAST RUN FOR ",
           Sys.Date(), "\n\n", brk, "\n", brk, "\n", brk))



## Get mean temperature from 2015-2019 for today for initial values
## Not included in this repo yet 
## SKIP TO LINE 87 to use default initial conditions

ncdf <- '../AEDmodels/output/ensemble_output.nc'
wtemp <- LakeEnsemblR::load_var(ncdf = ncdf, var = 'temp', return = 'list')

dtemp <- lapply(wtemp, function(x){
  xfilt <- dplyr::filter(x, month(datetime) == month(Sys.Date()), 
                         day(datetime) == day(Sys.Date()),
                         hour(datetime) == 0)
  xfilt <- tidyr::gather(xfilt, depth, temp_c, wtr_0:wtr_60)
})


## Pull most recent observations from JP vertical profiler
## Only works when connected to IBM VPN
strtdt <- Sys.Date() - 1
enddt <- Sys.Date() 
jpvp_baseurl <- "http://jp-was.watson.ibm.com/AssemblyPoint/Data/DeepThunderAssimilationReport.csv?sid=43"
strtend <- paste0("start=", stringr::str_replace_all(strtdt, "-", ""), "000000", "&end=", stringr::str_replace_all(enddt, "-", ""), "000000")
filename <- paste0("fileName=JP_VP_006A_ALL-Variables_", strtdt, "_", enddt)

obsvp <- read.csv(paste(jpvp_baseurl, strtend, filename, sep = "&"), skip = 10, na.strings = "-")

wt <- obsvp %>% select(MM.dd.yyyy, HH.mm.ss, PC.RFU.Water.EXO..RFU., Chla.RFU.Water.EXO..RFU.,
                       Temperature.Water.EXO..C., DO.Water.EXO..mg.L., RS232SensorDepth.Water.EXO..m.) 
colnames(wt) <- c("date", "time", "pc_rfu", "chla_rfu", "temp_c", "do_mgL", "depth_m")
wto <- filter(wt, mdy(date) == Sys.Date()-1, grepl("23:", time)) %>% select(date, time, depth_m, temp_c) %>% filter(!is.na(temp_c)) %>% 
  mutate(depth = plyr::round_any(depth_m, 0.5)) %>% group_by(depth) %>% summarize(temp_o = mean(temp_c))
wto <- rbind(data.frame(depth = 0, temp_o = wto$temp_o[1]), wto)


if(any(is.na(wto$temp_o))){
  do.call(rbind, dtemp) %>% filter(!is.na(temp_c)) %>% 
    tidyr::separate(depth, into = c("wtr", "depth"), sep = "_") %>% 
    mutate(depth = as.numeric(depth)) %>% 
    group_by(depth) %>% summarize(temp_c = mean(temp_c)) %>% 
    write.csv(paste0(here::here(), "/drivers/initial_4fc.csv"), row.names = FALSE)
}else{
  do.call(rbind, dtemp) %>% filter(!is.na(temp_c)) %>% 
    tidyr::separate(depth, into = c("wtr", "depth"), sep = "_") %>% 
    mutate(depth = as.numeric(depth)) %>% 
    group_by(depth) %>% summarize(temp_c = mean(temp_c)) %>%  
    left_join(wto, by = "depth") %>% mutate(temp_c = case_when(is.na(temp_o) ~ temp_c, TRUE ~ (temp_c + temp_o)/2)) %>% 
    select(depth, temp_c) %>% 
    write.csv(paste0(here::here(), "/drivers/initial_4fc.csv"), row.names = FALSE)
  
}


## RUNNING SIMULATIONS ####
source(paste0(here::here(), "/lake_fc.R"))
source(paste0(here::here(), "/GLMAED/run_glm-aed.R"))
source(paste0(here::here(), "/ATN/run_ATN.R"))
###########################

## workaround to get the markdown to render when automating 
Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/bin/quarto/bin")

## RENDER REPORT
rmarkdown::render(paste0(here::here(), "/forecast_report.Rmd"))

## SEND report via email
bodytext <- paste("The Lake George Forecast report for", Sys.Date(), "is attached")
email <- envelope() %>%
  from(Sys.getenv("mygmailID")) %>%
  to("borrejj@gmail.com") %>%
  subject("Weekly Lake George Forecast") %>%
  text(bodytext) %>%
  attachment(path = paste0(here::here(), "/forecast_report.pdf"))

smtp <- server(
  host = "smtp.gmail.com",
  port = 587,
  username = Sys.getenv("mygmailID"),
  password = Sys.getenv("mygmailAppPW")
)

smtp(email, verbose = FALSE)


cat(paste0(brk, "\n", brk, "\n", brk, "\n\n", "FORECAST COMPLETE",
           "\n\n", brk, "\n", brk, "\n", brk))
