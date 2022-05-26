
atn <- function(times, state, parms){
  with(as.list(c(state, parms)), {
    
    N <- state[1]
    Xi <- state[-1]
    
    #### FORCING FUNCTIONS #################################################
    
    # Nutrients
    Nflux <- ninfun(times) #- n.out
    Nlim <- N/(N + kN)
    # Carrying capacity
    # approx_Ks <- approxfun(Ks, rule = 2)
    # Ks <- approx_Ks(times)
    # Ks.t <- Ks
    #Klim <- Gi(state, cij, Ks.t, prod) 
    Klim <- 1
    # Light
    #approx_light <- approxfun(I, rule = 2)
    I.t <- approx_light(times)
    # Temperature
    #approx_temp <- approxfun(T.i, rule = 2)
    T.i <- approx_temp(times)
    
    #temperature
    c.t <- Q10^((T.i - T0)/10)
    x.t <- x.i * c.t
    yij.t <- yij #* temp_lim(T.i, 20, 40, 400)
    cp.t <- temp_lim(T.i, Topt = T0.phyt, Tmax = maxt, Tsteepness = tsteep)
    
    h <- approx_therm(times)
    # light lim
    prod2 <- prod
    prod2[1] <- 1
    sumphyt <- sum(Xi * prod2)
    I.z <- I.t * exp(-(a.phy * (sumphyt/1e9) + a.bg*h))
    c.l <- I.z/(I0 + I.z)
    r.l <- r * apply(cbind(Nlim, c.l), 1, min) * c(cp.t)
    #r.l <- r * apply(cbind(Nlim, c.l, cp.t), 1, min)
    
    #### DYNAMICS  #################################################
    
    w2 <- apply(w, 2, function(x){if(sum(x) > 0){x/sum(x)}else{rep(0, length(x))}})
    
    Nuptake <- sum(r.l * Xi * Klim * (1 - s) * n.frac)
    dN <- Nflux + Xi[1] * mineralization - Nuptake - N * n.out
    #dynamics
    dB <- r.l * Xi * Klim * (1 - s) +                                   
      # producer growth
      fa * x.t * Xi * colSums((yij.t * Fij(B0, q, d, p, Xi, w2))) -                           
      # gain from eating
      rowSums(t(apply(yij.t * Fij(B0, q, d, p, Xi, w2), 1, 
                      function(row){row * x.t * Xi}))/eij, na.rm = T) - 
      # loss from being eaten
      (fm * x.t * Xi) -                                                                        
      # maintenance loss
      m * Xi +                                                                                   
      # mortality (fish only)   
      x.in
    
    # increases in DOC
    DOCdelay <- 1#(1-exp(-.45 * ceiling(times)))^3500
    dB[1] <- (dB[1] + (sum(t(apply(yij.t * Fij(B0, q, d, p, Xi, w2), 1,
                                   function(row){row * x.t * Xi}))/eij  * (1-eij), na.rm = TRUE) +
                         # egestion
                         max(sum(r.l * Xi * Klim * s), 0))) * (1 - DOCloss) - (Xi[1] * mineralization)
    
    
    #chla <- (1/(6.4+45*Gi(state, cij, Ks.t, prod)*(1-c(0,c.l))))*prod
    
    return(list(c(dN, dB)))#, N.lim = Nlim, light = c.l[2], tlim = cp.t[2]))
  })
}

temp_lim <- function(temp, Topt = 25, Tmax = 40, Tsteepness = 150){
  temp_opt <- exp(-(((Topt - temp)^2)/Tsteepness))
  penalty <- 1 - (1/(1 + exp(-temp + Tmax)))
  
  return(temp_opt * penalty)
}

Fij <- function(B0, q, d, p, B, w){
  numer <-  w * B^q
  halfsat <- B0^q
  interfere <- 0#sapply(1:length(B), function(k){sum(d[k,] * p[,k] * B[k] * B0[k,]^q)}) # d * p * B * B0^q
  res <- colSums(w * B^q)
  
  fres <- numer/t(apply(halfsat, 1, function(x) x + interfere + res))
  fres[is.nan(fres)] <- 0
  return(fres)
}

goExtinct <- function(times, states, parms){
  with(as.list(states), {
    states[states < 10^-10] <- 0
    return(c(states))
  })
}



compute_errorstats <- function(dyn, phy, zoo, taxa = FALSE){
  mb <- filter(zoo) %>% group_by(date = ymd(date), Var2) %>% summarize(bio_l = mean(bio_l))
  mb2 <- filter(phy) %>% group_by(date = mdy(sample_date), Var2) %>% summarize(bio_l = mean(biovol))
  if(taxa){
    mse <- dyn %>% 
      mutate(date = Var1 + ymd("2015-05-01")) %>% 
      filter(Var2 %in% c("phyto", "zoop", "blue", "brown", "cal", "clad", "cyc", "flag", "green", "red", "rot"),
             date %in% ymd(rbind(mb, mb2)$date)) %>% 
      left_join(rbind(mb, mb2), by = c("Var2","date")) %>% 
      filter(!is.na(bio_l)) %>% 
      mutate(diff = value - bio_l) %>% 
      group_by(run, Var2) %>% 
      #group_by(run) %>% 
      summarize(rmse = sqrt(mean(diff^2, na.rm = TRUE)), 
                re = sum(abs(diff), na.rm = TRUE)/sum(bio_l))
    return(mse)
  }else{
    mse <- dyn %>% 
      mutate(date = Var1 + ymd("2015-05-01")) %>% 
      filter(Var2 %in% c("phyto", "zoop", "blue", "brown", "cal", "clad", "cyc", "flag", "green", "red", "rot"),
             #Var2 %in% c("clad"),
             date %in% ymd(rbind(mb, mb2)$date)) %>% 
      left_join(rbind(mb, mb2), by = c("Var2","date")) %>% 
      filter(!is.na(bio_l)) %>% 
      mutate(diff = value - bio_l) %>% 
      #group_by(run, Var2) %>% 
      group_by(run) %>% 
      summarize(rmse = sqrt(mean(diff^2, na.rm = TRUE)), 
                re = sum(abs(diff), na.rm = TRUE)/sum(bio_l))
    return(mse)
  }
}

## RUN #### 

spp <- c("detritus", "plant", "blue", "brown", "green", "red", "flag", 
         "prot", "rot", "clad", "cal", "cyc", "macroinv", "fish")

params <- readRDS("ATN/params_modelC.rds")

ltpath <- paste0(here::here(), "/temp_forecasts/fc_", initial_fc_date, ".csv")
laketemps <- data.table::fread(ltpath)

wlist <- readRDS("ATN/best_w_vals.rds")
statelist <- readRDS("ATN/best_state_vals.rds")
knlist <- readRDS("ATN/best_kN_vals.rds")

ens <- 1:31
t1 <- Sys.time()
atnOUT <- list()
for(j in 1:30){
  
  params$w <- apply(wlist[[j]], 2, function(x){if(sum(x) > 0){x/sum(x)}else{rep(0, length(x))}})
  params$kN[params$kN > 1] <- unlist(knlist[j,-c(1,8)])
  state2 <- unlist(statelist[j,-c(1,17)])
  names(state2) <- c("N", spp)
  
  
  
  outlist <- list()
  for(i in seq_along(ens)){
    meteo <- read.csv(paste0("drivers/met_files/LakeEnsemblr_meteo_gefs_ens", i,".csv"))
    pardf <- meteo %>% group_by(day = floor_date(ymd_hms(datetime), "1 day")) %>% 
      summarize(light = mean(Shortwave_Radiation_Downwelling_wattPerMeterSquared * 0.45)) %>% 
      mutate(day = 1:n())
    
    ltemp <- laketemps %>% filter(Var3 > -10.5, ens == i) %>% 
      group_by(day = floor_date(Var2, "1day")) %>% 
      summarize(value = mean(value)) %>% mutate(day = 1:n())
    
    params$approx_light <- approxfun(pardf) 
    params$approx_temp <- approxfun(ltemp)
    
    out <- ode(state2, times = 1:14, func = atn, parms = params,
               events = list(func = goExtinct, time = 1:14))
    
    outlist[[i]] <- reshape2::melt(out) %>% 
      filter(Var2 %in% c("blue", "green", "red", "brown", "rot", "clad", "cal", "cyc")) %>% 
      mutate(ens = i)
    
    #print(i)
  }
  
  atnOUT[[j]] <- data.table::rbindlist(outlist) %>%
    mutate(iter = j)
  print(j)
}
t2 <- Sys.time()
t2-t1

write.csv(data.table::rbindlist(atnOUT), "ATN_out/output/atn_fc.csv")

