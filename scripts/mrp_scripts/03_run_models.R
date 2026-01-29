#### MRP + Spatial Priors using CTIS Data ####

#Build a series of models representing vaccination probability (min. 1 dose)
#For June 12, 2021 using the COVID-19 Trends and Impacts Survey (CTIS)

#Post-stratification of age/sex/education at county level, according to what is in the model

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,here,SUMMER,data.table,beepr)

#if(!require("INLA")) install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) 
library(INLA)
options(scipen = 999)

here::i_am("scripts/mrp_scripts/03_run_models.R")

#### Functions ####

load_data = function(poststrat_vars){

  #load geo, dat, spatial graph, and post for analysis
  
  #post
  load(file = here("data/cleaned_input_data/clean_postrat_age_sex_county_with_edu.rda"))
  
  #dat
  load(file = here("data/cleaned_input_data/clean_data_county_with_edu.rda"))

 temp = dat %>% filter(sex == 1) 

  dat <<- dat
  post <<- post
  
  if (str_detect(poststrat_vars,"xfips")) {
    
    #geographic data (geo)
    load(file=here("data/cleaned_input_data/ca_geo_sf_county.rda"))
    
    #spatial graph
    load(file = here("data/cleaned_input_data/spatialgraph.rda"))
    geo <<- geo
    g.graph <<- g.graph
    
  }

  state_choice <<- "CA" #for pulling data from the US Census API with tidycensus

}

prepare_data = function(dat,post,poststrat_vars){
  
  #if there are poststratification categories missing, it means that
  #there isn't data for these in the Census (it's a zero/censored category)
  
  require(tidyverse,sf)
  
  if (poststrat_vars == "xfips,age,edu") {
  
    print("poststrat_vars == xfips,age,edu")
    dat = dat %>% 
      select(xfips,vax,sex,age,edu) 

    post = post %>% 
      filter(!is.na(edu)) %>% #drop total counts from census by age/sex, no edu
      group_by(xfips,sex,age,edu) %>%
      reframe(
        vax = NA,
        estimate = sum(estimate)) %>%
      distinct() %>%
      ungroup() %>%
      select(vax,everything())
  
    
    #add the empty prediction categories to the bottom of the real data
    pred_dat = post %>% select(colnames(dat))
    dat = rbind(dat, pred_dat)

    index = geo %>% select(xfips, id) %>% sf::st_drop_geometry()
    
    dat = left_join(dat, index, by =c("xfips"))
    
    dat <<- dat
    post <<- post
    pred_dat <<- pred_dat
    
  }
  
  if (poststrat_vars == "xfips,age") {
    
    print("poststrat_vars == xfips,age")
    
    dat = dat %>% 
      select(xfips,vax,sex,age) %>%
      distinct()
  
    post = post %>% 
      filter(is.na(edu)) %>%
      select(-edu) %>%
      group_by(xfips,sex,age) %>%
      reframe(
        vax = NA,
        estimate = sum(estimate)) %>%
      distinct() %>%
      ungroup() %>%
      select(xfips,sex,vax,age,estimate)
    
    pred_dat = post %>% select(colnames(dat))
    
    #add the empty prediction categories to the bottom of the real data
    dat = rbind(dat, pred_dat)

    index = geo %>% select(xfips, id) %>% sf::st_drop_geometry()
    
    dat = left_join(dat, index, by =c("xfips"))
    
    dat <<- dat
    post <<- post
    pred_dat <<- pred_dat
  }
  
  if (poststrat_vars == "xfips,edu") {
    
    print("poststrat_vars == xfips,edu")
   
    dat = dat %>% 
      select(xfips,vax,sex,edu) %>%
      distinct()
    
    post = post %>% 
      filter(!is.na(edu)) %>%
      select(-age) %>%
      group_by(xfips,sex,edu) %>%
      reframe(
        vax = NA,
        estimate = sum(estimate)) %>%
      distinct() %>%
      ungroup() %>%
      select(xfips,sex,vax,edu,estimate)
    
    pred_dat = post %>% select(colnames(dat))
    
    #add the empty prediction categories to the bottom of the real data
    dat = rbind(dat, pred_dat)
    
    index = geo %>% select(xfips, id) %>% sf::st_drop_geometry()
    
    dat = left_join(dat, index, by =c("xfips"))
    
    dat <<- dat
    post <<- post
    pred_dat <<- pred_dat
  }
  
  if (poststrat_vars == "age,edu") {
    
    print("poststrat_vars == age,edu")
    glimpse(dat)
    dat = dat %>% 
      select(vax,sex,age,edu) 
        
    post = post %>% 
      filter(!is.na(edu)) %>%
      select(-xfips) %>%
      group_by(sex,age,edu) %>%
      reframe(
        vax = NA,
        estimate = sum(estimate)) %>%
      distinct() %>%
      ungroup() %>%
      select(vax,everything())
   
    pred_dat = post %>% select(colnames(dat))
   
    #add the empty prediction categories to the bottom of the real data
    dat = rbind(dat, pred_dat)

    dat <<- dat
    post <<- post
    pred_dat <<- pred_dat
    
  }
  
  if (poststrat_vars == "edu") {
    
    print("poststrat_vars == edu")
    
    dat = dat %>% 
      select(vax,sex,edu) %>%
      distinct()
    
    post = post %>% 
      filter(!is.na(edu)) %>%
      select(-xfips, -age) %>%
      group_by(sex,edu) %>%
      reframe(
        vax = NA,
        estimate = sum(estimate)) %>%
      distinct() %>%
      ungroup() %>%
      select(vax,everything())
    
    pred_dat = post %>% select(colnames(dat))
    
    #add the empty prediction categories to the bottom of the real data
    dat = rbind(dat, pred_dat)
    
    dat <<- dat
    post <<- post
    pred_dat <<- pred_dat
  }
  
  if (poststrat_vars == "age") {
    
    print("poststrat_vars == age")
    
    dat = dat %>% 
      select(vax,sex,age) %>%
      distinct()
    
    post = post %>% 
      filter(is.na(edu)) %>%
      select(-xfips, -edu) %>%
      group_by(sex,age) %>%
      reframe(
        vax = NA,
        estimate = sum(estimate)) %>%
      distinct() %>%
      ungroup() %>%
      select(vax,everything())
    
    pred_dat = post %>% select(colnames(dat))
    
    #add the empty prediction categories to the bottom of the real data
    dat = rbind(dat, pred_dat)

    dat <<- dat
    post <<- post
    pred_dat <<- pred_dat
  }

}

run_save_model = function(mf,dat,sex_spec,specs) {

  require(tidyverse,INLA,data.table)
  
  #make copies of age and education for rw1 interaction on age and education

  if(str_detect(specs,"interaction")){
  dat$age_copy = dat$age
  dat$edu_copy = dat$edu
  }
  
  #number of covariates + number of outcome variables (1) + 1
  sim_start = length(mf)+2 
  
  #NB: For a spatial model to work, the formula MUST refer to the "id"
  #built into the spatial graph, NOT a code like xfips. This is because
  #it only works with integers -- and xfips aren't (can't be) integers. 
  #Later, we can join on xfips, but this will not work without referring 
  #to the graph id!
  
  model <<- inla(mf, 
               data = dat, 
               family="betabinomial", 
               Ntrials = 1,
               control.family=list(link='logit'),
               control.compute = list(config=TRUE,
                                      return.marginals.predictor=TRUE,
                                      dic = TRUE,
                                      cpo = TRUE),
               control.predictor = list(compute=TRUE), #need to expit it later
               verbose = FALSE)

  save(model, file=here(paste0('data/model_objects/',specs,"_sex_",sex_spec,"_model_object.RDS")))
  }

posterior_draw = function(model,post,poststrat_vars,pred_dat,dat,sex_spec) {

  #take a thousand draws of the posterior of the INLA model,
  #extract the joint posterior distribution draws for each
  #stratum in the poststratification matrix given as an row-bound
  #object affixed to the survey data in the model (used to predict Y
  #for all strata in all areas), determine the posterior median probabilities,
  #and credibility intervals, and calculate the small area estimates for each
  #stratum in each area based on these summary estimates of the 
  #posterior probabilities.
  
  #model == fitted INLA model
  #sim_start == the number of covariates + 1 (assures the draws only pertain
  #to the simulated strata, and skips the predictor variables' draws)
  #dat == the data given to the model as built in
  #2023-02_ctis_mrp_vaccination_with_sex_specific_models.R
  #post == the poststratification df as built in 
  #2023-02_ctis_mrp_vaccination_with_sex_specific_models.R 
  
  nsim <- 1000
  sim <- inla.posterior.sample(n=nsim, result=model)
  
  icar_index <- grep("Predictor", rownames(sim[[1]]$latent))
  
  #start with an empty list; 
  #make a list of lists of each #latent result, 
  #and expit() it (b/c model is on logit scale)
  sims <- list()
  for (i in 1:nsim){
    sims[[i]] <- as.vector(expit(sim[[i]]$latent))
  }
  sims <- t(do.call(rbind, sims))
  sims <- as.data.frame(sims)
  
  #filter to Predictors, which again are the observations 
  #(survey data + PS variables to use for prediction of outcome)
  sims <- sims[1:max(icar_index),]
  colnames(sims) <- NULL
  colnames(sims)[1:ncol(sims)] <- paste("sim_",1:nsim,sep="")
  sims <- cbind(dat, sims)
  
  #leave sims for predictions from post-strat data only -- 
  #get rid of sims for real survey data -- i.e, observations where vax is NA
  #we know it's right because it's the same length as the post and pred_dat objects
  sims = sims %>% filter(is.na(vax))

  # Now poststratify
  
  #sex is an implicit poststratification variable, so we add it here for ease:
  poststrat_vars = paste0("sex,",poststrat_vars)
  post <- post %>% select(-vax)
  sims <- left_join(sims, post, by = c(strsplit(poststrat_vars,split = ",")[[1]]))

  if("xfips" %in% colnames(sims)){
  sims = sims %>% 
    group_by(xfips) %>%
    distinct %>%
    mutate(N = sum(estimate, na.rm = TRUE),
           N= replace_na(N,0)) %>%
    ungroup() %>%  #ditto -- can't have zeroes for mrp; handles NAs introduced by estimate/zcta_pop; in this case, NA == zero for pop counts
    select(colnames(dat),N, estimate, everything())
  }

  `%nin%` = Negate(`%in%`)
  if("xfips" %nin% colnames(sims)){
    sims = sims %>% 
      distinct() %>%
      mutate(N = sum(estimate, na.rm = TRUE),
             N= replace_na(N,0)) %>%
      ungroup() %>%  #ditto -- can't have zeroes for mrp; handles NAs introduced by estimate/zcta_pop; in this case, NA == zero for pop counts
      select(colnames(dat), estimate, everything())
  }
 
  save(sims, file = here::here("data","raw_model_output", "sims", paste0(state_choice,"_vaccination_",specs,"sex",sex_spec,"_sims.rda")))

  sims <<- sims 

}

poststratify = function(sims,pred_dat,sex_spec){
  #begin post-stratification: multiply the xfips_pop by the sims across each row 
  #(remember: each row is an observation in the PS table)
  options(scipen=999)

  temp = pred_dat %>% select(-vax) 

  #mean sim estimation * estimate in stratum /estimate in county
  rowmedians <- sims %>%
    group_by(across(colnames(temp))) %>%
    pivot_longer(cols = sim_1:sim_1000, 
                 names_to = "sim", 
                 values_to = "sim_prob")  %>%
    ungroup()

  rowmedians <-
    rowmedians %>%
      mutate(weight = estimate/N,
             sex_specific_n_i = N) %>% # N = sex-specific n_i
      group_by(across(colnames(temp))) %>%
      mutate(y = sim_prob * estimate,
             median_y = round(median(y),2), #y_ijk med
             lw_y = round(quantile(y, probs = c(0.025)),2),#y_ijk lo
             up_y = round(quantile(y, probs = c(0.975))),2) %>% #y_ijk hi
      dplyr::select(colnames(temp), median_y, estimate,
                    N, lw_y, up_y,sex_specific_n_i) %>%
      distinct() %>%
      ungroup() %>%
      rename(strata_pop = estimate)
  
  rowmedians <- 
    rowmedians %>%
      mutate(strata_pop_est = median_y/N, #now calculate the proportions in the county
             strata_pop_lw.ci = lw_y/N,
             strata_pop_up.ci = up_y/N) %>%
      #group_by(across(colnames(temp))) %>%
      # mutate(sae_estimate = as.integer(sum(strata_pop_est)),
      #        sae_lw.ci = as.integer(sum(strata_pop_lw.ci)),
      #        sae_up.ci = as.integer(sum(strata_pop_up.ci))) %>%
      # ungroup() %>%
      distinct() %>%
      mutate(model = specs)
  
  fwrite(rowmedians, file = here(paste0("data/indirect_estimates/sex_specific_estimates/",state_choice,"_mrp_indirect_estimates_vax_",specs,"_sex",sex_spec,".csv")))
  
  rowmedians <<- rowmedians
  
}

run_sex_specific_steps = function(dat,post,mf,poststrat_vars,pred_dat,sex_spec,x,specs) {

  dat = dat %>%
    filter(sex == sex_spec)
  
  post = post %>%
    filter(sex == sex_spec)
  
  pred_dat = pred_dat %>%
    filter(sex == sex_spec)

  #warning about "identity link will be used to compute the fitted values 
  #for NA data" is fine. We leave warnings verbose to catch other issues.
  run_save_model(mf,dat,sex_spec,specs)
  gc()
  print(paste0("model ",x,": done running and saving model"))
  posterior_draw(model,post,poststrat_vars,pred_dat,dat,sex_spec)
  gc()
  print(paste0("model ",x,": done posterior draw"))
  poststratify(sims,pred_dat,sex_spec)
  gc()
  print(paste0("model ",x,": done poststratifying"))
  
}

run_analysis = function(model_formulae,mod,specs,state_choice,poststrat_vars) {

    prepare_data(dat,post,poststrat_vars)
    gc()
    print(paste0("model ",x,": done preparing data"))
    #build and run sex-specific models, run poststrat, and save the output
    sex_list = list(1,0)
    
      for (i in 1:length(sex_list)){
        
        print(paste0(specs," modeling sex: ", sex_list[i])) 
        
        sex_spec = as.integer(sex_list[i])
        run_sex_specific_steps(dat,post,mf,poststrat_vars,pred_dat,sex_spec,x,specs) 
        if (sex_spec == 1){
          rowmedians1 <- rowmedians #save the first sex-specific set of rowmedians to rbind below
          }
      }
    
    #bind both sex estimates together
    estimates = rbind(rowmedians1,rowmedians)
    
    fwrite(estimates,file=here(paste0("data/indirect_estimates/complete_estimates_both_sexes/",state_choice,"_mrp_indirect_estimates_vax_",specs,"_bothsexes.csv")))

  }


#### Run Full Analysis ####
#there are 42 uniquely described models; double it because we're running it for 
#both sexes (84 models to run total)

#because of the list of formulae and the df of model data, 
#it's simplest to for-loop

#contains two data objects that are ordered meaningfully to represent 
#the same model by list/row index (i.e., model_formulae[[1]] = models[1,]):
#1. a list of formulae (`model_formulae`), 
#and
#2. a dataframe of model information (`models`) 
load(file = here("data/cleaned_input_data/model_descriptions.RData"))

models = models[8,]
model_formulae = model_formulae[8]
#clean up
models$poststrat_vars = 
  models$poststrat_vars %>% 
  gsub("\"\"", "", ., fixed=TRUE)

#### Warning: this is a bit memory intensive as written. ####
#NB: each model runs !!twice!! because we run each model for each sex 
#independently.
#Runtime for the most complex models is less than 50 seconds on 
#a Macbook Pro 2.4 GHz 8-Core Intel Core i9 with 32GB of memory.

for (x in 1:nrow(models)) {

  print(paste("***","model",x,"***",sep=" "))
  print(Sys.time())
  start.time <- Sys.time()
  
  # mod = models[40,]
  # mf = as.formula(model_formulae[[40]])
  mod = models[x,] 
  mf = as.formula(model_formulae[[x]])
  
  specs = as.character(mod$specs)
  
  poststrat_vars =
    mod$poststrat_vars %>%
    gsub('\"',"",.)

  #reload the clean data
  load_data(poststrat_vars)
  
  run_analysis(mf,mod,specs,state_choice,poststrat_vars)
  
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time,2)
  print(paste0("execution time:",time.taken))
  rm(list=c())
  gc()
}

#done running when you hear:
beepr::beep(3)




