#### State-level model with poststratification ####

if (!require("pacman")) install.packages("pacman")
pacman::p_load(here,tidyverse,data.table,SUMMER,tidycensus)
library(INLA)


run_analysis = function(sex_spec,poststrat_vars){

  #### load data ####

#This is a simple state level GLM without spatial smoothing.

#load dat and post for analysis

  #post
  load(file = here("data/cleaned_input_data/clean_postrat_age_sex_county_with_edu.rda"))

  #aggregate to state level
  post = post %>%
    group_by(sex,age,edu) %>%
    reframe(vax=NA,
            estimate =sum(estimate)) %>%
    select(vax,everything()) %>%
    filter(!is.na(edu))

  #dat
  load(file = here("data/cleaned_input_data/clean_data_county_with_edu.rda"))

  dat = dat %>%
    select(vax,sex,age,edu)

  pred_dat = post %>% select(colnames(dat))

  #add the empty prediction categories to the bottom of the real data
  dat = rbind(dat, pred_dat)

  dat = dat %>%
    filter(sex == sex_spec)


  #formula
  mf = vax ~ 1 + age * edu # equivalent to age + education + age X education

  #number of covariates + number of outcome variables (1) + 1
  sim_start = length(mf)+2

#### run model ####
  model <- inla(mf,
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

  save(model, file=here(paste0('data/model_objects/state_level_ageedu_interaction_model_object_',sex_spec,'.RDS')))

#### posterior draw ####
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

       save(sims, file = here::here(paste0("data","raw_model_output", "sims","vaccination_state_level_ageedu_interaction_sims_",sex_spec,".rda")))

#### poststratification ####

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
     mutate(sex_specific_n_i = N) %>% # N = sex-specific n_i
     group_by(across(colnames(temp))) %>%
     mutate(y = sim_prob * estimate,
            median_prob = median(sim_prob),
            median_y = round(median(y),2), #y_ijk med
            lw_y = round(quantile(y, probs = c(0.025)),2),#y_ijk lo
            up_y = round(quantile(y, probs = c(0.975))),2) %>% #y_ijk hi
     dplyr::select(colnames(temp), median_y, median_prob,estimate,
                   N, lw_y, up_y,sex_specific_n_i) %>%
     distinct() %>%
     ungroup() %>%
     rename(strata_pop = estimate)

rowmedians <-
     rowmedians %>%
     mutate(strata_pop_est = median_y/strata_pop, #proportion in sex-age-edu stratum
            strata_pop_lw.ci = lw_y/strata_pop,
            strata_pop_up.ci = up_y/strata_pop,
            total_stratum_pop_est = median_y/N, #total population proportion by stratum
            total_stratum_pop_lw.ci = lw_y/N,
            total_stratum_pop_up.ci = up_y/N) %>%
  group_by(sex,age,edu) %>%
            mutate(total_pop_est = sum(median_y)/sum(strata_pop), #total population proportion by stratum
            total_pop_lw.ci = sum(lw_y)/sum(strata_pop),
            total_pop_up.ci = sum(up_y)/sum(strata_pop)) %>%
     #group_by(across(colnames(temp))) %>%
     # mutate(sae_estimate = as.integer(sum(strata_pop_est)),
     #        sae_lw.ci = as.integer(sum(strata_pop_lw.ci)),
     #        sae_up.ci = as.integer(sum(strata_pop_up.ci))) %>%
     # ungroup() %>%
     distinct()


   fwrite(rowmedians, file = here(paste0("data/indirect_estimates/sex_specific_estimates/mrp_indirect_estimates_vax_statelevel_ageedu_interaction_",sex_spec,".csv")))

   rm(model,post,pred_dat,rowmedians,sim,sims,temp,dat,icar_index)
   gc()
}

#poststrat_vars = "age,edu"

#sex_spec = 0 #1 = male, 0 = female

run_analysis(sex_spec = 1,poststrat_vars="age,edu")
run_analysis(sex_spec = 0,poststrat_vars="age,edu")




