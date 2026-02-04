#### Variation from Bradley et al.'s Benchmarking Code ####

# source: https://github.com/vcbradley/ddc-vaccine-US/blob/main/functions/functions_clean_CDC.R


pacman::p_load(tidyverse,here,data.table,haven,
               janitor,srvyr,tidycensus,tigris,
               labelled)

here::i_am("scripts/mrp_scripts/build_bradley_benchmark_estimates.R")


###National Estimates from Paper ####
#Bradley et al build national-level benchmarked estimates of vaccination like so:

getPctVaxxedUnder18 <- function() {
  
  cdcdata_age = fread(here("data/cdc_data/COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional.csv"))
  
  cdcdata_age = cdcdata_age %>%
    janitor::clean_names() %>%
    mutate(demographic_group = demographic_category,
           age_group = demographic_group,
           n_vaccinated = administered_dose1,
           n_fully_vaccinated = series_complete_yes)
  

  # limit to rows with age data
  cdcdata_age <- cdcdata_age %>%
    filter(!str_detect(age_group,"ale")) %>%
    mutate(age_flag = as.numeric(grepl("Ages", age_group)),
           #date = as.Date(date)
    ) %>%
    filter(age_flag == 1) %>%
    select(-age_flag)
  
  
  # calculate total number of vaccinations 1) in total 2) over 18 and 3) under 18
  # all TX reported as
  cdcdata_age <- cdcdata_age %>%
    group_by(date) %>%
    mutate(
      n_fully_vaccinated = as.numeric(n_fully_vaccinated),
      n_vaccinated = as.numeric(n_vaccinated),
      under_18_flag = as.numeric(age_group %in% c("Ages_<18yrs", "Ages_<12yrs", "Ages_12-15_yrs", "Ages_16-17_yrs"))
    ) %>%
    summarize(
      n_vax_notx = sum(n_vaccinated),
      n_vax_notx_fully = sum(n_fully_vaccinated),
      n_vax_notx_over18 = sum(n_vaccinated * (1 - under_18_flag)),
      n_vax_notx_fully_over18 = sum(n_fully_vaccinated * (1 - under_18_flag)),
      n_vax_notx_under18 = sum(n_vaccinated * under_18_flag),
      n_vax_notx_fully_under18 = sum(n_fully_vaccinated * under_18_flag)
    ) %>%
    mutate(
      pct_vax_under18 = n_vax_notx_under18 / n_vax_notx,
      pct_vax_fully_under18 = n_vax_notx_fully_under18 / n_vax_notx_fully
    )
  
  # we only really need
  # n_vax_notx: the total number of vaccinations outside of TX on each day
  # pct_vax_under18: the proportion of total one-dose vaccinations administered
  # to under-18s on each day
  cdcdata_age <- cdcdata_age %>%
    select(date, n_vax_notx, n_vax_notx_over18, pct_vax_under18)
  
  return(cdcdata_age)
}



getStatePopTotals <- function() {
  
  # "https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/state/detail/SCPRC-EST2019-18+POP-RES.csv"
  
  # read in pop totals
  state_pop_totals= fread(here('data/bradleyetal/SCPRC-EST2019-18+POP-RES.csv'))

  # get state abbreviations from state names
  state_pop_totals <- merge(state_pop_totals,
                            cbind(state.abb, state.name),
                            by.x = "NAME",
                            by.y = "state.name",
                            all.x = TRUE)
  
  # fix missing abbrevs
  state_pop_totals[NAME == "United States", state.abb := "US"]
  state_pop_totals[NAME == "District of Columbia", state.abb := "DC"]
  state_pop_totals[NAME == "Puerto Rico Commonwealth", state.abb := "PR"]
  
  # rename columns)
  setnames(state_pop_totals,
           old = c("NAME", "state.abb", "POPEST18PLUS2019"),
           new = c("state_name", "state", "pop_total"))
  
  return(state_pop_totals)
}

build_bradleyetal_benchmark_estimates = function(state_choice) {
  cdcdata = fread(here('data/bradleyetal/COVID-19_Vaccination_Trends_in_the_United_States_National_and_Jurisdictional.csv'))
  
  cdcdata = cdcdata %>%
    janitor::clean_names() %>%
    mutate(report_type = date_type,
           state = location,
           n_vaccinated_allages = admin_dose_1_cumulative,
           n_vaccinated_allages_fully = series_complete_cumulative) %>%
    filter(state == state_choice,
           report_type == "Admin")
  
  #above does below
  # rename columns
  # setnames(cdcdata, c("report_type", "date", "program", "n_vaccinated_allages", "n_vaccinated_allages_fully"))
  #
  # # filter to only administered doses at the national level (includes LTC doses)
  # cdcdata <- cdcdata %>%
  #   filter(report_type == "Admin" & program == "US") %>%
  #   mutate(
  #     date = as.Date(date),
  #     n_vaccinated_allages = as.numeric(n_vaccinated_allages),
  #     n_vaccinated_allages_fully = as.numeric(n_vaccinated_allages_fully)
  #   ) %>%
  #   rename(state = program)
  
  
  # we also need cdc doses by age data -- TX does not report doses by age, so need to impute
 cdc_by_age <- getPctVaxxedUnder18()

  # join age data onto rest of CDC data and use age data to impute the number of over-18 doses administered in TX
  cdcdata <- left_join(cdcdata, cdc_by_age, by = "date") %>%
    replace_na(list(n_vax_notx = 0, pct_vax_under18 = 0, n_vax_notx_over18 = 0)) %>%
    mutate(
      n_vax_tx = n_vaccinated_allages - n_vax_notx, # impute TX-only doses
      n_vax_tx_over18 = n_vax_tx * (1 - pct_vax_under18), # impute number of TX doses administered to adults
      n_pop_vaccinated = n_vax_tx_over18 + n_vax_notx_over18 # get total number of US adults with >=1 dose by adding TX and non-TX adult doses
    ) %>%
    select(date, state, n_pop_vaccinated)


  # get state pop totals to calculate % uptake
  state_pop_totals <- getStatePopTotals() %>%
    filter(state == state_choice) %>%
    select(-state_name,-STATE) 
  
  cdcdata <- left_join(cdcdata, state_pop_totals, by = c("state")) %>%
    mutate(pct_pop_vaccinated = n_pop_vaccinated / pop_total)
  
  return(cdcdata)
}


bradley_national = build_bradleyetal_benchmark_estimates("US")

#the date used for the benchmarks in the paper: 
#NB: this is based on the adult population, NOT the total (Bradley uses total population!)
bradley_national %>% filter(date=='05/25/2021')
  
#the date used in our paper currently
bradley_national %>% filter(date=='06/30/2021')
#0.7205854 of the adult population vaccinated on June 30, 2021

#### California estimates using same methodology ####
bradley_california = build_bradleyetal_benchmark_estimates("CA")
bradley_california %>% filter(date=='06/30/2021') #0.8122073 vaccinated in CA June 30, 2021

bradley_california %>% filter(date=='05/26/2021') 


#### Modification: estimates by age (as much as we can) ####

bradley_benchmark_estimates_byage <- function(state_choice) {
  
  #this can actually be calculated without the statelevel estimates -- they sum to the SAME denominators as what is in here as $census
  cdcdata_age = fread(here("data/cdc_data/COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional.csv"))

  cdcdata_age = cdcdata_age %>%
    janitor::clean_names() %>%
    filter(location == state_choice) %>%
    mutate(demographic_group = demographic_category,
           age_group = demographic_group,
           n_vaccinated = administered_dose1,
           n_fully_vaccinated = series_complete_yes,
           state= location)


  # calculate total number of vaccinations by age group, exclude non-adults
cdcdata_age <- 
  cdcdata_age %>%
    group_by(date) %>%
    filter(str_detect(age_group,"18-24|25-49|50-64|65+"),
           !str_detect(age_group,"ale|Sex|Unknown|65-74_yrs")) %>% 
    mutate(
      age_group = case_when(str_detect(age_group,"18-24") ~ 1,
                            str_detect(age_group,"25-49|50-64") ~ 2,
                            str_detect(age_group,"65+") ~ 3),
      total_census = sum(census)) %>% 
    group_by(date,age_group) %>%
    reframe(
      n_vaccinated = sum(n_vaccinated),
     # n_fully_vaccinated = sum(n_fully_vaccinated), #we are only concerned with min. 1 dose
      total_census,
      census = sum(census),
      proportion_agegroup_vaccinated = n_vaccinated/census) %>% #at least one dose, not fully
      group_by(date) %>%
      distinct() %>%
      mutate(proportion_ofadults_vaccinated = n_vaccinated/total_census) %>%
    select(date,age_group,census,n_vaccinated,census,
           proportion_agegroup_vaccinated,proportion_ofadults_vaccinated )
    
  

  return(cdcdata_age)
}


bradley_california_byage = bradley_benchmark_estimates_byage(state_choice = "CA")

#the date used for the benchmarks in the paper
bradley_california_byage %>% filter(str_detect(date,'06/30/2021'))%>% glimpse

#the date used in our paper currently
june30_bradley_byage = 
  bradley_california_byage %>% 
  filter(str_detect(date,'06/30/2021')) 

june30_bradley_byage = 
  june30_bradley_byage %>%
  mutate(se_agegroup = sqrt(proportion_agegroup_vaccinated * (1-proportion_agegroup_vaccinated)/sum(proportion_agegroup_vaccinated))) %>%
    select(age_group,proportion_agegroup_vaccinated,se_agegroup,n_vaccinated,census) %>%
  mutate(across(proportion_agegroup_vaccinated:se_agegroup, ~ round(.x, 3)),
         estimate = paste0(proportion_agegroup_vaccinated,", SE: ",se_agegroup)) %>%
  rename(census_count=census)

fwrite(june30_bradley_byage,here('results/tables/weights_comparison/bradleymethod_june30_cdc_3age_vax_proportion_estimates.csv'))


bradley_usnational_byage = bradley_benchmark_estimates_byage(state_choice = "US")

#the date used for the benchmarks in the paper
bradley_usnational_byage %>% filter(str_detect(date,'06/30/2021'))%>% glimpse

#the date used in our paper currently
june30_bradley_byage = 
  bradley_usnational_byage %>% 
  filter(str_detect(date,'06/30/2021')) 

june30_bradley_byage = 
  june30_bradley_byage %>%
  mutate(se_ofadults = sqrt(proportion_ofadults_vaccinated * (1-proportion_ofadults_vaccinated)/sum(proportion_ofadults_vaccinated))) %>%
  select(age_group,proportion_ofadults_vaccinated,se_ofadults,,n_vaccinated,census) %>%
  mutate(across(proportion_ofadults_vaccinated:se_ofadults, ~ round(.x, 3)),
         estimate = paste0(signif(proportion_ofadults_vaccinated,3),", SE: ",signif(se_ofadults,3))) %>%
  rename(census_count=census)

fwrite(june30_bradley_byage,here('results/tables/weights_comparison/bradleymethod_june30_usnational_cdc_3age_vax_proportion_estimates.csv'))


