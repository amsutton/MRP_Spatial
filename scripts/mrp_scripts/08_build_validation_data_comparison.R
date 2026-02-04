##### Build Axios-Ipsos BCC19 vaccination estimates by sex and age for CA #####

pacman::p_load(tidyverse,here,data.table,haven,
               janitor,srvyr,tidycensus,tigris,
               labelled)

here::i_am("scripts/mrp_scripts/build_bcc19_validation_data.R")

tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)
options(tigris_use_cache = TRUE)

dat = read_sav(here('data/original_ipsos_file/Covid-19 and Social Distancing_Wave3_Client_06242021.sav'))

#need to compare:
#USA CDC Estimates from Bradley 
#California CDC CA same method as Bradley for US
#Axios-Ipsos Coronavirus Tracker CA estimate 
#Zack's Ipsos CA estimate 


#q75a is about vaccination status, where 1 = received at least one dose.
#We keep either yes or no, no maybes or want-tos; sentiment not important, 
#estimating actual vaccine counts is the goal.
dat = 
    dat %>%
    clean_names() %>%
    select(ppgender,ppagecat,xfips,q75a,weight) %>% #could include education: ppeducat,ppeduc
    # filter(ppagecat != 99,
    #        q75a == 1|4) %>% #4 is not vaccinated.
    mutate(sex = ifelse(ppgender == 1, 1,0),
           age = case_when(ppagecat == 1 ~ 1,
                           ppagecat == 2 ~ 2,
                           ppagecat == 3 ~ 3,
                           ppagecat == 4 ~ 4,
                           ppagecat == 5 ~ 4,
                           ppagecat == 6 ~ 5,
                           ppagecat == 7 ~ 5),
           age2 = case_when(ppagecat == 1 ~ 1,
                            ppagecat == 2 ~ 2,
                            ppagecat == 3 ~ 2,
                            ppagecat == 4 ~ 2,
                            ppagecat == 5 ~ 2,
                            ppagecat == 6 ~ 3,
                            ppagecat == 7 ~ 3),
           xfips = sprintf("%05s",xfips)) %>%
    rename(county =xfips,
           vax=q75a) %>%
  select(county,vax,sex,age,age2,weight)

#get counties and states to match county in dat
geog = counties(year=2021)
geog = geog %>%
  clean_names() %>%
  select(statefp,geoid) %>%
  rename(state=statefp,
         county=geoid) %>%
  sf::st_drop_geometry()
  

#identify state using county geoids
dat = left_join(dat,geog, by = c("county")) 

design = as_survey_design(dat, 
                          ids = 1, 
                          weights = weight)

#these are coming out correct using srvyr
design  %>%
  filter(state == "06",
         age != 99) %>%
  group_by(interact(vax,age2)) %>%
  reframe(prop = survey_prop(vartype="se",na.rm=TRUE,level = 0.95)) %>%
  mutate(across(prop:prop_se, ~signif(.x))) %>% 
  filter(vax==1) %>%
  mutate(estimate = paste0(signif(prop,3),", SE: ",signif(prop_se,3))) %>%
  #mutate(across(prop:prop_upp, ~signif(.x))) %>% 
  fwrite(file=here('results/tables/weights_comparison/bcc19_3age_vax_proportion_estimates.csv'))

bcc19 = 
  design  %>%
  filter(state == "06") %>%
  group_by(interact(vax,age)) %>% #proportional total pop
  reframe(prop = survey_prop(vartype="se",na.rm=TRUE,level = 0.95)) %>%
  mutate(across(prop:prop_se, ~signif(.x))) %>%
  filter(vax==1)
  #mutate(across(prop:prop_upp, ~signif(.x,3)))
  
fwrite(bcc19,file=here('results/tables/weights_comparison/bcc19_5age_vax_proportion_estimates.csv'))

design  %>%
  filter(state == "06") %>%
  group_by(vax) %>%
  reframe(prop = survey_prop(vartype="se",na.rm=TRUE,level = 0.95))


#### build estimates for Axios-Ipsos Coronavirus Tracker data of same ####
dat = read_sav(here('data/axios_ipsos_coronavirus_index/Axios-Ipsos Wave 47.sav'))

#search for variables
# dat %>%  
#    labelled::look_for("CA") #CA == ppstaten 93


#dat$Q129: vaccination question
dat = dat %>%
  clean_names() %>%
  select(ppgender,ppage,ppstaten,q129,wt_final) %>% #could include education: ppeducat,ppeduc
 # filter(ppagecat != 99) %>% #4 is not vaccinated.
  mutate(sex = ifelse(ppgender == 1, 1,0),
         age = case_when(ppage %in% c(18:24) ~ 1,
                         ppage %in% c(25:64) ~ 2,
                         ppage %in% c(65:94) ~ 3#,
                         # ppage %in% c(45:64) ~ 4,
                         # ppage >= 65 ~ 5
                         ),
         # age2 = case_when(age == 1 ~ 1,
         #                  age == 2 ~ 2,
         #                  age == 3 ~ 2,
         #                  age == 4 ~ 2,
         #                  age == 5 ~ 3),
         ppstaten = sprintf("%02s",ppstaten)) %>%
  rename(state = ppstaten,
         vax=q129,
         weight=wt_final) %>%
  select(state,vax,sex,age,weight) %>%
  mutate(vax = ifelse(vax==2,1,vax))


design = as_survey_design(dat, 
                          ids = 1, 
                          weights = weight)

#these are coming out correct using srvyr
design  %>%
  filter(state == "93") %>%
  group_by(interact(vax,age)) %>%
  #group_by(interact(sex,age2)) %>% #proportional total pop
  reframe(prop = survey_prop(vartype="se",na.rm=TRUE,level = 0.95)) %>%
  mutate(across(prop:prop_se, ~signif(.x))) %>%
  filter(vax==1) %>%
  mutate(estimate = paste0(signif(prop,3),", SE: ",signif(prop_se,3))) %>%
  #mutate(across(prop:prop_upp, ~signif(.x))) %>% 
  fwrite(file=here('results/tables/weights_comparison/axiosipsos_coronavirustracker_3age_vax_proportion_estimates.csv'))
# 
# tracker = design  %>%
#   filter(state == "93",
#          vax==1|2) %>%
#   group_by(age) %>%
#   #group_by(interact(sex,age)) %>% #proportional total pop
#   reframe(prop = survey_prop(vartype="se",na.rm=TRUE,level = 0.95)) %>%
#   mutate(across(prop:prop_se, ~signif(.x,3)))
#  # mutate(across(prop:prop_upp, ~signif(.x,3)))
# 
# fwrite(tracker,file=here('results/tables/weights_comparison/axiosipsos_coronavirustracker_5age_vax_proportion_estimates.csv'))

tab = table(dat$vax,dat$sex,dat$age)
tab = ftable(tab) %>% as.data.frame()
tab %>%
  rename(vax=Var1,
         sex=Var2,
         age=Var3,
         count=Freq) %>%
  group_by(sex,age) %>%
  filter(vax!=-1) 

design  %>%
  filter(state == "93",
         vax==1|2) %>%
  group_by(interact(vax)) %>%
  reframe(prop = survey_prop(vartype="se",na.rm=TRUE,level = 0.95))


#compare tracker with bcc1compare
compare = left_join(bcc19,tracker,by=c("age"),suffix=c('.bcc19','.tracker'))
compare$prop_diff = compare$prop.bcc19 - compare$prop.tracker
# compare$prop_diff_low = compare$prop_low.bcc19 - compare$prop_low.tracker
# compare$prop_diff_upp = compare$prop_upp.bcc19 - compare$prop_upp.tracker
glimpse(compare)
compare %>%
  select(age,prop_diff) #:prop_diff_upp)

fwrite(compare,here('results/tables/weights_comparison/axiosipsoscoronavirustracker_vs_bcc19_5age_vax_proportion_estimates.csv'))



#read in bradley estimates
#national estimates taken from benchmark script
##the date used in our paper currently
#bradley_national %>% filter(date=='06/30/2021')
#0.7205854 of the adult population vaccinated on June 30, 2021






# 
# #Bradley et al CDC national estimates May 26 2021 (date says for May 25th in df?)
# brad = fread(here('/Users/Aja/Downloads/cdc_cleaned_2021-05-26.csv'))
# brad = brad %>%
#   filter(date =="2021-05-25") 
# 
# compare$bradleyetal_prop_national = brad$pct_pop_vaccinated
# 
# 
# #Bradley et al CDC national estimates May 26 2021 (date says for May 25th in df?)
# brad = fread(here('/Users/Aja/Downloads/benchmark_2021-05-26.csv'))
# table(brad$state)
# brad = brad %>%
#   filter(state=="CA",
#          date=="2021-05-25")
# brad
# 
# compare$bradleyetal_prop_ca = brad$pct_pop_vaccinated
# 
# 
# #maybe helpful? Separate estimates for sexes, and age groups for CA, May 30-June 26, 2021.
# nis = fread(here('/Users/Aja/Downloads/National_Immunization_Survey_Adult_COVID_Module__NIS-ACM___Vaccination_Status_and_Intent_by_Demographics.csv'))
# 
# nis = nis %>%
#   clean_names() %>%
#   filter(geography == "California",
#          year == "2021",
#          time_period == "May 30 - June 26",
#          str_detect(group_category,"18|50|65|Female|Male"),
#          str_detect(indicator_category,"1 dose|completed primary")) %>%
#   select(geography,group_name,group_category,indicator_category,time_period,year,
#          estimate_percent,x95_percent_ci_percent,sample_size)
# sum(compare$prop.bcc19)
# 
# glimpse(nis)
# 
# 
# #### covid vaccination counts by CA county, 18+ and 65+ ####
# 
# countycdc = fread(here('/Users/Aja/Documents/R Directories/MRP_Spatial/data/cdc_data/COVID-19_Vaccinations_in_the_United_States_County.csv'))
# 
# #from countycdc; won't show up in next call properly
# census2019_65_pop = countycdc %>%
#   clean_names() %>%
#   filter(recip_state == "CA",
#          fips !="UNK") %>%
#   group_by(fips) %>% 
#   select(census2019_65plus_pop) %>%
#   distinct() %>%
#   na.omit()
# 
# countycdc = countycdc %>%
#   clean_names() %>%
#   filter(recip_state == "CA",
#          fips != "UNK",
#          str_detect(date,"06/30/2021|06/12/2021|05/26/2021")) %>%
#   select(date,fips,recip_county,completeness_pct,
#          administered_dose1_recip_18plus,census2019_18plus_pop,
#          administered_dose1_recip_65plus)
# 
# countycdc = left_join(countycdc,census2019_65_pop,by="fips")


### CTIS estimates ####


#same denominator as Bradley
state_pop_totals= fread(here('data/bradleyetal/SCPRC-EST2019-18+POP-RES.csv'))

# get state abbreviations from state names
state_pop_totals <- merge(state_pop_totals,
                          cbind(state.abb, state.name),
                          by.x = "NAME",
                          by.y = "state.name",
                          all.x = TRUE)

# rename columns)
setnames(state_pop_totals,
         old = c("NAME", "state.abb", "POPEST18PLUS2019"),
         new = c("state_name", "state", "pop_total"))

state_pop_totals = state_pop_totals %>%
  filter(state=="CA")

filelist = list.files(
  path = here("data/indirect_estimates/complete_estimates_both_sexes/"),
  full.names = TRUE) # include the directory in the result)

dat = fs::dir_ls(here("data/indirect_estimates/complete_estimates_both_sexes/"), regexp = "\\.csv$")
dat = dat %>%
  map_dfr(read_csv)

dat %>%
  select(model_name,age,median_y,lw_y,up_y) %>%
  mutate(age = case_when(age == 1 ~ 1,
                         age == 2 ~ 2,
                         age == 3 ~ 2,
                         age == 4 ~ 2,
                         age == 5 ~ 3)) %>%
  group_by(model_name,age) %>%
  reframe(prop_sum = sum(median_y)/state_pop_totals$pop_total,
          lw_y = sum(lw_y)/state_pop_totals$pop_total,
          up_y = sum(up_y)/state_pop_totals$pop_total) %>% #denominator is for 18+ population
  group_by(model_name) %>%
  mutate(across(prop_sum:up_y, ~signif(.x,3)),
         estimate = paste0(prop_sum," (",lw_y,", ",up_y,")"),
         total = sum(prop_sum),
         total_lw = sum(lw_y),
         total_up  = sum(up_y))%>% 
  fwrite(here('/Users/Aja/Documents/R Directories/MRP_Spatial/results/tables/weights_comparison/ctis_5ages_aggregated_to_3ages_vax_estimates.csv'))





