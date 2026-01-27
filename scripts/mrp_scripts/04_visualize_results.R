#### MRP with Spatial Priors using CTIS Data ####

#Build a series of models representing vaccination probability (min. 1 dose)
#For June 30, 2021 using the COVID-19 Trends and Impacts Survey (CTIS)

pacman::p_load(tidyverse,here,fs,ggpubr,ggnewscale,
               tidycensus,tigris,janitor,ggbreak,
               sf,viridis,data.table)

here::i_am("scripts/mrp_scripts/04_visualize_results.R")

options(scipen=999)
options(tigris_use_cache = TRUE)
tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)

#Only have to run the one time (joins all the modeled results together) 
filelist <- list.files(
  path = here("data/indirect_estimates/complete_estimates_both_sexes/"),
  full.names = TRUE) # include the directory in the result)

dat = fs::dir_ls(here("data/indirect_estimates/complete_estimates_both_sexes/"), regexp = "\\.csv$")
dat = dat %>%
  map_dfr(read_csv)
save(dat, file = here("data/cleaned_input_data/clean_complete_model_indirect_estimates.rda"))

rm(filelist)
# 
# load("data/cleaned_input_data/clean_complete_model_indirect_estimates.rda")
# 

#load model information for intuitive labeling in figures

load(file = here("data/cleaned_input_data/model_descriptions.RData"))

#clean up
rm(model_formulae) #don't need it


models$poststrat_vars = 
  models$poststrat_vars %>% 
    gsub("\"\"", "", ., fixed=TRUE)

models = 
  models %>%
    select(model_name, specs)

write.csv(models,here("results/tables/Table2_model_list.csv"))

#model_name is now there for labels!
dat = right_join(dat,models, by = c("model"="specs"),relationship = "many-to-many")

#simplify
dat = dat %>%
  mutate(model_name=ifelse(model_name=="fixed effect on age, fixed effect on education, age-education RW1 interaction, IID by county","age-education RW1 interaction, IID by county",model_name))

#add geographic data for mapping
geo = counties(state = "CA", cb = FALSE, year = 2020)
geo = 
  geo %>% 
    clean_names() %>% 
    select(geoid,name) %>% 
    st_drop_geometry()

dat = left_join(dat, geo, by = c("xfips" = "geoid"),relationship = "many-to-many")

#### CDC STATE LEVEL estimates for CA by sex and age (3 categories match with CTIS data) ####
cdc = fread(here("data/cdc_data/COVID-19_Vaccination_Age_and_Sex_Trends_in_the_United_States__National_and_Jurisdictional.csv"))

#clean version for table in paper
cdc_bysexage_table = 
  cdc %>% 
    clean_names() %>%
    mutate(date = as.Date(date, "%m/%d/%Y")) %>%
    filter(location == "CA",
           date == "2021-06-30" | date == "2021-06-12",
           str_detect(demographic_category, "ale")) %>%
    mutate(sex = ifelse(str_detect(demographic_category, "Male"), "male", "female"),
           age = case_when(str_detect(demographic_category, "18") ~ "18-24 years", #18-24 yrs
                           str_detect(demographic_category, "25-4") ~ "25-64 years", #25-49
                           str_detect(demographic_category, "50") ~ "25-64 years", #50-64, aggregate with above to make 25-64 to match ctis data
                           str_detect(demographic_category, "65+") ~ "65+ years")) %>%
    select(date, sex, age, administered_dose1) %>%
    group_by(date, sex, age) %>%
    reframe(administered_dose1 = sum(administered_dose1, na.rm=TRUE)) %>%
    distinct() %>%
    na.omit() 

#population per xfips
pop =
  dat %>%
    select(xfips,N) %>%
    distinct()

#total state population
pop = sum(pop$N, na.rm = TRUE)


#taken from ACS 5 year 2023 -- Table B15001
m_1824_pop = 1830374
m_2564_pop= (2995085+2759655+4847781)
m_65plus_pop = 2702808

f_1824_pop = 1742201
f_2564_pop= (2839782+2639320+4865089)
f_65plus_pop = 3291678

cdc_bysexage_table = 
  cdc_bysexage_table %>%
    group_by(date,sex,age) %>%
    mutate(total_pop = pop,
           totalprop_estimate = administered_dose1/total_pop, #proportion of age/sex group with 1 dose, aged +18 years
           model = "CDC official count") %>%
    group_by(date) %>%
    mutate(totalpop_administered_doses = sum(administered_dose1), #total count, aged +18 years
           totalprop_estimate = round(totalprop_estimate,3),
           specpop = case_when(sex == "male" & age == "18-24 years" ~ m_1824_pop,
                               sex == "male" & age == "25-64 years" ~ m_2564_pop,
                               sex == "male" & age == "65+ years" ~ m_65plus_pop,
                               sex == "female" & age == "18-24 years" ~ f_1824_pop,
                               sex == "female" & age == "25-64 years" ~ f_2564_pop,
                               sex == "female" & age == "65+ years" ~ f_65plus_pop),
           stratumprop_estimate = administered_dose1/specpop) %>%
    ungroup() %>%
    select(date:administered_dose1,totalprop_estimate,totalpop_administered_doses,stratumprop_estimate)

fwrite(cdc_bysexage_table, file = "results/tables/table4_cdc_baseline_data_by_sexage.csv")


#both dates included because I noticed that June 12 looks a 
#lot like some of the earlier estimates
cdc = 
  cdc %>% 
    clean_names() %>%
    mutate(date = as.Date(date, "%m/%d/%Y")) %>%
    filter(location == "CA",
           date == "2021-06-30" | date == "2021-06-12",
           str_detect(demographic_category, "ale")) %>%
    mutate(sex = ifelse(str_detect(demographic_category, "Male"), 1, 0),
           age = case_when(str_detect(demographic_category, "18") ~ 1, #18-24 yrs
                           str_detect(demographic_category, "25-4") ~ 2, #25-49
                           str_detect(demographic_category, "50") ~ 2, #50-64, aggregate with above to make 25-64 to match ctis data
                           # str_detect(demographic_category, "65-74") ~ 3,
                           str_detect(demographic_category, "65+") ~ 3)) %>%
    select(date, sex, age, administered_dose1) %>%
    group_by(date, sex, age) %>%
    reframe(administered_dose1 = sum(administered_dose1, na.rm=TRUE)) %>%
    distinct() %>%
    na.omit() 

cdc = 
  cdc %>%
    group_by(date,sex,age) %>%
    mutate(total_pop = pop,
           estimate = administered_dose1/total_pop,
           model = "CDC official count") %>%
    group_by(date) %>%
    mutate(popest = sum(administered_dose1, na.rm=TRUE)) %>%
    ungroup

#rm(pop)

#### Read in direct estimates ####

dir = fread(here("data/direct_estimates/ctis_vax_direct_estimates_not_bootstrapped_with_edu_simplevax.csv"))

#### load state of CA county vaccination data from the CDC ####
countyvax = fread(here("data/cdc_data/COVID-19_Vaccinations_in_the_United_States_County.csv"))

countyvax = 
  countyvax %>%
    clean_names() %>%
    filter(recip_state == "CA",
           date == "06/30/2021" | date == "06/12/2021") %>%
    select(date, recip_county, fips,
           administered_dose1_recip_18plus,
           administered_dose1_recip_18plus_pop_pct)

#get acs counts of total population by county in CA
acs = get_acs(geography = "county",
              state = "CA",
              table = "B15001",
              survey = "acs5",
              geometry = FALSE,
              year = 2020)

#total population
acs = 
  acs %>% 
    clean_names() %>% 
    filter(variable == "B15001_001") %>% 
    select(geoid,estimate)

countyvax = left_join(countyvax,acs, by = c("fips"="geoid"))

countyvax = 
  countyvax %>% 
    rename(countypop = estimate)

#prop in original file is lower than CDC estimates
#because the denominator is different (we use 2021 5-year ACS rather
#than what they seem to have used: the 2019!)
countyvax = 
  countyvax %>%
    mutate(county = str_remove(recip_county," County")) %>% #make it match dat naming
    group_by(date, county) %>%
    mutate(prop = administered_dose1_recip_18plus/countypop) %>%
    ungroup %>%
    na.omit() %>% #CDC has a phantom unknown county with no data
    select(date,county,fips,administered_dose1_recip_18plus,countypop,prop)

# countyvax %>%
#   group_by(date) %>%
#   mutate(statedoses_fromcounty = sum(administered_dose1_recip_18plus),
#           stateprop = statedoses_fromcounty/sum(acs$estimate)) %>%
#   select(date,county,fips,prop,stateprop)

write.csv(countyvax, file = here("results/tables/table5_county_counts_cdc.csv"))

#for benchmark lines in plots at the state level

cdc_state = 
  cdc %>%
    #group_by(date, sex, age) %>%
    #mutate(estimate = sum(administered_dose1)) %>%
    group_by(date) %>%
    mutate(prop_est = sum(administered_dose1)/total_pop,
              model="CDC estimate") %>%
    ungroup() %>% 
    distinct

temp = 
  dat %>% 
    select(xfips,sex,N) %>% 
    distinct() %>%
    select(-sex) %>%
    group_by(xfips) %>%
    reframe(N_totalcounty = sum(N)) 


#### State Results ####

#simplify
models = 
  models %>%
  mutate(model_name=ifelse(model_name=="fixed effect on age, fixed effect on education, age-education RW1 interaction, IID by county","age-education RW1 interaction, IID by county",model_name))

# 
# `%nin%` = Negate(`%in%`)
# spatial_statelevel_order = models %>% 
#   filter(str_detect(model_name,"BYM2")) %>% 
#   select(model_name)
# spatial_statelevel_order =unlist(spatial_statelevel_order$model_name)
# spatial_statelevel_order = c("CDC reported count",
#                               "direct survey estimate",
#                               spatial_statelevel_order)

dat = left_join(dat,temp, by = 'xfips')
dat = left_join(dat,models, by ="model_name")

dat = 
  dat %>% 
    mutate(specs = str_replace_all(specs,"xfips","county"))

dat = 
  dat %>%
    select(xfips:strata_pop_up.ci,model_name:specs,N_totalcounty)


#I need:
# total number in each stratum in each county (==median_y, lw_y, up_y)
# total number in each stratum in the state (first mutate())
# proportion of each stratum in each county (already present)
# proportion of each stratum in the state (second mutate())

dat.state = 
  dat %>%
    group_by(sex,age,model_name,specs) %>%
    mutate(
      #N = sum(strata_pop),
           state_sae = sum(median_y),#total number each stratum in the state
           state_saelw = sum(lw_y),
           state_saeup = sum(up_y)) %>%
    distinct() %>% 
    mutate(state_prop = state_sae/sum(unique(strata_pop)), #total proportion in each stratum in the state (denominator is the sum of all unique values of population counts in the county)
            state_proplw = state_saelw/sum(unique(strata_pop)),
            state_propup = state_saeup/sum(unique(strata_pop))) %>% 
  distinct() %>%
    mutate(type = case_when(str_detect(specs, "bym2") ~ "Spatial",
                            str_detect(specs, "iid") ~ "Non-Spatial"),
           specs = str_replace_all(specs,"_"," ")) %>%
  select(xfips,sex:edu,specs,N,state_sae:state_propup,type,N_totalcounty,strata_pop) %>%
  ungroup()


#temp =dat.state %>% filter(xfips == "06001"|xfips == "06003")

# temp %>%
#  select(-xfips,-N) %>%
#   group_by(sex,model_name,specs) %>%
#   mutate(state_proportion_total = sum(state_sae)/pop) %>% glimpse



colors = wesanderson::wes_palette(3, name = "Darjeeling1", type = "discrete")
#colors = colors[2:3]

# 
# #if you wish to include the state-level model (age*education interaction): 
# #i.e., script 06_state_level_model.R
# f = fread(here("data/indirect_estimates/sex_specific_estimates/mrp_indirect_estimates_vax_statelevel_ageedu_interaction_0.csv"))
# m = fread(here("data/indirect_estimates/sex_specific_estimates/mrp_indirect_estimates_vax_statelevel_ageedu_interaction_1.csv"))
# 
# sm = rbind(f,m)
# rm(f,m)
# gc()
# 
# sm2 = sm %>%
#   group_by(age,sex) %>%
#   reframe(popest = sum(median_y)/sum(strata_pop),
#           popest_lw = sum(lw_y)/sum(strata_pop),
#           popest_up = sum(up_y)/sum(strata_pop),
#           median_y,lw_y,up_y) %>%
#   distinct() %>%
#   mutate(age = factor(case_when(age == 1 ~ "18-24 years",
#                                 age == 2 ~ "25-64 years",
#                                 age == 3 ~ "65+ years")),
#          sex = factor(ifelse(sex==1,"male","female")),
#          model_name = "state-level age-education interaction")
# 
# sm = sm %>%
#   group_by(age,sex,edu) %>%
#   reframe(popest = sum(median_y)/sum(strata_pop),
#           popest_lw = sum(lw_y)/sum(strata_pop),
#           popest_up = sum(up_y)/sum(strata_pop),
#          median_y,lw_y,up_y) %>%
#   ungroup() %>%
#   mutate(state_est = sum(median_y),
#          state_est_lw =sum(lw_y),
#          state_est_up =sum(up_y),
#          age = factor(case_when(age == 1 ~ "18-24 years",
#                          age == 2 ~ "25-64 years",
#                          age == 3 ~ "65+ years")),
#          sex = factor(ifelse(sex==1,"male","female")),
#          edu_order = factor(edu),
#          edu = case_when(edu == 1 ~ "less than\nhigh school",
#                          edu == 2 ~ "high school\nor equivalent",
#                          edu == 3 ~ "some college",
#                          edu == 4 ~ "associate's\ndegree",
#                          edu == 5 ~ "bachelor's\ndegree",
#                          edu == 6 ~ "professional or\ngraduate degree"),
#          model_name = "state-level age-education interaction")

statelevel_order = c("CDC reported count",
                     "CTIS survey estimate",
                  #   "state-level age-education interaction",
                     paste0(models[1,1]),
                     paste0(models[2,1]),
                     paste0(models[3,1]),
                     paste0(models[4,1]),
                     paste0(models[7,1]),
                     paste0(models[5,1]),
                     paste0(models[8,1]),
                     paste0(models[6,1]),
                     paste0(models[9,1])
)

# dat.state %>%
#   select(-xfips) %>%
#   distinct() %>%
#   na.omit() %>% 
#   group_by(sex,age,edu,model_name,specs) %>%
#   mutate(age = case_when(age == 1 ~ "18-24 years",
#                          age == 2 ~ "25-64 years",
#                          age == 3 ~ "65+ years"),
#          sex = ifelse(sex==1,"male","female"),
#          edu_order = factor(edu),
#          edu = case_when(edu == 1 ~ "less than\nhigh school",
#                                edu == 2 ~ "high school\nor equivalent",
#                                edu == 3 ~ "some college",
#                                edu == 4 ~ "associate's\ndegree",
#                                edu == 5 ~ "bachelor's\ndegree",
#                                edu == 6 ~ "professional or\ngraduate degree")) %>% 
#   arrange(model_name, sex,age,edu) %>% #filter(xfips == "06001", sex == "female", age == "18-24 years",str_detect(edu, "associate's")) %>% view 
#   ggplot(aes(x=factor(model_name), y= state_prop, color = factor(age), shape = factor(sex))) +
#   geom_pointrange(aes(ymin = state_proplw, ymax = state_propup),size = 0.3, position = position_dodge(0.5)) + #,position="jitter"
#   # geom_pointrange(data=sm, aes(y=popest,ymin = popest_lw, ymax = popest_up),size = 0.3, position = position_dodge(0.5)) +
#   scale_color_manual(values =colors,name = "age") +
#   theme_minimal() +
#   labs(x = "Model Parameters",
#        y="Population Proportion",
#        color = "Age",
#        shape = "Sex") +
#   scale_x_discrete(labels = ~ str_wrap(.x, 35)) +
#  # scale_y_continuous(limits = c(0, 0.08), breaks = seq(0, 0.08, by = 0.02)) +
#   scale_linetype_discrete(name = "") +
#   guides(color = guide_legend(override.aes = 
#                                 list(shape = 15,size=1,linetype=0,vjust=0.5),
#                               reverse=TRUE),
#          shape = guide_legend(override.aes = list(size=1,linetype=0)))+
#   theme(
#     legend.spacing = unit(1, 'cm'),
#     legend.position="right",
#     axis.text.x=element_text(size=8,angle = 90, hjust = 0),
#     axis.text.y = element_text(angle = 90,hjust=0),
#     legend.spacing.x = unit(0.5, 'cm'),
#     axis.title.x = element_text(angle=180),
#     legend.text = element_text(angle = 90,vjust=0.5,size = 9,margin = margin(t = 10),hjust = -1),
#     legend.title = element_text(angle = 90,vjust=0.5),
#       legend.title.position = "bottom",
#     legend.direction = "vertical",
#     strip.text = element_text(angle = 90)) +
#   facet_grid(~edu_order,
#              labeller = labeller(edu_order = c("1" = "less than\nhigh school",
#                                                "2" = "high school\nor equivalent",
#                                                "3" = "some college",
#                                                "4" = "associate's\ndegree",
#                                                "5" = "bachelor's\ndegree",
#                                                "6" = "professional or\ngraduate degree")),
#              switch ="y")
# 
# ggsave(file = here(paste0("plots/state_level/results_proportion_statelevel_age_education.pdf")),dpi = 320, bg = "white", width =25, height =7, unit = "in")

#sanity check: confirm all estimates are actually different, the differences are just very small!
dat.state %>%
  na.omit %>%
  group_by(sex,age,edu,model_name,specs) %>%
  mutate(age = case_when(age == 1 ~ "18-24 years",
                         age == 2 ~ "25-64 years",
                         age == 3 ~ "65+ years"),
         sex = ifelse(sex==1,"male","female"),
         edu = case_when(edu == 1 ~ "less than\nhigh school",
                         edu == 2 ~ "high school\nor equivalent",
                         edu == 3 ~ "some college",
                         edu == 4 ~ "associate's\ndegree",
                         edu == 5 ~ "bachelor's\ndegree",
                         edu == 6 ~ "professional or\ngraduate degree"),
  ) %>% 
  filter(sex == "female", age == "18-24 years",str_detect(edu,"associate's")) %>% 
  arrange(model_name, sex,age,edu) %>% #filter(xfips == "06001", sex == "female", age == "18-24 years",str_detect(edu, "associate's")) %>% view 
  ggplot(aes(x=factor(model_name), y= state_prop, color = factor(age), shape = factor(sex))) +
  geom_pointrange(aes(ymin = state_proplw, ymax = state_propup),linewidth = 0.5, position = position_dodge(0.6)) + #,position="jitter"
  # geom_hline(data =  cdc_state %>% filter(date == "2021-06-30"),
  #            aes(yintercept = sum(estimate)/total_pop,
  #                linetype = "CDC estimate: June 30, 2021"),
  #            color = "red",
  #            linewidth = 0.8) +
  # geom_hline(data =  dir %>% filter(vax == "Yes"),
  #            aes(yintercept = mean,
  #                linetype = "Direct Survey Estimate"),
  #            color = "darkorange",
  #            linewidth = 0.8) +
  scale_color_manual(values =colors,name = "age") +
  theme_minimal() +
  labs(x = "Model Parameters",
       y="Population Proportion",
       color = "Age",
       shape = "Sex") +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20)) +
  scale_linetype_discrete(name = "") +
  guides(color = guide_legend(override.aes = list(shape = 15,size=2) ))+
  theme(
    legend.spacing = unit(1, 'cm'),
    legend.position="bottom")

ggsave(file = here(paste0("plots/state_level/results_proportion_statelevel_age_education_SANITYCHECK.pdf")),dpi = 320, bg = "white", width =16, height =6, unit = "in")


modelnames = c(unique(dat.state$model_name))

temp = cdc_bysexage_table %>%
  filter(date == "2021-06-30") %>%
  select(sex,age,stratumprop_estimate)

temp$model_name = "CDC Estimate"#rep(modelnames,each=6)

temp2 = read_csv(here("/Users/Aja/Documents/R Directories/NSF_RAPID_COVID19_Survey/Aja/COVID_Behaviors/data/smrp/direct_estimates/ctis_vax_direct_estimates_not_boostrapped.csv"))
temp2$model_name = "CTIS Survey Estimate"


#temp = temp %>% rename(model_name = model)

colors = wesanderson::wes_palette(5, name = "Darjeeling1", type = "discrete")

dat.state = 
  dat %>%
  group_by(sex,age,model_name,specs) %>%
  mutate(
    #N = sum(strata_pop),
    state_sae = sum(median_y),#total number each stratum in the state
    state_saelw = sum(lw_y),
    state_saeup = sum(up_y)) %>%
  distinct() %>% 
  mutate(sexage_state_prop = state_sae/sum(unique(strata_pop)), #total proportion in each stratum in the state (denominator is the sum of all unique values of population counts in the county)
         sexage_state_proplw = state_saelw/sum(unique(strata_pop)),
         sexage_state_propup = state_saeup/sum(unique(strata_pop))) %>% 
  distinct() %>%
  mutate(type = case_when(str_detect(specs, "bym2") ~ "Spatial",
                          str_detect(specs, "iid") ~ "Non-Spatial"),
         specs = str_replace_all(specs,"_"," ")) %>%
  select(sex:edu,specs,type,sexage_state_prop:sexage_state_propup) %>%
  ungroup() %>%
  distinct()

# sm2 = sm2 %>%
#   mutate(model_name = factor(model_name,levels=statelevel_order))

dat.state %>%
  #select(model_name,sex,age,state_sae,state_saelw,state_saeup,strata_pop) %>%
  distinct() %>%
  na.omit() %>%
  group_by(model_name,sex,age) %>%
  mutate(
    model_name = factor(model_name,levels=statelevel_order),
    # specpop = case_when(sex == 1 & age == 1 ~ m_1824_pop,
    #                          sex == 1 & age == 2 ~ m_2564_pop,
    #                          sex == 1 & age == 3 ~ m_65plus_pop,
    #                          sex == 0 & age == 1 ~ f_1824_pop,
    #                          sex == 0 & age == 2 ~ f_2564_pop,
    #                          sex == 0 & age == 3 ~ f_65plus_pop),
         age = case_when(age == 1 ~ "18-24 years",
                         age == 2 ~ "25-64 years",
                         age == 3 ~ "65+ years"),
         sex = ifelse(sex==1,"male","female")) %>%
         # state_prop = sum(state_sae)/sum(unique(strata_pop)), #total proportion in each stratum in the state (denominator is the sum of all unique values of population counts in the county)
         #        state_proplw = sum(state_saelw)/sum(unique(strata_pop)),
         #        state_propup = sum(state_saeup)/sum(unique(strata_pop))) %>% 
         # total_sae  = sum(state_sae),
         # total_saelw  = sum(state_saelw),
         # total_saeup  = sum(state_saeup)) %>%  %>%
  ggplot(aes(x=factor(model_name), y= sexage_state_prop, color = factor(age), shape = factor(sex))) +
  geom_pointrange(aes(ymin = sexage_state_proplw, ymax = sexage_state_propup),position = position_dodge(0.8),size=0.5) + #,position="jitter"
  # geom_pointrange(data =  sm2 %>% 
  #                   select(model_name,age,sex,popest:popest_up) %>% 
  #                   distinct() %>%
  #                   mutate(model_name = factor(model_name)),
  #                 aes(y = popest, ymin = popest_lw, ymax=popest_up),
  #                 position = position_dodge2(0.8),size=0.5) +
  geom_point(data =  temp,
             aes(y = stratumprop_estimate),
             position = position_dodge2(0.8),size=2.5) +
  geom_pointrange(data =  temp2 %>% distinct(),
             aes(y = mean, ymin = mean_lw, ymax=mean_up),
             position = position_dodge2(0.8),size=0.5) +
  scale_color_manual(values =colors,name = "age") +
  theme_minimal() +
  coord_flip() +
  labs(x = "Model Parameters",
       y="Proportion of Population by Age and Sex",
       color = "age",
       shape = "sex") +
  scale_x_discrete(labels = ~ str_wrap(.x, 38)) +
  scale_y_continuous(limits = c(0.5, 1.42), breaks = seq(0.5, 1.5, by = 0.05)) +
  scale_linetype_discrete(name = "") +
  scale_y_break(c(1,1.35),scales=0.13,space = 0.3) +
  guides(color = guide_legend(override.aes = list(linetype=0,shape = 15,size=1)),
         shape= guide_legend(override.aes = list(linetype=0, size=1)))+
  theme(
    legend.spacing = unit(1, 'cm'),
    legend.position="bottom",
    panel.spacing.x =  unit(1, "cm"),
    axis.text.x = element_text(hjust=0.5),
    axis.text = element_text(colour = "black"),
    axis.text.x.top = element_blank())

ggsave(file = here(paste0("plots/state_level/results_proportion_statelevel_age.pdf")),dpi = 320, bg = "white", width =8, height =5, unit = "in")

signif(1112672 / 25296329,2) #0.04
signif(7610234 / 25296329,2) #0.3
signif(4480140  / 25296329,2)#0.18
signif(1044329 / 25296329,2) #0.4
7253978 / 25296329 #.29
3794976 / 25296329 #0.15
3794976 + 7253978 + 1044329 + 4480140 + 7610234 + 1112672 
0.04+0.3+0.18+0.4+0.29+0.15
  colors = wesanderson::wes_palette(3, name = "Darjeeling1", type = "discrete")
colors = colors[2:3]

dat %>%
  group_by(model_name,specs) %>%
  mutate(
    #N = sum(strata_pop),
    state_sae = sum(median_y),#total number each stratum in the state
    state_saelw = sum(lw_y),
    state_saeup = sum(up_y)) %>%
  distinct() %>% 
  mutate(state_prop = state_sae/sum(unique(strata_pop)), #total proportion in each stratum in the state (denominator is the sum of all unique values of population counts in the county)
         state_proplw = state_saelw/sum(unique(strata_pop)),
         state_propup = state_saeup/sum(unique(strata_pop))) %>% 
  distinct() %>%
  mutate(type = case_when(str_detect(specs, "bym2") ~ "Spatial",
                          str_detect(specs, "iid") ~ "Non-Spatial"),
         specs = str_replace_all(specs,"_"," ")) %>%
  select(specs,type,state_prop:state_propup) %>%
  ungroup() %>%
  distinct() %>%
  mutate(model_name = factor(model_name, levels =statelevel_order
                                  # c("fixed effect on age, fixed effect on education, IID by county",
                                  # "age-education interaction, IID by county" ,
                                  # "fixed effect on age, RW1 by education, IID by county",
                                  # "age-education RW1 interaction, IID by county",
                                  # "fixed effect on age, fixed effect on education, BYM2 by county",
                                  # "fixed effect on age, BYM2 by education",
                                  # "fixed effect on age, fixed effect on education, age-education BYM2 interaction",
                                  # "fixed effect on education, BYM2 by age",
                                  # "fixed effect on age, fixed effect on education, education-age BYM2 interaction"
                                  # )
                             )) %>% 
  #arrange(state_proportion_total) %>%
  ungroup() %>% #group_by(model_name,specs) %>%
    ggplot(aes(x=model_name, y= state_prop)) +
    geom_pointrange(aes(ymin =state_proplw, ymax = state_propup)) +
  # geom_pointrange(data=sm %>% 
  #                   mutate(model_name = "state-level age-education interaction") %>%
  #                   distinct(), 
  #                 aes(y = state_est/pop, ymin =state_est_lw/pop, ymax = state_est_up/pop)) +
      geom_hline(data =  cdc_state %>% 
                   filter(date == "2021-06-30") %>%
                   select(total_pop,popest) %>% 
                   distinct(),
                  aes(yintercept = popest/total_pop,
                      linetype = "CDC Estimate: June 30, 2021"),
                 color = "red",
                  linewidth = 0.8) +
      geom_hline(data =  dir %>% filter(vax == "Yes"),
               aes(yintercept = mean,
                   linetype = "CTIS Survey Estimate"),
               color = "darkorange",
               linewidth = 0.8) + #,position="jitter"
    #scale_color_manual(values =colors,name = "type") +
      theme_minimal() +
      guides(color="none") +
      labs(x = "Model Parameters",
           y="Population Proportion",
           color = "type") +
      scale_x_discrete(labels = ~ str_wrap(.x, 50)) +
      scale_y_continuous(limits = c(0.83, 0.915), breaks = seq(0.84, 0.90, by = 0.02)) +
      scale_linetype_discrete(name = "") +
      theme(
            legend.spacing = unit(1, 'cm'),
            legend.position="top",
            legend.direction = "horizontal",
            legend.key.justification = "left",
            legend.text = element_text(size = 8),
            axis.text = element_text(colour = "black")) +
  coord_flip() 

ggsave(file = here(paste0("plots/state_level/results_proportion_statelevel.pdf")),dpi = 320, bg = "white", width =7, height =4, unit = "in")
  

temp = dat %>%
  group_by(model_name,specs) %>%
  mutate(
    #N = sum(strata_pop),
    state_sae = sum(median_y),#total number each stratum in the state
    state_saelw = sum(lw_y),
    state_saeup = sum(up_y)) %>%
  distinct() %>% 
  mutate(state_prop = state_sae/sum(unique(strata_pop)), #total proportion in each stratum in the state (denominator is the sum of all unique values of population counts in the county)
         state_proplw = state_saelw/sum(unique(strata_pop)),
         state_propup = state_saeup/sum(unique(strata_pop))) %>% 
  distinct() %>%
  mutate(type = case_when(str_detect(specs, "bym2") ~ "Spatial",
                          str_detect(specs, "iid") ~ "Non-Spatial"),
         specs = str_replace_all(specs,"_"," ")) %>%
  select(specs,type,state_prop:state_propup) %>%
  ungroup() %>%
  distinct() %>%
  mutate(model_name = factor(model_name, levels =statelevel_order),
         ci = paste0("(",signif(state_proplw,3),", ",signif(state_propup,3),")"),
         state_prop = signif(state_prop,3)) %>%
  select(type,model_name,state_prop,ci)

fwrite(temp,file=here('/Users/Aja/Documents/R Directories/MRP_Spatial/results/tables/table8_supp_state_proportion_sums.csv'))

#direct estimate for plots
dirmean = dir %>% filter(vax == "Yes") %>% select(mean)
cdcest = cdc %>% filter(date == "2021-06-30") %>% reframe(est = sum(estimate))


#### County Results ####
ctdir = fread(here("data/direct_estimates/ctis_vax_direct_estimates_not_bootstrapped_with_edu_simplevax_bycounty.csv"))
ctdir = ctdir %>% filter(vax=="Yes")
#view(dat %>% filter(xfips == "06001"))
dat.county = 
dat %>%
  filter(!is.na(xfips)) %>%
  group_by(xfips,model_name,specs) %>%
  mutate(n_i = sum(unique(sex_specific_n_i))) %>% #sum the two sex estimates together to get total n_i
  group_by(xfips,model_name,specs) %>%
  mutate(
    county_sae = sum(median_y),#total number each stratum in the county
    county_saelw = sum(lw_y),
    county_saeup = sum(up_y)) %>%
  distinct() %>% 
  mutate(county_prop = county_sae/n_i, #total proportion in each stratum in the county (denominator is the sum of all unique values of population counts in the county)
         county_proplw = county_saelw/n_i,
         county_propup = county_saeup/n_i) %>% 
  distinct() %>%
  ungroup() %>%
  mutate(type = case_when(str_detect(specs, "bym2") ~ "Spatial",
                          !str_detect(specs, "bym2") ~ "Non-Spatial"),
         specs = str_replace_all(specs,"_"," ")) %>%
  select(xfips,type,model_name,county_prop:county_propup,specs) %>% #county_sae:county_saeup
  distinct()
# 
# dat.county %>%
#   group_by(model_name,specs,xfips) %>%
#   select(xfips,sex,specs,model_name,county_sae:county_saeup,type,n_i) %>%
#   mutate(county_proportion_total = sum(unique(county_sae))/unique(n_i), 
#          county_proportion_total_lw = sum(unique(county_saelw))/sum(unique(n_i)),
#          county_proportion_total_up = sum(unique(county_saeup))/sum(unique(n_i))) %>%
#   select(xfips,sex,model_name,specs,county_proportion_total:county_proportion_total_up,type) %>%
#   distinct %>%
#   ungroup() %>% glimpse

dat.county$xfips = as.integer(dat.county$xfips)

ctdir = ctdir %>%
    rename(county_prop = mean,
           county_proplw = mean_low,
           county_propup = mean_upp) %>%
    select(-V1,-ci) %>%
    mutate(model_name = "CTIS survey estimate",
           specs = "CTIS survey estimate") %>%
    select(xfips,model_name,specs,county_prop:county_propup)

ctdir1 = ctdir %>%
  mutate(type = "Spatial")
ctdir2 = ctdir %>%
  mutate(type = "Non-Spatial")

ctdir = rbind(ctdir1,ctdir2)
rm(ctdir1,ctdir2)

ctdir = ctdir %>% select(colnames(dat.county))

dat.county = rbind(dat.county,ctdir)

# glimpse(dat.county)
# glimpse(ctdir)
# # 
# # temp = dat %>% select(model_name,xfips) %>% na.omit %>% distinct() %>% mutate(xfips = as.integer(xfips))
# # glimpse(temp)
# #dat.county=left_join(dat.county, temp, by =c("xfips","model_name"),relationship = "many-to-many")
# 
# ctdir=left_join(ctdir, temp, by ="xfips",relationship = "many-to-many")
# glimpse(ctdir)
# dat.county = rbind(dat.county,ctdir)
# glimpse(dat.county)
# 
# spatial_colors = c("#7ad151","darkorange",
#                    "#22a884",
#                    "#2a788e",
#                    "#414487","#440154")
# 
# nonspatial_colors = c("darkorange","#7ad151",
#                       "#22a884",
#                       "#2a788e",
#                       "#414487","#440154")

#scales::show_col(viridis(5))

# for (i in types) {
# {
#   
#   lab = str_to_lower(types) # Converts the entire string to lower case
#   lab = str_replace(lab,' ', '_') 
#   lab = str_replace(lab,': ', '_') 
#   lab = str_replace(lab, " ","")
#   lab = str_replace(lab, "-","")
# # ?relevel
#  t2 =  dat.county %>%
#     mutate(specs = factor(specs),
#            specs = ifelse(type == "Spatial",
#                      forcats::fct_relevel(specs,"direct survey estimate",
#                                             "bym2 age",
#                                             "bym2 edu",
#                                             "rw1 age bym2 edu",
#                                             "rw1 edu bym2 age",
#                                             "rw2 edu bym2 age"),
#                      forcats::fct_relevel(specs,"direct survey estimate",
#                                                    "rw1 age iid xfips",
#                                                    "rw1 edu iid xfips",
#                                                    "rw2 edu iid xfips",
#                                                    "rw1 age rw1 edu iid xfips",
#                                                    "rw1 age rw2 edu iid xfips")))
#   glimpse(t2)


cdccounty = countyvax %>%
  filter(date == "06/30/2021") %>%
  rename(xfips = fips,
         name = county) %>%
  mutate(xfips = as.integer(xfips),
         model_name = "CDC reported count",
         proplw = prop, #for viz
         propup = prop,
         specs = "cdc reported count")

cdccounty2 = cdccounty
cdccounty$type = "Spatial"
cdccounty2$type = "Non-Spatial"
cdccounty = rbind(cdccounty,cdccounty2)
rm(cdccounty2)

cdccounty = cdccounty %>%
  rename(county_prop=prop,
         county_proplw = proplw,
         county_propup = propup)
cdccounty = cdccounty %>% select(name,colnames(dat.county))

temp = countyvax %>% select(county,fips) %>% rename(name = county) %>% mutate(fips = as.integer(fips))
dat.county = left_join(dat.county,temp, by=c("xfips"="fips"),relationship="many-to-many")
dat.county = dat.county %>% select(colnames(cdccounty))

dat.county =rbind(dat.county,cdccounty)
rm(temp)
gc()

dat.county = dat.county %>%
  mutate(specs = str_replace_all(specs,"xfips","county"))

nonspatial_countylevel_order = c("CDC reported count",
                              #  "CTIS survey estimate",
                                paste0(models[1,1]),
                                paste0(models[2,1]),
                                paste0(models[3,1]),
                                paste0(models[4,1])
                                )


dat.county = dat.county %>% distinct()


dat.county %>%
    filter(
      model_name != "CTIS survey estimate",
      type == "Non-Spatial") %>%
  distinct() %>% 
    ggplot(aes(x=name, y=county_prop, color=factor(model_name, levels = nonspatial_countylevel_order))) +
  geom_pointrange(data = . %>% filter(model_name == "CDC reported count"),
                  aes(ymin = county_proplw, ymax = county_propup)) +
  geom_pointrange(
    data = . %>% filter(model_name != "CDC reported count"),
    aes(ymin = county_proplw, ymax = county_propup),
    position = position_dodge2(0.8)) +
    theme_minimal() +
    labs(x = "County",
         y="Population Proportion",
         color = "Model") +
    scale_color_manual(values = c("red",
                                  #"darkorange",
                                  "#6DCD59FF",
                                  "#2a788e",
                                  "navy",
                                  "darkviolet"),
                       labels = label_wrap_gen(32)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    #scale_linetype_discrete(name = "") +
    guides(colour = guide_legend(override.aes = list(size=1,linetype=0),nrow=1)) +
    theme(text = element_text(size=12),
          title= element_text(size=14),
          legend.spacing = unit(1.5, 'cm'),
          legend.position="bottom",
          axis.text.x=element_text(size=12,angle = 90,hjust=1),
          axis.text.y = element_text(angle = 90,hjust=0),
          legend.spacing.x = unit(1, 'cm'),
         # axis.title.x = element_text(angle=180),
         legend.text = element_text(size =12),
         legend.title = element_text(size = 14),
          legend.title.position = "left",
          legend.direction = "horizontal",
          strip.text = element_blank()) 
    #scale_color_viridis(discrete = TRUE) +
  
  ggsave(file = here::here(paste0("plots/county_level/results_proportion_countylevel_nonspatial.pdf")),dpi = 320, bg = "white", width =23, height =8, unit = "in")

`%nin%` = Negate(`%in%`)
spatial_countylevel_order = c("CDC reported count",
                            #  "CTIS survey estimate",
                              statelevel_order[7:11])

#dat.county %>% filter(specs %nin% spatial_countylevel_order) %>% distinct(specs)
dat.county %>%
    filter(
      model_name != "CTIS survey estimate",
      type == "Spatial"
    ) %>%
    na.omit() %>%
    distinct() %>%
    ggplot(aes(x=name, y=county_prop,
               color=factor(model_name,levels=spatial_countylevel_order))) +
    geom_pointrange(data = . %>% filter(model_name == "CDC reported count"),
                                         aes(ymin = county_proplw, ymax = county_propup)) +
    geom_pointrange(
      data = . %>% filter(model_name != "CDC reported count"),
                    aes(ymin = county_proplw, ymax = county_propup),
                    position = position_dodge2(0.8)) +
    theme_minimal() +
    labs(x = "County",
         y="Population Proportion",
         color = "Model") +
  scale_color_manual(values = c("red",
                              #  "darkorange",
                                "darkorchid1",
                                "#440154",
                                "forestgreen",
                                "aquamarine2",
                                "deepskyblue2"),
                     labels = label_wrap_gen(35)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    #scale_linetype_discrete(name = "") +
  guides(colour = guide_legend(override.aes = list(size=1,linetype=0),nrow=1)) +
            theme(text = element_text(size=12),
                  title= element_text(size=14),
                  legend.spacing = unit(1.5, 'cm'),
                  legend.position="bottom",
                  axis.text.x=element_text(size=12,angle = 90,hjust=1),
                  axis.text.y = element_text(angle = 90,hjust=0),
                  legend.spacing.x = unit(1, 'cm'),
                  # axis.title.x = element_text(angle=180),
                  legend.text = element_text(size =12),
                  legend.title = element_text(size = 14),
                  legend.title.position = "left",
                  legend.direction = "horizontal",
                  strip.text = element_blank()) #+
  #scale_color_viridis(discrete = TRUE) +
  #facet_wrap(~type, nrow = 2)

  ggsave(file = here::here(paste0("plots/county_level/results_proportion_countylevel_spatial.pdf")),dpi = 320, bg = "white", width =23, height =8, unit = "in")

# joined_order = c(spatial_countylevel_order,nonspatial_countylevel_order)
# 
# joined_order = unique(joined_order)
# 
# dat.county %>%
#     filter(
#       !str_detect(specs,"bym2 county")
#     )%>%
#   na.omit() %>%
#     ggplot(aes(x=name, y=county_prop, 
#                color=factor(model_name,levels=joined_order),
#                group=type)) +
#     geom_pointrange(aes(ymin = county_proplw, ymax = county_propup),
#                     position = position_jitter(width = 0.3, height = NULL, seed = 2)) +
#     theme_minimal() +
#     labs(x = "County",
#          y="Proportion",
#          color = "Model") +
#     scale_color_manual(values = c("red",
#                                   "darkorange",
#                                   "#fde725",
#                                   "#b5de2b",
#                                   "#6ece58",
#                                   "#35b779",
#                                   "#1f9e89",
#                                   "#26828e",
#                                   "#31688e",
#                                   "#3e4989",
#                                   "#482878",
#                                   "#440154")) +
#     #scale_y_continuous(limits = c(0.45, 1), breaks = seq(0, 1, by = 0.1)) +
#     #scale_linetype_discrete(name = "") +
#     theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 45),
#           legend.spacing = unit(1, 'cm'),
#           legend.position="bottom") + #+
#   #scale_color_viridis(discrete = TRUE) +
#   facet_wrap(~type, nrow = 2)
#   
#     
# ggsave(file = here(paste0("plots/county_level/results_proportion_countylevel.pdf")),dpi = 320, bg = "white", width =20, height =14, unit = "in")

#save tables so that they can go into the paper
#tidy them up a little
dat.county = 
  dat.county %>%
  distinct() %>%
  group_by(xfips,type,model_name) %>%
  mutate(county_prop = signif(sum(county_prop),3),
         county_proplw = signif(sum(county_proplw),3),
         county_propup = signif(sum(county_propup),3)) %>%
  ungroup() %>%
  arrange(name,type,model_name) %>%
  select(type,name,model_name,county_prop:county_propup) %>%
  mutate(county_proplw = ifelse(model_name=="CDC reported count", NA,county_proplw),
         county_propup = ifelse(model_name=="CDC reported count", NA,county_propup)) %>%
  distinct()

#clean up duplicates needed for viz
dat.county = 
  dat.county %>%
  mutate(type = ifelse(model_name == "CDC reported count"|model_name=="CTIS survey estimate","Spatial",type)) %>%
  distinct() %>%
  mutate(type = ifelse(model_name == "CDC reported count"|model_name=="CTIS survey estimate","",type),
         county_proplw = ifelse(is.na(county_proplw),"-",as.character(county_proplw)),
         county_propup = ifelse(is.na(county_propup),"-",as.character(county_propup))) %>%
  arrange(name,type,model_name)

dat.state = dat.state %>%
  select(-specs) %>%
  group_by(type,model_name) %>%
  mutate(state_proportion_total = signif(sum(sexage_state_prop),3),
         state_proportion_total_lw = signif(sum(sexage_state_proplw),3),
         state_proportion_total_up = signif(sum(sexage_state_propup),3)) %>%
  ungroup() %>%
  arrange(type,model_name)  %>%
  select(type,model_name,state_proportion_total:state_proportion_total_up) %>%
  mutate(state_proportion_total_lw = ifelse(model_name=="CDC reported count", NA,state_proportion_total_lw),
         state_proportion_total_up = ifelse(model_name=="CDC reported count", NA,state_proportion_total_up)) %>%
  distinct()

table(dat.state$model_name)
write.csv(dat.county,here("results/tables/table6_county_regression_results.csv"))


#### Make maps ####
dat.county = 
dat.county %>%
  mutate(county_proplw = as.double(ifelse(county_proplw=="-", NA,county_proplw)),
         county_propup = as.double(ifelse(county_proplw=="-", NA,county_propup)))

geo = counties(state = "CA", cb = FALSE, year = 2020)
geo = 
  geo %>% 
  clean_names() %>% 
  select(geoid,name) 


dat.county = left_join(dat.county,geo, by = "name")
dat.county$geometry = st_simplify(dat.county$geometry)
table(dat.county$model_name)

dat.county=dat.county %>% mutate(type = ifelse(type=="","Validation Data",type),
                                 type = ifelse(is.na(type),"Non-Spatial",type))


maps = 
  dat.county %>%
  group_by(type,model_name) %>%
  distinct() %>%
  # filter(
  #   #type != "Spatial",
  #        model_name != "CTIS survey estimate",
  #        model_name != "CDC reported count") %>%
  pivot_longer(cols = county_prop:county_propup, names_to = 'estimates', values_to = 'prop') %>% 
  #filter(model_name == "fixed effect on age, fixed effect on education, IID by county") %>%
  mutate(estimates = factor(estimates, levels=c("county_proplw","county_prop","county_propup")),
         type = factor(type,levels=c("Validation Data","Non-Spatial","Spatial")),
         model_name = factor(model_name, levels = c(statelevel_order
                                                    # "CDC reported count",
                                                    # "CTIS survey estimate",
                                                    # "age-education interaction, IID by county",
                                                    # "fixed effect on age, fixed effect on education, IID by county",
                                                    # "age-education RW1 interaction, IID by county",
                                                    # "fixed effect on age, RW1 by education, IID by county",
                                                    # "fixed effect on age, fixed effect on education, BYM2 by county",
                                                    # "fixed effect on age, BYM2 by education",
                                                    # "fixed effect on age, fixed effect on education, age-education BYM2 interaction"
                                                    ))) %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(data = . %>% filter(prop < 1), aes(fill=prop)) +
  scale_fill_viridis(na.value = "grey",labels = ~ str_wrap(.x, 15),
                     limits = c(0, 1)) +
  labs(fill = "proportion") +
  new_scale("fill") +
  geom_sf(data = . %>% filter(prop > 1), aes(fill=prop))+#,aes(fill="lightpink")) +
  scale_fill_viridis(option = "magma") +
  labs(fill = "proportion > 1") +
  facet_grid(type + model_name~estimates,#factor(type,levels=c("Validation Data","Non-Spatial","Spatial")),
             switch = "y",
             labeller = labeller(estimates = c("county_prop" = "estimate",
                                             "county_proplw" = "95%CI: lower",
                                             "county_propup" = "95%CI: upper"),
                                model_name = label_wrap_gen(30))) +
  theme_pubclean() +
  theme(axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill= NA),
        panel.background = element_rect(fill= NA,color="black"),
        legend.background = element_rect(colour = NA),
        legend.title.position = "top"
        ) 

ggsave(maps,file=here('plots/maps/model_maps.pdf'),width = 8.5, height = 15, units = "in",dpi = 350)

maps = 
  dat.county %>%
  group_by(type,model_name) %>%
  filter(type != "Non-Spatial") %>%
  pivot_longer(cols = county_prop:county_propup, names_to = 'estimates', values_to = 'prop') %>% 
  #filter(model_name == "fixed effect on age, fixed effect on education, IID by county") %>%
  mutate(estimates = factor(estimates, levels=c("county_proplw","county_prop","county_propup")),
         model_name = factor(model_name, levels = c(statelevel_order
                                                    # "CDC reported count",
                                                    # "CTIS survey estimate",
                                                    # "fixed effect on age, BYM2 by education",
                                                    # "fixed effect on age, fixed effect on education, BYM2 by county",
                                                    # "fixed effect on age, fixed effect on education, age-education BYM2 interaction"
         ))) %>%
  st_as_sf() %>%
  ggplot() +
  geom_sf(data = . %>% filter(prop < 1), aes(fill=prop)) +
  scale_fill_viridis(na.value = "grey",labels = ~ str_wrap(.x, 15),
                     limits = c(0, 1), breaks = c(0,1)) +
  labs(fill = "proportion") +
  new_scale("fill") +
  geom_sf(data = . %>% filter(prop > 1), aes(fill=prop))+#,aes(fill="lightpink")) +
  scale_fill_viridis(option = "magma") +
  labs(fill = "proportion > 1") +
  facet_grid(model_name~estimates,
             switch = "y",
             labeller = labeller(estimates = c("county_prop" = "median estimate",
                                               "county_proplw" = "95%CI: lower",
                                               "county_propup" = "95%CI: upper"),
                                 model_name = label_wrap_gen(25))) +
  theme_pubclean() +
  theme(axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill= NA),
        panel.background = element_rect(fill= NA,color="black"),
        legend.background = element_rect(colour = NA),
        legend.title.position = "top"
  ) 

ggsave(file=here('plots/maps/spatial_model_maps.pdf'),width = 10, height = 20, units = "in",dpi = 500)



maps = dat.county %>%
  group_by(type,model_name) %>%
  filter(type != "Spatial") %>%
  pivot_longer(cols = county_prop:county_propup, names_to = 'estimates', values_to = 'prop') %>% 
  #filter(model_name == "fixed effect on age, fixed effect on education, IID by county") %>%
  mutate(estimates = factor(estimates, levels=c("county_proplw","county_prop","county_propup")),
         model_name = factor(model_name, levels = c("CDC reported count",
                                                     "CTIS survey estimate",
                                                     "age-education interaction, IID by county",
                                                     "fixed effect on age, fixed effect on education, IID by county",
                                                     "age-education RW1 interaction, IID by county",
                                                     "fixed effect on age, RW1 by education, IID by county"
         ))) %>% 
  st_as_sf() %>%
  ggplot() +
  geom_sf(data = . %>% filter(prop < 1), aes(fill=prop)) +
  scale_fill_viridis(na.value = "grey",labels = ~ str_wrap(.x, 15),
                     limits = c(0.8, 0.85), breaks = c(0.8,0.85, by = 0.01)) +
  labs(fill = "proportion") +
  new_scale("fill") +
  geom_sf(data = . %>% filter(prop > 1), aes(fill=prop))+#,aes(fill="lightpink")) +
  scale_fill_viridis(option = "magma") +
  labs(fill = "proportion > 1") +
  facet_grid(model_name~estimates,
             switch = "y",
             labeller = labeller(estimates = c("county_prop" = "median estimate",
                                               "county_proplw" = "95%CI: lower",
                                               "county_propup" = "95%CI: upper"),
                                 model_name = label_wrap_gen(25))) +
  theme_pubclean() +
  theme(axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill= NA),
        panel.background = element_rect(fill= NA,color="black"),
        legend.background = element_rect(colour = NA),
        legend.title.position = "top"
  ) 

ggsave(file=here('plots/maps/nonspatial_model_maps_SANITYCHECK.pdf'),width = 10, height = 16, units = "in",dpi = 500)

maps = 
  dat.county %>%
  group_by(type,model_name) %>%
  filter(type != "Non-Spatial",
         model_name!="fixed effect on age, fixed effect on education, BYM2 by county",
         model_name !="CDC reported count",
         model_name != "CTIS survey estimate") %>%
  pivot_longer(cols = county_prop:county_propup, names_to = 'estimates', values_to = 'prop') %>% 
  #filter(model_name == "fixed effect on age, fixed effect on education, IID by county") %>%
  mutate(estimates = factor(estimates, levels=c("county_proplw","county_prop","county_propup")),
         model_name = factor(model_name, levels = c("CDC reported count",
                                                    "CTIS survey estimate",
                                                    "fixed effect on age, BYM2 by education",
                                                    "fixed effect on age, fixed effect on education, BYM2 by county",
                                                    "fixed effect on age, fixed effect on education, age-education BYM2 interaction"
         ))) %>% 
  st_as_sf() %>%
  ggplot() +
  geom_sf(data = . %>% filter(prop < 1), aes(fill=prop)) +
  scale_fill_viridis(na.value = "grey",labels = ~ str_wrap(.x, 15))+#,
  #limits = c(0.6, 0.8), breaks = c(0.6,0.8, by = 0.1)) +
  labs(fill = "proportion") +
  new_scale("fill") +
  geom_sf(data = . %>% filter(prop > 1), aes(fill=prop))+#,aes(fill="lightpink")) +
  scale_fill_viridis(option = "magma") +
  labs(fill = "proportion > 1") +
  facet_grid(model_name~estimates,
             switch = "y",
             labeller = labeller(estimates = c("county_prop" = "median estimate",
                                               "county_proplw" = "95%CI: lower",
                                               "county_propup" = "95%CI: upper"),
                                 model_name = label_wrap_gen(25))) +
  theme_pubclean() +
  theme(axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text.y.left = element_text(angle = 0),
        axis.ticks = element_blank(),
        strip.background = element_rect(fill= NA),
        panel.background = element_rect(fill= NA,color="black"),
        legend.background = element_rect(colour = NA),
        legend.title.position = "top"
  ) 

ggsave(file=here('plots/maps/spatial_model_maps_SANITYCHECK.pdf'),width = 10, height = 20, units = "in",dpi = 500)



sanity = dat.county %>%
  group_by(type,model_name) %>%
  st_drop_geometry() %>%
  select(-geometry)
  #filter(type != "Non-Spatial") %>%
  
sanity = sanity %>%
  pivot_longer(cols = county_prop:county_propup, names_to = 'estimates', values_to = 'prop') %>% 
  #filter(model_name == "fixed effect on age, fixed effect on education, IID by county") %>%
  mutate(estimates = factor(estimates, levels=c("county_proplw","county_prop","county_propup")),
         model_name = factor(model_name, levels = c("CDC reported count",
                                                    "CTIS survey estimate",
                                                    "fixed effect on age, BYM2 by education",
                                                    "fixed effect on age, fixed effects on education, BYM2 by county",
                                                    "fixed effect on age, age-education BYM2 interaction",
                                                    "fixed effect on age, fixed effect on education, IID by county",
                                                    "age-education interaction, IID by county",
                                                    "fixed effect on age, RW1 by education, IID by county"
         )))

sanity %>%
  filter(str_detect(model_name,"BYM2"))


