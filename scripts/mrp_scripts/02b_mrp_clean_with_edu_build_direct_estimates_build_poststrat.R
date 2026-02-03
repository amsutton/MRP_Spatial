#### MRP + Spatial Priors using CTIS Data ####

#Script 2: Clean data, build survey direct estimates, and poststratification data 

#Build a series of models representing vaccination probability (min. 1 dose)
#For June 12, 2021 using the COVID-19 Trends and Impacts Survey (CTIS)

#Post-stratification of age/sex/education at county level

pacman::p_load(tidyverse, tidycensus, haven, here, janitor, srvyr, sf, spdep,data.table)


options(tigris_use_cache = TRUE)

tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)

here::i_am("scripts/mrp_scripts/02b_mrp_clean_with_edu_build_direct_estimates_build_poststrat.R")


#### set geographic and model string specifications for intuitive labeling ####
state_choice = "CA"
state = "California"
whole_time_period = FALSE


#### read/clean data: ctis 3, wave 11, june 2021 ####
ctis = fread(here("data/ctis_data_do_not_share/ctis_w11_2021-06.csv"),colClasses = c("fips"="character"))


ctis = ctis %>%
  filter(!is.na(V1)) %>%
  mutate(
    xfips = fips,
    xzip = A3,
    sex = D1,
    age = D2,
    edu = D8,
    mask = C14a,
    vax = V1) %>%
  select(StartDatetime, EndDatetime, weight, xzip,xfips, sex, age, edu, vax) %>%
  na.omit() %>%
  # mutate(state=NA) %>%  # we don't have state, but we do have zip code and 
  # fips code which we use to get state
  separate(StartDatetime, c("start_date", "start_time"), sep = " ") %>%
  separate(EndDatetime, c("end_date", "end_time"), sep = " ")

# if(whole_time_period == FALSE){
#   #only May 22 to June 12 to match the CDC data
#   ctis = ctis %>%
#     filter(start_date < as.Date("2021-06-13")
#   
# }

#make it binary for binomial models
#set 3 == i don't know to NA (nb:1=yes, 2=no)
ctis = ctis %>%
  mutate(vax= case_when(vax==1 ~ 1,
                        vax==2 ~ 0,
                        vax==3 ~ NA_real_))


#identify state by county fips codes; 
#more information to join geographically on later if needed
data("fips_codes")

fips_codes = fips_codes %>% 
  filter(state == state_choice) %>%
  mutate(xfips = paste0(state_code,county_code, sep=""))

#### match ctis data categories to ACS bins ####

#filter sex to ONLY male or female (male = 1, female = 2) 
#because eventually we'll use MRP w/ the ACS data as post-strat matrix
#and ACS only counts male and female, and they must match.
ctis = ctis %>%
  filter(sex == 1 | sex == 2) %>%
  mutate(sex = ifelse(sex == 1, 1, 0))


#to match the ACS categories listed later in this script under object "acs"
ctis = ctis %>%
  mutate(age = case_when(
    age == 1 ~ 1, #18-24
    age == 2 ~ 2, #25-34
    age == 3 ~ 3, #35-44
    age == 4 ~ 4, #45-54
    age == 5 ~ 4, #55-64
    age == 6 ~ 5, #65-74
    age == 7 ~ 5  #75+
  ))


#to match the ACS And the basic assumptions of the CDC baselines in file:
#01b_establish_baseline_proportions_cadoh_data.R

# ctis = ctis %>%
#   mutate(age = case_when(
#     age == 1 ~ 1, #18-24
#     age == 2 ~ 2, #25-34
#     age == 3 ~ 2, #35-44
#     age == 4 ~ 2, #45-54
#     age == 5 ~ 2, #55-64
#     age == 6 ~ 2, #65-74
#     age == 7 ~ 3  #75+
#   ))


#to match acs education values
ctis = ctis %>%
  mutate(edu = ifelse(edu==6 | edu == 7 | edu == 8, 6, edu))


#filter ctis such that only the fips codes present in fips_codes are left; 
#in this case, we are approximating California
dat = ctis %>% filter(str_starts(xfips, "06")) #I checked with the tigris object filter and it's still 57 counties; missing one!

#sanity check: it is a reasonable number of FIPS codes for the context
length(unique(dat$xfips)) # == missing a county, 57 out of 58

#### get county geographies using tigris ####
#county are only available by state for 2000 and 2010
geo = tigris::counties(cb=FALSE, year= 2010, state= state_choice) %>% clean_names()

dat = dat %>%
  select(weight, xfips, vax, sex, age, edu) 

#quick check: looks good
dat %>% group_by(vax) %>% summarize(n())

#### Build Post-stratification Matrix using ACS ####

#sex, age, education estimates for state of California
acs = get_acs(geography = "county",
               state = "CA",
               table = "B15001",
               survey = "acs5",
               geometry = FALSE,
               year = 2021,
               cache = TRUE)

vars = load_variables(2021, "acs5", cache=TRUE)
vars = vars %>% 
  filter(str_detect(name,"B15001"))
        #str_detect(label, "18 to 24 years:|25 to 34 years:|35 to 44 years:|45 to 64 years:|65 years and over:")) 

vars = vars %>%
  mutate(varstemp = str_remove_all(label,"Estimate!!Total:!!"))
demogs = str_split(vars$varstemp,":!!",simplify=TRUE) %>% as.data.frame()
demogs = demogs %>%
  rename(sex = V1,
         age = V2,
         edu = V3)
vars = cbind(demogs,vars)
rm(demogs)

vars = 
  vars %>% 
  select(-varstemp) %>%
  filter(edu != "",
         age != "")

#education using original table
vars =
  vars %>% 
  mutate(edu = case_when(str_detect(edu, "Less than 9th") ~ 1, #less than high school
                         str_detect(edu, "9th") ~ 1, #less than high school
                         str_detect(edu, "High school") ~ 2, 
                         str_detect(edu, "Some college, no degree") ~ 3, 
                         str_detect(edu, "Associate's degree") ~ 4,
                         str_detect(edu, "Bachelor's degree") ~ 5,
                         str_detect(edu, "professional degree") ~ 6)) #professional or graduate degree

acs = left_join(vars, acs,  by=c("name"="variable"))

# acs = acs %>%
#   mutate(sex = ifelse(str_detect(label, "Male"), 1, 0))  # 1= male, 0= female, like in ctis

# acs = acs %>%
#   mutate(age = case_when(str_detect(label, "18|20|21|22|24") ~ 1,
#                          str_detect(label, "29|30|35|40|45|50|55|60|62|65|67|70") ~ 2,
#                          str_detect(label, "75|80|85") ~ 3
#                          )) %>%
#       filter(!str_detect(concept, "ALONE|HISPANIC|RACES"))

acs = 
  acs %>% 
  mutate(age = case_when(str_detect(age, "18 to 24") ~ 1,
                         str_detect(age, "25 to 34") ~ 2,
                         str_detect(age, "35 to 44") ~ 3,
                         str_detect(age, "45 to 64") ~ 4,
                         str_detect(age, "65") ~ 5
                         
  )) %>%
  filter(!str_detect(concept, "ALONE|HISPANIC|RACES"),
         !is.na(age))

# 
# acs = acs %>%        
#   mutate(age = case_when(str_detect(label, "18") ~ 1,
#                          str_detect(label, "25|35|45") ~ 2,
#                          str_detect(label, "65") ~ 3 # a concession; doesn't match previous groups that are 75+
#   )) 

#this is the post-strat table; it has multiples in some categories 
#which we handle later when we build post and pred_strat
acs = 
  acs %>%
  clean_names() %>%
  rename(xfips = geoid) %>%
  select(xfips, sex, age, edu, estimate) %>%
  group_by(xfips, sex,age, edu) %>%
  summarize(estimate=sum(estimate)) %>%
  ungroup() 

post = acs


#### Build Adjacency Matrix ####

#compute adjacency matrix and make sure both column and row names correspond 
#to the xfips names
geo$xfips = as.factor(geo$geoid)

#we must assign the xfips an id 1:nrow(geo) because the adjacency matrix only works
#with integers 1:n -- and xfips aren't like that! NB: the 1:n order is artificial
#and needn't relate to the order of the xfips; it's simply important to have them
#in this format.
geo$id = 1:nrow(geo) 

#make an index of id-to-xfips so we can refer to and use it later
index = geo %>% 
  # select(id, zcta) %>% 
  mutate(xfips = geoid10) %>%
  select(id, xfips) %>%
  mutate(id = as.integer(factor(id)))

save(geo, file=here("data/cleaned_input_data/geo_index.rda"))

#tidy
dat = dat %>%
  select(xfips,vax,sex,age,edu) %>%
  mutate(age = as.integer(age))

dat = dat %>% na.omit()

# #write graph, acs, dat and post object to rda
save(geo, file=here("data/cleaned_input_data/ca_geo_sf_county.rda"))
save(dat, file = here("data/cleaned_input_data/clean_data_county_with_edu.rda"))
save(post, file = here("data/cleaned_input_data/clean_postrat_age_sex_county_with_edu.rda"))

