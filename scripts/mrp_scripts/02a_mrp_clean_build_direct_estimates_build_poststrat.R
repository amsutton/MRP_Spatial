#### MRP + Space using CTIS Data ####

#Script 02a: Clean data, build survey direct estimates, and poststratification data 

#Build a series of models representing vaccination probability (min. 1 dose)
#For June 12, 2021 using the COVID-19 Trends and Impacts Survey (CTIS)

#Post-stratification of age/sex/education at county level

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, tidycensus, haven, here, janitor, crosstable,
               srvyr, sf, spdep, data.table,'ipumsr')


options(tigris_use_cache = TRUE)

tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)

here::i_am("scripts/mrp_scripts/02a_mrp_clean_build_direct_estimates_build_poststrat.R")


#### set geographic and model string specifications for intuitive labeling ####
state_choice <- "CA"
state <- "California"
whole_time_period <- FALSE


 
#### read/clean survey data: ctis 3, wave 11, june 2021 ####
ctis = fread(here("data/ctis_data_do_not_share/ctis_w11_2021-06.csv"),colClasses = c("fips"="character"))

ctis <- ctis %>%
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
  mutate(state=NA) %>%  # we don't have state, but we do have zip code and fips code which we use to get state
  separate(StartDatetime, c("start_date", "start_time"), sep = " ") %>%
  separate(EndDatetime, c("end_date", "end_time"), sep = " ")

# if(whole_time_period == FALSE){
# #only May 22 to June 12 to match the CDC data
# ctis <- ctis %>%
#   filter(start_date < as.Date("2021-06-13"))
# 
# }

#make it binary for binomial models -- these are my outcomes!
#set 3 == i don't know to NA (nb:1=yes, 2=no)
ctis <- ctis %>%
  mutate(vax= case_when(vax==1 ~ 1,
                        vax==2 ~ 0,
                        vax==3 ~ NA_real_))


#identify the state using county fips codes; 
#for my purposes, this provides more information 
#to join geographically on later

data("fips_codes")

fips_codes <- fips_codes %>% 
  filter(state == state_choice) %>%
  mutate(xfips = paste0(state_code,county_code, sep=""))

#### match ctis data categories to ACS bins ####

#this is critical because your survey data and your census (or ACS, here) data
#need to have the same categories and same values within each -- they have to match
#to do MRP

#filter sex to ONLY male or female (male = 1, female = 2) 
#because eventually we'll use MRP w/ the ACS data as post-strat matrix
#and ACS only counts male and female, and they must match.
ctis <- ctis %>%
  filter(sex == 1 | sex == 2) %>%
  mutate(sex = ifelse(sex == 1, 1, 0))


#to match the ACS age categories
ctis <- ctis %>%
  mutate(age = case_when(
    age == 1 ~ 1, #18-24
    age == 2 ~ 2, #25-34
    age == 3 ~ 3, #35-44
    age == 4 ~ 4, #45-54
    age == 5 ~ 4, #55-64
    age == 6 ~ 5, #65-74
    age == 7 ~ 5  #75+
  ))


# ctis <- ctis %>%
#   mutate(age2 = case_when(
#     age == 1 ~ 1, #18-24
#     age == 2 ~ 2, #25-34
#     age == 3 ~ 2, #35-44
#     age == 4 ~ 2, #45-54
#     age == 5 ~ 2, #55-64
#     age == 6 ~ 2, #65-74
#     age == 7 ~ 3  #75+
#   ))

#to match acs education values
ctis <- ctis %>%
  mutate(edu = ifelse(edu==6 | edu == 7 | edu == 8, 6, edu))


#filter ctis such that only the fips codes present in fips_codes are left;
#in this case, we are approximating California
dat <- ctis %>%
  filter(xfips %in% fips_codes$xfips)

#number of respondents in California in June 2021 for whom we have all the data:
dat %>% select(-state) %>% na.omit() %>% nrow()

#sanity check: it is a reasonable number of ZIP codes for the context
length(unique(dat$xzip))

#### get county geographies using tigris ####
#county are only available by state for 2000 and 2010
geo <- tigris::counties(cb=FALSE, year= 2010, state= state_choice) %>% clean_names()

dat <- dat %>%
  select(weight, xfips, vax, sex, age, edu) %>%
  filter(xfips %in% as.list(geo$geoid10)) #counties in CA counties == CA

#quick check: looks good
dat %>% group_by(vax) %>% summarize(n())

#### assure geo object is filtered to area of interest using crswlk file ####

# crswlk <- crswlk %>% filter(state == "CA")
# geo <- geo %>% filter(zcta %in% as.list(crswlk$zcta))
#
# #much closer now that we can compare zcta from dat to zcta from geo (ie. survey vs. census data)
# length(unique(dat$zcta))
# length(unique(geo$zcta))
#
# (length(unique(dat$zcta))/length(unique(geo$zcta)))*100 #==81.9% so not great but tolerable

#### Build Post-stratification Table using ACS ####

#sex by age estimates for state of California
acs <- get_acs(geography = "county",
                state = "CA",
                table = "B01001",
                survey = "acs5",
                geometry = FALSE,
                year = 2021,
                cache_table = TRUE)

vars <- load_variables(2021, "acs5", cache=TRUE)

vars <- vars %>%
  filter(str_detect(label,"Estimate!!Total:!!"),
        str_detect(name, "B01001"),
         str_detect(label, "18|20|21|22|24|29|30|35|40|45|50|55|60|62|65|67|70|75|80|85"))

vars = vars %>%
  mutate(varstemp = str_remove_all(label,"Estimate!!Total:!!"))
demogs = str_split(vars$varstemp,":!!",simplify=TRUE) %>% as.data.frame()
demogs = demogs %>%
  rename(sex = V1,
         age = V2)
vars = cbind(demogs,vars)
rm(demogs)

vars = vars %>% select(-varstemp)

#education using original table
#vars <- vars %>%
#   mutate(edu = case_when(str_detect(label, "than 9th") ~ 1, #less than high school
#                          str_detect(label, "no diploma") ~ 1, #less than high school
#                          str_detect(label, "(includes equivalency)") ~ 2, #high school or equivalent
#                          str_detect(label, "no degree") ~ 3, #some college
#                          str_detect(label, "years:!!Associate's degree") ~ 4,
#                          str_detect(label, "years:!!Bachelor's degree") ~ 5,
#                          str_detect(label, "professional degree") ~ 6)) %>% #professional or graduate degree
#   na.omit() #to drop the variables we don't want

acs <- left_join(vars, acs,  by=c("name"="variable"))

# acs <- acs %>%
#   mutate(sex = ifelse(str_detect(label, "Male"), 1, 0))  # 1= male, 0= female, like in ctis

# acs <- acs %>%
#   mutate(age = case_when(str_detect(label, "18|20|21|22|24") ~ 1,
#                          str_detect(label, "29|30|35|40|45|50|55|60|62|65|67|70") ~ 2,
#                          str_detect(label, "75|80|85") ~ 3
#                          )) %>%
#       filter(!str_detect(concept, "ALONE|HISPANIC|RACES"))

keep_ages = c("18 and 19 years","20 to 24 years","25 to 29 years",
              "30 to 34 years","35 to 44 years","45 to 54 years",
              "55 to 64 years","65 to 74 years","75 to 79 years",
              "85 years and over")

acs <- 
  acs %>% #18-24, 25-34, 35-44, 45-54, 55-64, 65+
  filter(concept == "SEX BY AGE") %>% 
  mutate(age = case_when(str_detect(label, "18|20|21|22|24") ~ 1,
                         str_detect(label, "29|34") ~ 2,
                         str_detect(label, "35|44") ~ 3,
                         str_detect(label, "45|54") ~ 4,
                         str_detect(label, "55|60|62") ~ 5,
                         str_detect(label, "65|67|70|75|80|85") ~ 4
                        )) %>%
  filter(
         !is.na(age)) 

#this is the post-strat table; it has multiples in some categories
#which we handle later when we build post and pred_strat
acs <- acs %>%
  clean_names() %>%
  rename(fips = geoid) %>%
  select(fips, sex, age, estimate, moe)

acs <- acs %>%
  group_by(fips,sex,age) %>%
  reframe(age=age,
          sex=sex,
          estimate=sum(estimate,na.rm=TRUE),
          moe=sum(moe,na.rm=TRUE)) %>%
  distinct()

# tab = table(ctis$vax,ctis$sex,ctis$age,ctis$edu)
# ftable(tab) %>% 
  # as.data.frame() %>% 
  # fwrite(.,file="results/quick_demographics_crosstab_fullage_.csv")

#### Direct Estimates ####

tempsex = data.frame(sex = rep(c("Female","Male"),each=5))
tempage = data.frame(age = c("18-24 years","25-34 years","35-44 years", "45-64 years", "65+ years","18-24 years","25-34 years","35-44 years", "45-64 years", "65+ years"))

#we get the general CTIS, CPS, and ACS estimates first 
#before we get the direct vaccination estimates
ddi <- ipumsr::read_ipums_ddi(here("data/cpsmarch2017/cps_00001.xml"))
 
cps <- read_ipums_micro(ddi)

# cps=
#   cps%>%
#   clean_names() %>%
#   filter(statefip == 06, #California
#          age >= 18, #adults
#          sex != 9) %>% #no nulls
#   select(sex,age,wtfinl) %>%
#   mutate(sex = ifelse(sex==2,"Female","Male"),
#          age2 = case_when(age %in% c(18:24) ~ 1,
#                           age %in% c(25:64) ~ 2,
#                           age >= 65 ~ 3))


cps=
  cps%>%
  clean_names() %>%
  filter(statefip == 06, #California
         age >= 18, #adults
         sex != 9) %>% #no nulls
  select(sex,age,wtfinl) %>%
  mutate(sex = ifelse(sex==2,"Female","Male"),
         age2 = case_when(age %in% c(18:24) ~ 1,
                          age %in% c(25:34) ~ 2,
                          age %in% c(35:44) ~ 3,
                          age %in% c(45:64) ~ 4,
                          age >= 65 ~ 5))


# cps=
#   cps%>%
#   clean_names() %>%
#   filter(statefip == 06, #California
#          age >= 18, #adults
#          sex != 9) %>% #no nulls
#   select(sex,age,wtfinl) %>%
#   mutate(sex = ifelse(sex==2,"Female","Male"),
#          age2 = case_when(age %in% c(18:24) ~ 1,
#                           age %in% c(25:64) ~ 2,
#                           age >= 65 ~ 3))

#total population used for the CTIS (less than the ACS we use!)
sum(cps$wtfinl,na.rm = TRUE) #620133 people less

cps = cps %>%
  group_by(sex,age2) %>%
  reframe(estimate = sum(wtfinl,na.rm=TRUE)) 

cps$se = sqrt(cps$estimate * (1-cps$estimate / sum(cps$estimate)))
cps$me = cps$se * 1.96
cps$lower = cps$estimate - cps$me
cps$upper = cps$estimate + cps$me
cps$prop_estimate = signif(cps$estimate/sum(cps$estimate),3)
cps$prop_lower = signif(cps$lower/sum(cps$estimate),3)
cps$prop_upper = signif(cps$upper/sum(cps$estimate),3)
cps$prop_ci.95 = paste0("(",cps$prop_lower,", ",cps$prop_upper,")")
cps$age = tempage$age

fwrite(cps,here('results/tables/weights_comparison/currentpopulationsurvey_march2017_demographic_estimates_full_age_estimates.csv'))


#check Meta weights for CTIS
ctis$fips2 <- sprintf("%05s", ctis$xfips)
ctis$fips2 <- substr(ctis$fips2, 1, 2)

#to match CDC validation data
ctis <- ctis %>%
  mutate(age2 = case_when(
    age == 1 ~ 1, #18-24
    age == 2 ~ 2, #25-64
    age == 3 ~ 2, 
    age == 4 ~ 2, 
    age == 5 ~ 3 #65+
  ))

meta_design = as_survey_design(ctis, 
                        ids = 1, 
                        weights = weight, 
                        strata = fips2)

#these are coming out correct using srvyr
meta_proportions = meta_design  %>%
  filter(fips2 == "06") %>%
  group_by(sex,age) %>%
  reframe(prop = survey_prop(vartype="ci")) 


#the weights here are daily weights, meaning we have to divide the weights
#by the number of days in the survey month to get the correct weighting for
#survey totals: n_distinct(ctis$start_date), i.e., 30 days in June. 
#This isn't a concern for the proportions above because it's proportional!

meta_estimates = meta_design  %>%
  filter(fips2 == "06") %>%
  group_by(sex,age) %>%
  reframe(est = survey_total(vartype = "ci") / n_distinct(ctis$start_date))

#however, the Meta weights are still poorly estimating the population!
sum(meta_estimates$est) == sum(cps$estimate)

#it's underestimating by 8,038,947 people -- that's a problem.
sum(cps$estimate) - sum(meta_estimates$est)

direct = cbind(meta_estimates,meta_proportions[,3:5])

direct = direct %>%
  mutate(
    across(prop:prop_upp, ~ signif(.x,3)),
    across(est:est_upp, ~ round(.x)),
    est.ci95 = paste0("(",est_low,", ",est_upp,")"),
    prop.ci95 = paste0("(",prop_low,", ",prop_upp,")")) %>%
  select(matches("est"),matches("prop"))

direct = cbind(tempsex,tempage,direct)

write.csv(direct, here('results/tables/weights_comparison/ctis_direct_demographic_full_age_estimates.csv'))


#get direct estimates of vaccines by county w/ CI
#build CTIS weights and get estimates for vaccination

# ctis$vax = factor(ctis$vax)
# ctis$age2 = factor(ctis$age2)
# ctis$sex = factor(ctis$sex)
# ctis = ctis %>% filter(!is.na(vax))

#get vaccine estimates by sex and age
  meta_proportions = meta_design  %>%
    srvyr::filter(fips2 == "06",
                  vax == 1) %>%
    group_by(sex,age2) %>%
    reframe(prop = survey_prop(vartype="ci",na.rm=TRUE)) %>%
  pull(prop) %>%
  sum()


#again, recall:
#the weights here are daily weights, meaning we have to divide the weights
#by the number of days in the survey month to get the correct weighting for
#survey totals: n_distinct(ctis$start_date), i.e., 30 days in June.
#This isn't a concern for the proportions above because it's proportional!

meta_estimates = meta_design  %>%
  filter(fips2 == "06") %>%
  group_by(vax,sex,age) %>%
  reframe(est = survey_total(vartype = "ci",na.rm = TRUE) / n_distinct(ctis$start_date)) %>%
  filter(vax == 1)


direct = cbind(meta_estimates,meta_proportions[,4:6])

direct = direct %>%
  mutate(
    across(prop:prop_upp, ~ signif(.x,3)),
    across(est:est_upp, ~ round(.x)),
    est.ci95 = paste0("(",est_low,", ",est_upp,")"),
    prop.ci95 = paste0("(",prop_low,", ",prop_upp,")")) %>%
  select(matches("est"),matches("prop"))

direct = cbind(tempsex,tempage,direct)

write.csv(direct, here('data/direct_estimates/ctis_vax_direct_estimates_not_boostrapped_no_edu.csv'))

#compare with survey package estimates, to check
library(survey)
design  = svydesign(data= ctis,ids=~1,strata = ctis$fips2, weights= ctis$weight)

#california proportions vaccinated by sex and age
ca = subset(design,fips2=="06")
alt_vax_est = svyby(~vax,by = ~sex+age,ca,svymean,na.rm = TRUE) %>% as.data.frame()
rownames(alt_vax_est) = NULL
alt_vax_est$me = alt_vax_est$se * 1.96
alt_vax_est$lower = signif(alt_vax_est$vax - alt_vax_est$me,3)
alt_vax_est$upper = signif(alt_vax_est$vax + alt_vax_est$me,3)
alt_vax_est$vax = signif(alt_vax_est$vax,3)
alt_vax_est$prop_ci.95 = paste0("(",alt_vax_est$lower,", ",alt_vax_est$upper,")")

alt_vax_est = alt_vax_est %>%
  mutate(sex = ifelse(sex==0,"Female","Male"),
         age = case_when(age==1 ~ "18-24 years",
                         age==2 ~ "25-34 years",
                         age==3 ~ "35-44 years",
                         age==4 ~ "45-64 years",
                         age==5 ~ "65+ years")) %>%
  rename(prop=vax) %>%
  arrange(sex)

temp = alt_vax_est %>%
  select(-se,-me)

#population count estimates for same
alt_vax_est = svyby(~vax,by = ~sex+age,ca,svytotal,na.rm = TRUE) %>% as.data.frame()
rownames(alt_vax_est) = NULL
alt_vax_est$me = (alt_vax_est$se * 1.96)
alt_vax_est$lower = round((alt_vax_est$vax - alt_vax_est$me)/ 30) 
alt_vax_est$upper = round((alt_vax_est$vax + alt_vax_est$me)/ 30) 
alt_vax_est$vax = round(alt_vax_est$vax / 30)
alt_vax_est$prop_ci.95 = paste0("(",alt_vax_est$lower,", ",alt_vax_est$upper,")")

alt_vax_est = alt_vax_est %>%
  mutate(sex = ifelse(sex==0,"Female","Male"),
         age = case_when(age==1 ~ "18-24 years",
                         age==2 ~ "25-34 years",
                         age==3 ~ "35-44 years",
                         age==4 ~ "45-64 years",
                         age==5 ~ "65+ years")) %>%
  rename(est=vax,
         est_low = lower,
         est_upp  = upper,
         est_ci.95 = prop_ci.95) %>%
  arrange(sex) %>%
  select(-se,-me)

alt_vax_est = left_join(temp,alt_vax_est,by=c("sex","age"))

alt_vax_est

write.csv(alt_vax_est, here('data/direct_estimates/ctis_vax_direct_estimates_not_boostrapped_no_edu.csv'))

#acs estimates at the state level
acs_state = acs %>%
  group_by(sex,age) %>%
  reframe(estimate = sum(estimate,na.rm=TRUE),
          moe = sum(moe,na.rm=TRUE))

acs_state$lower = acs_state$estimate - acs_state$moe
acs_state$upper = acs_state$estimate + acs_state$moe
acs_state$prop_estimate = signif(acs_state$estimate/sum(acs_state$estimate),3)
acs_state$prop_lower = signif(acs_state$lower/sum(acs_state$estimate),3)
acs_state$prop_upper = signif(acs_state$upper/sum(acs_state$estimate),3)
acs_state$prop_ci.95 = paste0("(",acs_state$prop_lower,", ",acs_state$prop_upper,")")

acs_state = 
  acs_state %>%
  mutate(age = case_when(
            age == 1 ~ "18-24 years",
            age == 2 ~ "25-34 years",
            age == 3 ~ "35-44 years",
            age == 4 ~ "45-64 years",
            age == 5 ~ "65+ years"))


fwrite(acs_state,here('results/tables/weights_comparison/2021_acs_5yr_demographic_estimates.csv'))



# #### unweighted proportions in the CTIS for descriptive statistics####
# ctis_desc = ctis %>%
#   mutate(n = 1,
#          sex_total = sum(n),
#          age_total = sum(n),
#          edu_total = sum(n)) %>%
#   group_by(sex) %>%
#   mutate(sex_total_group = sum(n)) %>%
#   group_by(age) %>%
#   mutate(age_total_group = sum(n)) %>%
#   group_by(edu) %>%
#   mutate(edu_total_group = sum(n)) %>%
#   ungroup() %>%
#   mutate(sex_prop = sex_total_group/sex_total,
#          age_prop = age_total_group/age_total,
#          edu_prop = edu_total_group/edu_total) %>%
#   select(sex, age, edu, sex_prop, age_prop, edu_prop) %>%
#   distinct() %>%
#   group_by(sex, age, edu) %>%
#   arrange(sex, age, edu)
# 
# 
# write.csv(ctis_desc, here("data/direct_estimates/ctis_unweighted_descriptiveprops_age_sex_edu.csv"))
# 


#### Build Spatial Adjacency Matrix ####

#compute adjacency matrix and make sure both column and row names correspond 
#to the zcta names
# geo$zcta <- as.factor(geo$zcta)
geo$xfips = as.factor(geo$geoid)

#we must assign the zctas an id 1:nrow(geo) because the adjacency matrix only works
#with integers 1:n -- and zctas aren't like that! NB: the 1:n order is artificial
#and needn't relate to the order of the zctas; it's simply important to have them
#in this format.
geo$id <- 1:nrow(geo) 

#make an index of id-to-zcta so we can refer to and use it later
index <- geo %>% 
  # select(id, zcta) %>% 
  select(id, xfips) %>%
  mutate(id = as.integer(factor(id)))

#make the neighborhood structure
nb.r <- poly2nb(pl=geo, row.names=geo$id, queen = T)

#plot
coords = st_coordinates(st_centroid(st_geometry(geo)))
# plot(st_union(st_geometry(geo)), border="grey")
# plot(nb.r, coords, add=TRUE,pch=19)


#make graph readable by INLA; equivalent to nb2mat()
nb2INLA(here("map.adj"), nb.r)
g.graph <- INLA::inla.read.graph(filename = "map.adj")

save(g.graph,file=here('data/cleaned_input_data/spatialgraph.rda'))

#### Build Data (post-strat data and survey data): Case-Specific for Spatial MRP ####
# 
# # #make sure acs has id-to-zcta so can identify between later
# # acs <- left_join(acs, index, by = "zcta")
# 
# #poststratification data: combine variables where because of recoding, 
# #there are multiples (i.e. collapse strata that now are identical across id, 
# #sex, age, edu due to recode)
# 
# post <- acs %>%
#   group_by(id,sex,age,edu) %>%
#   summarize(estimate = sum(estimate)) %>%
#   select(id, sex, age, edu, estimate) %>%
#   ungroup() #double-check no grouping
# 
# 
# #make sure ctis has id-to-zcta info also
# dat <- left_join(dat, index, by = "zcta")
# 
# #Make sure nothing is a factor; you can assign in the formula later; tidy up data
# dat <- dat %>%
#   # select(weight, id, vax, sex, age, edu) %>% # don't need weights in full Bayesian!
#   select(id, vax, sex, age, edu) %>%
#   mutate(age = as.integer(age),
#          edu = as.integer(edu))


#### Build Data for Post Strat: No spatial id, no edu ####
# #make sure acs has id-to-zcta so can identify between later
# acs <- left_join(acs, index, by = "zcta")

#poststratification data: combine variables where because of recoding, 
#there are multiples (i.e. collapse strata that now are identical across id, 
#sex, age, edu due to recode)

post <- acs %>%
  group_by(fips, sex,age) %>%
  summarize(estimate = sum(estimate)) %>%
  select(fips, sex, age, estimate) %>%
  ungroup() #double-check no grouping

post = post %>% rename(xfips = fips)

#Make sure nothing is a factor; you can assign in the formula later; tidy up data
dat <- dat %>%
  # select(weight, id, vax, sex, age, edu) %>% # don't need weights in full Bayesian!
  select(xfips, vax, sex, age) %>%
  mutate(age = as.integer(age))


#write graph, acs, dat and post object to rda
save(geo, file=here("data/cleaned_input_data/ca_geo_sf_county.rda"))
save(dat, file = here("data/cleaned_input_data/clean_data_county.rda"))
save(post, file = here("data/cleaned_input_data/clean_postrat_age_sex_county.rda"))

