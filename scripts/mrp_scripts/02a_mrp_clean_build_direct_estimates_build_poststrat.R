#### MRP + Space using CTIS Data ####

#Script 1: Clean data, build survey direct estimates, and poststratification data 

#Build a series of models representing vaccination probability (min. 1 dose)
#For June 12, 2021 using the COVID-19 Trends and Impacts Survey (CTIS)

#Post-stratification of age/sex/education at county level

pacman::p_load(tidyverse, tidycensus, haven, here, janitor, survey,srvyr, sf, spdep)


options(tigris_use_cache = TRUE)

tidycensus::census_api_key(Sys.getenv("CENSUS_API_KEY"), overwrite = FALSE, install = FALSE)

here::i_am("scripts/mrp_scripts/01_mrp_clean_build_direct_estimates_build_poststrat.R")


#### set geographic and model string specifications for intuitive labeling ####
state_choice <- "CA"
state <- "California"
whole_time_period <- FALSE


 
#### read/clean survey data: ctis 3, wave 11, june 2021 ####
ctis <- read.csv(file=here::here("data/ctis_data_do_not_share/ctis_w11_2021-06.csv"),
                 header = TRUE, colClasses=c("fips"="character", "A3"="character"))

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


#to match the ACS categories listed later in this script under object "acs"
# ctis <- ctis %>%
#   mutate(age = case_when(
#     age == 1 ~ 1, #18-24
#     age == 2 ~ 2, #25-34
#     age == 3 ~ 3, #35-44
#     age == 4 ~ 4, #45-54
#     age == 5 ~ 4, #55-64
#     age == 6 ~ 5, #65-74
#     age == 7 ~ 5  #75+
#   ))


#to match the ACS And the basic assumptions of the CDC baselines in file.
#We also keep the original age category for the purposes of producing direct estimates
ctis <- ctis %>%
  mutate(age2 = case_when(
    age == 1 ~ 1, #18-24
    age == 2 ~ 2, #25-34
    age == 3 ~ 2, #35-44
    age == 4 ~ 2, #45-54
    age == 5 ~ 2, #55-64
    age == 6 ~ 2, #65-74
    age == 7 ~ 3  #75+
  ))


#to match acs education values
ctis <- ctis %>%
  mutate(edu = ifelse(edu==6 | edu == 7 | edu == 8, 6, edu))


#### RUCA data: identify urban/rural classification using ZIP code info ####
# # link 2010 RUCA (rural urban classification): https://www.ers.usda.gov/data-products/rural-urban-commuting-area-codes/
#
#   ruca <- read_csv(here::here("data/smrp/RUCA_2021_zipcode_ruralurban.csv"),
#                    col_types = cols(ZIP_CODE = col_factor()),
#                    col_names = T)
#
#     ruca <- ruca %>%
#       filter(RUCA1 != 99) %>%
#       mutate(urban = case_when(RUCA1 < 6 ~ 1,
#                                RUCA1 > 5 ~ 0)) %>%
#       clean_names() %>%
#       select(zip_code, state, urban)
#
#
#   #join the urban/rural classification information to our ctis data
#   ctis <- left_join(ctis, ruca, by = c("xzip" = "zip_code"))


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
#
# #tidy up so we have state, zcta, geometry
# geo <- geo %>%
#   select(statefp10, zcta5ce10, intptlat10, intptlon10, geometry) %>%
#   rename(state=statefp10, zcta=zcta5ce10)
#
# #### crosswalk needs check ####
# a crosswalk is when we convert one thing to another using a file provided by
# e.g. the government (as is the case here). ZIP codes and ZCTAs are not the same
# thing and they change over time, so you have to do some logical matching to make
# sure they're as close to matching as possible.
# #need to crosswalk -- there are fewer zip codes than zctas here, and they might not equate!
# length(unique(dat$xzip)) #1550
# length(unique(geo$zcta)) #1769

# complete crosswalk with dat (which is just the ctis data filtered down for CA)
#this might not be an ideal file because the zcta boundaries are from 2010
#in the tigris files and this is from 2021
# crswlk <- read.csv(file = here("Aja/COVID_Behaviors/data/ZiptoZcta_Crosswalk_2021.csv"),
#                    header = TRUE,
#                    colClasses = c("zip_code"="character")) %>%
#           clean_names()


# dat <- left_join(dat, crswlk, by=c("xzip"="zip_code")) %>%
#   select(weight, zcta, xzip, xfips, vax, sex, age, edu)

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
                year = 2020)

vars <- load_variables(2020, "acs5", cache=TRUE)

vars <- vars %>%
  filter(str_detect(name, "B01001"),
         str_detect(label, "18|20|21|22|24|29|30|35|40|45|50|55|60|62|65|67|70|75|80|85"))

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

acs <- acs %>%
  mutate(sex = ifelse(str_detect(label, "Male"), 1, 0))  # 1= male, 0= female, like in ctis

acs <- acs %>%
  mutate(age = case_when(str_detect(label, "18|20|21|22|24") ~ 1,
                         str_detect(label, "29|30|35|40|45|50|55|60|62|65|67|70") ~ 2,
                         str_detect(label, "75|80|85") ~ 3
                         )) %>%
      filter(!str_detect(concept, "ALONE|HISPANIC|RACES"))

#this is the post-strat table; it has multiples in some categories
#which we handle later when we build post and pred_strat
acs <- acs %>%
  clean_names() %>%
  rename(fips = geoid) %>%
  select(fips, sex, age, estimate)

acs <- acs %>%
  group_by(fips, sex,age) %>%
  reframe(age=age,
            sex=sex,
            estimate=sum(estimate)) %>%
  distinct()


#### Direct Estimates ####
#get direct estimates of masking in dat by county w/ CI

ctis$fips2 <- sprintf("%05s", ctis$xfips)
ctis$fips2 <- substr(ctis$fips2, 1, 2)

# #reweight to specifically CA population
# #uncertain about this:
# ctis$weight2 = 30513773 * (ctis$weight/sum(ctis$weight))

ctis_usd <- svydesign(ids = ~1, 
                      weights = ctis$weight, 
                      strata = ctis$fips2, 
                      data = ctis)


# 567289216/30513773*30
# 567289216/30
# 1830374+2995085+2759655+4847781+2702808+1742201+2839782+2639320+4865089+3291678
# 3794976/25296329


temp = c(
  as.data.frame(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 0 & age2 == 1), na.rm=T)),
  as.data.frame(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 0 & age2 == 2), na.rm=T)),
  as.data.frame(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 0 & age2 == 3), na.rm=T)),
  as.data.frame(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 1 & age2 == 1), na.rm=T)),
  as.data.frame(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 1 & age2 == 2), na.rm=T)),
  as.data.frame(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 1 & age2 == 3), na.rm=T))
  )

temp = as.data.frame(temp)
temp2 = temp[grepl("mean", names(temp))] %>% t()
rownames(temp2) = NULL

confint = c(
  as.data.frame(confint(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 0 & age2 == 1), na.rm=T),level = 0.95)),
  as.data.frame(confint(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 0 & age2 == 2), na.rm=T),level = 0.95)),
  as.data.frame(confint(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 0 & age2 == 3), na.rm=T),level = 0.95)),
  as.data.frame(confint(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 1 & age2 == 1), na.rm=T),level = 0.95)),
  as.data.frame(confint(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 1 & age2 == 2), na.rm=T),level = 0.95)),
  as.data.frame(confint(svymean(~vax, subset(ctis_usd,fips2=="06" & sex == 1 & age2 == 3), na.rm=T),level = 0.95))
)

confint = as.data.frame(confint)
confint_lw = confint[grepl("X2", names(confint))] %>% t()
confint_up = confint[grepl("X9", names(confint))] %>% t()
rownames(confint_lw) = NULL
rownames(confint_up) = NULL

tempsex = data.frame(sex = c("female","female","female","male","male","male"))
tempage = data.frame(age = c("18-24 years","25-64 years","65+ years","18-24 years","25-64 years","65+ years"))

temp = cbind(tempsex,tempage,temp2,confint_lw,confint_up)

direct = temp %>% rename("mean"="temp2","mean_lw"="confint_lw","mean_up"="confint_up")

direct = direct %>%
  mutate(
         mean = signif(mean, digits = 3),
         mean_lw = signif(mean_lw, digits = 3),
         mean_up = signif(mean_up, digits = 3),
         ci = paste0("(",mean_lw,", ",mean_up,")"))

write.csv(direct, here("/Users/Aja/Documents/R Directories/NSF_RAPID_COVID19_Survey/Aja/COVID_Behaviors/data/smrp/direct_estimates/ctis_vax_direct_estimates_not_boostrapped.csv"))


#### unweighted proportions in the CTIS for descriptive statistics####
ctis_desc = ctis %>%
  mutate(n = 1,
         sex_total = sum(n),
         age_total = sum(n),
         edu_total = sum(n)) %>%
  group_by(sex) %>%
  mutate(sex_total_group = sum(n)) %>%
  group_by(age) %>%
  mutate(age_total_group = sum(n)) %>%
  group_by(edu) %>%
  mutate(edu_total_group = sum(n)) %>%
  ungroup() %>%
  mutate(sex_prop = sex_total_group/sex_total,
         age_prop = age_total_group/age_total,
         edu_prop = edu_total_group/edu_total) %>%
  select(sex, age, edu, sex_prop, age_prop, edu_prop) %>%
  distinct() %>%
  group_by(sex, age, edu) %>%
  arrange(sex, age, edu)


write.csv(ctis_desc, here("/Users/Aja/Documents/R Directories/NSF_RAPID_COVID19_Survey/Aja/COVID_Behaviors/data/smrp/direct_estimates/ctis_unweighted_descriptiveprops_age_sex_edu.csv"))



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
plot(st_union(st_geometry(geo)), border="grey")
plot(nb.r, coords, add=TRUE,pch=19)


#make graph readable by INLA; equivalent to nb2mat()
nb2INLA(here("map.adj"), nb.r)
g.graph <- INLA::inla.read.graph(filename = "map.adj")

save(g.graph, file = here("Aja/COVID_Behaviors/data/smrp/cleaned_data/spatialgraph.rda"))

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

#add the post-strat data to the bottom of the dat data -- as the prediction covariates!
# pred_dat <- post %>%
#   mutate(vax = NA) %>%
#   relocate(vax, .after=id) %>%
#   select(-estimate)

# if(urban_include == TRUE){
#   #scratch code for joining pred_dat to ruca
#   index2 <- index %>% sf::st_drop_geometry()
#   index2 <- left_join(index2, crswlk, by=c("zcta"="zcta"))
#   index2 <- index2 %>% select(id, zcta, zip_code)
#   
#   pred_dat <- left_join(pred_dat, index2, by = "id")
#   
#   ruca <- left_join(ruca, crswlk, by=c("zip_code"="zip_code"))
#   ruca <- ruca %>% na.omit() %>% select(zcta, urban)
#   
#   pred_dat <- left_join(pred_dat, ruca, by=c("zcta"))
#   
#   pred_dat <- pred_dat %>% select(id, vax, sex, age, edu, urban)
# }

# #they match, so we can rbind
# colnames(pred_dat) == colnames(dat)

# dat <- rbind(dat, pred_dat)


#write graph, acs, dat and post object to rda
save(geo, file=here("Aja/COVID_Behaviors/data/smrp/cleaned_data/ca_geo_sf_county.rda"))
save(dat, file = here("Aja/COVID_Behaviors/data/smrp/clean_data_county.rda"))
save(post, file = here("Aja/COVID_Behaviors/data/smrp/clean_postrat_age_sex_county.rda"))

