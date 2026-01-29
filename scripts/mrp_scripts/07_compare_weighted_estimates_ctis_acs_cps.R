#### Compare design weighted estimates at the state level: CTIS vs CPS vs ACS ####

#The CTIS uses weights provided by Meta. 
#The third stage of these weights is a post-stratification step by sex and age 
#using the US Current Population Survey March 2017 estimates. 
#I have pulled these data from IPUMS for Sex and Age. 

#For our paper, as of Jan 29 2026, we are using the 
#2023 5-year American Community Survey Estimates. 
#If need be, we can easily change this, 
#or explore other years of the 5-year ACS (2019, perhaps?).

#Below is a comparison of the population estimates for sex and age between
#the CTIS, the CPS, and the ACS. These were build in script 02a.

#We are comparing: sex (Male, Female) and age (18-24, 25-64, 65+) because these
#are the strata made available in the CTIS.

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, data.table,scales)

here::i_am("scripts/mrp_scripts/07_compare_weighted_estimates_ctis_acs_cps.R")

options(scipen = 1)

#load our dataframes
acs = fread(here('results/tables/weights_comparison/2023_acs_5yr_demographic_estimates.csv'))
ctis = fread(here('results/tables/weights_comparison/ctis_direct_demographic_estimates.csv'))
cps = fread(here('results/tables/weights_comparison/currentpopulationsurvey_march2017_demographic_estimates.csv'))

acs = acs %>%
  select(-moe) %>%
  rename(estimate =estimate,
         estimate_low =lower,
         estimate_upp =upper,
         prop = prop_estimate,
         prop_low = prop_lower,
         prop_upp = prop_upper) %>%
  select(sex,age,estimate :estimate_upp,prop:prop_upp)

cps = cps %>%
  select(-me,-se) %>%
  rename(age=age2,
         estimate_low=lower,
         estimate_upp=upper,
         prop = prop_estimate,
         prop_low = prop_lower,
         prop_upp = prop_upper) %>%
  select(sex,age,estimate:estimate_upp,prop:prop_upp)

ctis = ctis %>%
  select(-V1) %>%
  rename(estimate=est,
         estimate_low = est_low,
         estimate_upp = est_upp) %>%
  select(sex,age,estimate:estimate_upp,prop:prop_upp)

ctis = ctis %>% pivot_longer(cols = c(estimate:prop_upp),names_to = "measure",values_to = "estimate")
acs = acs %>% pivot_longer(cols = c(estimate:prop_upp),names_to = "measure",values_to = "estimate")
cps = cps %>% pivot_longer(cols = c(estimate:prop_upp),names_to = "measure",values_to = "estimate")

test = rbind(ctis,cps,acs)
test$source = rep(c("CTIS","CPS March 2017","2023 5-year ACS"),each=36)

diffs = test %>%
  group_by(sex,age,measure) %>%
  reframe(diff = diff(estimate))
diffs$source = rep(c("CTIS vs CPS","CPS vs ACS"),times = 36)
diffs = diffs %>%
  select(source,sex,age,measure,diff)

test= pivot_wider(test,id_cols = c(source,sex:age),names_from = "measure",values_from = "estimate")

test = test %>%
  arrange(sex,age)

diffs = pivot_wider(diffs,id_cols = c(source,sex,age),names_from = "measure",values_from = "diff")

test = test %>%
  mutate(across(estimate:estimate_upp, ~ round(.x)),
         across(prop:prop_upp, ~ signif(.x,3)))

diffs = diffs %>%
  mutate(across(estimate:estimate_upp, ~ round(.x)),
         across(prop:prop_upp, ~ signif(.x,3))) %>%
  arrange(source,sex,age)


fwrite(test, here('/Users/Aja/Documents/R Directories/MRP_Spatial/results/tables/weights_comparison/full_weights_estimates_ctis_cps_acs.csv'))
fwrite(diffs, here('/Users/Aja/Documents/R Directories/MRP_Spatial/results/tables/weights_comparison/differences_in_full_weights_estimates_ctis_cps_acs.csv'))

table(test$source)
test %>%
  mutate(source= factor(source,levels=c("CTIS","CPS March 2017","2023 5-year ACS"))) %>%
  ggplot(aes(x= source,y=estimate,color=age,shape=sex)) +
  geom_pointrange(aes(ymin = estimate_low,ymax=estimate_upp),position=position_dodge(0.7))

test %>%
  mutate(source= factor(source,levels=c("CTIS","CPS March 2017","2023 5-year ACS"))) %>%
  ggplot(aes(x= source,y=prop,color=age,shape=sex)) +
  geom_pointrange(aes(ymin = prop_low,ymax=prop_upp),position=position_dodge(0.7))

diffs %>%
  mutate(source= factor(source,levels=c("CTIS vs CPS","CPS vs ACS"))) %>%
  ggplot(aes(x= source,y=prop,color=age,shape=sex)) +
  geom_pointrange(aes(ymin = prop_low,ymax=prop_upp),position=position_dodge(0.5))

diffs %>%
  mutate(source= factor(source,levels=c("CTIS vs CPS","CPS vs ACS"))) %>%
  ggplot(aes(x= source,y=estimate,color=age,shape=sex)) +
  geom_pointrange(aes(ymin = estimate_low,ymax=estimate_upp),position=position_dodge(0.5)) +
  scale_y_continuous(labels = label_comma())


