#### create large regression results table ####

pacman::p_load(tidyverse,fs,here,janitor,qdapRegex,viridis,data.table,stringr,sf,tigris)

options(scipen=999)
here::i_am("scripts/mrp_scripts/05_build_logcpo_results.R")

 dat = fs::dir_ls(here("/Users/Aja/Documents/R Directories/MRP_Spatial/data/raw_model_output/model_results"), regexp = "\\.csv$")
 

 dat = map_dfr(dat,fread, .id = 'id',colClasses = 'character')

 
 dat = dat %>%
   mutate(id = str_remove(id,"/Users/Aja/Documents/R Directories/MRP_Spatial/data/raw_model_output/model_results/"),
          id = str_remove(id,".csv"),
          sex = ifelse(substr(id,nchar(id),nchar(id))=="1","male","female"),
          id = str_remove(id,"_sex_0|_sex_1"))

 load(file = here("data/cleaned_input_data/model_descriptions.RData"))
 models =
   models %>% 
   select(model_name,specs)
 rm(model_formulae) #don't need it
models$specs
table(dat$id)
 #dat = dat %>% mutate(id = str_replace_all(id,", ","_"))
dat = left_join(dat,models, by = c("id"="specs"))

dat = dat %>%
  select(model_name,sex,everything())
# glimpse(dat)
# table(dat$model_name)
# models$model_name[3]
# 
# models$specs[4]
# 
# `%nin%` = Negate(`%in%`)
# 
# dat %>% filter(id %nin% models$specs) %>% distinct()
# n_distinct(dat$id)
# n_distinct(models$specs)
# 
# glimpse(dat)

#clean up table for pub

dat = 
  dat %>%
  mutate(sex = stringr::str_to_sentence(sex, locale = "en"),
        # name = str_replace_all(name,":",": "),
        # name = str_replace_all(name,"edu","education"),
        # name = str_replace_all(name,"id","county"),
         median = as.character(signif(as.double(median),3)),
         lower_ci = as.character(signif(as.double(lower_ci),3)),
         upper_ci = as.character(signif(as.double(upper_ci),3)),
        lcpo = as.character(signif(as.double(lcpo),4)),
        ci = paste0(lower_ci, ", ",upper_ci),
         mean = as.character(signif(as.double(mean),3)),
         sd = as.character(signif(as.double(sd),3)),
         ) %>%
  select(model_name,sex,name,variable,id_label,lcpo,median,ci,mean,sd) 

fwrite(dat,file=here('results/tables/inla_regression_results.csv'))





save(labs,file=here('/Users/Aja/Documents/R Directories/MRP_Spatial/data/model_objects/regression_results_simple_labels.Rdata'))


dat = 
  dat %>%
 filter(str_detect(model_name,"education-age|BYM2 by age")) 

#label the regression results







fwrite(dat,file=here('results/tables/inla_regression_results_extras.csv'))
