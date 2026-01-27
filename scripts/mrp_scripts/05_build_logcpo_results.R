#### Produce log-cpo summary ####

pacman::p_load(tidyverse,fs,here,janitor,qdapRegex,viridis,data.table,stringr,sf,tigris)
options(scipen=999)
here::i_am("scripts/mrp_scripts/05_build_logcpo_results.R")


#load list of models in final analysis
load(file = here("data/cleaned_input_data/model_descriptions.RData"))
models =
  models %>%
  select(model_name,specs)
rm(model_formulae) #don't need it


#only need to run once:
### tidy up regression results

files <- list.files(
    path = here("data/model_objects"),
    pattern = "*\\.RDS$",
    full.names = TRUE) # include the directory in the result

  #filter files that contain “rdata”
  files <- files %>%
    tibble() %>%
    rename(files = 1) %>%
    filter(stringr::str_detect(files, ".rds|.RDS")) %>%
    pull()

  list_of_files <- list() #create empty list


#loop through the files
  for (i in files) {
    print(i)
    list_of_files[[i]] <- get(load(paste0(i))) #add files to list position
  }

for (i in 1:length(list_of_files)) {

  model = list_of_files[[i]]

  model_specs = str_split(files[i],"model_objects/",simplify=TRUE)[2]
  model_specs = str_split( model_specs,"_model_object.RDS",simplify=TRUE)[1]

  model.summary = model$summary.fixed %>%
    janitor::clean_names()

  # model.summary = model.summary %>%
  # mutate(id_label = ifelse(variable == "fixed", paste(variable,name)),id_label)#,
          #  id_label = case_when(
          # variable == "random" & name == "age" & id == 1 ~ "18-24 years",
          #  variable == "random" & name == "age" & id == 2 ~ "25-64 years",
          #  variable == "random" & name == "age" & id == 3 ~ "65 years and over",
          #  variable == "random" & name == "edu" & id == 1 ~ "less than high school",
          #  variable == "random" & name == "edu" & id == 2 ~ "high school or equivalent",
          #  variable == "random" & name == "edu" & id == 3 ~ "some college",
          #  variable == "random" & name == "edu" & id == 4 ~ "associate's degree",
          #  variable == "random" & name == "edu" & id == 5 ~ "Bachelor's degree",
          #  variable == "random" & name == "edu" & id == 6 ~ "professional or graduate degree"))

  #
  # length_helper = model.summary %>% filter(!is.na(mean), name =="id") %>% n_distinct
  #
  # #if it's obviously a county random effect county
  # if (length_helper == 57 | length_helper == 58) { #there are models that don't estimate one county
  #
  #     #load id codes for counties: as "geo"
  #     load(here("data/cleaned_input_data/geo_index.rda"))
  #     geo = geo %>%
  #       select(id,namelsad10) %>%
  #       sf::st_drop_geometry() %>%
  #       rename(county_name = namelsad10)
  #
  #   not_county_model.summary =
  #     model.summary %>%
  #     filter(name !="id")
  #
  #   county_model.summary =
  #     model.summary %>%
  #     filter(name == "id") %>%
  #    # mutate(id = as.numeric(id)-1) %>%
  #     left_join(.,geo, by = "id") %>%
  #     mutate(id_label = county_name) %>%
  #     select(-county_name)
  #
  #   model.summary = rbind(not_county_model.summary,county_model.summary)
  # }

  model.summary$variable = rownames(model.summary)

  #get LCPO
  model.summary$lcpo = NA
  model.summary$lcpo = -mean(log(list_of_files[[i]]$cpo$cpo), na.rm = TRUE)

  failure_check = sum(model$cpo$failure, na.rm = TRUE)

  if (failure_check > 0) {
    stop("log-cpo failure")
  }

  if (failure_check == 0) {
    print(paste0("model ",files[i],": no failures in log-cpo"))
    rownames(model.summary) = NULL

    model.summary = model.summary %>%
      select(
        #variable:id,
             variable:lcpo,everything())

    fwrite(model.summary, file = here(paste0("data/raw_model_output/model_results/",model_specs,".csv")))

    #save proportions: tau, phi, etc.

    proportions = model$summary.hyperpar
    proportions$label = rownames(proportions)
    rownames(proportions) = NULL

    fwrite(proportions, file = here(paste0("data/raw_model_output/model_results/hypers/hyper_",model_specs,".csv")))


    rm(model.summary)
    gc()
    }


  }


#load("data/cleaned_input_data/clean_complete_model_indirect_estimates.rda")

rm(list_of_files,model,models,proportions)
gc()

#Only have to run the one time (joins all the modeled results together)
filelist <- list.files(
  path = here("data/raw_model_output/model_results"),
  full.names = TRUE) # include the directory in the result)

dat = fs::dir_ls(here("/Users/Aja/Documents/R Directories/MRP_Spatial/data/raw_model_output/model_results"), regexp = "\\.csv$")
 
dat = map_dfr(dat,fread, .id = 'id')

dat = 
  dat %>%
    filter(!str_detect(id,"state")) %>%
  mutate(sex = ifelse(ex_between(id,"sex_",".csv") == 1, "male","female"),
         specs = ex_between(id,"model_results/","_sex")) 

#save a legible version of the regression results
dat = dat %>%
  select(specs,sex,variable,lcpo,mean:mode) %>%
  distinct()

dat$specs = as.character(dat$specs)
dat = left_join(dat,models, by = c("specs"))

n_distinct(dat$model_name)


dat =
  dat %>%
  rename(lower_ci = x0_025quant,
         upper_ci = x0_975quant,
         median = x0_5quant) %>%
  #select(-kld) %>%
  mutate(lcpo = round(lcpo, 4),
         #kld = signif(kld, 3),
         #sex = ifelse(sex == 0, "female","male"),
         mean = signif(mean,2),
         sd = signif(sd,3),
         lower_ci = signif(lower_ci,2),
         median = signif(median,2),
         upper_ci = signif(upper_ci,2),
         mode = signif(mode,2),
         ci = as.character(paste0("(",lower_ci,", ",upper_ci,")"))) %>%
  select(model_name,sex:lcpo,median,ci,mean,sd)




#save regression results for paper before filtering to log-CPO table details
write.csv(dat,here("results/tables/table7_inla_regression_results.csv"))

dat1 =
  dat %>%
  select(model_name, sex, lcpo) %>%
  group_by(sex) %>%
  distinct() %>%
  arrange(sex,lcpo)

write.csv(dat1, file = here("results/tables/table3_model_logcpo_results.csv"))


#precision
dat = fs::dir_ls(here("/Users/Aja/Documents/R Directories/MRP_Spatial/data/raw_model_output/model_results/hypers"), regexp = "\\.csv$")

dat = map_dfr(dat,fread, .id = 'id')

dat = 
  dat %>%
  filter(!str_detect(id,"state")) %>%
  mutate(sex = ifelse(ex_between(id,"sex_",".csv") == 1, "Male","Female"),
         specs = ex_between(id,"hyper_","_sex")) 

dat %>% 
  select(specs,sex,label,mean:mode) %>%
 view
# #save a legible version of the regression results
# dat = dat %>%
#   select(specs,sex,variable,lcpo,mean:mode) %>%
#   distinct()

dat$specs = as.character(dat$specs)
dat = left_join(dat,models, by = c("specs"))

colnames(dat)
dat =
  dat %>%
  rename(lower_ci = `0.025quant`,
         upper_ci = `0.975quant`,
         median = `0.5quant`) %>%
  #select(-kld) %>%
  mutate(
        #lcpo = round(lcpo, 2),
         #kld = signif(kld, 3),
         #sex = ifelse(sex == 0, "female","male"),
         mean = signif(mean,2),
         sd = signif(sd,3),
         lower_ci = signif(lower_ci,2),
         median = signif(median,2),
         upper_ci = signif(upper_ci,2),
         mode = signif(mode,2),
         ci = as.character(paste0("(",lower_ci,", ",upper_ci,")"))) %>%
  select(model_name,sex,label,median,ci,mean,sd)

dat

view(dat)
#save regression results for paper before filtering to log-CPO table details
write.csv(dat,here("results/tables/table9_model_precision_results.csv"))



