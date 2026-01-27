#### Run from source #####
pacman::p_load(tidyverse, here, beepr,)

here::i_am("scripts/00_run_from_source.R")

#build models and data objects
{
  source(here("scripts/mrp_scripts/01_model_list.R"))
  rm(list=ls())
  source(here("scripts/mrp_scripts/02_mrp_clean_with_edu_build_direct_estimates_build_poststrat.R"))
  rm(list=ls())
  gc()
}

beep(2)

#run analyses
{
  #run models
  source(here("scripts/mrp_scripts/03_run_models.R"))
  #build
  rm(list=ls())
  gc()
  }

beep(2)

#build plots and tables
source(here("scripts/mrp_scripts/04_visualize_results.R"))
rm(list=ls())
gc()

source(here("scripts/mrp_scripts/05_build_logcpo_results.R"))
rm(list=ls())
gc()

beep(8)

