### simulating RCD in Once Bitten Data

library(lubridate)
library(tidyverse)



demo_file <- read.csv('from_amy/spat21_human_final_censored_data_for_dissertation_with_exposure_outcome_1MAR2020.csv', header = TRUE)
demo_file$malaria_symp = sapply(1:nrow(demo_file), function(x) ifelse(is.element('yes', unique(c(demo_file$fever[x], demo_file$Aches[x], demo_file$Vomiting[x], demo_file$Diarrhea[x], demo_file$Chills[x], demo_file$congestion[x], demo_file$Cough[x]))), 1, 0))

indiv_hh_village_file <- demo_file[,c('unq_memID', 'village_name', 'gender', 'age_cat_baseline', 'age_y_baseline', 'age_m_baseline', 'HH_ID')]
indiv_hh_village_file <- unique(indiv_hh_village_file)

sub_demo_file <- demo_file[,c('unq_memID', 'rdt_rst', 'sample_id_date', 'pf_pcr_infection_status', 'pfr364Q_std_combined', 'visit_type', 'village_name', 'gender', 'age_cat_baseline', 'age_y_baseline', 'age_m_baseline', 'malaria_symp', 'travelled', 'HH_ID')]

sub_demo_file$age <- ifelse(is.na(sub_demo_file$age_y_baseline), sub_demo_file$age_m_baseline/12, sub_demo_file$age_y_baseline)
sub_demo_file$date = ymd(sub_demo_file$sample_id_date)
#sub_demo_file <- filter(sub_demo_file, age < 10)
sub_demo_file$infect_type = ifelse(sub_demo_file$pf_pcr_infection_status == 'positive', ifelse(sub_demo_file$rdt_rst == 'positive' & sub_demo_file$malaria_symp == 1, 'symp', 'asymp'), 'neg')
length(unique(demo_file$unq_memID))
length(unique(demo_file$HH_ID))


### sample just symptomatic events and then sample household members from the most recent visit 


### human csp file
csp_haplo_file <- read.csv('from_amy/spat21_csp_haplotype_table_censored_final_version_with_moi_and_ids_CLEANVERSION_30SEPT2019.csv', header = TRUE)
human_csp_file <- filter(csp_haplo_file, sample_type == 'Human')
human_csp_file$date <- sapply(human_csp_file$sample_name_dbs, function(x) strsplit(as.character(x), split = '-')[[1]][2])
human_csp_file$date <- dmy(human_csp_file$date)
human_csp_file$unq_memID <- sapply(human_csp_file$sample_name_dbs, function(x) paste(strsplit(as.character(x), split = '-')[[1]][1], strsplit(as.character(x), split = '-')[[1]][3], sep = '_'))
human_csp_file$month_year <- paste(year(human_csp_file$date), month(human_csp_file$date), sep = '-')
human_csp_file$week <- week(human_csp_file$date)
human_csp_file$year <- year(human_csp_file$date)
human_csp_file$week_year <- paste(human_csp_file$week, human_csp_file$year, sep = '-')
#human_csp_file$biweek <- sapply(human_csp_file$week, function(x) find_biweek_func(x))
#human_csp_file$biweek_year <- paste(human_csp_file$biweek, human_csp_file$year, sep = '-')
human_csp_file$hh_id <- sapply(human_csp_file$unq_memID, function(x) strsplit(as.character(x), split = '_')[[1]][1])
human_csp_file <- left_join(human_csp_file, indiv_hh_village_file, by = 'unq_memID')
human_csp_file$unq_memID_date = paste(human_csp_file$unq_memID, human_csp_file$date, sep = '__')


## index events are only those with a sick visit where the individual was RDT positive 
sick_pos_visits <- filter(sub_demo_file, visit_type == 'sick visit' & rdt_rst == 'positive')
sick_pos_visits$unq_memID_date = paste(sick_pos_visits$unq_memID, sick_pos_visits$date, sep = '__')
## total 123 

### pull CSP haplotypes for index events and make a new ID that will make a label for that RCD event. reactve_event variable = yes if you are the index, no = if you are not the index
sick_pos_csp <- left_join(sick_pos_visits, human_csp_file, by = 'unq_memID_date') 
sick_pos_csp$reactive_id = seq(1, nrow(sick_pos_csp))
sick_pos_csp$reactive_id = paste('rcd', sick_pos_csp$reactive_id, sep = '_')
sick_pos_csp$reactive_event = 'yes'
sick_pos_csp$month_year_cfm <-  paste(year(sick_pos_csp$sample_id_date), month(sick_pos_csp$sample_id_date), sep = '-')

### you can set a day threshold to determine which household events you'll count as part of the reactive case event 
day_threshold = 45
final_csp_match_event <- list()
### mark all csp events as related to a rcd event 
for(ii in 1:nrow(sick_pos_csp)){
  ### some general code to look at matches based on household ids and samples that fall within the time window 
  date_event = ymd(sick_pos_csp$sample_id_date[ii])
  household_event = sick_pos_csp$HH_ID.x[ii]
  memID_event = sick_pos_csp$unq_memID.x[ii]
  possible_matches <- filter(human_csp_file, hh_id == household_event)
  possible_matches <- filter(possible_matches, abs(ymd(date) - date_event)<=day_threshold)
  possible_matches <- filter(possible_matches, memID_event != unq_memID)
  print(c(ii, nrow(possible_matches)))
  if(nrow(possible_matches)>0){
    ## remove duplicates, only one per person and choose the most recent one 
    possible_matches_indiv_ids <- sort(unique(possible_matches$unq_memID))
    possible_match_no_dup <- c()
    for(bb in 1:length(possible_matches_indiv_ids)){
      sub_data = filter(possible_matches, unq_memID == possible_matches_indiv_ids[bb])
      if(nrow(sub_data)>1){
        date_diff = abs(ymd(sub_data$date) - date_event)
        most_recent <- sub_data[which(date_diff == min(date_diff)),]
        possible_match_no_dup <- rbind(possible_match_no_dup, most_recent)
      }
      else{
        possible_match_no_dup <- rbind(possible_match_no_dup, sub_data)
      }
    }
    ## make an id for the data to be categorized as a RCD related event
    possible_match_no_dup$reactive_id = sick_pos_csp$reactive_id[ii]
    ## reactive_event variable: yes = index case, no = not index case
    possible_match_no_dup$reactive_event = 'no'
    final_csp_match_event[[sick_pos_csp$reactive_id[ii]]] = possible_match_no_dup
  }
}

## set a limit on the number of reactive events per month
monthly_limit <- 30
rcd_dates <- sick_pos_csp %>%
  select(reactive_id, month_year_cfm, sample_id_date) %>%
  filter(reactive_id %in% names(final_csp_match_event)) %>%
  arrange(sample_id_date) %>%
  group_by(month_year_cfm) %>%
  mutate(id = row_number()) %>%
  filter(id <= monthly_limit)

#restrict list to limited reactive events
final_csp_match_event_limited <- final_csp_match_event[rcd_dates$reactive_id]

#unlist into one df
hh_matches <- bind_rows(final_csp_match_event_limited, .id = "reactive_id")

#combine sick visits and hh matches
skinny_OB <- sick_pos_csp %>%
  filter(reactive_id %in% rcd_dates$reactive_id) %>%
  select(unq_memID.x, village_name.x, gender.x, age_cat_baseline.x, age_y_baseline.x, age_m_baseline.x, HH_ID.x, date.x, 18:319, sample_type, sample_name_dbs, haplotype_number, haplotype_reads, month_year_cfm, week, year, week_year, hh_id, reactive_id, reactive_event) %>%
  rename("unq_memID" = "unq_memID.x", "village_name" = "village_name.x", "gender" = "gender.x", "age_cat_baseline" = "age_cat_baseline.x", "age_m_baseline" = "age_m_baseline.x", "age_y_baseline" = "age_y_baseline.x", "HH_ID" = "HH_ID.x", "date" = "date.x" , "month_year" = "month_year_cfm") %>%
  rbind(hh_matches)


write.csv(skinny_OB, "data/generated_tables/skinny_OB.csv", row.names = FALSE)
write.csv(human_csp_file, "data/generated_tables/full_OB.csv", row.names = FALSE)

##I want to generate dfs for all possible cut offs, so here is an attempt at that
all_skinny_dfs <- list()
for(i in 1:24) {
  monthly_limit <- i
  rcd_dates <- sick_pos_csp %>%
    select(reactive_id, month_year_cfm, sample_id_date) %>%
    filter(reactive_id %in% names(final_csp_match_event)) %>%
    arrange(sample_id_date) %>%
    group_by(month_year_cfm) %>%
    mutate(id = row_number()) %>%
    filter(id <= monthly_limit)
  
  #restrict list to limited reactive events
  final_csp_match_event_limited <- final_csp_match_event[rcd_dates$reactive_id]
  
  #unlist into one df
  hh_matches <- bind_rows(final_csp_match_event_limited, .id = "reactive_id")
  
  #combine sick visits and hh matches
  skinny_OB <- sick_pos_csp %>%
    filter(reactive_id %in% rcd_dates$reactive_id) %>%
    select(unq_memID.x, village_name.x, gender.x, age_cat_baseline.x, age_y_baseline.x, age_m_baseline.x, HH_ID.x, date.x, 18:319, sample_type, sample_name_dbs, haplotype_number, haplotype_reads, month_year_cfm, week, year, week_year, hh_id, reactive_id, reactive_event) %>%
    rename("unq_memID" = "unq_memID.x", "village_name" = "village_name.x", "gender" = "gender.x", "age_cat_baseline" = "age_cat_baseline.x", "age_m_baseline" = "age_m_baseline.x", "age_y_baseline" = "age_y_baseline.x", "HH_ID" = "HH_ID.x", "date" = "date.x" , "month_year" = "month_year_cfm") %>%
    rbind(hh_matches)
  
  all_skinny_dfs[[i]] <- skinny_OB
  write.csv(skinny_OB, paste0("data/generated_tables/OB_dfs/skinny_OB-", i, ".csv"), row.names = FALSE)
}

library(rlist)
list.save(all_skinny_dfs, "data/generated_tables/OB_dfs/all_skinny_dfs.rds")

