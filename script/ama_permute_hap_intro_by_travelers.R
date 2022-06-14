# ----------------------------------------- #
#     permutation test for introductions    #
#             Turkana EMBATALK              #
#               ama1 target                 #
#              February 2021                #
#               C. Markwalter               #
# ----------------------------------------- #

#I am trying to find evidence to determine whether the haplotype introductions we observe in the dataset are simply artifacts or if they are "real." The criteria for a haplotype introduction are: (1) the haplotype has not appeared in individuals who report no travel in or before the month that the traveler haplotype was identified, and (2) the haplotype appears in individuals who report no travel in the months after the traveler haplotype is first seen. Here, I am defining an "introduction score" that is simply the proportion of possible traveler haplotype introductions that meet the criteria of an introduction. Using a permutation test to determine the "null" distribution of the haplotype score by permuting the assignement of whether someone traveled or not.

####----loading libraries----####
library(tidyverse)
library(lubridate)
options(dplyr.summarise.inform = FALSE)


####----reading in data----####
nodes_pcrpos <- read.csv("data/generated_tables/turkana_full_dataset_finaldates.csv") %>%
  filter(pf_pcr_infection_status == "positive", case_type == "facility" | case_type == "household" | (case_type == "bus/plane" & rdt == "negative"))
ama_merged <- readRDS("data/generated_tables/ama_merged.rds") %>%
  filter(study_id %in% nodes_pcrpos$study_id) %>%
  left_join(nodes_pcrpos %>% select(study_id, date_final, reported_travel), by = "study_id")

ama_merged$date_final <- ymd(ama_merged$date_final)

####----defining function that outputs introductions score----####

intro.score.hap <- function(ama_merged){
  hapdate <- ama_merged %>%
    select(date_final, reported_travel, haplotype) %>%
    distinct() %>%
    arrange(date_final) 
  first_date <- hapdate %>%
    filter(date_final == min(date_final))
  hapdate <- hapdate %>%
    filter(date_final > min(date_final))
  haplist <- unique(ama_merged$haplotype)
  introlist <- ""
  results <- data.frame(score = double(), introlist = character())
  count = 0
  sum = 0
  for (h in 1:length(haplist)) {
    count = count + 1
    intro <- hapdate %>%
      filter(haplotype == haplist[h]) %>%
      filter(date_final == min(date_final))
    
    later <- hapdate %>%
      filter(haplotype == haplist[h]) %>%
      filter(date_final > min(date_final))
    
    if("Yes" %in% intro$reported_travel & !"No" %in% intro$reported_travel & "No" %in% later$reported_travel & !"No" %in% first_date$reported_travel[first_date$haplotype == haplist[h]]) {
      sum = sum + 1
      introlist <- paste0(introlist, haplist[h], ", ")
    }
  }
  score <- sum/count
  results <- rbind(results, data.frame(score, introlist))
}



####----intro.score for our dataset----####



embatalk.score.hap <- intro.score.hap(ama_merged)

saveRDS(embatalk.score.hap, "data/generated_tables/ama.embatalk.score.hap.RDS")

####----permuting traveler assignments----####


perm.test.hap <- function(ama_merged) {
  temp <- ama_merged %>% group_by(study_id) %>% summarise(reported_travel = first(reported_travel))
  temp$reported_travel <- sample(temp$reported_travel)
  ama_temp <- merge(ama_merged %>% select(-reported_travel), temp)
  df <- intro.score.hap(ama_temp)
  return(df)
}

####----generating "null" distribution----####

set.seed(1)
many.perm.hap <- data.frame(t(replicate(1000, perm.test.hap(ama_merged), simplify = TRUE)))
percentile.hap <- ecdf(unlist(many.perm.hap$score))
percentile.hap(embatalk.score.hap$score)

saveRDS(many.perm.hap, "data/generated_tables/many.perm.hap_ama.RDS")

ggplot(many.perm.hap, aes(x = unlist(score))) +
  geom_histogram(binwidth = 0.0135, fill = "black", alpha = 0.5) +
  #geom_histogram(data = embatalk.dist, aes(x = embatalk.dist), fill = "darkred", alpha = 0.5, binwidth = 0.007) +
  geom_vline(aes(xintercept = embatalk.score.hap$score), color = "darkred") +
  #annotate(geom = "text", x = 11, y = 175, label = "embatalk", color = "red") +
  labs(x = "Proportion of haplotypes meeting importation criteria", y = "Frequency") +
  theme_bw()


ggplot(many.perm.hap, aes(x = unlist(score))) +
  geom_density(aes(y = ..scaled..),fill = "black", bw = 0.01, alpha = 0.5) +
  #geom_histogram(data = embatalk.dist, aes(x = embatalk.dist), fill = "darkred", alpha = 0.5, binwidth = 0.007) +
  geom_vline(aes(xintercept = embatalk.score.hap$score), color = "darkred") +
  #annotate(geom = "text", x = 11, y = 175, label = "embatalk", color = "red") +
  labs(x = "Proportion of haplotypes meeting importation criteria", y = "Frequency") +
  theme_bw()


ggsave("figures/ama_imports_permutation.png")



###---epi curves for imported haplotypes---###

imported <- ama_merged %>%
  filter(str_detect(embatalk.score.hap$introlist, haplotype))

ggplot(imported %>% filter(traveled %in% c("Yes", "No")), aes(x = ymd(date_final), fill = traveled)) +
  geom_dotplot(binwidth = 14, alpha = 0.5, dotsize = 1.5) +
  scale_y_continuous(NULL, breaks = NULL) +
  facet_wrap(~haplotype, ncol = 4) +
  labs(fill = "Individual reported travel?", x = "Date", y = "Cases") +
  theme(axis.text.x = element_text(angle = 90))

