## Script to create maps for Embatalk importation paper: map PfPR_2-10 and highlighting origins of travel
## Last updated by Hannah Meredith on March 26, 2022

library("raster")
library(sp)
library(rgeos)
library("rgdal")
library(gpclib)
library("ggplot2")
library("dplyr")
library(ggpubr)

# import data (setwd to Mapping)
Pf <- raster("Mapping/2020_Global_PfPR_KEN_2018.tiff") ## If you want the PfPR_1-99 instead, you can download it from https://malariaatlas.org/trends/country/KEN and import here
Kenya.counties.shp <- readOGR("Mapping/ken_admbnda_adm1_iebc_20191031.shp") 
trips <- readRDS("data/generated_tables/trips.RDS") %>%
  separate(address, into = c("county", "country"), sep = ", ") %>%
  filter(country == "KENYA")

trips$county[trips$county == "ELGEYO MARAKWET"] <- "ELGEYO-MARAKWET"
trips$county[trips$county == "MARASABIT"] <- "MARSABIT"
trips$county[trips$county == "NANDI COUNTY"] <- "NANDI"

trips <- trips %>%
  group_by(county) %>%
  tally()

trip.origins <- trips$county
trip.origins <- subset(Kenya.counties.shp@data, toupper(ADM1_EN) %in% trip.origins)[ , c("ADM1_PCODE", "ADM1_EN")] %>%
  mutate(county = toupper(ADM1_EN)) %>%
  left_join(trips) %>%
  dplyr::select(-county)

# turn raster and shapefiles into dataframes for plotting
Pf.points <- rasterToPoints(Pf, spatial = TRUE)
Pf.points.df <- data.frame(Pf.points)
Kenya.counties.f <- fortify(Kenya.counties.shp, region = "ADM1_PCODE", name = "ADM1_EN") # use fortify to create coordinates for borders, dataframe returns the info on the data portion of the shp file (i.e. adm name and id)
Kenya.counties.f$origin = ifelse(Kenya.counties.f$id %in% trip.origins$ADM1_PCODE, "origin", NA)  # ID origins for mapping
Kenya.counties.f <- Kenya.counties.f %>%
  left_join(trip.origins %>%
              rename("id" = "ADM1_PCODE") %>%
              dplyr::select(id, n))

Kenya.counties.f$n[is.na(Kenya.counties.f$n)] <- 0

# plot layers
prevalence_map <- ggplot()+
  geom_raster(data = Pf.points.df, aes(x = x, y = y, fill = X2020_Global_PfPR_KEN_2018 + 0.01))+  ## using log scale helps bring colors out more - would suggest changing the scalebar labels to be normal.
  scale_fill_viridis_c(option = "plasma", trans = "log10", limits = c(0.01,1.01), breaks = c(0.01,0.11,1.01), labels = c(0,.1,1))+  ## you may want to play around with the color scale - maybe make 0 data areas grey, no malaria = white, and the rest some color gradient
  geom_polygon(data = Kenya.counties.f, aes(x = long, y = lat, group = group), ## Plot shapefile to get county boundaries
               fill = NA, color = "black", size = 0.25, alpha = 0.5)+
  geom_polygon(data = Kenya.counties.f %>% filter(id == "KE023"), aes(x = long, y = lat, group = group), fill = NA, color = "orange", size = 1) +
  labs(fill = "PfPR 2-10")+
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

prevalence_map

saveRDS(prevalence_map, "Mapping/prevalence_map_fig.RDS")

## Note the white area in upper Turkana is due to the region being contested land. I think MAP must not have that area as part of Kenya. 

num_trips <- ggplot()+
    geom_polygon(data = Kenya.counties.f, aes(x = long, y = lat, group = group, fill = n + 1), color = "black", size = 0.25)+
  geom_polygon(data = Kenya.counties.f %>% filter(id == "KE023"), aes(x = long, y = lat, group = group), fill = NA, color = "orange", size = 1) +
    scale_fill_viridis_c(trans = "log10", breaks = c(1,11,101), labels = c(0,10,100)) +
    labs(fill = "Number\nof trips")+
    theme(panel.background = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "left")
num_trips


saveRDS(num_trips, "Mapping/num_trips_fig.RDS")
ggsave(plot = num_trips, filename = "manuscript/figures/num_trips_map.tiff", height = 6, width = 6)

#get mean PfPR per county
county.PfPR <- raster::extract(Pf, Kenya.counties.shp)
county.PfPR.mean <- unlist(lapply(county.PfPR, FUN = mean, na.rm = TRUE))

Kenya.counties.shp@data$id <- Kenya.counties.shp@data$ADM1_PCODE
Kenya.counties.shp@data <- data.frame(Kenya.counties.shp@data, mPfPR = county.PfPR.mean) %>%
  left_join(trip.origins %>%
              rename("id" = "ADM1_PCODE") %>%
              select(id, n))
Kenya.counties.shp@data$n[is.na(Kenya.counties.shp@data$n)] <- 0
Kenya.counties.f <- fortify(Kenya.counties.shp, region = "ADM1_PCODE", name = "ADM1_EN")
Kenya.counties.f <- left_join(Kenya.counties.f, Kenya.counties.shp@data, by = "id")


m_prevalence_map <- ggplot()+
  geom_polygon(data = Kenya.counties.f, aes(x = long, y = lat, group = group, fill = mPfPR), color = "black", size = 0.25)+
  scale_fill_viridis_c(option = "plasma")+
  geom_polygon(data = Kenya.counties.f %>% filter(id == "KE023"), aes(x = long, y = lat, group = group), fill = NA, color = "orange", size = 1) +
  labs(fill = "PfPR 2-10")+
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

m_prevalence_map

saveRDS(m_prevalence_map, "Mapping/m_prevalence_map_fig.RDS")


prev_trips_scatter <- Kenya.counties.f %>%
  select(id, ADM1_EN, mPfPR, n) %>%
  unique() %>%
  filter(n > 0) %>%
  ggplot(aes(x = n, y = mPfPR)) +
  scale_x_log10() +
  geom_point() +
  theme_bw()

prev_trips_scatter

ggarrange(num_trips,
          m_prevalence_map,
          nrow = 1)

summary <- Kenya.counties.f %>%
  select(id, ADM1_EN, mPfPR, n) %>%
  unique() %>%
  filter(n > 0)


saveRDS(summary, "data/generated_tables/county_mPfPR.RDS") 


