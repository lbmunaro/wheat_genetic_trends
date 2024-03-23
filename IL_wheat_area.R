# IL wheat area ----

# Objective ----
# - Plot the wheat cultivated area per IL county

rm(list=objects()) # clean workspace

# Packages ----

library(tidyverse) # R packages for data science
library(sf) # Simple Features for R
library(raster) # Geographic Data Analysis and Modeling
library(ggspatial) # Spatial Data Framework for ggplot2
library(RColorBrewer) # ColorBrewer Palettes

# Data ----

## Locations ----
loc_coord <- data.frame(loc=c("Brownstown", "Carmi", "Neoga", "Ridgway", "St. Jacob", "Urbana"),
                        lat=c(38.9948827,38.08553,39.23274,37.79930,38.74573,40.05833),
                        long=c(-88.95961,-88.19441,-88.38207,-88.27799,-89.78045,-88.22937)) |>
  glimpse()

## Cultivated area
# Data for the states
states <- c("illinois")
highlight_states <- c("illinois")

# Get the map data for the specified states
map_data <- map_data("state") %>%
  filter(region %in% states)
county <- map_data("county") |>
  filter(region %in% states) |>
  mutate(subregion=str_replace_all(subregion," ",""),
         subregion=str_replace_all(subregion,"[.]",""))

wheat_ac <- read.csv("data/2022_fsa_acres_web_082222.csv") |>
  mutate_if(is.character,~tolower(.)) |>
  dplyr::filter(Crop%in%"wheat") |>
  mutate(Planted.Acres=str_replace(Planted.Acres,",",""),
         Planted.Acres=as.numeric(Planted.Acres)) |>
  group_by(State,County,Crop) |>
  summarise(Planted.Acres=sum(Planted.Acres,na.rm=T)) |>
  mutate(region=State,
         subregion=County,
         subregion=str_replace_all(subregion,
                                   c("dewitt"="de witt","dupage"="du page","st. clair"="st clair","dekalb"="de kalb" ))) |>
  mutate(subregion=str_replace_all(subregion," ",""),
         subregion=str_replace_all(subregion,"[.]",""))|>
  glimpse()

# Define custom breaks for the classes
breaks <- c(1, 1000, 5000, 10000, 25000, 50000, 100000)
labels <- c("< 1000", "1001 to 5000", "5001 to 10000", "10001 to 25000", "25001 to 50000", "> 50001")

county_ac <- county |>
  left_join(wheat_ac) |>
  mutate(Planted.Acres = ifelse(Planted.Acres == 0, NA, Planted.Acres)) |>
  mutate(classes = cut(Planted.Acres, breaks=breaks, labels = labels, include.lowest = TRUE, right = FALSE)) |>
  glimpse()

# Plot ----
ggplot() +
  geom_polygon(data = county_ac, aes(x = long, y = lat, group = group, fill = classes), color = "gray70") +
  scale_fill_brewer(name = "Planted acres", palette = "Oranges", na.value = "white", drop= TRUE, 
                    labels = c("< 1000", "1001 to 5000", "5001 to 10000", "10001 to 25000", "25001 to 50000", "> 50001", "Not estimated")) +
  geom_polygon(data = map_data, aes(x = long, y = lat, group = group), fill = "transparent", color = "gray50") +
  geom_point(data = loc_coord, aes(x = long, y = lat, group = NA), color = "#13294B", size = 1) +
  geom_point(data = loc_coord, aes(x = long, y = lat), color = "#13294B", size = 2, shape = 21) + 
  geom_text(data = loc_coord, aes(x = long, y = lat, group = NA, label = loc),
            angle = 0, vjust = 0.5, hjust = -0.1, color = "#13294B", size = 4) +
  coord_fixed(1.3) +
  theme(panel.grid = element_blank()) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        panel.background=element_rect(fill="white"),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10),
        legend.key=element_rect(fill="white",color=NA),
        legend.position = "left",
        plot.margin = margin(0,0,0,0, "in")) +
  xlim(range(loc_coord$long) + c(-2,1)) +
  annotation_north_arrow(location="bl",pad_x=unit(0,"in"),pad_y=unit(0.2,"in"),
                         style = north_arrow_nautical,
                         height = unit(0.5,"in"),width = unit(0.5,"in"))
ggsave("figures/IL_wheat_map.png", width = 6, height = 4, units = "in", dpi = 320)

# Save
save.image("data/IL_wheat_area.RData")