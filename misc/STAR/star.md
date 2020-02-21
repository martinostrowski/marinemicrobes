
# State of Australia's Marine Microbial Communities:2019

#### Author: Martin Ostrowski
#### Contact: martin.ostrowski@uts.edu.au

This is an R notebook providing details of the data and code for the analyses presented in the State and Trends Report for Australia's Oceans and seas. The examples set out here provide a guide for accessing the processed Bacterial 16S data from the Australian Microbiome Initiative portal and observing the relative abundance of microbial taxa at different taxonomic levels.

1. Load the contextual data and plot a map and time series of the samples collected for 'pelagic' and 'coastal' datasets
2. Load the processed amplicon table 
3. Plot alpha and beta-diversity for the National Reference Station Time Series
4. Time series plots for specific taxonomic ranks

Load the package libraries

```r
library(tidyverse)
library(oce) #
```
1. Load the contextual metadata and streamline some formatting

```r
contextual.long <- read_csv('~/MarineMicrobes Dropbox/uniques/BRT_2019/input/contextual_marine_201907.csv')

colnames(contextual.long) <- gsub(" ", "_", colnames(contextual.long), fixed = TRUE) #remove spaces
colnames(contextual.long) <- gsub("[", "", colnames(contextual.long), fixed = TRUE) #remove brackets
colnames(contextual.long) <- gsub("]", "", colnames(contextual.long), fixed = TRUE)
colnames(contextual.long) <- gsub("/", "_per_", colnames(contextual.long), fixed = TRUE)
colnames(contextual.long) <- tolower(colnames(contextual.long)) #all lowercase

contextual.long <- contextual.long  %>%
  separate (sample_id, c("num", "code"), sep='/')
```

separate the YYYY-MM-DD format date into year, month and day. Add levels for months so they appear in the correct order and convert some flagged data points in the nutrient data to NAs

```r
contextual.long <- contextual.long %>%
    separate(date_sampled, c('year','month','day'), sep='-', remove=F) 

contextual.long$month.abb <- factor(month.abb[as.integer(contextual.long$month)], 
    levels=c(month.abb[seq(1,12,1)])) 

contextual.long$nitrate_nitrite_μmol_per_l[contextual.long$nitrate_nitrite_μmol_per_l == -999.000]<-NA;
contextual.long$phosphate_μmol_per_l[contextual.long$phosphate_μmol_per_l == -999.000]<-NA;
contextual.long$salinity_ctd_psu[contextual.long$salinity_ctd_psu == -999.000]<-NA;
contextual.long$silicate_μmol_per_l[contextual.long$silicate_μmol_per_l == -999.0000]<-NA;
contextual.long<- contextual.long[contextual.long$salinity_ctd_psu > 2,]
```

To examine the extent of the data in time and space plot the sites where the the samples have been collected from and a timeline of sample collection

```r
shape<-readRDS('~/MarineMicrobes Dropbox/uniques/winona/fortified_shoreline_ggplot_models.RDS')

 ggplot() +   geom_polygon(data = fortify(shape %>% filter(between(long, 100, 160), between(lat,-60,0))), 
               aes(x=long, y = lat, group = group), 
               color =  NA, fill = "grey80") + 
 
  geom_count(data=contextual.long %>% 
  filter (depth_m < 25, between(longitude_decimal_degrees,90,180)), 
  aes(x=longitude_decimal_degrees, y=latitude_decimal_degrees, color=nrs_location_code_voyage_code, alpha=0.8)) +
  scale_color_manual(values=mycols310, name="Voyage or Location") + 
  coord_fixed(1.1) + 
  theme_bw() +
  labs(y='Latitide (degrees North)', x='Longitude (degress East)')
  ```
