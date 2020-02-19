# Content coming

This is an R notebook describing the production of a cyano boosted regression tree model - eventually. The preamble has sosme exercises to warm up and visualise the extent of the data. There are some useful visualising routines including:

The polar_coords plot to view seasonality 

and the 

manysmall faceted dotplots

To Do: join the Martin-version of the voyage codes with the BPA-ID
get PAR
get Fe, for the IN16v03 and IN19v03 as well

```{r}

#install.packages('oce')

library(dada2)
library(oce)
library(tidyverse)
library(gridExtra)
library(tidyverse)

```
![Image description](link-to-image)

```{r}
petb<-read_tsv('~/MarineMicrobes Dropbox/Martin Ostrowski/Refined_provinces/global.syn.pro.june.2018.table')

petb.col<-read.table("~/MarineMicrobes Dropbox/Martin Ostrowski/Refined_provinces/petb2018/petb_col_codes2.txt", h=T, sep='\t')

shape<-readRDS('~/MarineMicrobes Dropbox/uniques/winona/fortified_shoreline_ggplot_models.RDS')

meta<-read_tsv('~/Dropbox/Refined_provinces/cyano.model.meta.csv/cyano.model.meta.csv')



#colnames(petb)[-1] %in% petb.col$clade

#jcols<-as.character(petb.col$col)
#names(jcols)<-petb.col$clade

petb.gather<- gather(petb, key='strain', value='abund', -c(Sampled))


petb$Sampled %in% contextual.long$code

goshship<-read_tsv('~/Dropbox/AMI/metadata/lat_lon/goship.txt')

in19v03<-read_tsv('~/DRopbox/AMI/metadata/IN19v03_latlonstn.txt')

indigo<-read_tsv('~/Dropbox/indigo_latlon.tsv')

amt<-read_tsv('~/DRopbox/Refined_provinces/combined/meta/fac_syn_A_MLD_v2.txt')

colnames(amt)<-c('sample', colnames(amt)[1:21])

PLot an overview of the data
```
```{r, fig.width=10}
ggplot(meta %>%  filter(depth < 200), aes(x=lon, y=lat)) + geom_point(aes(color=as.numeric(temp), size=depth)) + scale_color_gradientn(colors=oce.colorsTemperature(60)) + coord_fixed(1) + geom_polygon(data = fortify(shape ), #%>% filter(between(lat, -65,0), between(long, 90,180))
               aes(x=long, y = lat, group = group), 
               color =  NA, fill = "grey80", show.legend = F) +   geom_point(data=amt, aes(x=Longitude, y=Latitude, color=Temp))+ geom_point(data=in19v03, aes(x=lon, y=lat*-1))+  geom_point(data=goshship, aes(x=Longs+360, y=Lats))+ geom_point(data=indigo, aes(x=lon, y=lat))+
  coord_fixed(1.1) +  theme_classic(base_size=16)+
  theme(legend.key =  element_rect(fill = "NA"))  + labs(title='AMI + INDIGO + AMT15,18,19,22,23 + TARA Oceans petB Cyano. To Add, biosope, Ambition, boum, prosope99') + coord_map(projection = 'mollweide')

?mapproj::mapproject
ggsave('~/AMI.Tara.petB.temp.pdf', width=20, height=12)
```






    ggplot(cyano31.densities %>% filter(!is.na(subclade), depth_m < 100), 
    aes(x = `temperature_ctd_its-90,_deg_c`, y = subclade, height = scaled, fill=subclade, color=subclade)) +
    theme_bw(base_size=9) +
    geom_density_ridges( scale=2,stat = "identity", alpha=0.7, lwd=0.1) + 
    labs(y='sub clade', x='Temperature (˚C)') + theme(legend.position='none') + 
    scale_color_manual(values=jcols, name='subclade')+
    scale_fill_manual(values=jcols, name='subclade')  + 
  theme(strip.text.y = element_text(angle=0), panel.spacing=unit(0.0, 'cm')) + 
  scale_x_continuous(limits=c(5,35)) + 
  facet_grid(genus+depth_m ~.,scales='free', space='free')
