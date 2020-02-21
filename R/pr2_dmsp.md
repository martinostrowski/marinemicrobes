---
title: "R Notebook"
output: html_notebook
---

This notebook descsribes
plotting the 18S data from NRS

checking that all data are there
extracting the phytoplankton (according to PR2 definition)
taking a first look at the data
create a factor for high, low DMSP producers

Taxonomy was assigned on the unique zOTUs with > 9 sequences see ()
```
library(tidyverse)

pr2.tax<- read_csv('~/pelagic.PR2.taxonomy.noBS.csv')
```

```
phb<-read_csv('PHB.subset.csv')
mai<-read_csv('MAI.subset.csv')
nsi<-read_csv('NSI.subset.csv')


phb<-phb[, wanted]
nsi<-nsi[, wanted]
mai<-mai[, wanted]
```

mai<-read_csv('MAI.subset.csv')
mai<-read_csv('MAI.subset.csv')
mai<-read_csv('MAI.subset.csv')
mai<-read_csv('MAI.subset.csv')
yon<-yon[, wanted]
rot<-rot[, wanted]
kai<-kai[, wanted]
dar<-dar[, wanted]

 wanted<-c("zOTU","code", "abund", "year", "month", "date_sampled", "nrs_location_code_voyage_code", "month.abb", "nitrate_nitrite_μmol/l", "silicate_μmol/l", "temperature_ctd_its-90,_deg_c", "phosphate_μmol/l", "salinity_ctd_psu", "ammonium_μmol/l", "bottom_depth", "chlorophyll_a_μg/l", "date_sampled", "latitude_decimal_degrees", "longitude_decimal_degrees")

```
nrs<-bind_rows(phb, mai, nsi)

nrs.tax<- nrs %>%  left_join(pr2.tax, 'zOTU')

nrs.tax$month.abb <- factor(month.abb[as.integer(nrs.tax$month)], levels=c("Jan", "Feb","Mar", "Apr","May",  "Jun", "Jul","Aug","Sep","Oct","Nov", "Dec" ));

#nrs.tax <- nrs.tax %>% unite('FG', c(Family, Genus), sep='_')
#nrs.tax <- nrs.tax %>% unite('CF', c(Class, Family), sep='_')

write_csv(nrs.tax, "east.coast.nrs.pr2.18S.csv")

ami_treemap(nrs.tax)
```


```


ggplot(nrs.tax, aes(x=month.abb, y=abund, fill=fct_lump(FG, 49))) + geom_bar(stat='identity', position='fill') + scale_fill_manual(values=mycol50) + facet_wrap(. ~ nrs_location_code_voyage_code)
```

```
nrs.photo <- nrs.tax %>% filter((Division %in% c("Chlorophyta", "Dinophyta", "Cryptophyta", "Haptophyta", "Ochrophyta", "Cercozoa")) & !(Class %in% c("Syndiniales", "Sarcomonadea")))

```


```
ggplot(nrs.photo, aes(x=month.abb, y=abund, fill=fct_lump(FG, 49))) + geom_bar(stat='identity', position='fill') + scale_fill_manual(values=mycol50) + facet_grid(nrs_location_code_voyage_code ~ year) + theme_mo()
```




```
dmspHigh<-read_csv('dmsp_high_pr2.tax')
dmspLow<-read_csv('dmsp_low_pr2.tax')

dmsp<-rbind(dmspHigh, dmspLow)

colnames(dmsp)<-c("tax.Kingdom","tax.Supergroup","tax.Division","tax.Class","tax.Order","tax.Family","tax.Genus","tax.Species","boot.Kingdom","boot.Supergroup","boot.Division","boot.Class","boot.Order","boot.Family","boot.Genus","boot.Species","seq","dmsp")

write_csv(dmsp,'~/dmsp_factor_18s_pr2.csv')

nrs.photo<- nrs.photo %>% left_join(dmsp, c('Genus' ='tax.Genus'))

table(nrs.photo$Family %in% dmsp$tax.Genus)

table(nrs.photo$dmsp)
```


```
#{r, fig.width=8, fig.height=4}
ggplot(nrs.photo, aes(x=month.abb, y=abund, fill=dmsp)) + geom_bar(stat='identity', position='fill') + scale_fill_manual(values=c( 'red', 'blue')) + facet_grid(nrs_location_code_voyage_code ~ year) + theme_mo()
```
```

```
ami_treemap<-function(pr2) {
    # utilised with kind acknowldegement of PR2 and the taxomap team 
    # Define the levels
    level1 = "Class"
    level2 = "Genus"
    #color = 'color'
    # Group
    pr2_class <- pr2 %>% group_by_(level1, level2) %>% summarise(sequence_number = sum(abund))
    
    # Do a simple treemap
    treemap::treemap(pr2_class, index = c(level1, level2), vSize = "sequence_number", 
                      asp = 1, lowerbound.cex.labels = 0.2,
                       vColor = mycol40,format.legend = list(scientific = FALSE, big.mark = " "),fontsize.labels=c(18,12),                # size of labels. Give the size per level of aggregation: size for group, size for subgroup, sub-subgroups...
                     #palette = mycols200
                     fontcolor.labels=c("white","orange"),    # Color of labels
                     fontface.labels=c(2,1),                  # Font of labels: 1,2,3,4 for normal, bold, italic, bold-italic...
                     bg.labels=c("transparent"),              # Background color of labels
                     align.labels=list(
                         c("center", "center"), 
                         c("left", "top")
                     ),                                   # Where to place labels in the rectangle?
                     overlap.labels=1,                      # number between 0 and 1 that determines the tolerance of the overlap between labels. 0 means that labels of lower levels are not printed if higher level labels overlap, 1  means that labels are always printed. In-between values, for instance the default value .5, means that lower level labels are printed if other labels do not overlap with more than .5  times their area size.
                     inflate.labels=F,                        # If true, labels are bigger when rectangle is bigger.
                     aspRatio=3,
                     #title=paste(sample_type),
                     lwds=(0.5)
                    
                     
    )
}

mycol50<-c("#ff6943",
"#8c3918",
"#c65000",
"#ff9e61",
"#8f6500",
"#bd9c00",
"#ffe966",
"#fff29e",
"#aeb600",
"#c9c88e",
"#596800",
"#d3ff3c",
"#edffa8",
"#73b900",
"#bcff89",
"#5ed100",
"#2b6c00",
"#45ac00",
"#81ff5e",
"#019417",
"#5f895d",
"#97ffa1",
"#65ff8c",
"#01db68",
"#019247",
"#007b43",
"#78ffc7",
"#8ee4d2",
"#02e0d0",
"#38fff8",
"#6fc4ff",
"#0175ac",
"#0177cd",
"#70abff",
"#6073a6",
"#7c97ff",
"#5582ff",
"#014bc7",
"#2a5bfe",
"#7768ff",
"#623fe2",
"#5d4396",
"#d895ff",
"#8a2190",
"#94009b",
"#e943ea",
"#ffa3ed",
"#ff7ce9",
"#b10070",
"#e26391")
```

```
theme_mo<-function (base_size = 11, base_family = "") 
{
    blue <- "#2c3e50"
    green <- "#18BC9C"
    white <- "#FFFFFF"
    grey <- "grey80"
    theme_grey(base_size = base_size, base_family = base_family) %+replace% 
        theme(line = element_line(colour = blue, size = 0.5, 
            linetype = 1, lineend = "butt"), rect = element_rect(fill = white, 
            colour = blue, size = 0.5, linetype = 1), text = element_text(family = base_family, 
            face = "plain", colour = blue, size = base_size, 
            lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0, 
            margin = margin(), debug = FALSE), axis.line = element_blank(), 
            axis.text = element_text(size = rel(0.8)), axis.ticks = element_line(color = grey, 
                size = rel(1/3)), axis.title = element_text(size = rel(1)), 
            panel.background = element_rect(fill = white, color = NA), 
            panel.border = element_rect(fill = NA, size = rel(1/2), 
                color = blue), panel.grid.major = element_line(color = grey, 
                size = rel(1/3)), panel.grid.minor = element_line(color = NA, 
                size = rel(1/3)), panel.grid.minor.x = element_blank(), 
            panel.spacing = unit(0.02, "cm"), legend.key = element_rect(fill = white, 
                color = NA), legend.position = "bottom", strip.background = element_rect(fill = blue, 
                color = blue), strip.text = element_text(color = white, 
                size = rel(0.8)), plot.title = element_text(size = rel(1.2), 
                hjust = 0, margin = margin(t = 0, r = 0, b = 4, 
                  l = 0, unit = "pt")), plot.subtitle = element_text(size = rel(1.1), 
                hjust = 0, margin = margin(t = 0, r = 0, b = 3, 
                  l = 0, unit = "pt")), complete = TRUE)
}
```

