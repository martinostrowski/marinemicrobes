
#Data source: figshare

#example: Bacterial OTUs in the surface (<= 10m depth) samples at the seven IMOS National Reference Stations

library(plyr)

bacteria<-read.table("~/AMMBI_B16S_zotus_table_Silvav132_datarelease_20032018.txt", h=T, sep='\t', comment.char="")

metadata<-read.table("~/Table_S1_230318.txt", h=T, sep='\t', row.names=1, quote='\"') 

#extract data from the surface samples only

surface_samples<-row.names(metadata[metadata$Depth_m<=10,])
bact_surf<-bacteria[,surface_samples]

#aggregate on Silva Taxonomy and export

bact_surf$taxcode<-bacteria$Taxonomy_Silvav132
bact_surf_agg<-ddply(bact_surf, "taxcode", numcolwise(sum) )
write.table(bact_surf_agg, "~/Bacteria_aggregated_by_taxa.txt", sep='\t', quote=F)

#The 'unknown' row in causes some formatting issue, export, replace "" with unknown and reimport
#Transpose and aggregate on station

bact_surf_agg<-read.table("~/Bacteria_aggregated_by_taxa.txt", sep='\t', h=T, row.names=1)
bact_surf_aggt<-t(bact_surf_agg)
bact_surf_aggt<-as.data.frame(bact_surf_aggt)

stations<-metadata[surface_samples,10]
bact_surf_aggt$stations<-stations

bact_summary<-ddply(bact_surf_aggt, "stations", numcolwise(sum))

#Create colorblocks based on 5 groups
##Colormixing to obtain unique colors was carried out with [iwanthiue](http://tools.medialab.sciences-po.fr/iwanthue/)
#This can also be done within r using the iwanthue function or with RColorBrewer

row.names(bact_summary)<-bact_summary[,1]
plotting.matrix<-t(as.matrix(bact_summary[,-1]))

normalized.matrix<-plotting.matrix/rep(colSums(plotting.matrix), each=nrow(plotting.matrix))

barplot(normalized.matrix, col=rainbow(1320))

sites<-c("DAR","KAI","MAI","NSI","PHB","ROT", "YON")

#import colours prepared separately, one colour for each taxa
rolled_col<-read.table("~/rolled_col.txt", h=T, sep='\t')

cyanot<-normalized.matrix[c(323:353),]
alphat<-normalized.matrix[c(551:879),]
deltat<-normalized.matrix[c(880:930),]
gammat<-normalized.matrix[c(931:1265),]
othert<-normalized.matrix[c(1:322, 354:550, 1266:1319),]

rolled_data<-rbind(alphat,cyanot, gammat, deltat, othert)
rolled_data<-as.data.frame(rolled_data)

#rolled_col<-c(rev(alpha),cyano, gamma, delta, other, other2, other3)
#rolled_col<-as.data.frame(rolled_col)
#row.names(rolled_col)<-row.names(rolled_data)

#reorder the stations by latitude

rolled_data<-rolled_data[,c(3,2,5,6,4,7,1)]
sites<-sites[c(3,2,5,6,4,7,1)]

#Create a legend for the top 30 taxa by abundance across the dataset
top30<-as.data.frame(names(rev(sort(rowSums(normalized.matrix)))[1:30]))

names(top30)<-"taxa"
row.names(top30)<-top30[,1]
top30_col<-merge(top30, rolled_col, by=0)

#Use XPD and fine tune margins to make this plot clear. Refine the order of the top 30 taxa

par(mar=c(8.1,4.1,2.1,30), xpd=TRUE)
barplot(as.matrix(rolled_data), col=as.character(rolled_col$rolled_col), las=2, border=NA)
legend("right", legend=rev(top30_col$tax), fill=as.character(rev(top30_col$rolled_col)), cex=0.4, inset=c(-1.8,0))
