# Code used to produce metadata maps for Varkey, Mazard et al., "Spatial and temporal shifts in phytoplankton community structure in the Sydney Harbour estuary" submitted to Marine and Freshwater Research 
#
# The data used was first reported in Jeffries, T., Schmitz Fontes, M., Harrison, D., Van-Dongen-Vogels, V., Eyre, B., Ralph, P. and Seymour, J. (2016), 'Bacterioplankton dynamics within a large anthropogenically impacted urban estuary', Frontiers in Microbiology, vol 6 .
# Authors: Martin Ostrowski and Deepa varkey
# email: martin.ostrowski@mq.edu.au
# website: https://github.com/martinostrowski/marinemicrobes
# May 2018; Last revision: 09-May-2018
# high resolution shapefiles for NSW estuaries and ocean ecosystems were downloaded from http://data.environment.nsw.gov.au/dataset


install.packages('mapdata')
library(mapdata)
require(rgdal)
library(oce)

##prepare maps, read in the data and extract the data for the relevant estuaries


shape <- readOGR(dsn = "~/Dropbox/Syd_Harbour/Estuaries/Data/", layer = "Estuaries")
shapeocean <- readOGR(dsn = "~/Dropbox/Syd_Harbour/OceanEcosystems2002/Data/", layer = "OceanEcosystems2002")

my.estuaries<-c("PORT JACKSON HARBOUR","PARRAMATTA RIVER ESTUARY","LANE COVE RIVER ESTUARY","MIDDLE HARBOUR CREEK ESTUARY","PARRAMATTA RIVER FRESHWATER")
my.estuary.shapes<-shape[shape$NAMETYPE %in% my.estuaries,]

quartz(height=20, width=12) # set window dimensions on OSX

mycol<-"grey"
plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)

## load metadata

sepenv<-read.table("~/Dropbox/Syd_Harbour/HarbourdataforMartin.txt", h=T, sep='\t')


# prepare colour scales for each variable to be plotted


sepenv$temp <- oceColorsTemperature(100)[as.numeric(cut(sepenv$Temp.Feb,breaks = 100))]
sepenv$sal <- oceColorsSalinity(100)[as.numeric(cut(sepenv$Sal.Feb,breaks = 100))]
sepenv$peeash <- oceColorsPalette(100)[as.numeric(cut(sepenv$pH,breaks = 100))]
sepenv$teess <- oceColorsPalette(100)[as.numeric(cut(sepenv$TSS.Final,breaks = 100))]
sepenv$chl <- oceColorsChlorophyll(100)[as.numeric(cut(sepenv$Chl.Post.Cal,breaks = 100))]
sepenv$nocs <- oceColorsPalette(100)[as.numeric(cut(sepenv$NOx,breaks = 100))]
sepenv$amm <- oceColorsPalette(100)[as.numeric(cut(sepenv$NH4,breaks = 100))]
sepenv$totaln<- oceColorsPalette(100)[as.numeric(cut(sepenv$TN,breaks = 100))]
sepenv$phos <- oceColorsPalette(100)[as.numeric(cut(sepenv$PO4,breaks = 100))]
sepenv$sil <- oceColorsPalette(100)[as.numeric(cut(sepenv$Si,breaks = 100))]

# colour bars produced with the following function # thanks to JOhn Colby <http://www.colbyimaging.com/wiki/statistics/color-bars>

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
    scale = (length(lut)-1)/(max-min)

    #dev.new(width=1.75, height=5)
    plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab=title, main='', cex.lab=1.5, font=2)
    axis(2, ticks, las=1)
    for (i in 1:(length(lut)-1)) {
    	y = (i-1)/scale + min
    	rect(0,y,10,y+1/scale, col=lut[i], border=NA)
    }	
}

###Panel 1

layout(matrix(seq(1,15,1), nrow=5, ncol=3, byrow=T), widths=c(8, 32,32))

sepenv$temp <- oceColorsTemperature(100)[as.numeric(cut(sepenv$Temp.Feb,breaks = 100))]
par(mar=c(2,7,0,0), mgp=c(3.5,1,0))
color.bar(as.character(oceColorsTemperature(100)), min(round(sepenv$Temp.Feb,0)), max(round(sepenv$Temp.Feb,0)), title="Temperature (˚C)")

par(mar=c(2,2,1,2))

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(temp), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(temp), cex=2.5)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0))

sepenv$sal <- oceColorsSalinity(100)[as.numeric(cut(sepenv$Sal.Feb,breaks = 100))]
color.bar(as.character(oceColorsSalinity(100)), min(round(sepenv$Sal.Feb,0)), max(round(sepenv$Sal.Feb,0)), title="Salinity (psu)")

par(mar=c(2,2,1,2))

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(sal), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(sal), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0))

sepenv$peeash <- oceColorsPalette(100)[as.numeric(cut(sepenv$pH,breaks = 100))]
color.bar(as.character(oceColorsPalette(100)), min(round(sepenv$pH,0)), max(round(sepenv$pH,0)), title="pH")

par(mar=c(2,2,1,2))

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(peeash), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(peeash), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0))

sepenv$teess <- oceColorsPalette(100)[as.numeric(cut(sepenv$TSS.Final,breaks = 100))]
color.bar(as.character(oceColorsPalette(100)), min(round(sepenv$TSS.Final,0)), max(round(sepenv$TSS.Final,0)), title=expression(paste("TSS (mg.",L^-1,")")))

par(mar=c(2,2,1,2))

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(teess), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(teess), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0)
color.bar(as.character(oceColorsChlorophyll(100)), min(round(sepenv$Chl.Post.Cal,0)), max(round(sepenv$Chl.Post.Cal,0)), title=expression(paste("Chlorophyll (µg.",L^-1,")")))
 
par(mar=c(2,2,1,2))
plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(chl), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(1, at=c(151.01, 151.14, 151.28))

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(chl), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(1, at=c(151.010, 151.147, 151.284))
axis(4, at=c(-33.753, -33.878))


###Panel 2

layout(matrix(seq(1,15,1), nrow=5, ncol=3, byrow=T), widths=c(8, 32,32))


par(mar=c(2,7,0,0), mgp=c(3.5,1,0))
color.bar(as.character(oceColorsPalette(100)), min(round(sepenv$NOx,0)), max(round(sepenv$NOx,0)), title=expression(paste("Nitrate/Nitrite (µg.",L^-1,")")))

par(mar=c(2,2,1,2))
plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(nocs), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(nocs), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0))
color.bar(as.character(oceColorsPalette(100)), min(round(sepenv$NH4,0)), max(round(sepenv$NH4,0)), title=expression(paste("Ammonium (µg.",L^-1,")")))

par(mar=c(2,2,1,2))
plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(amm), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(amm), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0))
color.bar(as.character(oceColorsPalette(100)), min(round(sepenv$TN,0)), max(round(sepenv$TN,0)), title=expression(paste("Total N (µg.",L^-1,")")))

par(mar=c(2,2,1,2))
plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(totaln), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(totaln), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0))
color.bar(as.character(oceColorsPalette(100)), min(round(sepenv$PO4,0)), max(round(sepenv$PO4,0)), title=expression(paste("Phosphate (µg.",L^-1,")")))

par(mar=c(2,2,1,2))
plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(phos), cex=2.5, xaxt='n', yaxt='n', asp=1)

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(phos), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(4, at=c(-33.753, -33.878))

par(mar=c(2,7,0,0), mgp=c(3.5,1,0))
color.bar(as.character(oceColorsPalette(100)), min(round(sepenv$Si,0)), max(round(sepenv$Si,0)), title=expression(paste("Silicate (µg.",L^-1,")")))

par(mar=c(2,2,1,2))
plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv[1:27,], lat~lon, pch=16, cex=2.6)
points(data=sepenv[1:27,], lat~lon, pch=16, col=as.character(sil), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(1, at=c(151.010, 151.147, 151.284))

plot(my.estuary.shapes, col=mycol, xlim=c(151,151.3), lwd=0.25, border=NA)
plot(shapeocean, add=T, lwd=0.5, col="lightblue", border=NA)
points(data=sepenv, lat~lon, pch=16, cex=2.6)
points(data=sepenv[28:57,], lat~lon, pch=16, col=as.character(sil), cex=2.5, xaxt='n', yaxt='n', asp=1)
axis(1, at=c(151.010, 151.147, 151.284))
axis(4, at=c(-33.753, -33.878))



