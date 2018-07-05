
#Authors: Martin Ostrowski, Deepa Varkey and Mark Brown
#This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook Notebook that describes a boosted regression tree modelling workflow, including an integrated cross validation step, a batch processing step to display the repsonse variable contributions on a phylogenetic tree and and a clustering-based validation/correlation of the stackedd community predictions against the hold-back and input data.
#The input data format consists of a site by factor and species dataframe. Response variables are listed in the first x columns ("site", "depth","temp", "sal","sil", "no3", "po4", "N:P", "strat", "daylength", "zotu1"..."zoutux")

library(dismo)
library(gbm)
library(vegan)
library(clustsig)

#####read in the model input

zotus1000<-read.table("input.table.noRTOT.csv", h=T, sep=',')

mai<-read.table("MAI_Historical_env through_082017_no_silicate.csv", sep=',', h=T)
phb<-read.table("PH_Historical_Env_through_082017.csv", sep=',', h=T)
rot<-read.table("ROT_Historical_Env_through_082017.csv", sep=',', h=T)

colnames(mai)<-c("Sample", "Day" ,"Month" ,"Year","Lat","Lon" ,"depth","DL","strat","temp","sal","nox" )
colnames(phb)<-c("Sample", "Day" ,"Month" ,"Year","Bottom.depth","Lat","Lon" , "depth","DL","strat","temp","sal","sil", "nox", "po4")
colnames(rot)<-c("Sample", "Day" ,"Month" ,"Year","Lat","Lon" , "depth","DL","strat","temp","sal", "nox","sil")

#Now do the predictions

#"depth"                    "temp"                     "sal"                      "nox"                      "DL"                       "strat"  

zotus.batch<-sqrt(zotus.batch)
zotus1000<-cbind(zotus1000[,c(1:12)], zotus.batch)


nzotu<-ncol(zotus1000)-12

my.zotus<-vector('list',nzotu)
for (i in 1:(nzotu)){
  my.zotus[[i]]<-colnames(zotus1000)[i+12]
}

sp.names<-vector('double', nzotu)
preds.mai<-matrix(nrow=nrow(mai), ncol=nzotu)
preds.phb<-matrix(nrow=nrow(phb), ncol=nzotu)
preds.rot<-matrix(nrow=nrow(rot), ncol=nzotu)
my.mods<-vector('list', nzotu)

  eval.data.out<-vector('double', nzotu)  
  cv<-vector('double', nzotu)
  cvse<-vector('double', nzotu)
  ntrees<-vector('double', nzotu)
  contributions<-matrix(nrow=6, ncol=nzotu)
  rownames(contributions)<-c("temp","sal","depth","strat","nox", "DL")

for (i in 1:nzotu){

my.mod<-gbm.step(data=zotus1000,gbm.x = c(7:12), gbm.y = i+12, family = "gaussian", tree.complexity = 10, learning.rate = 0.001, bag.fraction = 0.5)


sp.names[i]<-my.mod[[28]][5][[1]]

preds.mai[,i]<-predict.gbm(my.mod, mai[,c(7,10:12,8:9)], n.trees=my.mod$gbm.call$best.trees, type="response")

preds.phb[,i]<-predict.gbm(my.mod, phb[,c(8,11:12,14,9:10)], n.trees=my.mod$gbm.call$best.trees, type="response")

preds.rot[,i]<-predict.gbm(my.mod, rot[,c(7,10:12,8:9)], n.trees=my.mod$gbm.call$best.trees, type="response")

#preds.goships[,i]<-predict.gbm(my.mod, rot[,c(7,10:12,8:9)], n.trees=my.mod$gbm.call$best.trees, type="response")

  eval.data.out<-unlist(sp.names)
  cv[i]<-my.mod$cv.statistics$correlation.mean
  cvse[i]<-my.mod$cv.statistics$correlation.se
  ntrees[i]<-my.mod$n.trees
  contributions[,i]<-my.mod[[32]][c("temp","sal","depth","strat","nox", "DL"),2]
  write.table(cbind(eval.data.out, cv, cvse, ntrees, t(contributions)), sep='\t', quote=F)
  
save.image(paste(my.mod[[28]][5][[1]], ".gaussian.RData", sep=""))
rm(my.mod)

}


out.data<-cbind(eval.data.out, cv, cvse, ntrees, t(contributions))
colnames(preds.mai)<-eval.data.out
colnames(preds.rot)<-eval.data.out
colnames(preds.phb)<-eval.data.out

write.table(out.data, file=paste(batch.number, "eval.g.data", sep=""), sep=',', quote=F)

write.table(preds.mai, file=paste(batch.number, "mai.g.preds", sep=""), sep=',', quote=F)

write.table(preds.phb, file=paste(batch.number, "phb.g.preds", sep=""), sep=',', quote=F)

write.table(preds.rot, file=paste(batch.number, "rot.g.preds", sep=""), sep=',', quote=F)

save.image(paste(batch.number, "gaussian.RData", sep=""))

