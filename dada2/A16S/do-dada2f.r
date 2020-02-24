#This is a pipeline produced for running an Illumina plate through Dada2
#In this version the plate name is: AHG6G, you can find and replace this code for your purposes.
#**When moving to the HPC replace multithread=16 with multithread=16

#For using this pipeline your Illumina reads must be paired reads, they must be demultiplexed already, they can still be zipped files (fastq.gz) or they can me unzipped files (fastq) either way is fine.

#This pipline is broken down into the following steps, if you wish to optimize them before running, please search for the following titles:
#1) Dada2 unzipfiles step
#2) Dada2 samplename step (double check your files are named the same way as this notebook assumes)
#3) Dada2 qualplots step
#4) Dada2 input primer sequences (change primer sequences according to your samples)
#5) Dada2 filtN step
#6) Dada2 cutadapt step
#7) Dada2 checkprimers step
#8) Dada2 trim step (CHANGE TRIM lengths in this step!) our filtration is truncLen=c(250,228), maxN=0, maxEE=c(2,6)
#9) Dada2 learn err rates step
#10) Dada2 derep step
#11) Dada2 mergers step

#some setup options for outputing markdown files; feel free to ignore these

#knitr::opts_chunk$set(eval = TRUE, 
#                      include = TRUE, 
#                      warning = FALSE, 
#                      message = FALSE,
#                      collapse = TRUE,
#                      dpi = 300,
#                      fig.dim = c(9, 9),
#                      out.width = '98%',
#                      out.height = '98%');
                      
library(R.utils);
library(dada2);
library(ShortRead); packageVersion("ShortRead"); # dada2 depends on this
library(tidyverse); packageVersion("dplyr"); # for manipulating data
library(Biostrings);  # for creating the final graph at the end of the pipeline
library(Hmisc); packageVersion("Hmisc"); # for creating the final graph at the end of the pipeline
#library(plotly); packageVersion("plotly"); # enables creation of interactive graphs, especially helpful for quality plots

#install the following make sure to do it twice. The first time, but after that you don't have to always rerun the libraries.

setwd("/shared/c3/bio_db/BPA/a16s/");
plates<-read_csv('A16SNRSplates.csv');
plates;
args = commandArgs(trailingOnly=TRUE);
p=args[1];
print(p);
plates[p,1];

#library("devtools");

#1) Dada2 unzipfiles step
#library(devtools);
library(dada2);
#unzip all downloaded files using MAC terminal and the following code: gunzip *.fastq.gz OR run the next bit of code to unzip your files from R.
#set home directory and then make a file for all graphs or saved RDS files

home.dir <- ("/shared/c3/bio_db/BPA/a16s/");
dir.create(paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports"));

exports <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports");

data.fp.gz <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/");

#List all files in shared folder to check path
#print(data.fp.gz);
files.fp.gz<-list.files(data.fp.gz,  pattern="fastq.gz");
list.files(files.fp.gz);
print(files.fp.gz);

for (i in seq_along (files.fp.gz)){
              gunzip(filename=paste0(data.fp.gz,files.fp.gz[i]), overwrite=T)};
              


#2) Dada2 samplename step
data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/"); 

#List all files in shared folder to check path
print(data.fp);
list.files(data.fp);

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(data.fp, pattern="_R1.fastq", full.names = TRUE));
fnRs <- sort(list.files(data.fp, pattern="_R2.fastq", full.names = TRUE));
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1);


#3) Dada2 qualplots step
#change the 20 to be the number of samples that you have or that you want to plot

fwd.plot.quals <- plotQualityProfile(fnFs[1:10]);
ggsave(file = paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports/fwd.qualplot.pdf"), fwd.plot.quals, device="pdf");

rev.plot.quals <- plotQualityProfile(fnRs[1:10]);
ggsave(file = paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports/rev.qualplot.pdf"), rev.plot.quals, device="pdf");



#4) Dada2 input primer sequences (change primer sequences according to your samples)
data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/");

FWD <- "TTCCGGTTGATCCYGCCGGA";  ## CHANGE ME to your forward primer sequence Archaea 16S: A2f/Arch21f Brown et al 
REV <- "GWATTACCGCGGCKGCTG"; ## CHANGE ME...Bacteria/Archaea 16S: 519r

allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
        RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD);
REV.orients <- allOrients(REV);


#5) Dada2 filtN step

data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/");
fnFs.filtN <- file.path(data.fp, "filtN", basename(fnFs)); # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(data.fp, "filtN", basename(fnRs));

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread=16, compress=FALSE, matchIDs=T);


#6) Dada2 cutadapt step (change cutaapt pathway to find it on your computer)

cutadapt <- "/shared/c3/apps/anaconda3/bin/cutadapt"; # CHANGE ME to the cutadapt path on your machine
#if you have issues finding the file path to cutadapt go into terminal and just type which cutadapt and paste in the whole file pathway above, don't use a relative pathway.
data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/filtN");
#system2(cutadapt.file, args = "--version"); # Run shell commands from R

path.cut <- file.path(data.fp, "cutadapt");
if(!dir.exists(path.cut)) dir.create(path.cut);
fnFs.cut <- file.path(path.cut, basename(fnFs));
fnRs.cut <- file.path(path.cut, basename(fnRs));

FWD.RC <- dada2:::rc(FWD);
REV.RC <- dada2:::rc(REV);
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC);
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC);
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}


#7) Dada2 checkprimers step

data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/filtN/cutadapt");

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0));
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
    FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
    REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
    REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]));

#there might be a few primers remaining here, it is ok if this happens, we will just check everything after merging and if we find any primer sequences we will remove those ASVs at that point.

#8) Dada2 trim step (CHANGE TRIM lengths in this step!)
data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/filtN/cutadapt/");
# Place filtered files in filtered/ subdirectory
trimFs <- file.path(data.fp, "trimmed", paste0(sample.names, "_R1_trim.fastq"));
trimRs <- file.path(data.fp, "trimmed", paste0(sample.names, "_R2_trim.fastq"));
names(trimFs) <- sample.names;
names(trimRs) <- sample.names;

out <- filterAndTrim(fnFs.cut, trimFs, fnRs.cut, trimRs, truncLen=c(250,228),
              maxN=0, maxEE=c(2,6), truncQ=2, rm.phix=TRUE,
              compress=TRUE, multithread=16, minLen = 50, matchIDs=TRUE) # On Windows set multithread=FALSE
head(out);

#now filter and trim your sequences, these parameters can be adjusted to fit your needs, truncLen=c(250,228) means that you are trimming the forward read or R1 to 250 bp (see your qual plot and notice that this is where the sequence quality begins to crash) AND you are trimming the reverse read or R2 to 228 bp.
#you are also allowing a maxEE of 6 for the R1 and 10 for the R2, this means that more errors are allowed for the reverse read, the standard error is (2,2). And you are allowing zero N's to make it through filtering which is why you don't need to prefilter N's
#match ID's = TRUE this is because if you lose one read R1 or R2 the other is automatically removed because you are dealing with paired end reads.
#minLen = 50 means that any reads with fewer than 50bp are automatically removed in this step, we will look at this again at the end and likely remove more sequences if any are shorter than around 350bp.

#9) Dada2 learn err rates step


data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/filtN/cutadapt/trimmed");

errF <- learnErrors(trimFs, verbose=TRUE, multithread=16);
errR <- learnErrors(trimRs, verbose=TRUE, multithread=16);

fwd.plot.errors <- plotErrors(errF, nominalQ=TRUE);
ggsave(file = paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports/fwd.errors.pdf"),fwd.plot.errors, width = 10, height = 10, device="pdf");

rev.plot.errors <- plotErrors(errR, nominalQ=TRUE);
ggsave(file = paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports/rev.errors.pdf"), rev.plot.errors, width = 10, height = 10, device="pdf");

#10) Dada2 derep step

data.fp <- paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], "/filtN/cutadapt/trimmed");
#run the following code to insure that no samples are missing from the trimming section

derepFs<-derepFastq(trimFs);
derepRs<-derepFastq(trimRs);

dadaFs <- dada(derepFs, err=errF, multithread=16);
dadaRs <- dada(derepRs, err=errR, multithread=16);

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE);
# Inspect the merger data.frame from the first sample
head(mergers[[1]]);

seqtab <- makeSequenceTable(mergers);
dim(seqtab);


#11) Dada2 mergers step
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)));
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=16, verbose=TRUE);
dim(seqtab.nochim);
sum(seqtab.nochim)/sum(seqtab);
getN <- function(x) sum(getUniques(x));
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim));
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim");
rownames(track) <- sample.names;
head(track);

plotLengthDistro <- function(st) {
  require(ggplot2)
  tot.svs <- table(nchar(colnames(st)))
  tot.reads <- tapply(colSums(st), nchar(colnames(st)), sum)
  df <- data.frame(Length=as.integer(c(names(tot.svs), names(tot.reads))),
                   Count=c(tot.svs, tot.reads),
                   Type=rep(c("SVs", "Reads"), times=c(length(tot.svs), length(tot.reads))))
  pp <- ggplot(data=df, aes(x=Length, y=Count, color=Type)) + geom_point() + facet_wrap(~Type, scales="free_y") + theme_bw() + xlab("Amplicon Length")
  pp
  }

plotLengthDistro(seqtab.nochim);
plotLengthDist.log10 <-plotLengthDistro(seqtab.nochim) + scale_y_log10();
ggsave(file = paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports/plotLengthDist.log10.pdf"), plotLengthDist.log10, width = 10, height = 5, device="pdf");

sample <- rownames(seqtab.nochim);
sequence <- colnames(seqtab.nochim);

#check the col names and check how many ASVs that you are losing in the pipeline
colnames(seqtab.nochim);
#what % had chimera's vs non-chimeras?
sum(seqtab.nochim)/sum(seqtab);
#this is the %
sum(rev(sort(colSums(seqtab.nochim)))[1:1000])/sum(colSums(seqtab.nochim));

# Flip table
seqtab.t <- as.data.frame(t(seqtab.nochim));
write.csv(seqtab.t, file = paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports/ASV.",plates[p,1], ".table.csv"), col.names=NA);

# tracking reads by percentage

track_pct <- track %>% 
  data.frame() %>%
  mutate(Sample = rownames(.),
         filtered_pct = ifelse(filtered == 0, 0, 100 * (filtered/input)),
         denoisedF_pct = ifelse(denoisedF == 0, 0, 100 * (denoisedF/filtered)),
         denoisedR_pct = ifelse(denoisedR == 0, 0, 100 * (denoisedR/filtered)),
         merged_pct = ifelse(merged == 0, 0, 100 * merged/((denoisedF + denoisedR)/2)),
         nonchim_pct = ifelse(nonchim == 0, 0, 100 * (nonchim/merged)),
         total_pct = ifelse(nonchim == 0, 0, 100 * nonchim/input)) %>%
  select(Sample, ends_with("_pct"));


track_pct_avg <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                            list(avg = mean));
head(track_pct_avg);

track_pct_med <- track_pct %>% summarize_at(vars(ends_with("_pct")), 
                                            list(avg = stats::median));
head(track_pct_avg);

head(track_pct_med);

track_plot <- track %>% 
  data.frame() %>%
  mutate(Sample = rownames(.)) %>%
  gather(key = "Step", value = "Reads", -Sample) %>%
  mutate(Step = factor(Step, 
                       levels = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"))) %>%
  ggplot(aes(x = Step, y = Reads)) +
  geom_line(aes(group = Sample), alpha = 0.2) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0)) + 
  stat_summary(fun.y = median, geom = "line", group = 1, color = "steelblue", size = 1, alpha = 0.5) +
  stat_summary(fun.y = median, geom = "point", group = 1, color = "steelblue", size = 2, alpha = 0.5) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5), 
               geom = "ribbon", group = 1, fill = "steelblue", alpha = 0.2) +
  geom_label(data = t(track_pct_avg[1:5]) %>% data.frame() %>% 
               rename(Percent = 1) %>%
               mutate(Step = c("filtered", "denoisedF", "denoisedR", "merged", "nonchim"),
                      Percent = paste(round(Percent, 2), "%")),
             aes(label = Percent), y = 1.1 * max(track[,2])) +
  geom_label(data = track_pct_avg[6] %>% data.frame() %>%
               rename(total = 1),
             aes(label = paste("Total\nRemaining:\n", round(track_pct_avg[1,6], 2), "%")), 
             y = mean(track[,6]), x = 6.5) +
  expand_limits(y = 1.1 * max(track[,2]), x = 7) +
  theme_classic();

ggsave(file = paste0("/shared/c3/bio_db/BPA/a16s/", plates[p,1], ".exports/track_plot.pdf"), track_plot, width = 6, height = 6, device="pdf");


