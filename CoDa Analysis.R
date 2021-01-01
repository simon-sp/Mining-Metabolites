### Set up environment
### Make sure to load and run the modified pairwise DA wrapper
options(stringsAsFactors = F)
#source("http://bioconductor.org/biocLite.R")
pacman::p_load(cowplot, vegan, ggplot2, dplyr, grid, ggrepel, reshape2, fossil, iNEXT, metagenomeSeq, coda.base, stringr, zCompositions, ALDEx2, RColorBrewer, metafolio, Tjazi, omixerRpm)
setwd('/myPath/')
counts   <- as.data.frame(read.csv("/myDir/genus_table.csv", header = T, row.names = 1))
metadata <- read.table("myDir/metadata.csv", header = T, sep = ',')

### Make sure metadata and counts table have samples in the correct order
metadata$Run == colnames(counts)
metadata = metadata[order(match(metadata$Run, colnames(counts))),]
metadata$Run == colnames(counts)

species   <- apply(counts, c(1,2),function(x) as.numeric(as.character(x)))
species   <- species[apply(species == 0, 1, sum) <= (ncol(counts) -2), ]    #remove rows with 2 or fewer hits
conds       <- c(rep("A", ncol(species)-10 ), rep("B", 10)) #If you have less than 12 animals, adjust!
species.clr <- aldex.clr(species, conds, mc.samples = 1000, denom="all", verbose=TRUE, useMC=FALSE) 
species.eff <- aldex.effect(species.clr, verbose = TRUE, include.sample.summary = TRUE)
colnames(species.eff) <- gsub(pattern = "rab.sample.", replacement = "", x = colnames(species.eff))
species.exp <- (species.eff[,c(4:(ncol(species.eff)-4))]) #remove the useless t-test-like results
res_posthoc = pairwise_DA_wrapper(species, 
                                   comparisons =combn(unique(metadata$Group),1),
                                   groups = metadata$Group,
                                   ignore.posthoc = F)
write.table(res_posthoc, quote = F, 'diff_microbes.csv', sep = ',')

### Generate GBM Counts
seqtab <- readRDS("myDir/seqtab_final.rds")
tax   <- readRDS("myDir/tax_final.rds")
samples.out <- rownames(seqtab)
samdf <- data.frame(ID = samples.out)
row.names(samdf) = samples.out
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))
pip = t(ps@otu_table@.Data)
pip_otu_table = pip
row.names(pip_otu_table) = c(paste("otu", 1:nrow(pip_otu_table), sep = ""))
write.table(data.frame("OTU"=rownames(pip_otu_table),pip_otu_table), 
            file = "Piphillin_otu_table.csv", row.names=FALSE, quote = F, sep = ",", )
repseqs <- c(paste(">otu", 1:nrow(pip_otu_table), "\n", row.names(pip), "\n", sep = ""))
write.table(unname(c(paste(">otu", 
                           1:nrow(pip_otu_table), 
                           "\n", row.names(pip), 
                           sep = ""))), 
            file = "Piphillin_representative.csv", quote = F, row.names = F, col.names = F)
### Submit these files into the Piphillin server

### Once piphilin is downloaded continue on
first <- read.table("myDir/ko_abund_table_unnorm.txt", sep = ,'\t', header=T)
colnames(first)[1] = "Pathway_number"
### Download all KEGG ortholog pathways and numbers (see KEGG.csv on GitHub)

second<- read.csv("/myPath/KEGG.csv",header=T)
third<-left_join(first,second,by="Pathway_number") #make sure the first column names are the same in two tables
write.csv(third,file="myPath/KEGG annotation results.csv")
kegg = read.csv('myPath/KEGG annotation results.csv', header = 1)
kegg = kegg[, colnames(kegg) != 'X' & colnames(kegg) != 'Pathway_name']
head(kegg)
library(devtools)
install_github('https://github.com/omixer/omixer-rpmR/releases/download/0.3.1/omixerRpm_0.3.1.tar.gz')
pacman::p_load(omixerRpm)
db <- loadDB("GBMs.v1.0")
mods <- rpm(kegg, minimum.coverage=0.3, annotation = 1, module.db = db)
mods@coverage[,1]
getNames(db, mods@annotation[1,])
mods
asDataFrame(mods, "coverage")
gbm = t(mods@abundance)
gbm 
colnames(gbm) = mods@annotation[1:nrow(mods@annotation),]
gbm   <- apply(gbm,c(1,2),function(x) as.numeric(as.character(x)))
gbm   <- gbm[apply(gbm == 0, 1, sum) <= (ncol(gbm) -2), ]    #remove rows with 2 or fewer hits
gbm = floor(gbm)
rownames(gbm) == metadata$Run
metadata = metadata[order(match(metadata$Run, rownames(gbm))),]
rownames(gbm) == metadata$Run

res_posthoc = pairwise_DA_wrapper(t(gbm), 
                                   comparisons = comparisons,
                                   groups = metadata$Group,
                                   ignore.posthoc = F)

