### Packages and Setting Directory
pacman::p_load(dada2, phyloseq, microbiome, tidyr)
setwd("/myDir/myFile/")

# Filename parsing
# If your reads aren't paired, you don't need to run any of the lines
# involving pathR, filtpathR, fastqRs etc.
pathF <- "F/" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
pathR <- "R/" # CHANGE ME to the directory containing your demultiplexed forward-read fastqs
filtpathF <- file.path(pathF, "Filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "Filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="fq.gz")) # Change to fastq.gz if necessary
fastqRs <- sort(list.files(pathR, pattern="fq.gz")) # Change to fastq.gz if necessary
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS

### Remove from the left and cut the rest until you reach desired length
filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              trimLeft = 37, truncLen=237, maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=FALSE)
# If your reads are single-end, delete the reverse file paths in the above function

filtFs <- list.files(filtpathF, pattern="fq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fq.gz", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
names(filtFs)
filtRs
# Check quality to see if it worked: 
# Use fastqc in R, check quality '/myPath/F/*.fastq.gz -o ...
# Check that quality scores for your bases are >28

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=FALSE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=FALSE)
# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Again when working with single-end reads, adjust the code above to get rid of derepR and ddR

# Construct sequence table and remove chimeras 
seqtab <- makeSequenceTable(mergers)
###Change the location
saveRDS(seqtab, "Processed/seqtab.rds") # CHANGE ME to where you want sequence table saved
# Remove chimeras
st1 <- readRDS("Processed/seqtab.rds")
seqtab <- removeBimeraDenovo(st1, method="consensus", multithread=TRUE)
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "refDataBasePath/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# Write to disk
saveRDS(seqtab, "Processed/seqtab_final.rds") # CHANGE ME to where you want sequence table saved
saveRDS(tax, "Processed/tax_final.rds") # CHANGE ME ...

seqtab <- readRDS("Procesed/seqtab_final.rds")
tax   <- readRDS("Processed/tax_final.rds")
samples.out <- rownames(seqtab)
samdf <- data.frame(ID = samples.out)
row.names(samdf) = samples.out

ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(tax))

count_table_tax  = t(rbind(ps@otu_table, t(ps@tax_table)))

count_table_one_word = data.frame(count_table_tax[,1:(nrow(ps@otu_table))])
View(head(count_table_one_word))
count_table_one_word$taxonomy = paste0(count_table_tax[,"Kingdom"], ";" , 
                                       count_table_tax[,"Phylum"] , ";" , 
                                       count_table_tax[,"Class"]  , ";" , 
                                       count_table_tax[,"Order"]  , ";" , 
                                       count_table_tax[,"Family"] , ";" , 
                                       count_table_tax[,"Genus"]  , ";")


write.csv(x = count_table_one_word, quote = F, file = "Processed/count_table.csv")


#Genus_count_table

agg <- aggregate_taxa(ps, level = "Genus", rm.na = F, verbose = TRUE)

colap <- unite(as.data.frame(agg@tax_table@.Data), newCol, -unique) 

genus_table <- agg@otu_table@.Data
View(head(genus_table))
row.names(genus_table) <- colap$newCol

write.csv(x = genus_table, quote = F, file = "Processed/genus_table.csv")
