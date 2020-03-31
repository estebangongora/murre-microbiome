#Loading required libraries
library(dada2)

#Setting paths
path <- "<path_to_output_folder>"
data_path <- "<path_to_fastq_files_folder"

#Getting fastq files and extracting sample names
fnFs <- sort(list.files(data_path, pattern = "_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(data_path, pattern = "_R2.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Checking data quality, filtering, and trimming
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "CCGYCAATTYMTTTRAGTTT"
FWD_LENGTH <- nchar(FWD)
REV_LENGTH <- nchar(REV)
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(245,240), trimLeft = c(FWD_LENGTH,REV_LENGTH), maxN = 0, maxEE = c(2,4),  multithread = T)
head(out)

#Checking if primers were removed correclty by filterAndTrim
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients
library("ShortRead")
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = filtRs[[1]]), #or fn = fnRs.filtN
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = filtFs[[1]]), #or fn = fnFs.filtN
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = filtRs[[1]])) #or fn = fnRs.filtN

#Cutadapt is only needed when primers are not removed using filterAndTrim (e.g. variable length amplicons/primers or overlapping amplicons)
cutadapt <- "<path_to_cutadapt.exe_file>"
system2(cutadapt, args = "--version")
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs.filtN[i], fnRs.filtN[i])) 
}
#Cutadapt over

#Leraning error rates
errF <- learnErrors(filtFs, nbases=1e+10, randomize=T, multithread=T)
errR <- learnErrors(filtRs, nbases=1e+10, randomize=T, multithread=T)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Dereplicating samples
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Run dada2 core algorithm
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[3]]
head(getUniques(dadaFs[[3]]))
head(getSequences(dadaFs[[3]]))
dadaFs_pool <- dada(derepFs, err=errF, pool = T, multithread=TRUE)
dadaRs_pool <- dada(derepRs, err=errR, pool = T, multithread=TRUE)
dadaFs_pool[[3]]
head(getUniques(dadaFs_pool[[3]]))
head(getSequences(dadaFs_pool[[3]]))

#Merging peired-end reads
mergers <- mergePairs(dadaFs_pool, derepFs, dadaRs_pool, derepRs, verbose=TRUE)
head(mergers[[3]])

#Creating ASV Table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(311,416)]
dim(seqtab2)
table(nchar(getSequences(seqtab2)))

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="pooled", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Tracking removed reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs_pool, getN), sapply(dadaRs_pool, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF_pooled", "denoisedR_pooled", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Assigning taxonomy
taxa <- assignTaxonomy(seqtab.nochim, "<path_to_silva_nr_v132_train_set.fa.gz_training_set_file>", multithread=TRUE)
taxa <- addSpecies(taxa, "<path_to_silva_species_assignment_v132.fa.gz_taxonomy_file>")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#Converting dada2 output for use in Qiime2
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#Creating fasta file of ASV sequences
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#Creating a count table
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.txt", sep="\t", quote=F)

#Creating a taxonomy table
asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "ASVs_taxonomy.txt", sep="\t", quote=F)

#Make sure you open the generated .txt files and add the ASV_ID header in your ASV counts and taxa table