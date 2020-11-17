#### IMPORT AND PROCESS HEALTH METADATA ####
cy.metadata <- read.csv("coyote_metadata.csv")

# Convert data to factors, where necessary.
for(i in c("PA.Empty", "PA.Prey", "PA.Veg", "PA.Anthro", "Analysis.Season")){
  cy.metadata[,i] <- as.factor(as.character(cy.metadata[,i]))
  }

# Identify and impute any missing data points using linear regressions on the remaining data.
ggplot(cy.metadata, aes(Mass)) + geom_histogram()
impute.model <- lm(Mass ~ Girth + Cementum.Age + KFI, data=cy.metadata) # Model for missing mass.
impute.data <- subset(cy.metadata, is.na(Mass)=="TRUE") 
for(i in 1:nrow(cy.metadata)){
  if(is.na(cy.metadata[i,"Mass"])=="TRUE"){
    cy.metadata[i,"Mass"] <- predict(impute.model, newdata=cy.metadata[i,]) } }

for(i in 1:nrow(cy.metadata)){ # Replace missing normalized spleen mass for coyote 73.
  if(is.na(cy.metadata[i,"Spleen.to.Body.Ratio"])=="TRUE"){
    cy.metadata[i,"Spleen.to.Body.Ratio"] <- cy.metadata[i,"Spleen.Mass"]/cy.metadata[i,"Mass"] } }

ggplot(cy.metadata, aes(Snout.to.Base)) + geom_histogram()
impute.model <- lm(Snout.to.Base ~ Girth + KFI + Mass, data=cy.metadata) # Missing body length
impute.data <- subset(cy.metadata, is.na(Snout.to.Base)=="TRUE") 
for(i in 1:nrow(cy.metadata)){
  if(is.na(cy.metadata[i,"Snout.to.Base"])=="TRUE"){
    cy.metadata[i,"Snout.to.Base"] <- predict(impute.model, newdata=cy.metadata[i,]) } }

ggplot(cy.metadata, aes(KFI)) + geom_histogram()
impute.model <- lm(KFI ~ Mass + Girth + Snout.to.Base, data=cy.metadata) # Missing KFI
impute.data <- subset(cy.metadata, is.na(KFI)=="TRUE") 
for(i in 1:nrow(cy.metadata)){
  if(is.na(cy.metadata[i,"KFI"])=="TRUE"){
    cy.metadata[i,"KFI"] <- predict(impute.model, newdata=cy.metadata[i,]) } }

ggplot(cy.metadata, aes(Cementum.Age)) + geom_histogram()
impute.model <- glm.nb(Cementum.Age ~ Mass + Girth + Snout.to.Base, data=cy.metadata)
impute.data <- subset(cy.metadata, is.na(Cementum.Age)=="TRUE") # Missing age.
for(i in 1:nrow(cy.metadata)){
  if(is.na(cy.metadata[i,"Cementum.Age"])=="TRUE"){
    cy.metadata[i,"Cementum.Age"] <- predict(impute.model, newdata=cy.metadata[i,]) } }

rm(impute.model, impute.data)

#### DADA2 - RAW SEQUENCE PROCESSING ####
# Initial sequence processing (based on the DADA2 tutorial: https://benjjneb.github.io/dada2/tutorial.html)
# Set a seed so that analyses are reproducible.
set.seed(100)

# Name the folder where the sequence files are stored.
input_path <- "~/coyote/sequences/"

# Identify forward (Fs) and reverse (Rs) reads.
fnFs <- sort(list.files(input_path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(input_path, pattern="_R2_001.fastq"))

# Identify sample names (should be the beginning of the filename).
sampleNames <- sapply(strsplit(fnFs, "_"), '[', 1)
fnFs <- file.path(input_path, fnFs)
fnRs <- file.path(input_path, fnRs)

# Create a new directory where filtered reads will be stored.
filt_path <- file.path(input_path, "filtered")
if(!file_test("-d", filt_path)) dir.create(filt_path)

# Name filtered sequences with the same names as input sequences with _filt added.
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sampleNames, "_R_filt.fastq.gz"))

# This filtering step allows a maximum of two expected errors per read, and zero ambiguous (N) bases.
# It also trims the forward reads to 240 bp and the reverse reads to 160 bp.
# The "out" table gives you a table with the number of input and output reads from the QC step.
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(240,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE)

# Dereplicate sequences.
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Assign sample names.
names(derepFs) <- sampleNames
names(derepRs) <- sampleNames

# Calculate error probabilities.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Use error probabilities to determine exact sequence variants.
dadaFs <- dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)

# Merged paired-end reads and construct an OTU table.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
seqtabAll <- makeSequenceTable(mergers)

# Examine the distribution of sequence lengths. Most sequences should be ~253bp.
table(nchar(getSequences(seqtabAll)))

# OPTIONAL: Trim sequence tables to remove amplicons that are not within 250-256bp.
seqtabAll <- seqtabAll[,nchar(colnames(seqtabAll)) %in% seq(250,256)]

# Remove chimeras.
seqtabNoC <- removeBimeraDenovo(seqtabAll, multithread=TRUE)

# This section produces a table that shows how many reads were lost during the processing steps.
getN <- function(x) sum(getUniques(x))
read_counts <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN),
               rowSums(seqtabAll), rowSums(seqtabNoC))
colnames(read_counts) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab", "nochimera")
rownames(read_counts) <- sampleNames

# Identify the reference database you will use for assigning taxonomy. This was downloaded from https://benjjneb.github.io/dada2/training.html.
fastaRef <- "~/coyote/dada2_rdp_train_set_16.fa.gz"

# Assign taxonomy to your sequence table.
taxtabNoC <- assignTaxonomy(seqtabNoC, refFasta=fastaRef, multithread=TRUE)

# Check to see that it worked.
unname(head(taxtabNoC))

# Import initial sample data
sample_data <- read.csv("~/coyote/initial_data.csv", row.names=SampleID)

# Import into phyloseq object and clean workspace
pseq.raw <- phyloseq(otu_table(seqtabNoC, taxa_are_rows=FALSE),
                     sample_data(sample_data),
                     tax_table(taxtabNoC))

rm(seqtabNoC, sample_data, taxtabNoC, fastaRef, seqtabAll, mergers, out, derepRs, derepFs, errF, errR,
   dadaFs, dadaRs, filt_path, input_path, fnFs, fnRs, sampleNames)

# Remove chloroplasts
temp <- as.data.frame(tax_table(pseq.raw))
temp$chloroplast <- temp$Class == "Chloroplast"
temp <- subset(temp, chloroplast=="TRUE")
badTaxa <- rownames(temp)
goodTaxa <- setdiff(taxa_names(complete.phyloseq.trimmed), badTaxa)
pseq.raw <- prune_taxa(goodTaxa, pseq.raw)
rm(temp, goodTaxa, badTaxa)

#### PROCESS MICROBIOME DATA ####
# Store taxa names as a 'refseq' object.
dna <- Biostrings::DNAStringSet(taxa_names(pseq.raw))
names(dna) <- taxa_names(pseq.raw)
pseq.raw <- merge_phyloseq(pseq.raw, dna)
rm(dna)
taxa_names(pseq.raw) <- paste0("ASV", seq(ntaxa(pseq.raw)))

# Determine which taxa are not assigned to the phylum level for use in a BLAST search.
temp.unclass <- subset_taxa(pseq.raw, Phylum !="NA")
temp.unclass <- as.data.frame(cbind(tax_table(temp.unclass)))
class.string <- rownames(temp.unclass)
temp.unclass <- subset_taxa(pseq.raw, !(taxa_names(pseq.raw) %in% class.string))

# Export the ASV sequences for these taxa. BLAST search to identify closest sequence identity.
seqtab.nochim <- as.matrix(otu_table(temp.unclass))
seqs <- refseq(temp.unclass)
ids <- taxa_names(temp.unclass)
db_out <- data.frame(ids=ids, seqs=seqs, count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
writeFasta(fasta, file = "unclassified.fna")
write.csv(db_out, "unclassified.csv")
rm(fasta, db_out, ids, seqs, temp.unclass, seqtab.nochim, class.string, i)

# Remove taxa that identify as archaea
pseq.filter <- subset_taxa(pseq.raw, !(Kingdom %in% c("Archaea")))

# Remove taxa that identify as mitochondria, chloroplasts, or other protists in a BLAST search of the file created above
pseq.filter <- subset_taxa(pseq.filter, !(taxa_names(pseq.filter) %in% 
                                                c("ASV44", "ASV67", "ASV227", "ASV310",
                                                  "ASV518", "ASV675", "ASV678", "ASV1053",
                                                  "ASV1170", "ASV1297", "ASV1908")))

# Identify contaminants using decontam.
sample_data(pseq.filter)$is.control <- sample_data(pseq.filter)$Study == "control" # Create logical identifying controls.
contamdf.prev <- isContaminant(pseq.filter, method="prevalence", neg="is.control", threshold=0.25)
table(contamdf.prev$contaminant)

# Plot presence/absence in negative controls vs. true samples
contam.data <- list()

ps.pa <- transform_sample_counts(pseq.filter, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Study == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Study != "control", ps.pa)
contam.data$prevalence <- data.frame(pa.pos=taxa_sums(ps.pa.pos),
                              pa.neg=taxa_sums(ps.pa.neg),
                              contaminant=contamdf.prev$contaminant)

ggplot(data=contam.data$prevalence, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

rm(ps.pa, ps.pa.neg, ps.pa.pos)

# Check relative abundances of contaminants
contamdf.prev <- subset(contamdf.prev, contaminant==TRUE) # Subset to true contaminants
contam.data$abundance <- transform_sample_counts(pseq.filter, function(x) 100*x/sum(x))
contam.data$abundance <- subset_taxa(contam.data$abundance, 
                                     taxa_names(contam.data$abundance) %in% as.character(rownames(contamdf.prev))) # Select only contaminants
contam.data$abundance <- cbind(as.data.frame(cbind(tax_table(contam.data$abundance))),
                     as.data.frame(t(otu_table(contam.data$abundance))))
contam.data[["abundance"]]$means <- rowMeans(contam.data$abundance[,c(7:99)])
contam.data[["abundance"]]$ASV <- rownames(contam.data[["abundance"]]) 

ggplot(contam.data$abundance, aes(x=ASV, y=means)) + geom_bar(stat="identity", position="dodge") +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

# Remove contaminants from samples
pseq.filter <- subset_taxa(pseq.filter, !(taxa_names(pseq.filter) %in% as.character(rownames(contamdf.prev)))) # Remove 22 contaminants.
sample_data(pseq.filter)$is.neg <- NULL
rm(contamdf.prev)

# Remove negative controls from phyloseq object.
pseq.filter <- subset_samples(pseq.filter, Study !="control") # Remove control samples from further analysis.
pseq.filter <- prune_taxa(taxa_sums(pseq.filter) > 0, pseq.filter) # 

# Examine read count distribution to determine appropriate cutoff threshold.
temp <- as.data.frame(sample_sums(pseq.filter))
ggplot(temp, aes(x=temp[,1])) + geom_histogram() # Looks like a clear divide, with five samples far below 4,500 reads
rm(temp)

# Evaluate sample completeness to confirm appropriate cutoff threshold.
temp <- as.data.frame(t(otu_table(pseq.filter)))
iNext.result <- iNEXT(temp, q=0, datatype="abundance") # All samples, except the five low-read samples, have completion >98%
rm(temp)

# Remove singletons and samples with fewer than 4,500 reads prior to rarefaction.
pseq.filter <- subset_samples(pseq.filter, sample_sums(pseq.filter) > 4500) # Remove low-read samples
pseq.filter <- prune_taxa(taxa_sums(pseq.filter) > 1, pseq.filter) # Remove singleton taxa

# Rarefy samples to the minimum remaining library size (approx. 4500 reads). Average across 1,000 rarefactions.
rarefaction.average <- list()
for(i in 1:1000){
  temp.rarefy <- rarefy_even_depth(pseq.rarefied, 
                                   sample.size = min(sample_sums(pseq.rarefied)),
                                   replace = FALSE,
                                   trimOTUs = FALSE,
                                   rngseed=i)
  rarefaction.average[[i]] <- as.data.frame(otu_table(temp.rarefy))
}

dfAvg <- Reduce("+", rarefaction.average)/length(rarefaction.average)
dfAvg <- round(dfAvg, 0)
dfAvg <- otu_table(dfAvg, taxa_are_rows=FALSE)

# Replace feature table in phyloseq object with new, rarefied feature table.
otu_table(pseq.rarefied) <- dfAvg
pseq.rarefied <- prune_taxa(taxa_sums(pseq.rarefied) > 0, pseq.rarefied)
rm(temp.rarefy, rarefaction.average, dfAvg)

# Construct a phylogenetic tree for both filtered and rarefied data, for use in a phyloseq object.
pseq.list <- list(pseq.filter, pseq.rarefied)

for(i in c(1:2)){
  vector <- taxa_names(pseq.list[[i]])
  taxa_names(pseq.list[[i]]) <- refseq(pseq.list[[i]])
  seqtabNoC <- as.matrix(as.data.frame(otu_table(pseq.list[[i]])))
  seqs <- getSequences(seqtabNoC)

  names(seqs) <- seqs
  alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA, verbose=FALSE)
  phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
  dm <- dist.ml(phangAlign)
  treeNJ <- NJ(dm)
  fit = pml(treeNJ, data=phangAlign)
  fitGTR <- update(fit, k=4, inv=0.2)
  fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE, rearrangement="stochastic", control=pml.control(trace=0))

  phy_tree(pseq.list[[i]]) <- fitGTR$tree
  taxa_names(pseq.list[[i]]) <- vector
  
  rm(fitGTR, fit, treeNJ, dm, phangAlign, alignment, seqs, seqtabNoC, vector)
}

pseq.filter <- pseq.list[[1]]
pseq.rarefied <- pseq.list[[2]]
rm(pseq.list)

# Specify one phyloseq object for use in CLR-based analyses #
pseq.clr <- pseq.filter

# Add sample metadata to each phyloseq object.
rownames(cy.metadata) <- cy.metadata$FecesID
sample_data(pseq.rarefied) <- subset(cy.metadata, FecesID %in% sample_names(pseq.rarefied))
sample_data(pseq.clr) <- subset(cy.metadata, FecesID %in% sample_names(pseq.clr))

# Replace 'NA' taxa with 'Uncl' to prevent errors during agglomeration.
pseq.rarefied.sub <- pseq.rarefied
for(i in 1:nrow(tax_table(pseq.rarefied.sub))){
  for(j in 2:6){
    if(is.na(tax_table(pseq.rarefied.sub)[i,j])==TRUE){
      if(substr(tax_table(pseq.rarefied.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(pseq.rarefied.sub)[i,j] <- tax_table(pseq.rarefied.sub)[i,j-1]}
      else {
        tax_table(pseq.rarefied.sub)[i,j] <- paste0("Uncl_", tax_table(pseq.rarefied.sub)[i,j-1])}}
  }}

pseq.clr.sub <- pseq.clr
for(i in 1:nrow(tax_table(pseq.clr.sub))){
  for(j in 2:6){
    if(is.na(tax_table(pseq.clr.sub)[i,j])==TRUE){
      if(substr(tax_table(pseq.clr.sub)[i,j-1], 1, 4)=="Uncl"){
        tax_table(pseq.clr.sub)[i,j] <- tax_table(pseq.clr.sub)[i,j-1]}
      else {
        tax_table(pseq.clr.sub)[i,j] <- paste0("Uncl_", tax_table(pseq.clr.sub)[i,j-1])}}
  }}

# Calculate alpha diversity metrics from rarefied data.
alpha.data <- estimate_richness(pseq.rarefied) # Alpha diversity.
alpha.data$FecesID <- sample_names(pseq.rarefied)
cy.metadata <- merge(cy.metadata, alpha.data, by="FecesID", all=TRUE)

# Calculate PD and NTI from rarefied data using picante.
temp.phylo.class = phy_tree(pseq.rarefied)
temp.vegan.otu.table = as.data.frame(otu_table(pseq.rarefied))
picante.pd.result <- pd(temp.vegan.otu.table, temp.phylo.class, include.root=FALSE)

picante.phydist <- cophenetic(temp.phylo.class)
picante.ses.mntd <- ses.mntd(temp.vegan.otu.table, picante.phydist, null.model="taxa.labels")
NTI <- picante.ses.mntd$mntd.obs.z*(-1)

picante.results <- data.frame(FecesID=rownames(picante.pd.result),
                              PD = picante.pd.result$PD,
                              NTI = NTI)

cy.metadata <- merge(cy.metadata, picante.results, by="FecesID", all=TRUE)

rm(temp.phylo.class, temp.vegan.otu.table, picante.pd.result, picante.phydist,
   picante.ses.mntd, NTI, picante.results, alpha.data)

# Calculate richness and diversity using rarefaction and extrapolation from unrarefied data (iNext)
temp <- as.data.frame(t(otu_table(pseq.filter)))
iNext.result <- iNEXT(temp, q=0, datatype="abundance") # All samples, except the five low-read samples, have completion >98%
rm(temp)

temp <- iNext.result$AsyEst
temp.richness <- subset(temp, Diversity=="Species richness")
temp.richness$FecesID <- temp.richness$Site
temp.richness$Observed_Extrap <- temp.richness$Estimator
temp.richness$Observed_Extrap <- round(temp.richness$Observed_Extrap, digits = 0)

temp.diversity <- subset(temp, Diversity=="Shannon diversity")
temp.diversity$FecesID <- temp.diversity$Site
temp.diversity$Shannon_Extrap <- log(temp.diversity$Estimator)

cy.metadata <- merge(cy.metadata, temp.richness[,c("FecesID", "Observed_Extrap")], by="FecesID", all=TRUE)
cy.metadata <- merge(cy.metadata, temp.diversity[,c("FecesID", "Shannon_Extrap")], by="FecesID", all=TRUE)

rm(temp.diversity, temp.richness, temp)

# Calculate PD using rarefaction and extrapolation from unrarefied data
library(phylobase)
temp <- as.data.frame(t(otu_table(pseq.filter)))
labels <- rownames(temp)
tree <- as(phy_tree(pseq.filter), "phylo4")
tree <- as(tree, "phylog")

iNextPD.result <- iNextPD::iNextPD(temp, labels, q=0, tree, datatype="abundance")

temp <- as.data.frame(iNextPD.result$AsyPDEst)
temp <- subset(temp, Var3=="Estimator" & Var2=="q = 0")
temp$FecesID <- temp$Var1
temp$PD_Extrap <- temp$Freq

cy.metadata <- merge(cy.metadata, temp[,c("FecesID","PD_Extrap")], by="FecesID", all=TRUE)
rm(temp, tree, labels)

# Calculate NTI from unrarefied data
temp.phylo.class = phy_tree(pseq.clr)
temp.vegan.otu.table = as.data.frame(otu_table(pseq.clr))

picante.phydist <- cophenetic(temp.phylo.class)
picante.ses.mntd <- ses.mntd(temp.vegan.otu.table, picante.phydist, null.model="taxa.labels")
NTI <- picante.ses.mntd$mntd.obs.z*(-1)

picante.results <- data.frame(FecesID=sample_names(pseq.clr),
                              NTI_Raw = NTI)

cy.metadata <- merge(cy.metadata, picante.results[,c("FecesID","NTI_Raw")], by="FecesID", all=TRUE)

rm(temp.phylo.class, temp.vegan.otu.table, picante.phydist, picante.ses.mntd, NTI, picante.results)

# Calculate distance matrices (BC, Jac, Aitchison, wUF, and uwUF).
dist.cy.list <- list()
dist.cy.list[["Bray-Curtis"]] <- vegdist(as.data.frame(otu_table(pseq.rarefied)), method="bray")
dist.cy.list[["Jaccard"]] <- vegdist(as.data.frame(otu_table(pseq.rarefied)), method="jaccard")
dist.cy.list[["Aitchison"]] <- vegdist(as.data.frame(otu_table(microbiome::transform(pseq.clr, "clr"))), method="euclidean")
dist.cy.list[["wUF"]] <- UniFrac(pseq.rarefied, weighted=TRUE, normalized=TRUE)
dist.cy.list[["uwUF"]] <- UniFrac(pseq.rarefied, weighted=FALSE, normalized=TRUE)

# Calculate general microbiome dispersion.
dispersion <- data.frame(FecesID=sample_names(pseq.rarefied),
                         BC.Dispr=betadisper(dist.cy.list[["Bray-Curtis"]], group=rep("feces", 88))[["distances"]],
                         AIT.Dispr=betadisper(dist.cy.list[["Aitchison"]], group=rep("feces", 88))[["distances"]])
cy.metadata <- merge(cy.metadata, dispersion, by="FecesID", all=TRUE)

rm(alpha.data, dispersion)

# Create a summary of all taxon abundances and prevalences, for subsetting purposes later.
prev.data <- list()
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")
for(i in c(1:6)){
  temp.physeq <- pseq.clr.sub
  if(i<6){temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]])}
  if(i<6){taxa_names(temp.physeq) <- as.data.frame(cbind(tax_table(temp.physeq)))[,i+1]}
  prev <- as.data.frame(microbiome::prevalence(temp.physeq, detection=0, count=TRUE))
  colnames(prev) <- "Prevalence"
  prev$Taxon <- rownames(prev)
  temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  temp <- otu_table(temp.physeq)
  prev$Abundance <- colMeans(temp)
  prev.data[[i]] <- prev
  rm(prev, temp.physeq, temp)
}
names(prev.data) <- c(taxranks, "ASV")
                                       
#### CALCULATE HEALTH INDICES ####
rownames(cy.metadata) <- cy.metadata$FecesID
health.index.data <- cy.metadata[, c("Mass","Snout.to.Base", "Girth", "Cementum.Age",
                                     "Spleen.to.Body.Ratio", "KFI")]

# Choose variables for dimension reduction:
index1 <- c("Mass", "Snout.to.Base", "Girth", "KFI")

index.summary <- data.frame(matrix(ncol=0, nrow=3))
rownames(index.summary) <- c("cases", "PC1%", "PC2%")

# Index 1: Mass, Size, Girth, KFI
health.data.subset <- health.index.data[,index1]
health.data.subset <- health.data.subset[complete.cases(health.data.subset), ]
pca <- rda(health.data.subset)
scores <- scores(pca, display="sites", c(1:3))
colnames(scores) <- c("index1.PC1", "index1.PC2", "index1.PC3")
scores <- cbind(scores, health.data.subset)
scores$FecesID <- rownames(scores)
health.index.cors <- Hmisc::rcorr(as.matrix(scores[,sapply(scores, is.numeric)]), type="pearson")
for(i in 1:length(health.index.cors)){
  health.index.cors[[i]] <- subset(health.index.cors[[i]], rownames(health.index.cors[[i]]) %in% index1)
  health.index.cors[[i]] <- health.index.cors[[i]][,c(1:3)]}
openxlsx::write.xlsx(rbind(health.index.cors[["r"]], vector(length=3), health.index.cors[["P"]]),
           "~/health/index1.cors.xlsx", row.names=TRUE)
cy.metadata <- merge(cy.metadata, scores[, c("FecesID", "index1.PC1")], by="FecesID", all.x=TRUE)
rownames(cy.metadata) <- cy.metadata$FecesID

index.summary$index1 <- c(nrow(health.data.subset),
                          100*pca[["CA"]][["eig"]][[1]]/pca[["CA"]]$tot.chi, 
                          100*pca[["CA"]][["eig"]][[2]]/pca[["CA"]]$tot.chi)

rm(health.data.subset, pca, scores, health.index.cors, health.index.data)

# Calculate z-scores for the new health indices
cy.metadata.zscore <- cy.metadata
for(i in 1:ncol(cy.metadata.zscore)){
  if(is.numeric(cy.metadata.zscore[,i])=="TRUE"){
    cy.metadata.zscore[,i] <- (cy.metadata.zscore[,i] - mean(cy.metadata.zscore[,i], na.rm=TRUE))/(sd(cy.metadata.zscore[,i], na.rm=TRUE))
  }
}

#### SUMMARIZE TAXONOMIC DATA AND HEALTH DATA ####
# Add z-scores to cy.metadata
cy.metadata.zscore.temp <- as.data.frame(cy.metadata.zscore[,c("index1.PC1")])
colnames(cy.metadata.zscore.temp) <- "index1.PC1.z"
cy.metadata.zscore.temp$FecesID <- cy.metadata.zscore$FecesID
cy.metadata <- merge(cy.metadata, cy.metadata.zscore.temp, by="FecesID", all=TRUE)
rownames(cy.metadata) <- cy.metadata$FecesID
rm(cy.metadata.zscore.temp)

# Create bins for health scores 
hist(cy.metadata.zscore$index1.PC1, breaks=10)
cy.metadata.zscore <- dplyr::as_tibble(cy.metadata.zscore) %>% 
  dplyr::mutate(index1.z.3.bins = dplyr::case_when(
    index1.PC1 < -0.5 ~ "unhealthy",
    index1.PC1 >= -0.5 & index1.PC1 <= 0.5 ~ "average",
    index1.PC1 > 0.5 ~ "healthy"))
cy.metadata$index1.z.3.bins <- cy.metadata.zscore$index1.z.3.bins

hist(cy.metadata.zscore$index1.PC1, breaks=10)
cy.metadata.zscore <- dplyr::as_tibble(cy.metadata.zscore) %>% 
  dplyr::mutate(index1.z.5.bins = dplyr::case_when(
    index1.PC1 < -1.5 ~ "v.unhealthy",
    index1.PC1 >= -1.5 & index1.PC1 < -0.5 ~ "unhealthy",
    index1.PC1 >= -0.5 & index1.PC1 <= 0.5 ~ "average",
    index1.PC1 >= 0.5 & index1.PC1 < 1.5 ~ "healthy",
    index1.PC1 >= 1.5 ~ "v.healthy"))
cy.metadata$index1.z.5.bins <- cy.metadata.zscore$index1.z.5.bins

#### STABLE ISOTOPE MIXING MODELS ####
# Identify consumer stable isotope values.
mix <- as.matrix(cy.metadata[,c("d13C", "d15N")])
colnames(mix) <- c("d13C", "d15N")

# Identify source stable isotope values.
s_names <- as.vector(mixing.model.data[["sources"]]$Source)
s_means <- matrix(c(mixing.model.data[["sources"]]$Meand13C,
                    mixing.model.data[["sources"]]$Meand15N),
                  nrow=3, ncol=2)
s_sds <- matrix(c(mixing.model.data[["sources"]]$SDd13C,
                  mixing.model.data[["sources"]]$SDd15N),
                ncol=2, nrow=3)

# Identify discrimination factors for each source.
c_means <- matrix(c(mixing.model.data[["discrimination.factors"]]$Meand13C,
                    mixing.model.data[["discrimination.factors"]]$Meand15N),
                  ncol=2, nrow=3)
c_sds <- matrix(c(mixing.model.data[["discrimination.factors"]]$SDd13C,
                  mixing.model.data[["discrimination.factors"]]$SDd15N),
                ncol=2, nrow=3)

# Identify groups.
grp <- as.character(cy.metadata$Location)

# Run mixing model.
simmr_groups_indv = simmr_load(mixtures=mix,
                               source_names=s_names,
                               source_means=s_means,
                               source_sds=s_sds,
                               correction_means=c_means,
                               correction_sds=c_sds,
                               group=grp)

plot(simmr_groups_indv, group=1:2, xlab=expression(paste(delta^13, "C (\u2030)",sep="")), 
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")))

simmr_groups_indv_out = simmr_mcmc(simmr_groups_indv)

mixing.model.results <- summary(simmr_groups_indv_out, type="statistics", group=c(1:2))$statistics
mixing.model.results <- cbind(mixing.model.results[[1]],
                              mixing.model.results[[2]])
colnames(mixing.model.results) <- c("rural_mean", "rural_sd", "urban_mean", "urban_sd")
mixing.model.results <- subset(mixing.model.results, rownames(mixing.model.results) %in% c("anthro","prey","fruit"))
mixing.model.results <- 100*mixing.model.results
mixing.model.results <- as.data.frame(mixing.model.results)
mixing.model.results$Source <- rownames(mixing.model.results)

rm(c_means, c_sds, s_names, s_means, s_sds, mix, grp, simmr_groups_indv_out, simmr_groups_indv)

#### DIET, HEALTH, AND DIVERSITY MEASURES: CONTROLLED SIGNIFICANCE TESTS ####
health.test.data <- cy.metadata

# Log transform stomach volumes to achieve a normal distribution.
health.test.data$Vol.Prey <- log(health.test.data$Vol.Prey+0.01)
health.test.data$Vol.Anthro <- log(health.test.data$Vol.Anthro+0.01)

# Significance tests for continuous measures (linear regression).
health.signif.tests <- data.frame(matrix(ncol=4, nrow=0))

for(i in c("Mass", "Snout.to.Base", "Girth", "KFI", "Spleen.to.Body.Ratio", 
           "index1.PC1", "d13C","d15N","Vol.Prey","Vol.Anthro", "Shannon", 
           "Shannon_Extrap", "PD", "PD_Extrap", "NTI", "NTI_Raw")){
  p <- Anova(lm(cy.metadata[,i] ~ Sex+Cementum.Age+Location, cy.metadata))
  health.signif.tests[i,1] <- i
  health.signif.tests[i,2] <- p$`F value`[1]
  health.signif.tests[i,3] <- p$Df[1]
  health.signif.tests[i,4] <- p$`Pr(>F)`[1]
  health.signif.tests[i,5] <- p$`F value`[2]
  health.signif.tests[i,6] <- p$Df[2]
  health.signif.tests[i,7] <- p$`Pr(>F)`[2]
  health.signif.tests[i,8] <- p$`F value`[3]
  health.signif.tests[i,9] <- p$Df[3]
  health.signif.tests[i,10] <- p$`Pr(>F)`[3]
  rm(p)
}

colnames(health.signif.tests) <- c("Variable", "Sex_F", "Sex_df", "Sex_p", "Age_F", "Age_df", "Age_p",
                                   "Location_F", "Location_df", "Location_p")
health.summary <- cy.metadata %>% 
  dplyr::group_by(Location) %>% 
  dplyr::summarise_if(is.numeric, .funs=list(mean), na.rm=TRUE)
rownames(health.summary) <- health.summary$Location
health.summary$Location <- NULL
health.summary <- as.data.frame(t(health.summary))

health.signif.tests <- merge(health.summary, health.signif.tests, by=0, all.x=FALSE)
rm(health.summary)

# Test age separately.
Anova(lm(Cementum.Age ~ Sex + Location, cy.metadata)) 

# Test richness using a negative binomial distribution.
Anova(glm.nb(Observed ~ Sex+Cementum.Age+Location, health.test.data))
Anova(glm.nb(Observed_Extrap ~ Sex+Cementum.Age+Location, health.test.data))

# Test for differences in the prevalence of dietary items or E. multilocularis (logistic regressions).
logit.signif.tests <- data.frame(matrix(ncol=4, nrow=0))
for(i in c("PA.Empty", "PA.Anthro", "PA.Prey", "PA.Veg", "Em.PCR.Overall")){
  p <- Anova(glm(health.test.data[,i] ~ Sex+Cementum.Age+Location, 
                 family="binomial"(link="logit"), health.test.data))
  logit.signif.tests[i,1] <- i
  logit.signif.tests[i,2] <- p$`LR Chisq`[1]
  logit.signif.tests[i,3] <- p$Df[1]
  logit.signif.tests[i,4] <- p$`Pr(>Chisq)`[1]
  logit.signif.tests[i,5] <- p$`LR Chisq`[2]
  logit.signif.tests[i,6] <- p$Df[2]
  logit.signif.tests[i,7] <- p$`Pr(>Chisq)`[2]
  logit.signif.tests[i,8] <- p$`LR Chisq`[3]
  logit.signif.tests[i,9] <- p$Df[3]
  logit.signif.tests[i,10] <- p$`Pr(>Chisq)`[3]
  rm(p)
}
colnames(logit.signif.tests) <- c("Variable", "Sex_LRChisq", "Sex_df", "Sex_p", "Age_LRChisq", "Age_df", "Age_p",
                                   "Location_LRChisq", "Location_df", "Location_p")

# Test for the effects on and of E. multilocularis infection on continuous variables.
health.test.data <- subset(health.test.data, is.na(Observed_Extrap)=="FALSE")

emulti.signif.tests <- data.frame(matrix(ncol=4, nrow=0))
for(i in c("Mass", "Snout.to.Base", "Girth", "KFI", "Spleen.to.Body.Ratio", "index1.PC1",
           "d13C","d15N","Vol.Prey","Vol.Anthro", "Shannon", "Shannon_Extrap", "PD", "PD_Extrap", "NTI", "NTI_Raw")){
  p <- Anova(lm(health.test.data[,i] ~ Location+Em.PCR.Overall, health.test.data))
  emulti.signif.tests[i,1] <- i
  emulti.signif.tests[i,2] <- p$`F value`[1]
  emulti.signif.tests[i,3] <- p$Df[1]
  emulti.signif.tests[i,4] <- p$`Pr(>F)`[1]
  emulti.signif.tests[i,5] <- p$`F value`[2]
  emulti.signif.tests[i,6] <- p$Df[2]
  emulti.signif.tests[i,7] <- p$`Pr(>F)`[2]
  rm(p)
}
colnames(emulti.signif.tests) <- c("Variable", "Location_F", "Location_df", "Location_p",
                                  "Emulti_F", "Emulti_df", "Emulti_p")

# Test richness using a negative binomial distribution.
Anova(glm.nb(Observed ~ Location+Em.PCR.Overall, health.test.data))
Anova(glm.nb(Observed_Extrap ~ Location+Em.PCR.Overall, health.test.data))

rm(health.test.data)

#### RICHNESS, DIVERSITY, PD, and NTI: MODELS EVALUATING ALL PREDICTORS ####
# This section uses model-averaged coefficients to compare the importance of predictors of alpha diversity.
# Prepare the data.
adiv.model.results <- list()
model.data <- cy.metadata[,c("Mass", "Snout.to.Base", "Girth", "KFI", "Spleen.to.Body.Ratio", "Cementum.Age",
                             "index1.PC1", "d13C","d15N","Vol.Prey","Vol.Anthro",
                             "Observed", "Observed_Extrap", "Shannon", "Shannon_Extrap", 
                             "PD", "PD_Extrap", "NTI", "NTI_Raw", "Location", "Em.PCR.Overall", "Sex")]
model.data <- subset(model.data, is.na(Observed)=="FALSE")
model.data[,c(1:11)] <- scale(model.data[,c(1:11)]) # Scale data.

# Use a negative binomial model to model species richness.
for(p in c("Observed", "Observed_Extrap")){
  # Construct an initial model and evaluate all subsets.
  model <- glm.nb(model.data[,p] ~ Sex+Cementum.Age+Location+index1.PC1+d13C+d15N+Vol.Prey+Vol.Anthro+Spleen.to.Body.Ratio+Em.PCR.Overall,
                  model.data, na.action="na.fail")
  model.dredge <- dredge(model)
  num.models <- nrow(subset(model.dredge, delta < 2))
  
  # Calculate model-averaged coefficients standardized by partial standard deviation.
  temp2 <- model.avg(model.dredge, subset=delta < 2, beta="partial.sd")
  temp2 <- cbind(confint(temp2, level=0.95, full=TRUE)[,1],
                 confint(temp2, level=0.5, full=TRUE)[,1],
                 temp2[["coefficients"]][1,],
                 confint(temp2, level=0.5, full=TRUE)[,2],
                 confint(temp2, level=0.95, full=TRUE)[,2])
  temp2 <- as.data.frame(temp2)
  colnames(temp2) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")
  
  # Calculate summed weights for each variable, across all models.
  model.dredge <- as.data.frame(model.dredge)
  weights <- data.frame(matrix(nrow=0, ncol=1))
  for(j in c(2:11)){
    temp <- model.dredge
    temp <- subset(temp, is.na(temp[,j])=="FALSE")
    weights[j-1,1] <- sum(temp[,ncol(temp)])
    rm(temp)
  }
  rownames(weights) <- colnames(model.dredge[,c(2:11)])
  colnames(weights) <- "cum_weight"
  
  temp2$merge <- rownames(temp2)
  weights$merge <- as.factor(rownames(weights))
  levels(weights$merge)[levels(weights$merge)=="Em.PCR.Overall"] <- "Em.PCR.OverallY"
  levels(weights$merge)[levels(weights$merge)=="Location"] <- "Locationurban"
  levels(weights$merge)[levels(weights$merge)=="Sex"] <- "SexM"
  
  results <- merge(temp2, weights, by="merge", all=TRUE)
  results <- subset(results, merge !="(Intercept)")
  results <- results[order(-results$cum_weight),]
  results$merge <- forcats::fct_inorder(results$merge)
  adiv.model.results[[p]] <- results
  rm(results, temp2, model.dredge, weights, num.models)
}

# Use a linear model to model diversity, PD, and NTI.
for(p in c("Shannon", "Shannon_Extrap", "PD", "PD_Extrap", "NTI", "NTI_Raw")){
  # Construct an initial model and evaluate all subsets.
  model <- lm(model.data[,p] ~ Sex+Cementum.Age+Location+index1.PC1+d13C+d15N+Vol.Prey+Vol.Anthro+Spleen.to.Body.Ratio+Em.PCR.Overall,
              model.data, na.action="na.fail")
  model.dredge <- dredge(model)
  num.models <- nrow(subset(model.dredge, delta < 2))
  
  # Calculate model-averaged coefficients standardized by partial standard deviation.
  temp2 <- model.avg(model.dredge, subset=delta < 2, beta="partial.sd")
  temp2 <- cbind(confint(temp2, level=0.95, full=TRUE)[,1],
                 confint(temp2, level=0.5, full=TRUE)[,1],
                 temp2[["coefficients"]][1,],
                 confint(temp2, level=0.5, full=TRUE)[,2],
                 confint(temp2, level=0.95, full=TRUE)[,2])
  temp2 <- as.data.frame(temp2)
  colnames(temp2) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")
  
  # Calculate summed weights across all models.
  model.dredge <- as.data.frame(model.dredge)
  weights <- data.frame(matrix(nrow=0, ncol=1))
  for(j in c(2:11)){
    temp <- model.dredge
    temp <- subset(temp, is.na(temp[,j])=="FALSE")
    weights[j-1,1] <- sum(temp[,ncol(temp)])
    rm(temp)
  }
  rownames(weights) <- colnames(model.dredge[,c(2:11)])
  colnames(weights) <- "cum_weight"
  
  temp2$merge <- rownames(temp2)
  weights$merge <- as.factor(rownames(weights))
  levels(weights$merge)[levels(weights$merge)=="Em.PCR.Overall"] <- "Em.PCR.OverallY"
  levels(weights$merge)[levels(weights$merge)=="Location"] <- "Locationurban"
  levels(weights$merge)[levels(weights$merge)=="Sex"] <- "SexM"
  
  results <- merge(temp2, weights, by="merge", all=TRUE)
  results <- subset(results, merge !="(Intercept)")
  results <- results[order(-results$cum_weight),]
  results$merge <- forcats::fct_inorder(results$merge)
  adiv.model.results[[p]] <- results
  rm(results, temp2, model.dredge, weights, num.models)
}

rm(model, model.data)

#### CONTROLLED DIFFERENTIAL ABUNDANCE TESTS ####
# (1) ~ Sex + Age ####
diff.abund.just.sex.and.age<- list()
for(k in c("Sex")){
  for(i in 1:5){
    # Agglomerate taxa to the studied rank.
    temp.physeq <- tax_glom(pseq.clr.sub, taxrank=taxranks[[i]], NArm=FALSE)
    temp.data <- subset(cy.metadata, FecesID %in% sample_names(temp.physeq))
    temp.otu.table <- as.data.frame(t(otu_table(temp.physeq)))
    temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
    rownames(temp.otu.table) <- temp.tax[,i+1]
    
    covariates <- data.frame("Sex" = temp.data$Sex,
                             "Age" = temp.data$Cementum.Age)
    mm <- model.matrix(~ Sex + Age, covariates)
    x <- aldex.clr(temp.otu.table, mm, denom="all", verbose=F)
    glm.test <- aldex.glm(x, mm)
    
    # Calculate summary statistics for each taxon.
    # Transform to relative abundance.
    temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
    
    # Extract OTU table and rename features with taxon names.
    temp.transform.otu <- as.data.frame(otu_table(temp.transform))
    if(i<6){colnames(temp.transform.otu) <- temp.tax[,i+1]}
    
    # Append intestinal site.
    temp.transform.otu[,k] <- sample_data(temp.transform)[,k]
    
    # Define a data frame where calculation outputs can be stored. Number of rows is 7 intestinal sites * 3 statistics per site (mean, sd, median).
    summary.data <- data.frame(matrix(nrow=7)) 
    
    # For each taxon...
    for(j in 1:(ncol(temp.transform.otu)-1)){
      # Calculate mean relative abundance.
      means <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), mean)
      means$Group.1 <- NULL
      rownames(means) <- levels(temp.transform.otu[,k])
      rownames(means) <- paste("mean", rownames(means))
      
      # Calculate standard deviation of relative abundance.
      sds <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), sd)
      sds$Group.1 <- NULL
      rownames(sds) <- levels(temp.transform.otu[,k])
      rownames(sds) <- paste("stdev", rownames(sds))
      
      # Calculate median relative abundance.
      medians <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), median)
      medians$Group.1 <- NULL
      rownames(medians) <- levels(temp.transform.otu[,k])
      rownames(medians) <- paste("median", rownames(medians))
      
      # Effect size
      eff <- effsize::cohen.d(temp.transform.otu[,j] ~ temp.transform.otu[,k], hedges.correction=TRUE)$estimate
      
      # Combine mean, sd, and medians into one data frame.
      temp.summary <- rbind(means, sds, medians, eff)
      colnames(temp.summary) <- colnames(temp.transform.otu)[j]
      
      # Append results for each taxon into one data frame.
      summary.data <- cbind(summary.data, temp.summary)
    }
    
    # Remove artefactual column.
    summary.data[,1] <- NULL
    
    # Add prevalence information to summary data.
    prev <- microbiome::prevalence(temp.physeq, detection=0, count=TRUE)
    
    # Add taxonomic information to summary data.
    summary.data <- as.data.frame(t(summary.data))
    summary.data <- cbind(temp.tax, prev, summary.data)
    if(i<6){rownames(summary.data) <- summary.data[,i+1]}
    
    # Combine differential abundance results and summary data.
    diff.abund.just.sex.and.age[[i]] <- transform(merge(glm.test, summary.data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    rm(means, sds, medians, sites, temp.transform, temp.transform.otu, summary.data, prev,
       temp.physeq, temp.otu.table, temp.summary, temp.tax)
  }
  names(diff.abund.just.sex.and.age) <- taxranks
}


# (2) ~ Sex + Age + Location ####
diff.abund.location <- list()
for(k in c("Location")){
  for(i in 1:6){
    # Agglomerate taxa to the studied rank.
    temp.physeq <- pseq.clr.sub
    if(i<6){temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]], NArm=FALSE)}
    temp.data <- subset(cy.metadata, FecesID %in% sample_names(temp.physeq))
    temp.otu.table <- as.data.frame(t(otu_table(temp.physeq)))
    temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
    if(i<6){rownames(temp.otu.table) <- temp.tax[,i+1]}
    
    covariates <- data.frame("Sex" = temp.data$Sex,
                             "Age" = temp.data$Cementum.Age,
                             "Location" = temp.data$Location)
    mm <- model.matrix(~ Sex + Age + Location, covariates)
    x <- aldex.clr(temp.otu.table, mm, denom="all", verbose=F)
    glm.test <- aldex.glm(x, mm)
    
    # Calculate summary statistics for each taxon.
    # Transform to relative abundance.
    temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
    
    # Extract OTU table and rename features with taxon names.
    temp.transform.otu <- as.data.frame(otu_table(temp.transform))
    if(i<6){colnames(temp.transform.otu) <- temp.tax[,i+1]}
    
    # Append intestinal site.
    temp.transform.otu[,k] <- sample_data(temp.transform)[,k]
    
    # Define a data frame where calculation outputs can be stored. Number of rows is 7 intestinal sites * 3 statistics per site (mean, sd, median).
    summary.data <- data.frame(matrix(nrow=7)) 
    
    # For each taxon...
    for(j in 1:(ncol(temp.transform.otu)-1)){
      # Calculate mean relative abundance.
      means <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), mean)
      means$Group.1 <- NULL
      rownames(means) <- levels(temp.transform.otu[,k])
      rownames(means) <- paste("mean", rownames(means))
      
      # Calculate standard deviation of relative abundance.
      sds <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), sd)
      sds$Group.1 <- NULL
      rownames(sds) <- levels(temp.transform.otu[,k])
      rownames(sds) <- paste("stdev", rownames(sds))
      
      # Calculate median relative abundance.
      medians <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), median)
      medians$Group.1 <- NULL
      rownames(medians) <- levels(temp.transform.otu[,k])
      rownames(medians) <- paste("median", rownames(medians))
      
      # Effect size
      eff <- effsize::cohen.d(temp.transform.otu[,j] ~ temp.transform.otu[,k], hedges.correction=TRUE)$estimate
      
      # Combine mean, sd, and medians into one data frame.
      temp.summary <- rbind(means, sds, medians, eff)
      
      colnames(temp.summary) <- colnames(temp.transform.otu)[j]
      
      # Append results for each taxon into one data frame.
      summary.data <- cbind(summary.data, temp.summary)
    }
    
    # Remove artefactual column.
    summary.data[,1] <- NULL
    
    # Add prevalence information to summary data.
    prev <- microbiome::prevalence(temp.physeq, detection=0, count=TRUE)
    
    # Add taxonomic information to summary data.
    summary.data <- as.data.frame(t(summary.data))
    summary.data <- cbind(temp.tax, prev, summary.data)
    if(i<6){rownames(summary.data) <- summary.data[,i+1]}
    
    # Combine differential abundance results and summary data.
    diff.abund.location[[i]] <- transform(merge(glm.test, summary.data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    rm(means, sds, medians, sites, temp.transform, temp.transform.otu, summary.data, prev,
       temp.physeq, temp.otu.table, temp.summary, temp.tax)
  }
  names(diff.abund.location) <- c(taxranks, "ASV")
}

# (3) ~ Location + Em ####
diff.abund.loc.Em <- list()
for(k in c("Em.PCR.Overall")){
  for(i in 1:5){
    # Agglomerate taxa to the studied rank.
    temp.physeq <- tax_glom(pseq.clr.sub, taxrank=taxranks[[i]], NArm=FALSE)
    temp.data <- subset(cy.metadata, FecesID %in% sample_names(temp.physeq))
    temp.otu.table <- as.data.frame(t(otu_table(temp.physeq)))
    temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
    rownames(temp.otu.table) <- temp.tax[,i+1]
    
    covariates <- data.frame("Location" = temp.data$Location,
                             "Em.PCR.Overall" = temp.data$Em.PCR.Overall)
    mm <- model.matrix(~ Location + Em.PCR.Overall, covariates)
    x <- aldex.clr(temp.otu.table, mm, denom="all", verbose=F)
    glm.test <- aldex.glm(x, mm)
    
    # Calculate summary statistics for each taxon.
    # Transform to relative abundance.
    temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
    
    # Extract OTU table and rename features with taxon names.
    temp.transform.otu <- as.data.frame(otu_table(temp.transform))
    if(i<6){colnames(temp.transform.otu) <- temp.tax[,i+1]}
    
    # Append intestinal site.
    temp.transform.otu[,k] <- sample_data(temp.transform)[,k]
    
    # Define a data frame where calculation outputs can be stored. Number of rows is 7 intestinal sites * 3 statistics per site (mean, sd, median).
    summary.data <- data.frame(matrix(nrow=7)) 
    
    # For each taxon...
    for(j in 1:(ncol(temp.transform.otu)-1)){
      # Calculate mean relative abundance.
      means <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), mean)
      means$Group.1 <- NULL
      rownames(means) <- levels(temp.transform.otu[,k])
      rownames(means) <- paste("mean", rownames(means))
      
      # Calculate standard deviation of relative abundance.
      sds <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), sd)
      sds$Group.1 <- NULL
      rownames(sds) <- levels(temp.transform.otu[,k])
      rownames(sds) <- paste("stdev", rownames(sds))
      
      # Calculate median relative abundance.
      medians <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), median)
      medians$Group.1 <- NULL
      rownames(medians) <- levels(temp.transform.otu[,k])
      rownames(medians) <- paste("median", rownames(medians))
      
      # Effect size
      eff <- effsize::cohen.d(temp.transform.otu[,j] ~ temp.transform.otu[,k], hedges.correction=TRUE)$estimate
      
      # Combine mean, sd, and medians into one data frame.
      temp.summary <- rbind(means, sds, medians, eff)
      colnames(temp.summary) <- colnames(temp.transform.otu)[j]
      
      # Append results for each taxon into one data frame.
      summary.data <- cbind(summary.data, temp.summary)
    }
    
    # Remove artefactual column.
    summary.data[,1] <- NULL
    
    # Add prevalence information to summary data.
    prev <- microbiome::prevalence(temp.physeq, detection=0, count=TRUE)
    
    # Add taxonomic information to summary data.
    summary.data <- as.data.frame(t(summary.data))
    summary.data <- cbind(temp.tax, prev, summary.data)
    if(i<6){rownames(summary.data) <- summary.data[,i+1]}
    
    # Combine differential abundance results and summary data.
    diff.abund.loc.Em[[i]] <- transform(merge(glm.test, summary.data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    rm(means, sds, medians, sites, temp.transform, temp.transform.otu, summary.data, prev,
       temp.physeq, temp.otu.table, temp.summary, temp.tax)
  }
  names(diff.abund.loc.Em) <- taxranks
}


#### RANDOM FOREST MODELS: PREDICTING LOCATION AND TESTING MISCLASSIFIED SAMPLES ####
# This code is adapted from https://rpubs.com/michberr/randomforestmicrobe
# (1) Run the predictive model and evaluate accuracy. ####
# Prepare data.
temp <- as.data.frame(otu_table(transform_sample_counts(pseq.clr, function(x) 100*x/sum(x))))
temp <- as.data.frame(t(temp))
temp$Mean <- rowMeans(temp)
temp <- subset(temp, Mean > 0.01)
predictors <- as.data.frame(otu_table(microbiome::transform(subset_taxa(pseq.clr, taxa_names(pseq.clr) %in% rownames(temp)), "clr")))

response <- subset(cy.metadata, FecesID %in% rownames(predictors))$Location
rf.data <- data.frame(response, predictors)

# Run the model.
set.seed(100)
rf.raw.results <- randomForest(response~., data = rf.data, ntree = 1000)
print(rf.raw.results)

# Extract the 40 most important ASVs, based on their mean decrease in the Gini coefficient.
rf.importance <- randomForest::importance(rf.raw.results)
rf.importance <- data.frame(predictors = rownames(rf.importance), rf.importance)
rf.importance <- dplyr::arrange(rf.importance, desc(MeanDecreaseGini))
rf.importance$predictors <- factor(rf.importance$predictors,
                                   levels = rf.importance$predictors)
rf.importance <- rf.importance[1:40, ]
rownames(rf.importance) <- rf.importance$predictors
rf.importance$predictors <- NULL

# Add taxonomy information to the predictors.
temp <- tax_table(pseq.clr.sub)
temp <- subset(temp, rownames(temp) %in% rownames(rf.importance))
rf.importance <- merge(rf.importance, temp, by="row.names", all=TRUE)

# Organize by Gini score.  
rf.importance <- dplyr::arrange(rf.importance, desc(MeanDecreaseGini))

# Add information about taxon abundances for each predictor.
temp.pseq <- transform_sample_counts(pseq.clr, function(x) 100*x/sum(x))
temp.pseq <- subset_taxa(temp.pseq, taxa_names(temp.pseq) %in% rf.importance$Row.names)
temp.pseq.otu <- as.data.frame(otu_table(temp.pseq))

# Calculate mean taxon abundances in urban and rural coyotes and overall mean abundance in the sample.
temp.pseq.otu$Location <- subset(cy.metadata, is.na(Observed)=="FALSE")$Location
temp.pseq.tab <- temp.pseq.otu %>% dplyr::group_by(Location) %>% dplyr::summarise_all(list(mean=mean))
temp.pseq.tab2 <- temp.pseq.otu %>% dplyr::summarise_if(is.numeric, .funs=list(mean=mean))
temp.pseq.tab <- as.data.frame(temp.pseq.tab)
temp.pseq.tab2 <- as.data.frame(temp.pseq.tab2)
rownames(temp.pseq.tab) <- temp.pseq.tab$Location
rownames(temp.pseq.tab2) <- "Overall"

temp.pseq.tab[,1] <- NULL

colnames(temp.pseq.tab) <- taxa_names(temp.pseq)
colnames(temp.pseq.tab2) <- taxa_names(temp.pseq)

temp.pseq.tab <- as.data.frame(t(rbind(temp.pseq.tab, temp.pseq.tab2)))

rownames(rf.importance) <- rf.importance$Row.names
rf.importance <- merge(rf.importance, temp.pseq.tab, by="row.names", all=TRUE)
rownames(rf.importance) <- rf.importance$Row.names
rf.importance$Row.names <- NULL
rf.importance$Row.names <- NULL
rm(temp.pseq, temp.pseq.otu, temp.pseq.tab, temp.pseq.tab2)

# Store data and clean workspace.
rm(predictors, response, temp)

# (2) Diet/condition/alpha diversity differences between correctly and incorrectly classified urban coyotes ####
misclassified <- data.frame(Predicted = rf.raw.results$predicted, 
                            Actual = subset(cy.metadata, is.na(Observed)=="FALSE")$Location)
misclassified$Equal <- misclassified$Predicted == misclassified$Actual

# Create a data frame with coyote metadata and their random forest classification accuracy ('Equal')
rf.test <- cbind(subset(cy.metadata, FecesID %in% sample_names(pseq.clr)), misclassified$Equal)
colnames(rf.test) <- c(colnames(rf.test)[1:(ncol(rf.test)-1)], "Equal")
rownames(rf.test) <- rf.test$FecesID

rural.data <- subset(rf.test, Microbiome=="Y" & Location=="rural")
urban.data <- subset(rf.test, Microbiome=="Y" & Location=="urban")

rural.data <- rural.data[,c("Equal","Sex","Em.PCR.Overall","Mass","Snout.to.Base", "Girth", "KFI","index1.PC1",
                            "Spleen.to.Body.Ratio", "Cementum.Age","d13C","d15N","Vol.Prey","Vol.Anthro",
                            "Observed", "Observed_Extrap", "Shannon", "Shannon_Extrap", "PD", "PD_Extrap", "NTI", "NTI_Raw")]
urban.data <- urban.data[,c("Equal","Sex","Em.PCR.Overall","Mass","Snout.to.Base", "Girth", "KFI","index1.PC1",
                            "Spleen.to.Body.Ratio", "Cementum.Age","d13C","d15N","Vol.Prey","Vol.Anthro",
                            "Observed", "Observed_Extrap", "Shannon", "Shannon_Extrap", "PD", "PD_Extrap", "NTI", "NTI_Raw")]

test.data <- urban.data
test.data$Vol.Prey <- log(test.data$Vol.Prey+0.01)
test.data$Vol.Anthro <- log(test.data$Vol.Anthro+0.01)

# Test for differences in condition, diet, and microbiome measures.
rf.miscl.signif.tests <- data.frame(matrix(ncol=4, nrow=0))
for(i in c("Mass", "Snout.to.Base", "Girth", "KFI", "Spleen.to.Body.Ratio", "Cementum.Age", "index1.PC1",
           "d13C","d15N","Vol.Prey","Vol.Anthro","Observed", "Observed_Extrap", "Shannon", "Shannon_Extrap", "PD",
           "PD_Extrap", "NTI", "NTI_Raw")){
  p <- Anova(lm(urban.data[,i] ~ Equal, test.data))
  rf.miscl.signif.tests[i,1] <- mean(subset(test.data, Equal=="TRUE")[,i], na.rm=TRUE)
  rf.miscl.signif.tests[i,2] <- mean(subset(test.data, Equal=="FALSE")[,i], na.rm=TRUE)
  rf.miscl.signif.tests[i,3] <- i
  rf.miscl.signif.tests[i,4] <- p$`F value`[1]
  rf.miscl.signif.tests[i,5] <- p$Df[1]
  rf.miscl.signif.tests[i,6] <- p$`Pr(>F)`[1]
  rm(p)
}
colnames(rf.miscl.signif.tests) <- c("Urban_Correct", "Urban_Incorrect", "Measure", "F", "df", "p")

rural.data[,c(4:ncol(rural.data))] <- scale(rural.data[,c(4:ncol(rural.data))])
urban.data[,c(4:ncol(urban.data))] <- scale(urban.data[,c(4:ncol(urban.data))])

rm(rf.test)

# Evaluate the most important differences using information criteria.
# Note that PD and NTI cannot go into the model because they cause perfect separation.
model <- glm(Equal ~ Sex+Cementum.Age+Spleen.to.Body.Ratio+index1.PC1+d13C+d15N+Vol.Anthro+Vol.Prey+Observed_Extrap+Shannon_Extrap,
             family="binomial"(link="logit"), data=urban.data, na.action="na.fail", maxit=100)
model.dredge <- dredge(model)

# Calculate model-averaged coefficients standardized by partial standard deviation.
temp2 <- model.avg(model.dredge, subset=delta < 2, beta="partial.sd")
temp2 <- cbind(confint(temp2, level=0.95, full=TRUE)[,1],
               confint(temp2, level=0.5, full=TRUE)[,1],
               temp2[["coefficients"]][1,],
               confint(temp2, level=0.5, full=TRUE)[,2],
               confint(temp2, level=0.95, full=TRUE)[,2])
temp2 <- as.data.frame(temp2)
colnames(temp2) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")

# Calculate summed weights for each variable, across all models.
model.dredge <- as.data.frame(model.dredge)
weights <- data.frame(matrix(nrow=0, ncol=1))
for(j in c(2:11)){
  temp <- model.dredge
  temp <- subset(temp, is.na(temp[,j])=="FALSE")
  weights[j-1,1] <- sum(temp[,ncol(temp)])
  rm(temp)
}
rownames(weights) <- colnames(model.dredge[,c(2:11)])
colnames(weights) <- "cum_weight"

temp2$merge <- rownames(temp2)
weights$merge <- as.factor(rownames(weights))
levels(weights$merge)[levels(weights$merge)=="Location"] <- "Locationurban"
levels(weights$merge)[levels(weights$merge)=="Sex"] <- "SexM"

rf.model.results <- merge(temp2, weights, by="merge", all=TRUE)
rf.model.results <- subset(rf.model.results, merge !="(Intercept)")
rf.model.results <- rf.model.results[order(-rf.model.results$cum_weight),]
rf.model.results$merge <- forcats::fct_inorder(rf.model.results$merge)
rm(model, temp2, model.dredge, weights)

# (3) Differential abundance between correctly and incorrectly classified coyotes. ####
rf.test.diff.urban <- list()
for(k in c("Equal")){
  for(i in 1:6){
    temp.physeq <- subset_samples(pseq.clr.sub, sample_names(pseq.clr.sub) %in% rownames(urban.data))
    
    # Agglomerate taxa to the studied rank.
    if(i<6){temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]], NArm=FALSE)}
    temp.data <- subset(cy.metadata, FecesID %in% sample_names(temp.physeq))
    temp.otu.table <- as.data.frame(t(otu_table(temp.physeq)))
    temp.tax <- as.data.frame(cbind(tax_table(temp.physeq)))
    if(i<6){rownames(temp.otu.table) <- temp.tax[,i+1]}
    
    covariates <- data.frame("Sex" = temp.data$Sex,
                             "Age" = temp.data$Cementum.Age,
                             "Equal" = urban.data$Equal)
    mm <- model.matrix(~ Sex + Age + Equal, covariates)
    x <- aldex.clr(temp.otu.table, mm, denom="all", verbose=F)
    glm.test <- aldex.glm(x, mm)
    
    # Calculate summary statistics for each taxon.
    # Transform to relative abundance.
    temp.transform <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
    
    # Extract OTU table and rename features with taxon names.
    temp.transform.otu <- as.data.frame(otu_table(temp.transform))
    if(i<6){colnames(temp.transform.otu) <- temp.tax[,i+1]}
    
    # Append intestinal site.
    temp.transform.otu[,k] <- urban.data$Equal
    temp.transform.otu[,k] <- as.factor(temp.transform.otu[,k])
    
    # Define a data frame where calculation outputs can be stored. Number of rows is 7 intestinal sites * 3 statistics per site (mean, sd, median).
    summary.data <- data.frame(matrix(nrow=7)) 
    
    # For each taxon...
    for(j in 1:(ncol(temp.transform.otu)-1)){
      # Calculate mean relative abundance.
      means <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), mean)
      means$Group.1 <- NULL
      rownames(means) <- levels(temp.transform.otu[,k])
      rownames(means) <- paste("mean", rownames(means))
      
      # Calculate standard deviation of relative abundance.
      sds <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), sd)
      sds$Group.1 <- NULL
      rownames(sds) <- levels(temp.transform.otu[,k])
      rownames(sds) <- paste("stdev", rownames(sds))
      
      # Calculate median relative abundance.
      medians <- aggregate(temp.transform.otu[,j], list(temp.transform.otu[,k]), median)
      medians$Group.1 <- NULL
      rownames(medians) <- levels(temp.transform.otu[,k])
      rownames(medians) <- paste("median", rownames(medians))
      
      # Effect size
      eff <- effsize::cohen.d(temp.transform.otu[,j] ~ temp.transform.otu[,k], hedges.correction=TRUE)$estimate
      
      # Combine mean, sd, and medians into one data frame.
      temp.summary <- rbind(means, sds, medians, eff)
      
      colnames(temp.summary) <- colnames(temp.transform.otu)[j]
      
      # Append results for each taxon into one data frame.
      summary.data <- cbind(summary.data, temp.summary)
    }
    
    # Remove artefactual column.
    summary.data[,1] <- NULL
    
    # Add prevalence information to summary data.
    prev <- microbiome::prevalence(temp.physeq, detection=0, count=TRUE)
    
    # Add taxonomic information to summary data.
    summary.data <- as.data.frame(t(summary.data))
    summary.data <- cbind(temp.tax, prev, summary.data)
    if(i<6){rownames(summary.data) <- summary.data[,i+1]}
    
    # Combine differential abundance results and summary data.
    rf.test.diff.urban[[i]] <- transform(merge(glm.test, summary.data, by=0, all=TRUE), row.names=Row.names, Row.names=NULL)
    rm(means, sds, medians, sites, temp.transform, temp.transform.otu, summary.data, prev,
       temp.physeq, temp.otu.table, temp.summary, temp.tax)
  }
  names(rf.test.diff.urban) <- c(taxranks, "ASV")
}

#### DIMENSIONAL REDUCTION APPROACHES: PERMANOVA, ENVFIT ####
cy.metadata.pseq <- subset(cy.metadata, FecesID %in% sample_names(pseq.clr))

# PERMANOVA for sex, age, location.
permanova.sex.age.loc <- list()
for(i in 1:length(dist.cy.list)){
  permanova.sex.age.loc[[i]] <- adonis(dist.cy.list[[i]] ~ Sex+Cementum.Age+Location, 
                                       cy.metadata.pseq, permutations=1000)$aov.tab}
names(permanova.sex.age.loc) <- names(dist.cy.list)

# PERMANOVA for E. multilocularis.
permanova.em <- list()
for(i in 1:length(dist.cy.list)){
  permanova.em[[i]] <- adonis(dist.cy.list[[i]] ~ Em.PCR.Overall, 
                              cy.metadata.pseq, permutations=1000)$aov.tab}
names(permanova.em) <- names(dist.cy.list)

# Calculate ordinations (PCoA for BC, Jac, w/uwUF distances; PCA for Aitchison distance).
community.ords <- list()
for(i in 1:length(dist.cy.list)){
  if(i %in% c(1,2,4,5)){
    community.ords[[i]] <- capscale(dist.cy.list[[i]] ~ 1)
  }
  if(i==3){
    temp <- as.data.frame(otu_table(microbiome::transform(pseq.clr, "clr")))
    community.ords[[i]] <- rda(temp)
    rm(temp)
  }
}
names(community.ords) <- names(dist.cy.list)

# Fit vectors onto the ordinations.
temp.data <- cy.metadata.pseq[,c("Spleen.to.Body.Ratio","index1.PC1","Vol.Anthro","Vol.Prey","d13C","d15N",
                                 "Cementum.Age","Sex","Em.PCR.Overall","Location")]
temp.data$Location <- factor(temp.data$Location)
community.envfit <- list()
for(i in 1:length(dist.cy.list)){
  community.envfit[[i]] <- envfit(community.ords[[i]], temp.data, permutations=1000, na.rm=TRUE)
}
names(community.envfit) <- names(dist.cy.list)

# Ordination scores
community.scores <- list()
for(i in 1:length(dist.cy.list)){
  community.scores[[i]] <- as.data.frame(scores(community.ords[[i]], display = "sites", choices=c(1:3)))
  community.scores[[i]] <- cbind(temp.data, community.scores[[i]])
}
names(community.scores) <- names(dist.cy.list)

# Calculate error ellipses and plot to visualize differences.
community.ellipses <- list()
for(i in 1:length(community.ords)){
  plot.new()
  temp <- ordiellipse(community.ords[[i]], temp.data$Location,
                      display="sites", kind="sd", label=T, conf=0.95)
  community.ellipses[[i]] <- data.frame()
  for(g in levels(temp.data$Location)){
    community.ellipses[[i]] <- rbind(community.ellipses[[i]], 
                                     cbind(as.data.frame(with(temp.data[temp.data$Location==g,],
                                                                 veganCovEllipse(temp[[g]]$cov,
                                                                                 temp[[g]]$center,
                                                                                 temp[[g]]$scale)))
                                              ,Type=g))
  }
  rm(temp)
}
names(community.ellipses) <- names(dist.cy.list)

community.plots <- list()
for(i in 1:length(community.scores)){
  plot <- ggplot(community.scores[[i]], aes_string(x=community.scores[[i]][,11], 
                                                            y=community.scores[[i]][,12])) +
    geom_point(aes(color=Location, shape=Location), size=2) +
    geom_path(data=community.ellipses[[i]], aes(x=community.ellipses[[i]][,1],
                                                y=community.ellipses[[i]][,2], color=Type),
              size=0.5, linetype=1) +
    scale_color_manual(values=c("forestgreen", "purple")) + theme_bw()
  community.plots[[i]] <- plot
  rm(plot)
}

#### UNIVARIATE CORRELATIONS: TAXON ABUNDANCE vs. METADATA VARIABLES ####
# (1) Subset to highly abundant taxa for this analysis. ####
taxranks <- c("Phylum", "Class", "Order", "Family", "Genus")

taxnames <- list()
for(i in 1:6){
  if(i==1){
    taxnames[[1]] <- c("Actinobacteria","Bacteroidetes","Firmicutes","Fusobacteria","Proteobacteria")
  }
  if(i==2){
    temp <- prev.data[[2]]
    remove <- subset(temp, grepl('Uncl_', Taxon) )$Taxon
    temp <- subset(temp, !(Taxon %in% remove) & Abundance > 0.1)
    temp <- temp[order(-temp$Prevalence), ]
    temp$Taxon <- as.character(temp$Taxon)
    taxnames[[2]] <- temp$Taxon
  }
  if(i==3){
    temp <- prev.data[[3]]
    remove <- subset(temp, grepl('Uncl_', Taxon) )$Taxon
    temp <- subset(temp, !(Taxon %in% remove) & Abundance > 0.1)
    temp <- temp[order(-temp$Prevalence), ]
    temp$Taxon <- as.character(temp$Taxon)
    taxnames[[3]] <- temp$Taxon
  }
  if(i==4){
    temp <- prev.data[[4]]
    remove <- subset(temp, grepl('Uncl_', Taxon) )$Taxon
    temp <- subset(temp, !(Taxon %in% remove) & Abundance > 0.1)
    temp <- temp[order(-temp$Prevalence), ]
    temp <- temp[c(1:19), ]
    temp$Taxon <- as.character(temp$Taxon)
    taxnames[[4]] <- temp$Taxon
  }
  if(i==5){
    temp <- prev.data[[5]]
    remove <- subset(temp, grepl('Uncl_', Taxon) )$Taxon
    temp <- subset(temp, !(Taxon %in% remove) & Abundance > 0.1)
    temp <- temp[order(-temp$Prevalence), ]
    temp <- temp[c(1:25), ]
    temp$Taxon <- as.character(temp$Taxon)
    taxnames[[5]] <- temp$Taxon
  }
  if(i==6){
    temp <- prev.data[[6]]
    remove <- subset(temp, grepl('Uncl_', Taxon) )$Taxon
    temp <- subset(temp, !(Taxon %in% remove) & Abundance > 0.5)
    temp <- temp[order(-temp$Prevalence), ]
    temp <- temp[c(1:36), ]
    temp$Taxon <- as.character(temp$Taxon)
    taxnames[[6]] <- temp$Taxon
  }
}

rm(temp, remove)

names(taxnames) <- c(taxranks, "ASV")

feces.strict.sum.rel <- list()
for(i in 1:6){
  temp.physeq <- pseq.clr.sub
  
  # Agglomerate to each taxonomic rank and calculate relative abundance.
  if(i<6){temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]], NArm=FALSE)}
  temp.physeq <- transform_sample_counts(temp.physeq, function(x) 100*x/sum(x))
  
  # Extract and label OTU table.
  temp.otu.tab <- as.data.frame(otu_table(temp.physeq))
  if(i<6){colnames(temp.otu.tab) <- as.data.frame(cbind(tax_table(temp.physeq)))[,i+1]}
  
  # Only consider the top 50% of taxa... (i.e. the most abundant 50%)
  temp.otu.tab <- as.data.frame(t(temp.otu.tab))
  temp.otu.tab <- subset(temp.otu.tab, rownames(temp.otu.tab) %in% taxnames[[i]])
  temp.otu.tab <- as.data.frame(t(temp.otu.tab))
  
  # Add to list.
  feces.strict.sum.rel[[i]] <- temp.otu.tab
  rm(temp.physeq, temp.otu.tab)
}
names(feces.strict.sum.rel) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV") 

feces.strict.sum.clr <- list()
for(i in 1:6){
  temp.physeq <- pseq.clr.sub
  
  # Agglomerate to each taxonomic rank and calculate CLR-transformed abundance.
  if(i<6){temp.physeq <- tax_glom(temp.physeq, taxrank=taxranks[[i]], NArm=FALSE)}
  if(i<6){taxa_names(temp.physeq) <- as.data.frame(cbind(tax_table(temp.physeq)))[,i+1]}
  temp.physeq <- subset_taxa(temp.physeq, taxa_names(temp.physeq) %in% colnames(feces.strict.sum.rel[[i]]))
  temp.physeq <- microbiome::transform(temp.physeq, "clr")
  temp.otu.tab <- as.data.frame(otu_table(temp.physeq))
  
  temp.otu.tab <- temp.otu.tab[,colnames(feces.strict.sum.rel[[i]])]
  
  # Add to list.
  feces.strict.sum.clr[[i]] <- temp.otu.tab
  rm(temp.physeq, temp.otu.tab)
}
names(feces.strict.sum.clr) <- c("Phylum", "Class", "Order", "Family", "Genus", "ASV") 

# CALCULATE NEW SUMMARIES FOR PLOTTING AND MODELING....
comp.strict.rel <- list()
comp.strict.clr <- list()

for(i in 1:length(feces.strict.sum.rel)){
  feces.strict.sum.rel[[i]]$FecesID <- rownames(feces.strict.sum.rel[[i]])
  comp.strict.rel[[i]] <- merge(cy.metadata, feces.strict.sum.rel[[i]], by="FecesID", all=TRUE)
  feces.strict.sum.rel[[i]]$FecesID <- NULL
}

for(i in 1:length(feces.strict.sum.clr)){
  feces.strict.sum.clr[[i]]$FecesID <- rownames(feces.strict.sum.clr[[i]])
  comp.strict.clr[[i]] <- merge(cy.metadata, feces.strict.sum.clr[[i]], by="FecesID", all=TRUE)
  feces.strict.sum.clr[[i]]$FecesID <- NULL
}

names(comp.strict.rel) <- c(taxranks, "ASV")
names(comp.strict.clr) <- c(taxranks, "ASV")

# (2) Calculate univariate correlations ####
corr.results <- list()
for(i in 1:length(comp.strict.clr)){
  temp <- comp.strict.clr[[i]]

  # Reduce the data frame to taxa preset in at least 20% of samples.
  # (They were already filtered at a 0.01% abundance threshold).
  temp <- temp[,sapply(temp, is.numeric)]
  temp <- as.data.frame(t(temp))
  temp1 <- subset(temp, rownames(temp) %in% rownames(subset(prev.data[[i]], Prevalence > 17)))
  
  # Subset to the variables of interest.
  temp2 <- subset(temp, rownames(temp) %in% c("Mass","Snout.to.Base","Girth",
                                            "Cementum.Age","Spleen.to.Body.Ratio","KFI", "index1.PC1",
                                            "d13C", "d15N", "Vol.Prey", "Vol.Anthro", "Vol.Veg"))

  # Create a data frame for calculating correlations.
  temp <- rbind(temp1, temp2)
  rm(temp1, temp2)
  temp <- as.data.frame(t(temp))

  # Calculate correlations between each taxa and metadata variable.
  temp.spearman <- as.matrix(temp)
  temp.spearman <- Hmisc::rcorr(temp.spearman, type="spearman")
  for(j in 1:3){
    temp.spearman[[j]] <- subset(temp.spearman[[j]],
                                 rownames(temp.spearman[[j]]) %in% rownames(subset(prev.data[[i]], Prevalence > 17)))
    temp.spearman[[j]] <- temp.spearman[[j]][,!(colnames(temp.spearman[[j]]) %in% rownames(temp.spearman[[j]]))]
    temp.spearman[[j]] <- as.data.frame(temp.spearman[[j]])
  }
  
  # Save data
  corr.results[[i]] <- temp.spearman
  rm(temp.spearman, temp)
}
names(corr.results) <- names(comp.strict.clr)
  
#### REGRESSIONS: HEALTH/SPLEEN ~ SEX + ABUNDANCE ####
# (1) Health score regression ####
coef.health.regression.ranks <- list()
for(i in 1:length(feces.strict.sum.clr)){
  model.data <- subset(comp.strict.clr[[i]], is.na(Observed)=="FALSE")
  model.data <- model.data[order(model.data$index1.PC1),]
  model.data$index1.PC1 <- c(1:nrow(model.data))
  coefficients <- data.frame()
  
  abundance.data <- subset(comp.strict.rel[[i]], Microbiome=="Y")
  
  for(j in 54:ncol(model.data)){
    model.data <- model.data[order(model.data[,j]),]
    model.data[,j] <- c(1:nrow(model.data))
    temp.model <- lm(index1.PC1 ~ model.data[,j] + Sex,
                     data=model.data)

    coefficients[(j-53),1] <- colnames(model.data)[j] # Taxon name
    coefficients[(j-53),2] <- temp.model[["coefficients"]][[2]] - sqrt(diag(vcov(temp.model)))[[2]]
    coefficients[(j-53),3] <- temp.model[["coefficients"]][[2]]
    coefficients[(j-53),4] <- temp.model[["coefficients"]][[2]] + sqrt(diag(vcov(temp.model)))[[2]]
    coefficients[(j-53),5] <- AIC(temp.model)
    coefficients[(j-53),6] <- AIC(lm(index1.PC1 ~ Sex,
                                     data=model.data))-AIC(temp.model)
    
  }
  colnames(coefficients) <- c("Taxon", "Lower", "Coefficient", "Upper", "AIC", "dAIC")
  coefficients <- merge(coefficients, prev.data[[i]], by="Taxon", all=FALSE)
  
  coef.health.regression.ranks[[i]] <- coefficients
  rm(coefficients, model.data)
}
names(coef.health.regression.ranks) <- c(taxranks, "ASV")

# (2) Spleen mass regression ####
coef.spleen.regression.ranks <- list()
for(i in 1:length(feces.strict.sum.clr)){
  model.data <- subset(comp.strict.rel[[i]], is.na(Observed)=="FALSE")
  model.data <- model.data[order(model.data$Spleen.to.Body.Ratio),]
  model.data$Spleen.to.Body.Ratio <- c(1:nrow(model.data))
  coefficients <- data.frame()
  
  abundance.data <- subset(comp.strict.rel[[i]], Microbiome=="Y")
  
  for(j in 54:ncol(model.data)){
    model.data <- model.data[order(model.data[,j]),]
    model.data[,j] <- c(1:nrow(model.data))
    temp.model <- lm(Spleen.to.Body.Ratio ~ model.data[,j] + Sex,
                     data=model.data)
    
    coefficients[(j-53),1] <- colnames(model.data)[j] # Taxon name
    coefficients[(j-53),2] <- temp.model[["coefficients"]][[2]] - sqrt(diag(vcov(temp.model)))[[2]]
    coefficients[(j-53),3] <- temp.model[["coefficients"]][[2]]
    coefficients[(j-53),4] <- temp.model[["coefficients"]][[2]] + sqrt(diag(vcov(temp.model)))[[2]]
    coefficients[(j-53),5] <- AIC(temp.model)
    coefficients[(j-53),6] <- AIC(lm(Spleen.to.Body.Ratio ~ Sex,
                                     data=model.data))-AIC(temp.model)
  }
  colnames(coefficients) <- c("Taxon", "Lower", "Coefficient", "Upper", "AIC", "dAIC")
  coefficients <- merge(coefficients, prev.data[[i]], by="Taxon", all=FALSE)
  
  coef.spleen.regression.ranks[[i]] <- coefficients
  rm(coefficients, model.data)
}
names(coef.spleen.regression.ranks) <- c(taxranks, "ASV")

#### RCCA ANALYSIS ####
# Create a subset of important environmental metadata.
temp.env <- subset(cy.metadata, FecesID %in% sample_names(pseq.clr))
rownames(temp.env) <- temp.env$FecesID
temp.env <- temp.env[,c("d13C","d15N","Vol.Prey","Vol.Anthro","index1.PC1","Spleen.to.Body.Ratio","Cementum.Age")]

# Subset ASV abundance data.
temp.otu <- feces.strict.sum.clr[["Genus"]]

# Perform rCCA; examine scree plot and network scores.
rgcca.diet <- rcc(temp.otu, temp.env, ncomp = 3, method = 'shrinkage')
plot(rgcca.diet, scree.type = "barplot")  

rgcca.cim <- cim(rgcca.diet, margins=c(10,10), dist.method=c("correlation","correlation"))

# Identify the taxa that have the strongest influences in the network, based on the sum of the absolute values
# of their relevance scores.
rgcca.scores <- as.data.frame(rgcca.cim$mat)
rgcca.scores$Important <- rowSums(abs(rgcca.scores))

# Validate relationships for taxa that are important.
rgcca.taxon.models <- list()

model.data <- cbind(comp.strict.clr[["Genus"]], comp.strict.clr[["Family"]])
model.data <- subset(model.data, is.na(Observed)=="FALSE")

model.data[,c("Cementum.Age","index1.PC1","d13C","d15N","Vol.Prey","Vol.Anthro","Spleen.to.Body.Ratio")] <-
  scale(model.data[,c("Cementum.Age","index1.PC1","d13C","d15N","Vol.Prey","Vol.Anthro","Spleen.to.Body.Ratio")])

for(p in c("Erysipelotrichaceae", "Coriobacteriaceae", "Lachnospiraceae", # Spleen families
           "Lactobacillaceae", "Enterococcaceae", "Streptococcaceae", # Anthro families
           "Lactobacillus", "Enterococcus", "Streptococcus", # Anthro genera
           "Fusobacteriaceae",  "Succinivibrionaceae", "Sutterellaceae", # Protein families
           "Alloprevotella", "Bacteroides", "Sutterella", "Anaerobiospirillum", "Fusobacterium", # Protein genera
           "Peptostreptococcaceae", "Clostridium_XlVa")){
  # Construct an initial model and evaluate all subsets.
  model <- lm(model.data[,p] ~ 
                Sex+Cementum.Age+index1.PC1+d13C+d15N+Vol.Prey+Vol.Anthro+Spleen.to.Body.Ratio+Em.PCR.Overall,
              model.data, na.action="na.fail")
  model.dredge <- dredge(model)
  num.models <- nrow(subset(model.dredge, delta < 2))
  
  # Calculate model-averaged coefficients standardized by partial standard deviation.
  temp2 <- model.avg(model.dredge, subset=delta < 2, beta="partial.sd")
  temp2 <- cbind(confint(temp2, level=0.95, full=TRUE)[,1],
                 confint(temp2, level=0.5, full=TRUE)[,1],
                 temp2[["coefficients"]][1,],
                 confint(temp2, level=0.5, full=TRUE)[,2],
                 confint(temp2, level=0.95, full=TRUE)[,2])
  temp2 <- as.data.frame(temp2)
  colnames(temp2) <- c("psd_2.5", "psd_25", "psd_Coef", "psd_75", "psd_97.5")
  
  # Calculate summed weights across all models.
  model.dredge <- as.data.frame(model.dredge)
  weights <- data.frame(matrix(nrow=0, ncol=1))
  for(j in c(2:10)){
    temp <- model.dredge
    temp <- subset(temp, is.na(temp[,j])=="FALSE")
    weights[j-1,1] <- sum(temp[,ncol(temp)])
    rm(temp)
  }
  rownames(weights) <- colnames(model.dredge[,c(2:10)])
  colnames(weights) <- "cum_weight"
  
  temp2$merge <- rownames(temp2)
  weights$merge <- as.factor(rownames(weights))
  levels(weights$merge)[levels(weights$merge)=="Em.PCR.Overall"] <- "Em.PCR.OverallY"
  levels(weights$merge)[levels(weights$merge)=="Sex"] <- "SexM"
  
  results <- merge(temp2, weights, by="merge", all=TRUE)
  results <- subset(results, merge !="(Intercept)")
  results <- results[order(-results$cum_weight),]
  results$merge <- forcats::fct_inorder(results$merge)
  rgcca.taxon.models[[p]] <- results
  rm(results, temp2, model.dredge, weights, num.models)
}

#### STRUCTURAL EQUATION MODELS ####
# (1) Prepare data for SEM ####
feces.SEM.clr <- cbind(subset(comp.strict.clr[["Genus"]], is.na(Observed)=="FALSE"), feces.strict.sum.clr[["Family"]])

# Convert factors to dummy variables for lavaan.
levels(feces.SEM.clr$Em.PCR.Overall)[levels(feces.SEM.clr$Em.PCR.Overall)=="Y"] <- 1
levels(feces.SEM.clr$Em.PCR.Overall)[levels(feces.SEM.clr$Em.PCR.Overall)=="N"] <- 0
feces.SEM.clr$Em.PCR.Overall <- as.numeric(as.character(feces.SEM.clr$Em.PCR.Overall))

levels(feces.SEM.clr$Location)[levels(feces.SEM.clr$Location)=="urban"] <- 1
levels(feces.SEM.clr$Location)[levels(feces.SEM.clr$Location)=="rural"] <- 0
feces.SEM.clr$Location <- as.numeric(as.character(feces.SEM.clr$Location))

levels(feces.SEM.clr$Sex)[levels(feces.SEM.clr$Sex)=="M"] <- 1
levels(feces.SEM.clr$Sex)[levels(feces.SEM.clr$Sex)=="F"] <- 0
feces.SEM.clr$Sex <- as.numeric(as.character(feces.SEM.clr$Sex))

# Scale stomach content variables to meet a normal distribution.
feces.SEM.clr$Vol.Prey <- log(feces.SEM.clr$Vol.Prey+0.01)
feces.SEM.clr$Vol.Anthro <- log(feces.SEM.clr$Vol.Anthro+0.01)

# Reduce the magnitude of some variables so that variances are equal.
feces.SEM.clr$Vol.Anthro <- feces.SEM.clr$Vol.Anthro/10
feces.SEM.clr$Vol.Prey <- feces.SEM.clr$Vol.Prey/10
feces.SEM.clr$Observed_Extrap <- feces.SEM.clr$Observed_Extrap/100
feces.SEM.clr$PD_Extrap <- feces.SEM.clr$PD_Extrap/10

# (2) SEM FOR ANTHROPOGENIC FOOD / POOR HEALTH CONNECTIONS ####
# For simplicity, only the initial model and the final model are shown.
# From the null model, additional paths were added as recommended by modification indices.
# Non-significant (p>0.1) paths were removed.
# All models were evaluated based on AIC, NNFI, SRMR, RMSEA, and CFI, as shown.
# Model selection stopped when the addition or removal of a path caused model AIC to increase,
# or other fit parameters (NNFI, SRMR, RMSEA, CFI) to get worse, 
# even if the path being removed was not significant.
# Criteria for NNFI, SRMR, RMSEA, CFI, Chi-sq: https://www.cscu.cornell.edu/news/Handouts/SEM_fit.pdf 

# Create storage areas for model texts and model fits.
div.gen.models <- list()
div.gen.fits <- list()

# Construct the initial model.
div.gen.models[[1]] <- '
# regressions
index1.PC1 ~ Location + d13C + Vol.Anthro + Vol.Anthro + Em.PCR.Overall + Enterococcus + Streptococcus + PD_Extrap
Spleen.to.Body.Ratio ~ Location + d13C + Vol.Anthro + Em.PCR.Overall + PD_Extrap + Enterococcus + Streptococcus
d13C ~ Location
Vol.Anthro ~ Location
PD_Extrap ~ Vol.Anthro + d13C + Location
Enterococcus ~ Location + d13C + Vol.Anthro
Streptococcus ~ Location + d13C + Vol.Anthro
Em.PCR.Overall ~ Location + d13C + Vol.Anthro + PD_Extrap
Enterococcus ~~ Streptococcus
PD_Extrap ~~ Streptococcus
PD_Extrap ~~ Enterococcus'

# Calculate and evaluate the SEM fit.
div.gen.fits[[1]] <- sem(div.gen.models[[1]], data=feces.SEM.clr)
summary(div.gen.fits[[1]], standardized=TRUE)

# Determine if any paths should be added.
mi <- modindices(div.gen.fits[[1]])
print(mi[mi$mi > 3.0,])

# Run the next model after adding or removing paths, based on the previous model results.
div.gen.models[[2]] <- '
# regressions
index1.PC1 ~ Streptococcus + PD_Extrap
Spleen.to.Body.Ratio ~ Location + Em.PCR.Overall + Enterococcus
d13C ~ Location
PD_Extrap ~ d13C + Location + Vol.Anthro
Enterococcus ~ Vol.Anthro
Streptococcus ~ Vol.Anthro
Em.PCR.Overall ~ PD_Extrap
Enterococcus ~~ Streptococcus
PD_Extrap ~~ Streptococcus
PD_Extrap ~~ Enterococcus'

div.gen.fits[[2]] <- sem(div.gen.models[[2]], data=feces.SEM.clr)
summary(div.gen.fits[[2]], standardized=TRUE)
mi <- modindices(div.gen.fits[[2]])
print(mi[mi$mi > 3.0,])

# Compare the two models.
anova(div.gen.fits[[2]], div.gen.fits[[1]])

# Evaluate model AIC.
aictab.lavaan(div.gen.fits,
              paste0("Model", seq(length(div.gen.fits))))

# Evaluate additional fit statistics.
div.gen.fits.stats <- data.frame(matrix(nrow=42))
for(i in 1:length(div.gen.fits)){
  vector <- fitmeasures(div.gen.fits[[i]])
  div.gen.fits.stats[,i] <- vector
  rownames(div.gen.fits.stats) <- names(vector)
}
div.gen.fits.stats <- subset(div.gen.fits.stats, rownames(div.gen.fits.stats) %in% c("npar", "chisq", "df", "pvalue",
                                                                                     "cfi", "rfi", "aic", "rmsea", "srmr", "gfi", "nnfi"))
colnames(div.gen.fits.stats) <- paste0("Model", seq(length(div.gen.fits)))
for(i in 1:ncol(div.gen.fits.stats)){div.gen.fits.stats[,i] <- as.numeric(div.gen.fits.stats[,i])}
div.gen.fits.stats <- as.data.frame(t(div.gen.fits.stats))

# Determine the R-sq for each variable.
lavInspect(div.gen.fits[[2]], "rsquare")

# Export to CytoScape for figure design.
t <- semPlot::semPaths(div.gen.fits[[16]], "std", title=FALSE, layout="spring", 
                       residuals=FALSE, intercepts=FALSE, thresholds=FALSE, nCharNodes=0, sizeMan=10)
t <- igraph::as.igraph(t)
igraph::write.graph(t, file = "~/anthro.SEM.gen.gml", format = "gml")

# (3) SEM FOR PROTEIN / GOOD HEALTH CONNECTIONS ####
# Create storage areas for model texts and model fits.
protein.gen.models <- list()
protein.gen.fits <- list()

# Propose the initial model.
protein.gen.models[[1]] <- '
# regressions
index1.PC1 ~ Location + d15N + Vol.Prey + Vol.Prey + Em.PCR.Overall + Anaerobiospirillum + Fusobacterium + Sutterella
Spleen.to.Body.Ratio ~ Location + d15N + Em.PCR.Overall + Anaerobiospirillum + Fusobacterium + Sutterella
d15N ~ Location
Vol.Prey ~ Location
Anaerobiospirillum ~ Location + d15N + Vol.Prey
Fusobacterium ~ Location + d15N + Vol.Prey
Sutterella ~ Location + d15N + Vol.Prey
Em.PCR.Overall ~ Location + d15N + Vol.Prey
Anaerobiospirillum ~~ Fusobacterium
Anaerobiospirillum ~~ Sutterella
Sutterella ~~ Fusobacterium'

# Calculate and evaluate the SEM fit.
protein.gen.fits[[1]] <- sem(protein.gen.models[[1]], data=feces.SEM.clr)
summary(protein.gen.fits[[1]], standardized=TRUE)

# Determine if any paths should be added.
mi <- modindices(protein.gen.fits[[1]])
print(mi[mi$mi > 3.0,])

# Run the next model.
# (Accelerated here to show the final model).
protein.gen.models[[2]] <- '
# regressions
index1.PC1 ~ Location + d15N
Spleen.to.Body.Ratio ~ Location + Em.PCR.Overall + Anaerobiospirillum + Sutterella
d15N ~ Location
Vol.Prey ~ Spleen.to.Body.Ratio
Anaerobiospirillum ~ d15N
Fusobacterium ~ Location + d15N
Sutterella ~ d15N
Em.PCR.Overall ~ Location + Vol.Prey
Anaerobiospirillum ~~ Fusobacterium
Anaerobiospirillum ~~ Sutterella
Sutterella ~~ Fusobacterium'

protein.gen.fits[[2]] <- sem(protein.gen.models[[2]], data=feces.SEM.clr)
summary(protein.gen.fits[[2]], standardized=TRUE)
mi <- modindices(protein.gen.fits[[2]])
print(mi[mi$mi > 3.0,])

# Compare to the previous model.
anova(protein.gen.fits[[6]], protein.gen.fits[[5]])

# Evaluate model AICs.
aictab.lavaan(protein.gen.fits,
              paste0("Model", seq(length(protein.gen.fits))))

# Evaluate additional fit statistics.
protein.gen.fits.stats <- data.frame(matrix(nrow=42))
for(i in 1:length(protein.gen.fits)){
  vector <- fitmeasures(protein.gen.fits[[i]])
  protein.gen.fits.stats[,i] <- vector
  rownames(protein.gen.fits.stats) <- names(vector)
}
protein.gen.fits.stats <- subset(protein.gen.fits.stats, rownames(protein.gen.fits.stats) %in% c("npar", "chisq", "df", "pvalue",
                                                                                                 "cfi", "nnfi", "aic", "rmsea", "srmr", "gfi"))
colnames(protein.gen.fits.stats) <- paste0("Model", seq(length(protein.gen.fits)))
for(i in 1:ncol(protein.gen.fits.stats)){protein.gen.fits.stats[,i] <- as.numeric(protein.gen.fits.stats[,i])}
protein.gen.fits.stats <- as.data.frame(t(protein.gen.fits.stats))

# Determine the R-sq value for each variable.
lavInspect(protein.gen.fits[[6]], "rsquare")

# Export to CytoScape for figure design.
t <- semPlot::semPaths(protein.gen.fits[[6]], "std", title=FALSE, layout="spring", 
                       residuals=FALSE, intercepts=FALSE, thresholds=FALSE, nCharNodes=0, sizeMan=10)
t <- igraph::as.igraph(t)
igraph::write.graph(t, file = "~/protein.gen.SEM.gml", format = "gml")

# (4) SEM FOR SPLEEN MASS CONNECTIONS ####
spleen.models <- list()
spleen.fits <- list()

spleen.models[[1]] <- '
# regressions
index1.PC1 ~ Location + d13C + Vol.Anthro + Em.PCR.Overall + Erysipelotrichaceae + Coriobacteriaceae + Lachnospiraceae
Spleen.to.Body.Ratio ~ Location + d13C + Vol.Anthro + Em.PCR.Overall + Erysipelotrichaceae + Coriobacteriaceae + Lachnospiraceae
d13C ~ Location
Vol.Anthro ~ Location
Erysipelotrichaceae ~ Location + d13C + Vol.Anthro
Coriobacteriaceae ~ Location + d13C + Vol.Anthro
Lachnospiraceae ~ Location + d13C + Vol.Anthro
Em.PCR.Overall ~ Location + d13C + Vol.Anthro
Erysipelotrichaceae ~~ Coriobacteriaceae
Coriobacteriaceae ~~ Lachnospiraceae
Erysipelotrichaceae ~~ Lachnospiraceae'

spleen.fits[[1]] <- sem(spleen.models[[1]], data=feces.SEM.clr)
summary(spleen.fits[[1]], standardized=TRUE)
mi <- modindices(spleen.fits[[1]])
print(mi[mi$mi > 3.0,])

spleen.models[[9]] <- '
# regressions
index1.PC1 ~ Location + Lachnospiraceae
Spleen.to.Body.Ratio ~ Location + Em.PCR.Overall + Erysipelotrichaceae
d13C ~ Location
Erysipelotrichaceae ~ Location + Vol.Anthro
Coriobacteriaceae ~ Location + Vol.Anthro
Lachnospiraceae ~ Location
Em.PCR.Overall ~ Location + d13C
Erysipelotrichaceae ~~ Coriobacteriaceae
Coriobacteriaceae ~~ Lachnospiraceae
Erysipelotrichaceae ~~ Lachnospiraceae'

spleen.fits[[9]] <- sem(spleen.models[[9]], data=feces.SEM.clr)
summary(spleen.fits[[9]], standardized=TRUE)
mi <- modindices(spleen.fits[[9]])
print(mi[mi$mi > 3.0,])

anova(spleen.fits[[9]], spleen.fits[[1]])

aictab.lavaan(spleen.fits,
              paste0("Model", seq(length(spleen.fits))))

spleen.fits.stats <- data.frame(matrix(nrow=42))
for(i in 1:length(spleen.fits)){
  vector <- fitmeasures(spleen.fits[[i]])
  spleen.fits.stats[,i] <- vector
  rownames(spleen.fits.stats) <- names(vector)
}
spleen.fits.stats <- subset(spleen.fits.stats, rownames(spleen.fits.stats) %in% c("npar", "chisq", "df", "pvalue",
                                                                                  "cfi", "nnfi", "aic", "rmsea", "srmr", "gfi"))
colnames(spleen.fits.stats) <- paste0("Model", seq(length(spleen.fits)))
for(i in 1:ncol(spleen.fits.stats)){spleen.fits.stats[,i] <- as.numeric(spleen.fits.stats[,i])}
spleen.fits.stats <- as.data.frame(t(spleen.fits.stats))

lavInspect(spleen.fits[[9]], "rsquare")

t <- semPlot::semPaths(spleen.fits[[9]], "std", title=FALSE, layout="spring", 
                       residuals=FALSE, intercepts=FALSE, thresholds=FALSE, nCharNodes=0, sizeMan=10)
t <- igraph::as.igraph(t)
igraph::write.graph(t, file = "~/spleen.SEM.gml", format = "gml")

# FUNCTIONS ####
AICc.lavaan<-function(object, second.ord=TRUE, c.hat = 1, return.K = FALSE){
  object <- as.list(fitMeasures(object))
  npar <- object$baseline.df - object$df + 2 # number of parameters, including coefficients of control variables (if there are any), intercept and residual variance (thus the +2)
  if(return.K == TRUE) return(npar)
  if(second.ord == FALSE && c.hat>1) return(-2*object$logl/c.hat+2*npar)
  if(second.ord == FALSE) return(-2*object$logl + 2*npar)
  if(c.hat>1) return( -2*object$logl/c.hat+2*npar + 2*( npar*(object$npar+1))/(object$ntotal-npar-1))
  -2*object$logl + 2*npar + 2*( npar*(npar+1))/(object$ntotal-npar-1)
}
