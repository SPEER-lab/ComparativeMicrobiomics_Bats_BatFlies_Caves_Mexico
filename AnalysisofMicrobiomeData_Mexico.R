#Analysis of Microbiome Data
#This script takes filtered 16S data, imports it into phyloseq and creates balances between closely related taxa within a sample
#Not only does this allow for normalization of data that does not rely on a minimum sequencing cutoff
#It also takes into account the compositional nature of microbiome data
#CITATION:Silverman et al. 2017
#CITATION: Gloor et al. 2017

##############################Load the necessary libraries#######################################################
library(phyloseq)
library(philr)
library(ape)
library(tidyverse)
library(ggtree)
library(vegan)
library(Matrix)
library(glmnet)
library(decontam)
library(cluster)    
library(factoextra)
library(SpiecEasi)
library(huge)
library(microbiome)
library(viridis)
library(igraph)
library(scales)
library(reshape2)
library(lme4)
library('pez')
library('caret')
library('ranger')
library('e1071')

#If revisiting this script, search for "START HERE" to load filtered data into the phyloseq object BrazilMicrob


#############################Format QIIME2-exported Filtered Data################################################
#If data was exported from QIIME2, the biom file does not have the necessary taxonomy information associated with it
#I edited the instructions in the pdf posted to this GitHub query (https://github.com/joey711/phyloseq/issues/821)
#to format the taxonomy file correctly. The feature table file is only needed to eliminate some rows from the taxonomy file

#Open dowloaded feature-table in Text Wrangelr and remove any rows before header row and remove # sign before OTUID
#Read in the .txt version of the feature table, which should now have a column header "OTUID", not "#OTUID"
features <- read.table(file="QIIME2/exported-feature-table/feature-table.txt", header=TRUE, row.names=1)
head(features)

# features_loose <- read.table(file="QIIME2/exported-feature-table-loose/feature-table.txt", header=TRUE, row.names=1)
# head(features_loose)

#Read in the .tsv version of the taxonomy table, which should also have a column header "OTUID", not "Feature ID"
tax <- read.table(file="QIIME2/exported-feature-table/taxonomy.tsv", sep='\t', header=TRUE, row.names=1)
head(tax)

#eliminate OTUs present in tax, but not features
#The remaining table should have the same number of rows as in the features data frame
tax_filtered <- tax[row.names(tax) %in% row.names(features),]
head(tax_filtered)

# tax_filtered_loose <- tax[row.names(tax) %in% row.names(features_loose),]
# head(tax_filtered_loose)

#Separate the "Taxon" column in the tax_filtered data frame by semicolon so that each step of the taxonomy (e.g., kingdom, phylum, class, etc.) is its own column
tax_filtered <- separate(tax_filtered, Taxon, c("Kingdom","Phylum","Class","Order", "Family", "Genus","Species"), sep= ";", remove=TRUE)
head(tax_filtered)

# tax_filtered_loose <- separate(tax_filtered_loose, Taxon, c("Kingdom","Phylum","Class","Order", "Family", "Genus","Species"), sep= ";", remove=TRUE)
# head(tax_filtered_loose)

#Remove the prefixes from the taxonomic information
tax_filtered$Kingdom <- sub("D_0__", "", tax_filtered$Kingdom)
tax_filtered$Phylum <- sub("D_1__", "", tax_filtered$Phylum)
tax_filtered$Class <- sub("D_2__", "", tax_filtered$Class)
tax_filtered$Order <- sub("D_3__", "", tax_filtered$Order)
tax_filtered$Family <- sub("D_4__", "", tax_filtered$Family)
tax_filtered$Genus <- sub("D_5__", "", tax_filtered$Genus)
tax_filtered$Species <- sub("D_6__", "", tax_filtered$Species)
head(tax_filtered)

# tax_filtered_loose$Kingdom <- sub("D_0__", "", tax_filtered_loose$Kingdom)
# tax_filtered_loose$Phylum <- sub("D_1__", "", tax_filtered_loose$Phylum)
# tax_filtered_loose$Class <- sub("D_2__", "", tax_filtered_loose$Class)
# tax_filtered_loose$Order <- sub("D_3__", "", tax_filtered_loose$Order)
# tax_filtered_loose$Family <- sub("D_4__", "", tax_filtered_loose$Family)
# tax_filtered_loose$Genus <- sub("D_5__", "", tax_filtered_loose$Genus)
# tax_filtered_loose$Species <- sub("D_6__", "", tax_filtered_loose$Species)
# head(tax_filtered_loose)

#remove "Confidence" column
tax_filtered$Confidence<-NULL
head(tax_filtered)

# tax_filtered_loose$Confidence<-NULL
# head(tax_filtered_loose)

#write one outfile containing the OTUID and taxonomic info
write.csv(tax_filtered, file="phyloseq/taxonomy_phyloseq.csv", row.names=TRUE)
# write.csv(tax_filtered_loose, file="phyloseq/taxonomy_loose_phyloseq.csv", row.names=TRUE)

#Re-format the sample name in the feature data frame to remove the X at the beginning of each sample name 
#and replace the . with an _
colnames(features)<- sub("X", "Sample", colnames(features))
colnames(features)<- sub("^", "Sample", colnames(features))
colnames(features)<- sub("SampleSample", "Sample", colnames(features))
colnames(features)<- sub(".", "_", colnames(features), fixed=TRUE)
head(features)

# colnames(features_loose)<- sub("X", "Sample", colnames(features_loose))
# colnames(features_loose)<- sub("^", "Sample", colnames(features_loose))
# colnames(features_loose)<-sub("SampleSample", "Sample", colnames(features_loose))
# colnames(features_loose)<- sub(".", "_", colnames(features_loose), fixed=TRUE)
# colnames(features_loose)

#write out the reformatted feature table
write.csv(features, file="phyloseq/features_phyloseq.csv", row.names=TRUE)
# write.csv(features_loose, file="phyloseq/features_loose_phyloseq.csv", row.names=TRUE)


#####################Format the metadata#########################################################################
#Read in the metadata for the microbiome study
metadata <- read.table("Metadata_Mexico_filtered.tsv", sep='\t', row.names = 1, header =TRUE)
head(metadata)

# metadata_loose <- read.table("Metadata_Mexico.tsv", sep='\t', row.names = 1, header =TRUE)
# head(metadata_loose)

#Format the Metadata sample IDs to match the sample IDs in the features matrix
row.names(metadata) <- sub("-", "_", row.names(metadata))
row.names(metadata) <- sub("^", "Sample", row.names(metadata))
row.names(metadata)

# row.names(metadata_loose) <- sub("-", "_", row.names(metadata_loose))
# row.names(metadata_loose) <- sub("^", "Sample", row.names(metadata_loose))
# head(metadata_loose)

#Make sure sample names in features and metadata match
length(colnames(features))
length(row.names(metadata))
ToRemove<-setdiff(row.names(metadata), colnames(features))
setdiff(colnames(features), row.names(metadata))

# length(colnames(features_loose))
# length(row.names(metadata_loose))
# ToRemove<-setdiff(row.names(metadata_loose), colnames(features_loose))
# setdiff(colnames(features_loose), row.names(metadata_loose))

metadata<- metadata[!(row.names(metadata) %in% ToRemove),]
setdiff(row.names(metadata), colnames(features))

# metadata_loose<- metadata_loose[!(row.names(metadata_loose) %in% ToRemove),]
# setdiff(row.names(metadata_loose), colnames(features_loose))

#Add Library Concetration Information and Bat ID as separate columns to filtered metadata
LibraryConc<- read.csv("LibraryConcentration_Decontam.csv", header=TRUE, row.names=1)

#Check that Samples match between metadata and LibraryConc
setdiff(row.names(metadata), row.names(LibraryConc))
setdiff(row.names(LibraryConc), row.names(metadata))

#Merge LibConc with metadata
metadata <- merge(metadata, LibraryConc, by="row.names")
row.names(metadata)<-metadata$Row.names
metadata$Row.names<-NULL
head(metadata)

#Add Baja versus Mainland as a broad locality identifier
metadata$BroadLocality <- as.character(metadata$CollectionLocality)
metadata$BroadLocality[metadata$BroadLocality == "Carmen, Loreto, Baja California Sur" | 
                         metadata$BroadLocality == "Chivato, Sierra Cacachilas, Baja California Sur" |
                         metadata$BroadLocality == "La Gitana, Sierra Cacachilas, Baja California Sur"]  <- "Baja"
metadata$BroadLocality[metadata$BroadLocality == "Cueva de la Fabrica, Coquimatlan, Colima, MX" | 
                         metadata$BroadLocality == "Isla Panchito, near Chamela, Jalisco, MX"]  <- "Mainland"

#Add a column that says the sample type "Bat, Cave, Trichobius, Nycterophilia"
metadata$SampleSource[metadata$BatFlyGenus == "Trichobius"] <- "Trichobius"
metadata$SampleSource[metadata$BatFlyGenus == "Nycterophilia"] <- "Nycterophilia"
metadata$SampleSource[metadata$BatSwab == "Y"] <- "Bat"
metadata$SampleSource[metadata$CaveSwab == "Y"] <- "Cave"

# #Add Sample or control column to metadata_loose
# metadata_loose$Control[row.names(metadata_loose) == "SampleExt_" | row.names(metadata_loose) == "SamplePCR_" |
#                                  row.names(metadata_loose) == "SampleMicrobialplus" | row.names(metadata_loose) == "SampleMicrobDNAplus"]  <- "TRUE"
# metadata_loose$Control[is.na(metadata_loose$Control)]  <- "FALSE"

#write the metdata to an output file
write.csv(metadata, "phyloseq/metadata_phyloseq.csv", row.names=TRUE)
# write.csv(metadata_loose, "phyloseq/metadata_loose_phyloseq.csv", row.names=TRUE)


###########################################Import QIIME2-filtered Data into Phyloseq#############################
#Now we can load the data into phyloseq
#If you are starting here and don't have the tax_filtered and features data frames loaded, do that now
features <- read.csv("phyloseq/features_phyloseq.csv", sep=',', row.names=1, header=TRUE)
features <- as.matrix(features)
head(features)

tax_filtered <- read.csv("phyloseq/taxonomy_phyloseq.csv", sep=',', row.names=1, header=TRUE)
tax_filtered <- as.matrix(tax_filtered)
head(tax_filtered)

metadata<- read.csv("phyloseq/metadata_phyloseq.csv", sep=',', row.names=1, header=TRUE)
head(metadata)

# features_loose <- read.csv("phyloseq/features_loose_phyloseq.csv", sep=',', row.names=1, header=TRUE)
# features_loose <- as.matrix(features_loose)
# head(features_loose)
# 
# tax_filtered_loose <- read.csv("phyloseq/taxonomy_loose_phyloseq.csv", sep=',', row.names=1, header=TRUE)
# tax_filtered_loose <- as.matrix(tax_filtered_loose)
# head(tax_filtered_loose)
# 
# metadata_loose<- read.csv("phyloseq/metadata_loose_phyloseq.csv", sep=',', row.names=1, header=TRUE)
# head(metadata_loose)

#Read in rooted tree as phyloseq object
phy_tree <- read_tree("QIIME2/exported-tree-rooted/tree.nwk")
# phy_tree_loose <- read_tree("QIIME2/exported-tree-rooted-loose/tree.nwk")

#Import as phyloseq objects
OTU <- otu_table(features, taxa_are_rows=TRUE)
TAX <- tax_table(tax_filtered)
META <- sample_data(metadata)

# OTU_loose <- otu_table(features_loose, taxa_are_rows=TRUE)
# TAX_loose <- tax_table(tax_filtered_loose)
# META_loose <- sample_data(metadata_loose)

#Are OTU names consistent across objects?
taxa_names(OTU)
taxa_names(TAX)
taxa_names(phy_tree)

# taxa_names(OTU_loose)
# taxa_names(TAX_loose)
# taxa_names(phy_tree_loose)

#Do files have the same sample names?
sample_names(OTU)
sample_names(META)
setdiff(sample_names(META), sample_names(OTU))

# sample_names(OTU_loose)
# sample_names(META_loose)
# setdiff(sample_names(META_loose), sample_names(OTU_loose))

#Merge into one phyloseq object
MexicoMicrob.unfiltered <- phyloseq(OTU, TAX, META, phy_tree)
MexicoMicrob.unfiltered
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7341 taxa and 177 samples ]
# sample_data() Sample Data:       [ 177 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 7341 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7341 tips and 7339 internal nodes ]

# MexicoMicrob.unfiltered.loose <- phyloseq(OTU_loose, TAX_loose, META_loose, phy_tree_loose)
# MexicoMicrob.unfiltered.loose
# 
# # phyloseq-class experiment-level object
# # otu_table()   OTU Table:         [ 9758 taxa and 181 samples ]
# # sample_data() Sample Data:       [ 181 samples by 32 sample variables ]
# # tax_table()   Taxonomy Table:    [ 9758 taxa by 7 taxonomic ranks ]
# # phy_tree()    Phylogenetic Tree: [ 9758 tips and 9757 internal nodes ]


######################################Filter phyloseq data using decontam########################################
#https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

#Filter using abundance with the strictly filtered data (BA24 is assumed to be a negative control)
#Estimate which ASVs might be contaminants based on their occurence across samples with different concentrations
#Contaminant frequency is expected to be inversely proportional to input DNA concentration.
contamdf.freq <- isContaminant(MexicoMicrob.unfiltered, method="frequency", conc="Concentration")
head(contamdf.freq)

#Summarize how many contaminant features there are in the data set and their rank in abundance in the overall dataset
table(contamdf.freq$contaminant)
# FALSE  TRUE 
# 7330    11 
(which(contamdf.freq$contaminant))
# [1] 2514 2662 2699 3816 4052 5262 6034 6078 6842 6844 6846

set.seed(100)
plot_frequency(MexicoMicrob.unfiltered, taxa_names(MexicoMicrob.unfiltered)[sample(which(contamdf.freq$contaminant))], conc="Concentration") +
  xlab("DNA Concentration (PicoGreen fluorescent intensity)")

MexicoMicrob.filtered <- prune_taxa(!contamdf.freq$contaminant, MexicoMicrob.unfiltered)
MexicoMicrob.filtered
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7330 taxa and 177 samples ]
# sample_data() Sample Data:       [ 177 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 7330 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7330 tips and 7328 internal nodes ]

#We must eliminate 7 additional bat flies that were likely contaminated, evidenced by the top 4 symbionts not
#comprising at least 40% of their microbiomes. These samples are:

ToRemove<- c("Sample249595_1","Sample249715_2","Sample249716","Sample249725_3","Sample249726_2","Sample249743_2","Sample249753_2","Sample249785_1","Sample249746_2","Sample249733_2","Sample256054_1")

#Remove the suspicious bat flies from the phyloseq object to create a filtered object
MexicoMicrob<- prune_samples(!(sample_names(MexicoMicrob.filtered) %in% ToRemove), MexicoMicrob.filtered)
MexicoMicrob
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7330 taxa and 166 samples ]
# sample_data() Sample Data:       [ 166 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 7330 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7330 tips and 7328 internal nodes ]

# #Filter using the negative and positive controls (BA24 is not assumed to be a control)
# #add "SampleorControl" column to metadata
# sample_data(MexicoMicrob.unfiltered.loose)$Control<- as.logical(sample_data(MexicoMicrob.unfiltered.loose)$Control)
# contamdf.prev <- isContaminant(MexicoMicrob.unfiltered.loose, method="prevalence", neg="Control")
# table(contamdf.prev$contaminant)
# # FALSE  TRUE 
# # 9719    39 
# 
# which(contamdf.prev$contaminant)
# # [1]  209  342  455  471  489  501  503  517  520  524  542  642  648 1002 1005 1020 1112 1149 1541 1560 1568 1573
# # [23] 1816 1883 2628 2649 2706 2743 3440 3976 4580 4633 6482 7260 7753 8033 8341 8439 8611
# 
# #Using a lower threshold to determine a contaminant
# contamdf.prev05 <- isContaminant(MexicoMicrob.unfiltered.loose, method="prevalence", neg="Control", threshold=0.5)
# table(contamdf.prev05$contaminant)
# # FALSE  TRUE 
# # 9688    70

# #Plot the prevalence of taxa in True Samples by the prevalence in taxa in negative controls using both thresholds
# MexicoMicrob.unfiltered.loose.pa <- transform_sample_counts(MexicoMicrob.unfiltered.loose, function(abund) 1*(abund>0))
# MexicoMicrob.unfiltered.loose.pa.neg <- prune_samples(sample_data(MexicoMicrob.unfiltered.loose.pa)$Control == TRUE, MexicoMicrob.unfiltered.loose.pa)
# MexicoMicrob.unfiltered.loose.pa.pos <- prune_samples(sample_data(MexicoMicrob.unfiltered.loose.pa)$Control == FALSE, MexicoMicrob.unfiltered.loose.pa)
# 
# # Make data.frame of prevalence in positive and negative samples
# df.pa <- data.frame(MexicoMicrob.unfiltered.loose.pa.pos=taxa_sums(MexicoMicrob.unfiltered.loose.pa.pos), MexicoMicrob.unfiltered.loose.pa.neg=taxa_sums(MexicoMicrob.unfiltered.loose.pa.neg),
#                     contaminant=contamdf.prev$contaminant)
# ggplot(data=df.pa, aes(x=MexicoMicrob.unfiltered.loose.pa.neg, y=MexicoMicrob.unfiltered.loose.pa.pos, color=contaminant)) + geom_point() +
#   xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
# 
# df.pa.2 <- data.frame(MexicoMicrob.unfiltered.loose.pa.pos=taxa_sums(MexicoMicrob.unfiltered.loose.pa.pos), MexicoMicrob.unfiltered.loose.pa.neg=taxa_sums(MexicoMicrob.unfiltered.loose.pa.neg),
#                     contaminant=contamdf.prev05$contaminant)
# ggplot(data=df.pa.2, aes(x=MexicoMicrob.unfiltered.loose.pa.neg, y=MexicoMicrob.unfiltered.loose.pa.pos, color=contaminant)) + geom_point() +
#   xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")


##########################Filter data to remove low abundance and low variation OTUs############################################
#We differ in the filtering recommended in the phyloseq tutorial
#We eliminate any ASV from a sample that was not detected at least 5 or 25 times 
MexicoMicrob
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7330 taxa and 166 samples ]
# sample_data() Sample Data:       [ 166 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 7330 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7330 tips and 7328 internal nodes ]

#Remove samples that fall below a minimum sequencing depth of 1000
MM<-prune_samples(sample_sums(MexicoMicrob)>=1000, MexicoMicrob)
setdiff(sample_names(MexicoMicrob), sample_names(MM))
MM
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7330 taxa and 163 samples ]
# sample_data() Sample Data:       [ 163 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 7330 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7330 tips and 7328 internal nodes ]

#Take the otu_table from the phyloseq object filtered to min read depth
MM.otu<- as.matrix(otu_table(MM))
MM.otu25 <- as.matrix(otu_table(MM))

#If an ASV has fewer than 25 reads or 5 reads in a sample, replace the read copy number with 0, making it so that ASV is no longer detected in that sample
MM.otu[MM.otu < 5] <- 0
MM.otu25[MM.otu25 < 25] <- 0

#Recreate the phyloseq object using the new filtered otu_table
MM<-phyloseq(otu_table(MM.otu), TAX, META, phy_tree)
MM.25 <- phyloseq(otu_table(MM.otu25), TAX, META, phy_tree)

#Remove ASVs in the new phyloseq object that are not present in any samples
MM<-prune_taxa(taxa_sums(MM)>0, MM)
MM.25<-prune_taxa(taxa_sums(MM.25)>0, MM.25)

MM
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6027 taxa and 163 samples ]
# sample_data() Sample Data:       [ 163 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 6027 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6027 tips and 6025 internal nodes ]

MM.25
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2706 taxa and 163 samples ]
# sample_data() Sample Data:       [ 163 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 2706 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2706 tips and 2704 internal nodes ]

#Split the phyloseq objects by Sample Type (bat swab, cave swab, bat fly, and by bat fly species)
MM.Bat <- subset_samples(MM, BatSwab=="Y")
MM.Bat <-prune_taxa(taxa_sums(MM.Bat)>0, MM.Bat)
MM.Bat
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4050 taxa and 42 samples ]
# sample_data() Sample Data:       [ 42 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 4050 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 4050 tips and 4048 internal nodes ]

MM.BatFly <- subset_samples(MM, BatFly=="Y")
MM.BatFly <-prune_taxa(taxa_sums(MM.BatFly)>0, MM.BatFly)
MM.BatFly
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1746 taxa and 116 samples ]
# sample_data() Sample Data:       [ 116 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 1746 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1746 tips and 1745 internal nodes ]

MM.Trichobius <- subset_samples(MM, BatFlyGenus=="Trichobius")
MM.Trichobius <-prune_taxa(taxa_sums(MM.Trichobius)>0, MM.Trichobius)
MM.Trichobius
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1042 taxa and 57 samples ]
# sample_data() Sample Data:       [ 57 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 1042 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1042 tips and 1041 internal nodes ]

MM.Nycterophilia <- subset_samples(MM, BatFlyGenus=="Nycterophilia")
MM.Nycterophilia <-prune_taxa(taxa_sums(MM.Nycterophilia)>0, MM.Nycterophilia)
MM.Nycterophilia
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 756 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 756 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 756 tips and 755 internal nodes ]

MM.Cave <- subset_samples (MM, CaveSwab =="Y")
MM.Cave <-prune_taxa(taxa_sums(MM.Cave)>0, MM.Cave)
MM.Cave
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 609 taxa and 5 samples ]
# sample_data() Sample Data:       [ 5 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 609 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 609 tips and 608 internal nodes ]

#read out the filtered otu_table, tax_table, and phylogeny
write.csv(otu_table(MM), "phyloseq/features_phyloseq_filtered.csv", row.names=TRUE)
write.csv(tax_table(MM), "phyloseq/taxonomy_phyloseq_filtered.csv", row.names=TRUE)
write.tree(phy_tree(MM), "phyloseq/tree_phyloseq_filtered.nwk")
write.csv(sample_data(MM), "phyloseq/metadata_phyloseq_filtered.csv", row.names=TRUE)

write.csv(otu_table(MM.25), "phyloseq/features_phyloseq_filtered25.csv", row.names=TRUE)
write.csv(tax_table(MM.25), "phyloseq/taxonomy_phyloseq_filtered25.csv", row.names=TRUE)
write.tree(phy_tree(MM.25), "phyloseq/tree_phyloseq_filtered25.nwk")
write.csv(sample_data(MM.25), "phyloseq/metadata_phyloseq_filtered25.csv", row.names=TRUE)

write.csv(otu_table(MM.BatFly), "phyloseq/features_phyloseq_BatFly.csv", row.names=TRUE)
write.csv(tax_table(MM.BatFly), "phyloseq/taxonomy_phyloseq_BatFly.csv", row.names=TRUE)
write.tree(phy_tree(MM.BatFly), "phyloseq/tree_phyloseq_BatFly.nwk")
write.csv(sample_data(MM.BatFly), "phyloseq/metadata_phyloseq_BatFly.csv", row.names=TRUE)

write.csv(otu_table(MM.Trichobius), "phyloseq/features_phyloseq_Trichobius.csv", row.names=TRUE)
write.csv(tax_table(MM.Trichobius), "phyloseq/taxonomy_phyloseq_Trichobius.csv", row.names=TRUE)
write.tree(phy_tree(MM.Trichobius), "phyloseq/tree_phyloseq_Trichobius.nwk")
write.csv(sample_data(MM.Trichobius), "phyloseq/metadata_phyloseq_Trichobius.csv", row.names=TRUE)

write.csv(otu_table(MM.Nycterophilia), "phyloseq/features_phyloseq_Nycterophilia.csv", row.names=TRUE)
write.csv(tax_table(MM.Nycterophilia), "phyloseq/taxonomy_phyloseq_Nycterophilia.csv", row.names=TRUE)
write.tree(phy_tree(MM.Nycterophilia), "phyloseq/tree_phyloseq_Nycterophilia.nwk")
write.csv(sample_data(MM.Nycterophilia), "phyloseq/metadata_phyloseq_Nycterophilia.csv", row.names=TRUE)

write.csv(otu_table(MM.Bat), "phyloseq/features_phyloseq_Bat.csv", row.names=TRUE)
write.csv(tax_table(MM.Bat), "phyloseq/taxonomy_phyloseq_Bat.csv", row.names=TRUE)
write.tree(phy_tree(MM.Bat), "phyloseq/tree_phyloseq_Bat.nwk")
write.csv(sample_data(MM.Bat), "phyloseq/metadata_phyloseq_Bat.csv", row.names=TRUE)

write.csv(otu_table(MM.Cave), "phyloseq/features_phyloseq_Cave.csv", row.names=TRUE)
write.csv(tax_table(MM.Cave), "phyloseq/taxonomy_phyloseq_Cave.csv", row.names=TRUE)
write.tree(phy_tree(MM.Cave), "phyloseq/tree_phyloseq_Cave.nwk")
write.csv(sample_data(MM.Cave), "phyloseq/metadata_phyloseq_Cave.csv", row.names=TRUE)

#Re-create phyloseq object from filtered data
##############START HERE IF REVISITING THIS SCRIPT!###################################################################
#For decontaminated data with ASVs removed if they have fewer than 5 reads per sample
features <- read.csv("phyloseq/features_phyloseq_filtered.csv", sep=',', row.names=1, header=TRUE)
features <- as.matrix(features)
head(features)

tax <- read.csv("phyloseq/taxonomy_phyloseq_filtered.csv", sep=',', row.names=1, header=TRUE)
tax <- as.matrix(tax)
head(tax)

metadata<- read.csv("phyloseq/metadata_phyloseq_filtered.csv", sep=',', row.names=1, header=TRUE)
head(metadata)

#Read in rooted tree as phyloseq object
phy_tree <- read_tree("phyloseq/tree_phyloseq_filtered.nwk")

#Import as phyloseq objects
OTU <- otu_table(features, taxa_are_rows=TRUE)
TAX <- tax_table(tax)
META <- sample_data(metadata)

#Are OTU names consistent across objects?
taxa_names(OTU)
taxa_names(TAX)
taxa_names(phy_tree)

#Do files have the same sample names?
sample_names(OTU)
sample_names(META)

#Merge into one phyloseq object
MM <- phyloseq(OTU, TAX, META, phy_tree)
MM
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 6027 taxa and 163 samples ]
# sample_data() Sample Data:       [ 163 samples by 34 sample variables ]
# tax_table()   Taxonomy Table:    [ 6027 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 6027 tips and 6025 internal nodes ]

max(sample_sums(MM))

#For decontaminated data with ASVs removed that have fewer than 25 reads in a sample
features.25 <- read.csv("phyloseq/features_phyloseq_filtered25.csv", sep=',', row.names=1, header=TRUE)
features.25 <- as.matrix(features.25)
head(features.25)

tax.25 <- read.csv("phyloseq/taxonomy_phyloseq_filtered25.csv", sep=',', row.names=1, header=TRUE)
tax.25 <- as.matrix(tax.25)
head(tax.25)

#Read in rooted tree as phyloseq object
phy_tree.25 <- read_tree("phyloseq/tree_phyloseq_filtered25.nwk")

#Import as phyloseq objects
OTU.25 <- otu_table(features.25, taxa_are_rows=TRUE)
TAX.25 <- tax_table(tax.25)

#Are OTU names consistent across objects?
taxa_names(OTU.25)
taxa_names(TAX.25)
taxa_names(phy_tree.25)

#Do files have the same sample names?
sample_names(OTU.25)
sample_names(META)

#Merge into one phyloseq object
MM.25 <- phyloseq(OTU.25, TAX.25, META, phy_tree.25)
MM.25
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2706 taxa and 163 samples ]
# sample_data() Sample Data:       [ 163 samples by 34 sample variables ]
# tax_table()   Taxonomy Table:    [ 2706 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 2706 tips and 2704 internal nodes ]

#For Bat Fly Samples and Individual Bat Fly Species
features_BF <- read.csv("phyloseq/features_phyloseq_BatFly.csv", sep=',', row.names=1, header=TRUE)
features_BF <- as.matrix(features_BF)
head(features_BF)

tax_BF <- read.csv("phyloseq/taxonomy_phyloseq_BatFly.csv", sep=',', row.names=1, header=TRUE)
tax_BF <- as.matrix(tax_BF)
head(tax_BF)

metadata_BF<- read.csv("phyloseq/metadata_phyloseq_BatFly.csv", sep=',', row.names=1, header=TRUE)
head(metadata_BF)

phy_tree_BF <- read_tree("phyloseq/tree_phyloseq_BatFly.nwk")

OTU_BF <- otu_table(features_BF, taxa_are_rows=TRUE)
TAX_BF <- tax_table(tax_BF)
META_BF <- sample_data(metadata_BF)

taxa_names(OTU_BF)
taxa_names(TAX_BF)
taxa_names(phy_tree_BF)

sample_names(OTU_BF)
sample_names(META_BF)

MM.BatFly <- phyloseq(OTU_BF, TAX_BF, META_BF, phy_tree_BF)
MM.BatFly
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1746 taxa and 116 samples ]
# sample_data() Sample Data:       [ 116 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 1746 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1746 tips and 1745 internal nodes ]

features_Trich <- read.csv("phyloseq/features_phyloseq_Trichobius.csv", sep=',', row.names=1, header=TRUE)
features_Trich <- as.matrix(features_Trich)
head(features_Trich)

tax_Trich <- read.csv("phyloseq/taxonomy_phyloseq_Trichobius.csv", sep=',', row.names=1, header=TRUE)
tax_Trich <- as.matrix(tax_Trich)
head(tax_Trich)

metadata_Trich<- read.csv("phyloseq/metadata_phyloseq_Trichobius.csv", sep=',', row.names=1, header=TRUE)
head(metadata_Trich)

phy_tree_Trich <- read_tree("phyloseq/tree_phyloseq_Trichobius.nwk")

OTU_Trich <- otu_table(features_Trich, taxa_are_rows=TRUE)
TAX_Trich <- tax_table(tax_Trich)
META_Trich <- sample_data(metadata_Trich)

taxa_names(OTU_Trich)
taxa_names(TAX_Trich)
taxa_names(phy_tree_Trich)

sample_names(OTU_Trich)
sample_names(META_Trich)

MM.Trichobius <- phyloseq(OTU_Trich, TAX_Trich, META_Trich, phy_tree_Trich)
MM.Trichobius
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1042 taxa and 57 samples ]
# sample_data() Sample Data:       [ 57 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 1042 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 1042 tips and 1041 internal nodes ]

features_Nyct <- read.csv("phyloseq/features_phyloseq_Nycterophilia.csv", sep=',', row.names=1, header=TRUE)
features_Nyct <- as.matrix(features_Nyct)
head(features_Nyct)

tax_Nyct <- read.csv("phyloseq/taxonomy_phyloseq_Nycterophilia.csv", sep=',', row.names=1, header=TRUE)
tax_Nyct <- as.matrix(tax_Nyct)
head(tax_Nyct)

metadata_Nyct<- read.csv("phyloseq/metadata_phyloseq_Nycterophilia.csv", sep=',', row.names=1, header=TRUE)
head(metadata_Nyct)

phy_tree_Nyct <- read_tree("phyloseq/tree_phyloseq_Nycterophilia.nwk")

OTU_Nyct <- otu_table(features_Nyct, taxa_are_rows=TRUE)
TAX_Nyct <- tax_table(tax_Nyct)
META_Nyct <- sample_data(metadata_Nyct)

taxa_names(OTU_Nyct)
taxa_names(TAX_Nyct)
taxa_names(phy_tree_Nyct)

sample_names(OTU_Nyct)
sample_names(META_Nyct)

MM.Nycterophilia <- phyloseq(OTU_Nyct, TAX_Nyct, META_Nyct, phy_tree_Nyct)
MM.Nycterophilia
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 756 taxa and 59 samples ]
# sample_data() Sample Data:       [ 59 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 756 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 756 tips and 755 internal nodes ]

#For Bat Swabs
features_Bat <- read.csv("phyloseq/features_phyloseq_Bat.csv", sep=',', row.names=1, header=TRUE)
features_Bat <- as.matrix(features_Bat)
head(features_Bat)

tax_Bat <- read.csv("phyloseq/taxonomy_phyloseq_Bat.csv", sep=',', row.names=1, header=TRUE)
tax_Bat <- as.matrix(tax_Bat)
head(tax_Bat)

metadata_Bat<- read.csv("phyloseq/metadata_phyloseq_Bat.csv", sep=',', row.names=1, header=TRUE)
head(metadata_Bat)

phy_tree_Bat <- read_tree("phyloseq/tree_phyloseq_Bat.nwk")

OTU_Bat <- otu_table(features_Bat, taxa_are_rows=TRUE)
TAX_Bat <- tax_table(tax_Bat)
META_Bat <- sample_data(metadata_Bat)

taxa_names(OTU_Bat)
taxa_names(TAX_Bat)
taxa_names(phy_tree_Bat)

sample_names(OTU_Bat)
sample_names(META_Bat)

MM.Bat <- phyloseq(OTU_Bat, TAX_Bat, META_Bat, phy_tree_Bat)
MM.Bat
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 4050 taxa and 42 samples ]
# sample_data() Sample Data:       [ 42 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 4050 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 4050 tips and 4048 internal nodes ]

#For Cave Swabs
features_Cave <- read.csv("phyloseq/features_phyloseq_Cave.csv", sep=',', row.names=1, header=TRUE)
features_Cave <- as.matrix(features_Cave)
head(features_Cave)

tax_Cave <- read.csv("phyloseq/taxonomy_phyloseq_Cave.csv", sep=',', row.names=1, header=TRUE)
tax_Cave <- as.matrix(tax_Cave)
head(tax_Cave)

metadata_Cave<- read.csv("phyloseq/metadata_phyloseq_Cave.csv", sep=',', row.names=1, header=TRUE)
head(metadata_Cave)

phy_tree_Cave <- read_tree("phyloseq/tree_phyloseq_Cave.nwk")

OTU_Cave <- otu_table(features_Cave, taxa_are_rows=TRUE)
TAX_Cave <- tax_table(tax_Cave)
META_Cave <- sample_data(metadata_Cave)

taxa_names(OTU_Cave)
taxa_names(TAX_Cave)
taxa_names(phy_tree_Cave)

sample_names(OTU_Cave)
sample_names(META_Cave)

MM.Cave <- phyloseq(OTU_Cave, TAX_Cave, META_Cave, phy_tree_Cave)
MM.Cave
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 609 taxa and 5 samples ]
# sample_data() Sample Data:       [ 5 samples by 33 sample variables ]
# tax_table()   Taxonomy Table:    [ 609 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 609 tips and 608 internal nodes ]

#philr transformed data
mm.philr<- read.csv(file="philr/MexicoMicrobBalances_philr.csv", header= TRUE, row.names= 1)
mm.25.philr<-read.csv(file="philr/MexicoMicrobBalances_25_philr.csv", header= TRUE, row.names= 1)
mm.Bat.philr<-read.csv(file="philr/MexicoMicrobBalances_Bat_philr.csv", header= TRUE, row.names= 1)
mm.Cave.philr<-read.csv( file="philr/MexicoMicrobBalances_Cave_philr.csv", header= TRUE, row.names= 1)
mm.BatFly.philr<-read.csv(file="philr/MexicoMicrobBalances_BatFly_philr.csv", header= TRUE, row.names= 1)
mm.Trichobius.philr<-read.csv( file="philr/MexicoMicrobBalances_Trichobius_philr.csv", header= TRUE, row.names= 1)
mm.Nycterophilia.philr<-read.csv(file="philr/MexicoMicrobBalances_Nycterophilia_philr.csv", header= TRUE, row.names= 1)
mm.NoEndo.philr<-read.csv(file="philr/MexicoMicrobBalances_NoEndo_philr.csv", header= TRUE, row.names= 1)

MM.philr<-readRDS("philr/MM.philr.RData")
MM.25.philr<-readRDS("philr/MM25.philr.RData")
MM.Cave.philr<-readRDS("philr/Cave.philr.RData")
MM.Bat.philr<-readRDS("philr/Bat.philr.RData")
MM.BatFly.philr<-readRDS("philr/Fly.philr.RData")
MM.Trichobius.philr<-readRDS( "philr/Trichobius.philr.RData")
MM.Nycterophilia.philr<-readRDS("philr/Nycterophilia.philr.RData")
MM.BatCave.philr<- readRDS("philr/BatCave.philr.RData")
MM.BatCaveTrichobius.philr <- readRDS("philr/BatCaveTrichobius.philr.RData")
MM.BatCaveNycterophilia.philr<- readRDS("philr/BatCaveNycterophilia.philr.Rdata")
MM.NoEndo.philr<- readRDS("philr/NoEndo.philr.Rdata")

#########Transform the data using philr##########################################################################
#https://bioconductor.org/packages/3.7/bioc/vignettes/philr/inst/doc/philr-intro.html

#Add a pseudocount of 1 to the remaining OTUs to avoid calculating log-ratios involving zeros
MM.philr <- transform_sample_counts(MM, function(x) x+1)
MM.25.philr <- transform_sample_counts(MM.25, function(x) x+1)
MM.Bat.philr <- transform_sample_counts(MM.Bat, function(x) x+1)
MM.Cave.philr <- transform_sample_counts(MM.Cave, function(x) x+1)
MM.BatFly.philr <- transform_sample_counts(MM.BatFly, function(x) x+1)
MM.Trichobius.philr <- transform_sample_counts(MM.Trichobius, function(x) x+1)
MM.Nycterophilia.philr <- transform_sample_counts(MM.Nycterophilia, function(x) x+1)
#A phyloseq object with no Arsenophonus or endosymbiont3
MM.NoEndo.philr <- subset_taxa(MM, !(Genus %in% c("endosymbionts3","Arsenophonus")))
MM.NoEndo.philr <- transform_sample_counts(MM.NoEndo.philr, function(x) x+1)

#Second: Process the Tree
#Check that the tree is rooted
is.rooted(phy_tree(MM.philr))
is.rooted(phy_tree(MM.25.philr))
is.rooted(phy_tree(MM.Bat.philr))
is.rooted(phy_tree(MM.Cave.philr))
is.rooted(phy_tree(MM.BatFly.philr))
is.rooted(phy_tree(MM.Trichobius.philr))
is.rooted(phy_tree(MM.Nycterophilia.philr))
is.rooted(phy_tree(MM.NoEndo.philr))

#Check that there are not multichotomies
is.binary.tree(phy_tree(MM.philr)) #False
is.binary.tree(phy_tree(MM.25.philr)) #False
is.binary.tree(phy_tree(MM.Bat.philr)) #False
is.binary.tree(phy_tree(MM.Cave.philr)) #True
is.binary.tree(phy_tree(MM.BatFly.philr)) #True
is.binary.tree(phy_tree(MM.Trichobius.philr)) #True
is.binary.tree(phy_tree(MM.Nycterophilia.philr)) #True
is.binary.tree(phy_tree(MM.NoEndo.philr)) #False

#If False, you the function multi2di from the ape package to replace multichotomies with a series
#dichotomies with one or several branches of zero length
phy_tree(MM.philr) <- multi2di(phy_tree(MM.philr))
phy_tree(MM.25.philr) <- multi2di(phy_tree(MM.25.philr))
phy_tree(MM.Bat.philr) <- multi2di(phy_tree(MM.Bat.philr))
phy_tree(MM.NoEndo.philr) <- multi2di(phy_tree(MM.NoEndo.philr))

is.binary.tree(phy_tree(MM.philr)) #True
is.binary.tree(phy_tree(MM.25.philr)) #True
is.binary.tree(phy_tree(MM.Bat.philr)) #True
is.binary.tree(phy_tree(MM.Cave.philr)) #True
is.binary.tree(phy_tree(MM.BatFly.philr)) #True
is.binary.tree(phy_tree(MM.Trichobius.philr)) #True
is.binary.tree(phy_tree(MM.Nycterophilia.philr)) #True
is.binary.tree(phy_tree(MM.NoEndo.philr)) #True

#Name the internal nodes fo the tree (n1, n2, n3...)
phy_tree(MM.philr) <- makeNodeLabel(phy_tree(MM.philr), method="number", prefix='n')
phy_tree(MM.25.philr) <- makeNodeLabel(phy_tree(MM.25.philr), method="number", prefix='n')
phy_tree(MM.Bat.philr) <- makeNodeLabel(phy_tree(MM.Bat.philr), method="number", prefix='n')
phy_tree(MM.Cave.philr) <- makeNodeLabel(phy_tree(MM.Cave.philr), method="number", prefix='n')
phy_tree(MM.BatFly.philr) <- makeNodeLabel(phy_tree(MM.BatFly.philr), method="number", prefix='n')
phy_tree(MM.Trichobius.philr) <- makeNodeLabel(phy_tree(MM.Trichobius.philr), method="number", prefix='n')
phy_tree(MM.Nycterophilia.philr) <- makeNodeLabel(phy_tree(MM.Nycterophilia.philr), method="number", prefix='n')
phy_tree(MM.NoEndo.philr) <- makeNodeLabel(phy_tree(MM.NoEndo.philr), method="number", prefix='n')

#To find out how the tree is rooted, use name.balance. The tree imported from QIIME2 is midpoint
#rooted, but if Archaea is used as the outgroup, than the balance will say "Kingdom_Bacteria/Kingdom_Bacteria"
name.balance(phy_tree(MM.philr), tax_table(MM.philr), 'n1')
name.balance(phy_tree(MM.25.philr), tax_table(MM.25.philr), 'n1')
name.balance(phy_tree(MM.Bat.philr), tax_table(MM.Bat.philr), 'n1')
name.balance(phy_tree(MM.Cave.philr), tax_table(MM.Cave.philr), 'n1')
name.balance(phy_tree(MM.BatFly.philr), tax_table(MM.BatFly.philr), 'n1')
name.balance(phy_tree(MM.Trichobius.philr), tax_table(MM.Trichobius.philr), 'n1')
name.balance(phy_tree(MM.Nycterophilia.philr), tax_table(MM.Nycterophilia.philr), 'n1')
name.balance(phy_tree(MM.NoEndo.philr), tax_table(MM.NoEndo.philr), 'n1')

saveRDS(MM.philr, "philr/MM.philr.RData")
saveRDS(MM.25.philr, "philr/MM25.philr.RData")
saveRDS(MM.Cave.philr, "philr/Cave.philr.RData")
saveRDS(MM.Bat.philr, "philr/Bat.philr.RData")
saveRDS(MM.BatFly.philr, "philr/Fly.philr.RData")
saveRDS(MM.Trichobius.philr, "philr/Trichobius.philr.RData")
saveRDS(MM.Nycterophilia.philr, "philr/Nycterophilia.philr.RData")
saveRDS(MM.NoEndo.philr, "philr/NoEndo.philr.RData")

#Third: Investigate the Dataset Components
#Because philr uses the compositions package, taxa must be columns and samples must be rows
#We have to transpose our data to get it in this format, because it is exported from QIIME2
#with the taxa as rows and samples as columns

otu.table<- t(otu_table(MM.philr))
write.csv(otu.table, "phyloseq/philr/MM.philr.otu.csv", row.names = TRUE)
tree <- phy_tree(MM.philr)
write.tree(tree, "phyloseq/philr/MM.philr.tree.nwk")
metadata<- sample_data(MM.philr)
write.csv(metadata, "phyloseq/philr/MM.philr.metadata.csv", row.names =TRUE)
tax <- tax_table(MM.philr)
write.csv(tax, "phyloseq/philr/MM.philr.tax.csv", row.names = TRUE)

otu.table.25<- t(otu_table(MM.25.philr))
write.csv(otu.table.25, "phyloseq/philr/MM25.philr.otu.csv", row.names = TRUE)
tree.25 <- phy_tree(MM.25.philr)
write.tree(tree.25, "phyloseq/philr/MM25.philr.tree.nwk")
metadata.25<- sample_data(MM.25.philr)
write.csv(metadata.25, "phyloseq/philr/MM25.philr.metadata.csv", row.names =TRUE)
tax.25 <- tax_table(MM.25.philr)
write.csv(tax.25, "phyloseq/philr/MM25.philr.tax.csv", row.names = TRUE)

otu.table.Cave<- t(otu_table(MM.Cave.philr))
write.csv(otu.table.Cave, "phyloseq/philr/MMCave.philr.otu.csv", row.names = TRUE)
tree.Cave <- phy_tree(MM.Cave.philr)
write.tree(tree.Cave, "phyloseq/philr/MMCave.philr.tree.nwk")
metadata.Cave<- sample_data(MM.Cave.philr)
write.csv(metadata.Cave, "phyloseq/philr/MMCave.philr.metadata.csv", row.names =TRUE)
tax.Cave <- tax_table(MM.Cave.philr)
write.csv(tax.Cave, "phyloseq/philr/MMCave.philr.tax.csv", row.names = TRUE)

otu.table.Bat<- t(otu_table(MM.Bat.philr))
write.csv(otu.table.Bat, "phyloseq/philr/MMBat.philr.otu.csv", row.names = TRUE)
tree.Bat <- phy_tree(MM.Bat.philr)
write.tree(tree.Bat, "phyloseq/philr/MMBat.philr.tree.nwk")
metadata.Bat<- sample_data(MM.Bat.philr)
write.csv(metadata.Bat, "phyloseq/philr/MMBat.philr.metadata.csv", row.names =TRUE)
tax.Bat <- tax_table(MM.Bat.philr)
write.csv(tax.Bat, "phyloseq/philr/MMBat.philr.tax.csv", row.names = TRUE)

otu.table.BatFly<- t(otu_table(MM.BatFly.philr))
write.csv(otu.table.BatFly, "phyloseq/philr/MMBatFly.philr.otu.csv", row.names = TRUE)
tree.BatFly <- phy_tree(MM.BatFly.philr)
write.tree(tree.BatFly, "phyloseq/philr/MMBatFly.philr.tree.nwk")
metadata.BatFly<- sample_data(MM.BatFly.philr)
write.csv(metadata.BatFly, "phyloseq/philr/MMBatFly.philr.metadata.csv", row.names =TRUE)
tax.BatFly <- tax_table(MM.BatFly.philr)
write.csv(tax.BatFly, "phyloseq/philr/MMBatFly.philr.tax.csv", row.names = TRUE)

otu.table.Trichobius<- t(otu_table(MM.Trichobius.philr))
write.csv(otu.table.Trichobius, "phyloseq/philr/MMTrichobius.philr.otu.csv", row.names = TRUE)
tree.Trichobius <- phy_tree(MM.Trichobius.philr)
write.tree(tree.Trichobius, "phyloseq/philr/MMTrichobius.philr.tree.nwk")
metadata.Trichobius<- sample_data(MM.Trichobius.philr)
write.csv(metadata.Trichobius, "phyloseq/philr/MMTrichobius.philr.metadata.csv", row.names =TRUE)
tax.Trichobius <- tax_table(MM.Trichobius.philr)
write.csv(tax.Trichobius, "phyloseq/philr/MMTrichobius.philr.tax.csv", row.names = TRUE)

otu.table.Nycterophilia<- t(otu_table(MM.Nycterophilia.philr))
write.csv(otu.table.Nycterophilia, "phyloseq/philr/MMNycterophilia.philr.otu.csv", row.names = TRUE)
tree.Nycterophilia <- phy_tree(MM.Nycterophilia.philr)
write.tree(tree.Nycterophilia, "phyloseq/philr/MMNycterophilia.philr.tree.nwk")
metadata.Nycterophilia<- sample_data(MM.Nycterophilia.philr)
write.csv(metadata.Nycterophilia, "phyloseq/philr/MMNycterophilia.philr.metadata.csv", row.names =TRUE)
tax.Nycterophilia <- tax_table(MM.Nycterophilia.philr)
write.csv(tax.Nycterophilia, "phyloseq/philr/MMNycterophilia.philr.tax.csv", row.names = TRUE)

otu.table.NoEndo<- t(otu_table(MM.NoEndo.philr))
write.csv(otu.table.NoEndo, "philr/MMNoEndo.philr.otu.csv", row.names = TRUE)
tree.NoEndo <- phy_tree(MM.NoEndo.philr)
write.tree(tree.NoEndo, "philr/MMNoEndo.philr.tree.nwk")
metadata.NoEndo<- sample_data(MM.NoEndo.philr)
write.csv(metadata.NoEndo, "philr/MMNoEndo.philr.metadata.csv", row.names =TRUE)
tax.NoEndo <- tax_table(MM.NoEndo.philr)
write.csv(tax.NoEndo, "philr/MMNoEndo.philr.tax.csv", row.names = TRUE)

#Check that the samples are the rows and the taxa are the columns
otu.table[1:2, 1:2]
otu.table.25[1:2, 1:2]
otu.table.Bat[1:2, 1:2]
otu.table.Cave[1:2, 1:2]
otu.table.BatFly[1:2, 1:2]
otu.table.Trichobius[1:2, 1:2]
otu.table.Nycterophilia[1:2, 1:2]
otu.table.NoEndo[1:2, 1:2]

#Fourth: Transform the data using PhILR
#1)The function philr::philr() implements a user friendly wrapper for the key steps in the philr transform.
#2)Convert the phylogenetic tree to its sequential binary partition (SBP) representation using the function philr::phylo2sbp()
#3)Calculate the weighting of the taxa (aka parts) or use the user specified weights
#4)Built the contrast matrix from the SBP and taxa weights using the function philr::buildilrBasep()
#5)Convert OTU table to relative abundance (using philr::miniclo()) and ‘shift’ dataset using the weightings [@egozcue2016] using the function philr::shiftp().
#6)Transform the data to PhILR space using the function philr::ilrp()
#7)Weight the resulting PhILR space using phylogenetic distance. These weights are either provided by the user or can be calculated by the function philr::calculate.blw().

mm.philr <- philr(otu.table, tree, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.philr[1:5,1:5]
# n1        n2          n3          n4          n5
# Sample249594_1 0.03422963 0.1172247  0.02277608 -0.05674337 -0.03829157
# Sample249594_2 0.04812641 0.1171539  0.02277608 -0.05674337 -0.03829157
# Sample249595_2 0.08976235 0.1164949  0.02277608 -0.05674337 -0.03829157
# Sample249600_1 0.06440504 0.1157579  0.02277608 -0.05674337 -0.03829157
# Sample249600_2 0.02844889 0.2225585 -0.12995969 -0.05674337 -0.03829157

mm.25.philr <- philr(otu.table.25, tree.25, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.25.philr[1:5,1:5]
# n1         n2         n3         n4         n5
# Sample249594_1 -0.060085148 0.06776035 0.02427531 0.02492734 0.04489876
# Sample249594_2 -0.055334446 0.06776035 0.02427531 0.02492734 0.04489876
# Sample249595_2  0.004576144 0.06776035 0.02427531 0.02492734 0.04489876
# Sample249600_1 -0.029592265 0.06776035 0.02427531 0.02492734 0.04489876
# Sample249600_2 -0.080224068 0.06707494 0.02427531 0.02629122 0.04393488

mm.Bat.philr <- philr(otu.table.Bat, tree.Bat, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.Bat.philr[1:5,1:5]
# n1          n2          n3          n4          n5
# Sample475V 0.557287249 -0.02364489 0.003778679 -0.09700996 -0.09029978
# Sample478V 0.001030591  0.09864435 0.003778679 -0.09700996 -0.09029978
# Sample482V 0.163286711  0.07453813 0.003778679 -0.09700996 -0.09029978
# Sample483V 0.035974058  0.12312931 0.003778679 -0.09700996 -0.09029978
# Sample488V 0.189507101  0.06619169 0.003778679 -0.09700996 -0.09029978

mm.Cave.philr <- philr(otu.table.Cave, tree.Cave, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.Cave.philr[1:5,1:5]
# n1          n2          n3          n4          n5
# SampleCarmen   -3.047482  1.01506360  1.64988675 -0.68521526  1.37083384
# SampleChivato  -1.078135 -0.28794263  0.08668882 -0.28137164  0.08310712
# SampleFabrica  -1.323376 -0.13358925  0.08668882 -0.10834760  0.16809749
# SampleLaGitana -1.256973 -0.08649807  0.08668882 -0.05556025  0.19402693
# SamplePanchito  5.756402 -0.68837614 -1.79365469  0.96581927 -1.44945951

mm.BatFly.philr <- philr(otu.table.BatFly, tree.BatFly, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.BatFly.philr[1:5,1:5]
# n1         n2         n3          n4           n5
# Sample249594_1 0.05348690 0.01612632  0.1712401 -0.00690442  0.006978459
# Sample249594_2 0.08181427 0.01589888  0.1712401 -0.00690442  0.007385305
# Sample249595_2 0.32369802 0.01380532  0.1712401 -0.00690442  0.011130317
# Sample249600_1 0.26096703 0.01152694  0.1712401 -0.00690442 -0.026984240
# Sample249600_2 0.03416105 0.23688325 -0.1502843 -0.09089419  0.011440982

mm.Trichobius.philr <- philr(otu.table.Trichobius, tree.Trichobius, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.Trichobius.philr[1:5,1:5]
# n1        n2          n3        n4        n5
# Sample249595_2 -0.1044917 0.2126220 -0.04768257 0.1259911 0.1603523
# Sample249600_1 -0.5960342 0.2068742 -0.11197051 0.1221232 0.1507724
# Sample249601_2 -0.8580391 0.2205042 -0.05686135 0.1334551 0.1462591
# Sample249601_3 -0.9231108 0.2219568 -0.05855292 0.1348307 0.1269344
# Sample249614_1 -0.7302412 0.2042480 -0.41448594 0.1264161 0.1471593

mm.Nycterophilia.philr <- philr(otu.table.Nycterophilia, tree.Nycterophilia, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.Nycterophilia.philr[1:5,1:5]
# n1          n2           n3         n4        n5
# Sample249594_1 -0.2285950 -0.08141277 -0.008192069 0.07753846 0.1952403
# Sample249594_2 -0.1086718 -0.08273446 -0.008192069 0.07846243 0.1939959
# Sample249600_2 -0.4266147  0.52130995 -0.116305201 0.08725881 0.1821494
# Sample249600_3 -0.8455820 -0.07836922 -0.008192069 0.07541079 0.1981057
# Sample249601_1 -0.3337377 -0.07880901 -0.008192069 0.07571824 0.1976917

mm.NoEndo.philr <- philr(otu.table.NoEndo, tree.NoEndo, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.NoEndo.philr[1:5,1:5]
# n1        n2          n3          n4          n5
# Sample249594_1 0.10467874 0.1172247  0.02277608 -0.05674337 -0.03829157
# Sample249594_2 0.11423158 0.1171539  0.02277608 -0.05674337 -0.03829157
# Sample249595_2 0.07881886 0.1164949  0.02277608 -0.05674337 -0.03829157
# Sample249600_1 0.07737474 0.1157579  0.02277608 -0.05674337 -0.03829157
# Sample249600_2 0.10838155 0.2225585 -0.12995969 -0.05674337 -0.03829157

#Now the transformed data is represented in terms of balances and since each balance is associated with a single internal node of the tree, we denote the balances using the same names we assigned to the internal nodes (e.g., n1).
#Write the balance matrix out
write.csv(mm.philr, file="MexicoMicrobBalances_philr.csv")
write.csv(mm.25.philr, file="MexicoMicrobBalances_25_philr.csv")
write.csv(mm.Bat.philr, file="MexicoMicrobBalances_Bat_philr.csv")
write.csv(mm.Cave.philr, file="MexicoMicrobBalances_Cave_philr.csv")
write.csv(mm.BatFly.philr, file="MexicoMicrobBalances_BatFly_philr.csv")
write.csv(mm.Trichobius.philr, file="MexicoMicrobBalances_Trichobius_philr.csv")
write.csv(mm.Nycterophilia.philr, file="MexicoMicrobBalances_Nycterophilia_philr.csv")
write.csv(mm.NoEndo.philr, file="philr/MexicoMicrobBalances_NoEndo_philr.csv")

###Do PCoA on distance between philr balances###########################################
#https://bioconductor.org/packages/3.7/bioc/vignettes/philr/inst/doc/philr-intro.html
mm.philr <- read.csv("MexicoMicrobBalances_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.philr <- as.matrix(mm.philr)

mm.25.philr <- read.csv("MexicoMicrobBalances_25_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.25.philr <- as.matrix(mm.25.philr)

mm.Bat.philr <- read.csv("MexicoMicrobBalances_Bat_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.Bat.philr <- as.matrix(mm.Bat.philr)

mm.Cave.philr <- read.csv("MexicoMicrobBalances_Cave_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.Cave.philr <- as.matrix(mm.Cave.philr)

mm.BatFly.philr <- read.csv("MexicoMicrobBalances_BatFly_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.BatFly.philr <- as.matrix(mm.BatFly.philr)

mm.Trichobius.philr <- read.csv("MexicoMicrobBalances_Trichobius_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.Trichobius.philr <- as.matrix(mm.Trichobius.philr)

mm.Nycterophilia.philr <- read.csv("MexicoMicrobBalances_Nycterophilia_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.Nycterophilia.philr <- as.matrix(mm.Nycterophilia.philr)

mm.NoEndo.philr <- read.csv("philr/MexicoMicrobBalances_NoEndo_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.NoEndo.philr <- as.matrix(mm.NoEndo.philr)

#create distance matrix using euclidean distance
mm.dist <- dist(mm.philr, method="euclidean")
mm.25.dist <- dist(mm.25.philr, method="euclidean")
mm.Bat.dist <- dist(mm.Bat.philr, method="euclidean")
mm.Cave.dist <- dist(mm.Cave.philr, method="euclidean")
mm.BatFly.dist <- dist(mm.BatFly.philr, method="euclidean")
mm.Trichobius.dist <- dist(mm.Trichobius.philr, method="euclidean")
mm.Nycterophilia.dist <- dist(mm.Nycterophilia.philr, method="euclidean")
mm.NoEndo.dist <- dist(mm.NoEndo.philr, method="euclidean")

#Do PCoA on distance matrix
mm.pcoa <- ordinate(MM.philr, 'PCoA', distance=mm.dist)
mm.25.pcoa <- ordinate(MM.25.philr, 'PCoA', distance=mm.25.dist)
mm.Bat.pcoa <- ordinate(MM.Bat.philr, 'PCoA', distance=mm.Bat.dist)
mm.Cave.pcoa <- ordinate(MM.Cave.philr, 'PCoA', distance=mm.Cave.dist)
mm.BatFly.pcoa <- ordinate(MM.BatFly.philr, 'PCoA', distance=mm.BatFly.dist)
mm.Trichobius.pcoa <- ordinate(MM.Trichobius.philr, 'PCoA', distance=mm.Trichobius.dist)
mm.Nycterophilia.pcoa <- ordinate(MM.Nycterophilia.philr, 'PCoA', distance=mm.Nycterophilia.dist)
mm.NoEndo.pcoa <- ordinate(MM.NoEndo.philr, 'PCoA', distance=mm.NoEndo.dist)

#Plot PCoA of MM.philr with Sample Source as the color and Mainland/Baja as the shape
sample_data(MM.philr)$SampleSource <- factor(sample_data(MM.philr)$SampleSource, levels=c("Cave","Bat","Nycterophilia","Trichobius"))

pdf(file = "PCoAPlots/MM_SampleSource_BroadLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.philr, mm.pcoa, color='SampleSource', shape='BroadLocality') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

#Plot PCoA of MM.philr with Collection Locality as Color and Sample Source as Shape
sample_data(MM.philr)$CollectionLocality<- factor(sample_data(MM.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur",
                                                                                                             "Chivato, Sierra Cacachilas, Baja California Sur",
                                                                                                             "La Gitana, Sierra Cacachilas, Baja California Sur",
                                                                                                             "Isla Panchito, near Chamela, Jalisco, MX",
                                                                                                             "Cueva de la Fabrica, Coquimatlan, Colima, MX"))

pdf(file = "PCoAPlots/MM_CollectionLoc_SampleSource.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.philr, mm.pcoa, color='CollectionLocality', shape='SampleSource') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

#Plot PCoA of MM.25.philr with Sample Source as the color and Mainland/Baja as the shape
sample_data(MM.25.philr)$SampleSource <- factor(sample_data(MM.25.philr)$SampleSource, levels=c("Cave","Bat","Nycterophilia","Trichobius"))

pdf(file = "PCoAPlots/MM25_SampleSource_BroadLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.25.philr, mm.25.pcoa, color='SampleSource', shape='BroadLocality') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

#Plot PCoA of MM.25.philr with Collection Locality as Color and Sample Source as Shape
sample_data(MM.25.philr)$CollectionLocality<- factor(sample_data(MM.25.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur",
                                                                                                     "Chivato, Sierra Cacachilas, Baja California Sur",
                                                                                                     "La Gitana, Sierra Cacachilas, Baja California Sur",
                                                                                                     "Isla Panchito, near Chamela, Jalisco, MX",
                                                                                                     "Cueva de la Fabrica, Coquimatlan, Colima, MX"))

pdf(file = "PCoAPlots/MM25_CollectionLoc_SampleSource.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.25.philr, mm.25.pcoa, color='CollectionLocality', shape='SampleSource') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

#Plot PCoA of MM.Cave.philr with Collection Locality as the shape
pdf(file = "PCoAPlots/Cave_CollectionLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Cave.philr, mm.Cave.pcoa, shape ='CollectionLocality') 
p + geom_point(size=2, stroke=0) #+ scale_color_viridis_d(option = "magma")
dev.off()

#Plot PCoA of MM.Bat.philr with Collection Locality as the color and Mainland/Baja as shape
#Order the Collection Localities by Mainland or Baja
sample_data(MM.Bat.philr)$CollectionLocality<- factor(sample_data(MM.Bat.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur",
                                                                                                             "Chivato, Sierra Cacachilas, Baja California Sur",
                                                                                                             "La Gitana, Sierra Cacachilas, Baja California Sur",
                                                                                                             "Isla Panchito, near Chamela, Jalisco, MX",
                                                                                                             "Cueva de la Fabrica, Coquimatlan, Colima, MX"))

pdf(file = "PCoAPlots/Bat_CollectionLoc_BroadLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Bat.philr, mm.Bat.pcoa, color='CollectionLocality',shape='BroadLocality') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#Plot PCoA of MM.Bat.philr and MM.Cave.philr with Collection Locality as the color and Sample type as the shape
#Subset MM.philr, do philr transform, calculate euclidean distance and pcoa
MM.BatCave.philr <- subset_samples(MM.philr, SampleSource == "Bat" | SampleSource == "Cave")
MM.BatCave.philr <-prune_taxa(taxa_sums(MM.BatCave.philr)>0, MM.BatCave.philr)
saveRDS(MM.BatCave.philr, "phyloseq/philr/BatCave.philr.RData")

otu.table.BatCave<- t(otu_table(MM.BatCave.philr))
tree.BatCave <- phy_tree(MM.BatCave.philr)
metadata.BatCave<- sample_data(MM.BatCave.philr)
tax.BatCave <- tax_table(MM.BatCave.philr)
str(metadata.BatCave)

mm.BatCave.philr <- philr(otu.table.BatCave, tree.BatCave, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.BatCave.philr[1:5,1:5]
# n1         n2          n3          n4          n5
# Sample475V  0.457248142 0.02606351 -0.09719067 -0.09467946 -0.08573297
# Sample478V -0.026846070 0.13106907 -0.09719067 -0.09467946 -0.08573297
# Sample482V  0.112955119 0.11047566 -0.09719067 -0.09467946 -0.08573297
# Sample483V  0.002036003 0.15193166 -0.09719067 -0.09467946 -0.08573297
# Sample488V  0.138488394 0.10317084 -0.09719067 -0.09467946 -0.08573297

write.csv(mm.BatCave.philr, file="MexicoMicrobBalances_BatCave_philr.csv")
mm.BatCave.philr <- read.csv("philr/MexicoMicrobBalances_BatCave_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.BatCave.philr <- as.matrix(mm.BatCave.philr)

mm.BatCave.dist <- dist(mm.BatCave.philr, method="euclidean")
mm.BatCave.pcoa <- ordinate(MM.BatCave.philr, 'PCoA', distance=mm.BatCave.dist)

sample_data(MM.BatCave.philr)$CollectionLocality<- factor(sample_data(MM.BatCave.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur",
                                                                                                                     "Chivato, Sierra Cacachilas, Baja California Sur",
                                                                                                                     "La Gitana, Sierra Cacachilas, Baja California Sur",
                                                                                                                     "Isla Panchito, near Chamela, Jalisco, MX",
                                                                                                                     "Cueva de la Fabrica, Coquimatlan, Colima, MX"))

pdf(file = "PCoAPlots/BatCave_CollectionLoc_SampleSource.pdf", width = 8, height = 5)
p<-plot_ordination(MM.BatCave.philr, mm.BatCave.pcoa, color='CollectionLocality', shape='SampleSource') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma") +
  stat_ellipse(aes(data=BroadLocality), type = "euclid")
dev.off()

#Plot PCoA of MM.BatFly.philr with Parasite Species as color and Collection Locality as shape
sample_data(MM.BatFly.philr)$CollectionLocality<- factor(sample_data(MM.BatFly.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur",
                                                                                                                     "Chivato, Sierra Cacachilas, Baja California Sur",
                                                                                                                     "La Gitana, Sierra Cacachilas, Baja California Sur",
                                                                                                                     "Isla Panchito, near Chamela, Jalisco, MX",
                                                                                                                     "Cueva de la Fabrica, Coquimatlan, Colima, MX"))


pdf(file = "PCoAPlots/Fly_Parasite_ColLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.BatFly.philr, mm.BatFly.pcoa, shape='CollectionLocality', color='BatFlyGenus') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of MM.BatFly.philr, color = Parasite Sex, shape = Parasite Species
pdf(file = "PCoAPlots/Fly_ParasiteSex_Parasite.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.BatFly.philr, mm.BatFly.pcoa, color='BatFlySex', shape='BatFlyGenus') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of MM.Trichobius, shape = Bat Sex, color = Collection Locality
pdf(file = "PCoAPlots/Trichobius_CollectionLoc_BatSex.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Trichobius.philr, mm.Trichobius.pcoa, shape='BatSex', color='CollectionLocality') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of MM.Trichobius, color = Parasite Sex, shape = Broad Locality
pdf(file = "PCoAPlots/Trichobius_ParasiteSex_BroadLocality.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Trichobius.philr, mm.Trichobius.pcoa, color='BatFlySex', shape = "BroadLocality") 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of MM.Trichobius, color = Bat Sex, shape = Parasite Sex
pdf(file = "PCoAPlots/Trichobius_BatSex_ParasiteSex.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Trichobius.philr, mm.Trichobius.pcoa, color='BatSex', shape = "BatFlySex") 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of Combined Cave, Bat, and Trichobius Samples, color = Sample Source
#Create a new phyloseq object and redo philr transform
MM.BatCaveTrichobius.philr <- subset_samples(MM.philr, SampleSource == "Bat" | SampleSource == "Cave" | SampleSource == "Trichobius")
MM.BatCaveTrichobius.philr <-prune_taxa(taxa_sums(MM.BatCaveTrichobius.philr)>0, MM.BatCaveTrichobius.philr)
saveRDS(MM.BatCaveTrichobius.philr, "phyloseq/philr/BatCaveTrichobius.philr.RData")

otu.table.BatCaveTrichobius<- t(otu_table(MM.BatCaveTrichobius.philr))
tree.BatCaveTrichobius <- phy_tree(MM.BatCaveTrichobius.philr)
metadata.BatCaveTrichobius<- sample_data(MM.BatCaveTrichobius.philr)
tax.BatCaveTrichobius <- tax_table(MM.BatCaveTrichobius.philr)
str(metadata.BatCaveTrichobius)

mm.BatCaveTrichobius.philr <- philr(otu.table.BatCaveTrichobius, tree.BatCaveTrichobius, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.BatCaveTrichobius.philr[1:5,1:5]

write.csv(mm.BatCaveTrichobius.philr, file="MexicoMicrobBalances_BatCaveTrichobius_philr.csv")
mm.BatCaveTrichobius.philr <- read.csv("MexicoMicrobBalances_BatCaveTrichobius_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.BatCaveTrichobius.philr <- as.matrix(mm.BatCaveTrichobius.philr)

mm.BatCaveTrichobius.dist <- dist(mm.BatCaveTrichobius.philr, method="euclidean")
mm.BatCaveTrichobius.pcoa <- ordinate(MM.BatCaveTrichobius.philr, 'PCoA', distance=mm.BatCaveTrichobius.dist)

sample_data(MM.BatCaveTrichobius.philr)$CollectionLocality<- factor(sample_data(MM.BatCaveTrichobius.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur",
                                                                                                                     "Chivato, Sierra Cacachilas, Baja California Sur",
                                                                                                                     "La Gitana, Sierra Cacachilas, Baja California Sur",
                                                                                                                     "Isla Panchito, near Chamela, Jalisco, MX",
                                                                                                                     "Cueva de la Fabrica, Coquimatlan, Colima, MX"))

pdf(file = "PCoAPlots/BatCaveTrichobius_SampleSource.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.BatCaveTrichobius.philr, mm.BatCaveTrichobius.pcoa, color='SampleSource') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d()
dev.off()

#PCoA of MM.Nycterophilia, color = Bat Sex, shape = Collection Locality
pdf(file = "PCoAPlots/Nycterophilia_BatSex_CollectionLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Nycterophilia.philr, mm.Nycterophilia.pcoa, color='BatSex', shape='CollectionLocality') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of MM.Nycterophilia, color = Parasite Sex, shape = Broad Locality
pdf(file = "PCoAPlots/Nycterophilia_ParasiteSex_BroadLocality.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Nycterophilia.philr, mm.Nycterophilia.pcoa, color='BatFlySex', shape = "BroadLocality") 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of MM.Nycterophilia, color = Bat Sex, shape = Parasite Sex
pdf(file = "PCoAPlots/Nycterophilia_BatSex_ParasiteSex.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Nycterophilia.philr, mm.Nycterophilia.pcoa, color='BatSex', shape = "BatFlySex") 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option = "magma")
dev.off()

#PCoA of Combined Cave, Bat, and Nycterophilia Samples, color = Sample Source
#Create a new phyloseq object and redo philr transform
MM.BatCaveNycterophilia.philr <- subset_samples(MM.philr, SampleSource == "Bat" | SampleSource == "Cave" | SampleSource == "Nycterophilia")
MM.BatCaveNycterophilia.philr <-prune_taxa(taxa_sums(MM.BatCaveNycterophilia.philr)>0, MM.BatCaveNycterophilia.philr)
saveRDS(MM.BatCaveNycterophilia.philr, "phyloseq/philr/BatCaveNycterophilia.philr.RData")

otu.table.BatCaveNycterophilia<- t(otu_table(MM.BatCaveNycterophilia.philr))
tree.BatCaveNycterophilia <- phy_tree(MM.BatCaveNycterophilia.philr)
metadata.BatCaveNycterophilia<- sample_data(MM.BatCaveNycterophilia.philr)
tax.BatCaveNycterophilia <- tax_table(MM.BatCaveNycterophilia.philr)
str(metadata.BatCaveNycterophilia)

mm.BatCaveNycterophilia.philr <- philr(otu.table.BatCaveNycterophilia, tree.BatCaveNycterophilia, part.weights='enorm.x.gm.counts', ilr.weights='blw.sqrt')
mm.BatCaveNycterophilia.philr[1:5,1:5]

write.csv(mm.BatCaveNycterophilia.philr, file="MexicoMicrobBalances_BatCaveNycterophilia_philr.csv")
mm.BatCaveNycterophilia.philr <- read.csv("MexicoMicrobBalances_BatCaveNycterophilia_philr.csv", header = TRUE, stringsAsFactors = F, check.names=F, row.names=1)
mm.BatCaveNycterophilia.philr <- as.matrix(mm.BatCaveNycterophilia.philr)

mm.BatCaveNycterophilia.dist <- dist(mm.BatCaveNycterophilia.philr, method="euclidean")
mm.BatCaveNycterophilia.pcoa <- ordinate(MM.BatCaveNycterophilia.philr, 'PCoA', distance=mm.BatCaveNycterophilia.dist)

sample_data(MM.BatCaveNycterophilia.philr)$CollectionLocality<- factor(sample_data(MM.BatCaveNycterophilia.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur",
                                                                                                                                         "Chivato, Sierra Cacachilas, Baja California Sur",
                                                                                                                                         "La Gitana, Sierra Cacachilas, Baja California Sur",
                                                                                                                                         "Isla Panchito, near Chamela, Jalisco, MX",
                                                                                                                                         "Cueva de la Fabrica, Coquimatlan, Colima, MX"))

pdf(file = "PCoAPlots/BatCaveNycterophilia_SampleSource.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.BatCaveNycterophilia.philr, mm.BatCaveNycterophilia.pcoa, color='SampleSource') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d()
dev.off()

#PCoA of Bartonella and Wolbachia presence/absence
#First subset the Trichobius and Nycterophilia phyloseq objects by whether the samples have Bartonella
MM.Trichobius.Bart = subset_taxa(MM.Trichobius, Genus=="Bartonella")
MM.Trichobius.Bart = prune_samples(sample_sums(MM.Trichobius.Bart)> 0, MM.Trichobius.Bart)

MM.Nycterophilia.Bart = subset_taxa(MM.Nycterophilia, Genus=="Bartonella")
MM.Nycterophilia.Bart = prune_samples(sample_sums(MM.Nycterophilia.Bart)> 0, MM.Nycterophilia.Bart)

#Next subset the Trichobius and Nycterophilia phyloseq objects by whether the samples have Wolbachia
MM.Trichobius.Wol = subset_taxa(MM.Trichobius, Genus=="Wolbachia")
MM.Trichobius.Wol = prune_samples(sample_sums(MM.Trichobius.Wol)> 0, MM.Trichobius.Wol)

MM.Nycterophilia.Wol = subset_taxa(MM.Nycterophilia, Genus=="Wolbachia")
MM.Nycterophilia.Wol = prune_samples(sample_sums(MM.Nycterophilia.Wol)> 0, MM.Nycterophilia.Wol)

#Create a list of sample names in each subset
MM.Trichobius.Bart.l <- sample_names(MM.Trichobius.Bart)
MM.Nycterophilia.Bart.l <- sample_names(MM.Nycterophilia.Bart)
MM.Trichobius.Wol.l <- sample_names(MM.Trichobius.Wol)
MM.Nycterophilia.Wol.l <- sample_names(MM.Nycterophilia.Wol)

#Create a list of the samples that have both Wolbachia and Bartonella
MM.Trichobius.BartWol.l <- intersect(MM.Trichobius.Bart.l, MM.Trichobius.Wol.l)
MM.Nycterophilia.BartWol.l <- intersect(MM.Nycterophilia.Bart.l, MM.Nycterophilia.Wol.l) #No shared infections

#Create a column in the phyloseq sample data for Trichobius and Nycterophilia that indicates whether they have
#Bartonella, Wolbachia, or both
#Nycterophilia has no shared infections
sample_data(MM.Trichobius.philr)$BartWol[sample_names(MM.Trichobius.philr) %in% MM.Trichobius.Bart.l]<-"Bartonella"
sample_data(MM.Trichobius.philr)$BartWol[sample_names(MM.Trichobius.philr) %in% MM.Trichobius.Wol.l]<-"Wolbachia"
sample_data(MM.Trichobius.philr)$BartWol[sample_names(MM.Trichobius.philr) %in% MM.Trichobius.BartWol.l]<-"Both"
sample_data(MM.Trichobius.philr)$BartWol[is.na(sample_data(MM.Trichobius.philr)$BartWol)]<-"Not Detected"

sample_data(MM.Nycterophilia.philr)$BartWol[sample_names(MM.Nycterophilia.philr) %in% MM.Nycterophilia.Bart.l]<-"Bartonella"
sample_data(MM.Nycterophilia.philr)$BartWol[sample_names(MM.Nycterophilia.philr) %in% MM.Nycterophilia.Wol.l]<-"Wolbachia"
sample_data(MM.Nycterophilia.philr)$BartWol[is.na(sample_data(MM.Nycterophilia.philr)$BartWol)]<-"Not Detected"

#Now plot the PCoA 
sample_data(MM.Trichobius.philr)$BartWol <- factor(sample_data(MM.Trichobius.philr)$BartWol, levels=c("Bartonella", "Wolbachia", "Both", "Not Detected"))

pdf(file = "PCoAPlots/Trichobius_BartWol_FlySex.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Trichobius.philr, mm.Trichobius.pcoa, color='BartWol', shape='BatFlySex') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

pdf(file = "PCoAPlots/Trichobius_BartWol_ColLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Trichobius.philr, mm.Trichobius.pcoa, color='BartWol', shape='CollectionLocality') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

sample_data(MM.Nycterophilia.philr)$BartWol <- factor(sample_data(MM.Nycterophilia.philr)$BartWol, levels=c("Bartonella", "Wolbachia", "Not Detected"))

pdf(file = "PCoAPlots/Nycterophilia_BartWol_FlySex.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Nycterophilia.philr, mm.Nycterophilia.pcoa, color='BartWol', shape='BatFlySex') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

pdf(file = "PCoAPlots/Nycterophilia_BartWol_ColLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.Nycterophilia.philr, mm.Nycterophilia.pcoa, color='BartWol', shape='CollectionLocality') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

#Plot a PCoA with no endosymbionts in either bat fly
sample_data(MM.NoEndo.philr)$SampleSource <- factor(sample_data(MM.NoEndo.philr)$SampleSource, levels=c("Cave", "Bat", "Nycterophilia", "Trichobius"))
sample_data(MM.NoEndo.philr)$CollectionLocality <- factor(sample_data(MM.NoEndo.philr)$CollectionLocality, levels=c("Carmen, Loreto, Baja California Sur", 
                                                                                                                    "Chivato, Sierra Cacachilas, Baja California Sur", 
                                                                                                                    "La Gitana, Sierra Cacachilas, Baja California Sur", 
                                                                                                                    "Isla Panchito, near Chamela, Jalisco, MX", 
                                                                                                                    "Cueva de la Fabrica, Coquimatlan, Colima, MX"))

pdf(file = "PCoAPlots/NoEndo_SampleSource_ColLoc.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.NoEndo.philr, mm.NoEndo.pcoa, color='SampleSource', shape="CollectionLocality") 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1)
dev.off()

pdf(file = "PCoAPlots/NoEndo_SampleSource_Ellipse.pdf", width = 6.5, height = 5)
p<-plot_ordination(MM.NoEndo.philr, mm.NoEndo.pcoa, color='SampleSource') 
p + geom_point(size=3, stroke=0) + scale_colour_viridis_d(option="magma", direction = -1) + stat_ellipse(type="t") + stat_ellipse(type="euclid", linetype = 2)
dev.off()

##################################Permanova on philr PCoAs using the vegan command adonis########################################################
#PERMANOVA can only be done on datasets of one sample type, because PERMANOVA can't accept missing data
#So Bat Fly Sex can't be examined in a dataset that contains sample types other than bat flies
metadata.df<-data.frame(sample_data(MM))
mm.dist.mat<-as.matrix(mm.dist)

metadata.25.df<-data.frame(sample_data(MM.25))
mm.25.dist.mat<-as.matrix(mm.25.dist)

metadata.Cave.df<-data.frame(sample_data(MM.Cave))
mm.Cave.dist.mat<-as.matrix(mm.Cave.dist)

metadata.Bat.df<-data.frame(sample_data(MM.Bat))
mm.Bat.dist.mat<-as.matrix(mm.Bat.dist)

metadata.BatFly.df<-data.frame(sample_data(MM.BatFly))
mm.BatFly.dist.mat<-as.matrix(mm.BatFly.dist)

metadata.Trichobius.df<-data.frame(sample_data(MM.Trichobius.philr))
mm.Trichobius.dist.mat<-as.matrix(mm.Trichobius.dist)

metadata.Nycterophilia.df<-data.frame(sample_data(MM.Nycterophilia.philr))
mm.Nycterophilia.dist.mat<-as.matrix(mm.Nycterophilia.dist)

metadata.BatCave.df<-data.frame(sample_data(MM.BatCave.philr))
mm.BatCave.dist.mat<-as.matrix(mm.BatCave.dist)

metadata.BatCaveTrichobius.df<-data.frame(sample_data(MM.BatCaveTrichobius.philr))
mm.BatCaveTrichobius.dist.mat<-as.matrix(mm.BatCaveTrichobius.dist)

metadata.BatCaveNycterophilia.df<-data.frame(sample_data(MM.BatCaveNycterophilia.philr))
mm.BatCaveNycterophilia.dist.mat<-as.matrix(mm.BatCaveNycterophilia.dist)

metadata.NoEndo.df<-data.frame(sample_data(MM.NoEndo.philr))
mm.NoEndo.dist.mat<-as.matrix(mm.NoEndo.dist)

#There are some missing values in the ParasiteSexcolumns that need to be assigned "NA"
#Only do this for the datsets that are entirely composed of bat flies
metadata.BatFly.df$BatFlySex[metadata.BatFly.df$BatFlySex==""]<- NA
metadata.Trichobius.df$BatFlySex[metadata.Trichobius.df$BatFlySex==""]<- NA
metadata.Nycterophilia.df$BatFlySex[metadata.Nycterophilia.df$BatFlySex==""]<- NA

#Remove the NAs from the ParasiteSex Column of the datasets composed entirely of bat flies
metadata.BatFly.df <- metadata.BatFly.df[!is.na(metadata.BatFly.df$BatFlySex),]
table(is.na(metadata.BatFly.df$BatFlySex))

metadata.Trichobius.df <- metadata.Trichobius.df[!is.na(metadata.Trichobius.df$BatFlySex),]
table(is.na(metadata.Trichobius.df$BatFlySex))

metadata.Nycterophilia.df <- metadata.Nycterophilia.df[!is.na(metadata.Nycterophilia.df$BatFlySex),]
table(is.na(metadata.Nycterophilia.df$BatFlySex))

#Now make sure the samples in the distance matrix and in the metadata are the same
#for those where missing values were eliminated
mm.BatFly.dist.mat <- mm.BatFly.dist.mat[row.names(mm.BatFly.dist.mat) %in% row.names(metadata.BatFly.df), colnames(mm.BatFly.dist.mat) %in% row.names(metadata.BatFly.df)]
ncol(mm.BatFly.dist.mat) == nrow(metadata.BatFly.df)
nrow(mm.BatFly.dist.mat) == nrow(metadata.BatFly.df)

mm.Trichobius.dist.mat <- mm.Trichobius.dist.mat[row.names(mm.Trichobius.dist.mat) %in% row.names(metadata.Trichobius.df), colnames(mm.Trichobius.dist.mat) %in% row.names(metadata.Trichobius.df)]
ncol(mm.Trichobius.dist.mat) == nrow(metadata.Trichobius.df)
nrow(mm.Trichobius.dist.mat) == nrow(metadata.Trichobius.df)

mm.Nycterophilia.dist.mat <- mm.Nycterophilia.dist.mat[row.names(mm.Nycterophilia.dist.mat) %in% row.names(metadata.Nycterophilia.df), colnames(mm.Nycterophilia.dist.mat) %in% row.names(metadata.Nycterophilia.df)]
ncol(mm.Nycterophilia.dist.mat) == nrow(metadata.Nycterophilia.df)
nrow(mm.Nycterophilia.dist.mat) == nrow(metadata.Nycterophilia.df)

#MM PERMANOVAs
adonis(mm.dist.mat ~ SampleSource + CollectionLocality, data=metadata.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
#                     Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleSource         3   2444.69  814.90 213.688 0.79575 0.0001 ***
# CollectionLocality   4     36.42    9.10   2.387 0.01185 0.0217 *  
# Residuals          155    591.09    3.81         0.19240           
# Total              162   3072.20                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.dist.mat), metadata.df$SampleSource)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups      3  362.8 120.920 2.8128    999  0.045 *
#   Residuals 159 6835.3  42.989                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.dist.mat), metadata.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      4   6253 1563.27 6.3257    999  0.001 ***
#   Residuals 158  39047  247.13                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.dist.mat ~ BroadLocality, data=metadata.df, permutations=9999, method="euclidean", strata=metadata.df$SampleSource)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# BroadLocality   1    116.14 116.144  6.3258 0.03781 0.0474 *
#   Residuals     161   2956.05  18.361         0.96219         
# Total         162   3072.20                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.dist.mat), metadata.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      1   4783  4783.1 16.796    999  0.001 ***
#   Residuals 161  45849   284.8                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.dist.mat ~ BroadLocality, data=metadata.df, permutations=9999, method="euclidean")
beta <- betadisper(dist(mm.dist.mat), metadata.df$BroadLocality)
permutest(beta)

#MM.25 PERMANOVAs
adonis(mm.25.dist.mat ~ SampleSource + CollectionLocality, data=metadata.25.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleSource         3    2957.5  985.82 229.058 0.80687 0.0001 ***
#   CollectionLocality   4      40.8   10.20   2.371 0.01113 0.0261 *  
#   Residuals          155     667.1    4.30         0.18200           
# Total              162    3665.4                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.25.dist.mat), metadata.25.df$SampleSource)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups      3  537.6 179.194 3.4403    999   0.02 *
#   Residuals 159 8281.7  52.086                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.25.dist.mat), metadata.25.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      4   7984 1995.96 6.0066    999  0.001 ***
#   Residuals 158  52502  332.29                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.25.dist.mat ~ BroadLocality, data=metadata.25.df, permutations=9999, method="euclidean", strata=metadata.df$SampleSource)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# BroadLocality   1     139.8 139.799  6.3842 0.03814 0.0649 .
# Residuals     161    3525.6  21.898         0.96186         
# Total         162    3665.4                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.25.dist.mat), metadata.25.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups      1   6239  6239.3 16.555    999  0.001 ***
#   Residuals 161  60680   376.9                         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#MM.Cave PERMANOVAs
adonis(mm.Cave.dist.mat ~ BroadLocality, data=metadata.Cave.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 119
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# BroadLocality  1    211.52  211.52  1.9179 0.38999    0.1
# Residuals      3    330.86  110.28         0.61001       
# Total          4    542.38                 1.00000

beta <- betadisper(dist(mm.Cave.dist.mat), metadata.Cave.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 119
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1   7.07   7.070 0.1029    119 0.9083
# Residuals  3 206.20  68.735 

#MM.Bat PERMANOVAs
adonis(mm.Bat.dist.mat ~ CollectionLocality, data=metadata.Bat.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# CollectionLocality  4    155.15  38.789  2.8148 0.23331  1e-04 ***
#   Residuals          37    509.87  13.780         0.76669           
# Total              41    665.02                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.Bat.dist.mat), metadata.Bat.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     4  88.02  22.005 1.8025    999  0.139
# Residuals 37 451.70  12.208

adonis(mm.Bat.dist.mat ~ BatSex, data=metadata.Bat.df, permutations=9999, method="euclidean")
beta <- betadisper(dist(mm.Bat.dist.mat), metadata.Bat.df$BatSex)
permutest(beta)

#MM.BatCave PERMANOVAs
adonis(mm.BatCave.dist.mat ~ CollectionLocality, data=metadata.BatCave.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# CollectionLocality  4    133.31  33.329  2.9133 0.21719  1e-04 ***
#   Residuals          42    480.49  11.440         0.78281           
# Total              46    613.80                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatCave.dist.mat ~ SampleSource + CollectionLocality, data=metadata.BatCave.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleSource        1     16.92  16.917  1.4943 0.02756 0.0684 .  
# CollectionLocality  4    132.72  33.180  2.9308 0.21623 0.0001 ***
#   Residuals          41    464.16  11.321         0.75621           
# Total              46    613.80                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatCave.dist.mat), metadata.BatCave.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq     F N.Perm Pr(>F)  
# Groups     4 110.65  27.664 2.173    999  0.082 .
# Residuals 42 534.69  12.731                      
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatCave.dist.mat), metadata.BatCave.df$SampleSource)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1  31.37  31.368 2.0472    999  0.171
# Residuals 45 689.52  15.323

adonis(mm.BatCave.dist.mat ~ CollectionLocality, data=metadata.BatCave.df, permutations=9999, method="euclidean", strata=metadata.BatCave.df$BroadLocality)
beta <- betadisper(dist(mm.BatCave.dist.mat), metadata.BatCave.df$BroadLocality)
permutest(beta)

#MM.BatFly PERMANOVAs
adonis(mm.BatFly.dist.mat ~ SampleSource + CollectionLocality, data=metadata.BatFly.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleSource         1    5875.0  5875.0  738.61 0.86470 0.0001 ***
#   CollectionLocality   4      44.3    11.1    1.39 0.00652 0.2278    
# Residuals          110     875.0     8.0         0.12878           
# Total              115    6794.3                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatFly.dist.mat ~ CollectionLocality, data=metadata.BatFly.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# CollectionLocality   4     776.5 194.123  3.5806 0.11429 0.0043 **
#   Residuals          111    6017.8  54.214         0.88571          
# Total              115    6794.3                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatFly.dist.mat ~ BroadLocality, data=metadata.BatFly.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# BroadLocality   1     266.1 266.149  4.6477 0.03917 0.0285 *
#   Residuals     114    6528.1  57.264         0.96083         
# Total         115    6794.3                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatFly.dist.mat ~ SampleSource, data=metadata.BatFly.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)    
# SampleSource   1    5875.0  5875.0  728.56 0.8647  1e-04 ***
#   Residuals    114     919.3     8.1         0.1353           
# Total        115    6794.3                 1.0000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatFly.dist.mat ~ CollectionLocality, data=metadata.BatFly.df, permutations=9999, method="euclidean", strata=metadata.BatFly.df$SampleSource)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# CollectionLocality   4     776.5 194.123  3.5806 0.11429  0.808
# Residuals          111    6017.8  54.214         0.88571       
# Total              115    6794.3                 1.00000 

adonis(mm.BatFly.dist.mat ~ BroadLocality, data=metadata.BatFly.df, permutations=9999, method="euclidean", strata = metadata.BatFly.df$SampleSource)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# BroadLocality   1     266.1 266.149  4.6477 0.03917 0.4853
# Residuals     114    6528.1  57.264         0.96083       
# Total         115    6794.3                 1.00000 

beta <- betadisper(dist(mm.BatFly.dist.mat), metadata.BatFly.df$SampleSource)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
# Groups      1  1452.1 1452.11 8.8045    999  0.004 **
#   Residuals 114 18801.8  164.93                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatFly.dist.mat), metadata.BatFly.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups      4  22805  5701.3 2.8085    999  0.038 *
#   Residuals 111 225329  2030.0                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatFly.dist.mat), metadata.BatFly.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      1   3091  3090.7 1.3317    999  0.242
# Residuals 114 264570  2320.8

#MM.Trichobius
adonis(mm.Trichobius.dist.mat ~ BroadLocality, data=metadata.Trichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# BroadLocality  1     42.85  42.846 0.83213 0.0149 0.4079
# Residuals     55   2831.91  51.489         0.9851       
# Total         56   2874.76                 1.0000

adonis(mm.Trichobius.dist.mat ~ CollectionLocality, data=metadata.Trichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# CollectionLocality  4    163.77  40.942 0.78531 0.05697 0.5935
# Residuals          52   2710.99  52.134         0.94303       
# Total              56   2874.76                 1.00000

adonis(mm.Trichobius.dist.mat ~ BatFlySex, data=metadata.Trichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# BatFlySex  1    456.04  456.04   10.37 0.15863  6e-04 ***
#   Residuals 55   2418.72   43.98         0.84137           
# Total     56   2874.76                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Trichobius.dist.mat ~ BatSex, data=metadata.Trichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# BatSex     1     25.83  25.831 0.49868 0.00899 0.5953
# Residuals 55   2848.93  51.799         0.99101       
# Total     56   2874.76                 1.00000

adonis(mm.Trichobius.dist.mat ~ AMCC, data=metadata.Trichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# AMCC       1     91.97  91.972  1.8178 0.03199 0.1577
# Residuals 55   2782.79  50.596         0.96801       
# Total     56   2874.76                 1.00000 

adonis(mm.Trichobius.dist.mat ~ BatFlySex, data=metadata.Trichobius.df, permutations=9999, method="euclidean", strata=metadata.Trichobius.df$CollectionLocality)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# BatFlySex  1    456.04  456.04   10.37 0.15863  3e-04 ***
#   Residuals 55   2418.72   43.98         0.84137           
# Total     56   2874.76                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Trichobius.dist.mat ~ BartWol, data=metadata.Trichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# BartWol    3    809.08 269.694  6.9197 0.28144  2e-04 ***
#   Residuals 53   2065.68  38.975         0.71856           
# Total     56   2874.76                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Trichobius.dist.mat ~ BatFlySex + BartWol, data=metadata.Trichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# BatFlySex  1    456.04  456.04 13.1448 0.15863  1e-04 ***
#   BartWol    3    614.68  204.89  5.9058 0.21382  6e-04 ***
#   Residuals 52   1804.05   34.69         0.62755           
# Total     56   2874.76                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.Trichobius.dist.mat), metadata.Trichobius.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1     45   44.69 0.0589    999  0.816
# Residuals 55  41718  758.52

beta <- betadisper(dist(mm.Trichobius.dist.mat), metadata.Trichobius.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     4   1267  316.87 0.3914    999  0.815
# Residuals 52  42097  809.55 

beta <- betadisper(dist(mm.Trichobius.dist.mat), metadata.Trichobius.df$BatFlySex)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq     F N.Perm Pr(>F)   
# Groups     1  5680.3  5680.3 11.43    999  0.002 **
#   Residuals 55 27332.1   496.9                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.Trichobius.dist.mat), metadata.Trichobius.df$BatSex)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1    993  992.82 1.3916    999  0.258
# Residuals 55  39238  713.42

beta <- betadisper(dist(mm.Trichobius.dist.mat), metadata.Trichobius.df$AMCC)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups    40 9720.1 243.002 6.7163    999  0.048 *
#   Residuals 16  578.9  36.181                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.Trichobius.dist.mat), metadata.Trichobius.df$BartWol)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq    F N.Perm Pr(>F)    
# Groups     3  10908  3635.9 7.97    999  0.001 ***
#   Residuals 53  24179   456.2                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#MM.BatCaveTrichobius PERMANOVAs
adonis(mm.BatCaveTrichobius.dist.mat ~ CollectionLocality, data=metadata.BatCaveTrichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# CollectionLocality   4    130.75  32.688  1.6883 0.06386 0.0934 .
# Residuals           99   1916.80  19.362         0.93614         
# Total              103   2047.55                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatCaveTrichobius.dist.mat ~ BroadLocality, data=metadata.BatCaveTrichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# BroadLocality   1     39.89  39.893  2.0268 0.01948 0.1179
# Residuals     102   2007.66  19.683         0.98052       
# Total         103   2047.55                 1.00000 

adonis(mm.BatCaveTrichobius.dist.mat ~ SampleSource, data=metadata.BatCaveTrichobius.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# SampleSource   2    1026.3  513.16  50.751 0.50124  1e-04 ***
#   Residuals    101    1021.2   10.11         0.49876           
# Total        103    2047.5                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatCaveTrichobius.dist.mat), metadata.BatCaveTrichobius.df$SampleSource)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)   
# Groups      2  2018.5 1009.25 7.0978    999  0.003 **
#   Residuals 101 14361.5  142.19                        
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatCaveTrichobius.dist.mat), metadata.BatCaveTrichobius.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups     4  1429.8  357.44 2.0156    999  0.076 .
# Residuals 99 17556.3  177.34                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatCaveTrichobius.dist.mat), metadata.BatCaveTrichobius.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      1   381.8  381.84 1.9971    999  0.166
# Residuals 102 19501.6  191.19 

#MM.Nycterophilia PERMANOVAs
adonis(mm.Nycterophilia.dist.mat ~ BroadLocality, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# BroadLocality  1     15.72 15.7224  2.5241 0.04241 0.0635 .
# Residuals     57    355.04  6.2289         0.95759         
# Total         58    370.77                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Nycterophilia.dist.mat ~ CollectionLocality, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# CollectionLocality  3     37.39 12.4636  2.0562 0.10085 0.0626 .
# Residuals          55    333.38  6.0614         0.89915         
# Total              58    370.77                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Nycterophilia.dist.mat ~ BatFlySex, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# BatFlySex  1     12.80  12.805   2.039 0.03454 0.1039
# Residuals 57    357.96   6.280         0.96546       
# Total     58    370.77                 1.00000

adonis(mm.Nycterophilia.dist.mat ~ BatSex, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# BatSex     1      4.48  4.4817 0.69743 0.01209 0.4761
# Residuals 57    366.29  6.4261         0.98791       
# Total     58    370.77                 1.00000

adonis(mm.Nycterophilia.dist.mat ~ AMCC, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# AMCC       1     10.42 10.4190  1.6481 0.0281 0.1558
# Residuals 57    360.35  6.3219         0.9719       
# Total     58    370.77                 1.0000 

adonis(mm.Nycterophilia.dist.mat ~ BatFlySex, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean", strata=metadata.Nycterophilia.df$CollectionLocality)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# BatFlySex  1     12.80  12.805   2.039 0.03454 0.1378
# Residuals 57    357.96   6.280         0.96546       
# Total     58    370.77                 1.00000 

adonis(mm.Nycterophilia.dist.mat ~ BartWol, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)  
# BartWol    2     33.70  16.851  2.7996 0.0909 0.0283 *
#   Residuals 56    337.07   6.019         0.9091         
# Total     58    370.77                 1.0000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Nycterophilia.dist.mat ~ BatFlySex + BartWol, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# BatFlySex  1     12.80 12.8048  2.1969 0.03454 0.0990 .
# BartWol    2     37.39 18.6932  3.2071 0.10084 0.0195 *
#   Residuals 55    320.58  5.8286         0.86463         
# Total     58    370.77                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Nycterophilia.dist.mat ~ BartWol + CollectionLocality, data=metadata.Nycterophilia.df, permutations=9999)
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# BartWol             2     33.70 16.8506  2.9611 0.09090 0.0263 *
#   CollectionLocality  3     35.45 11.8182  2.0767 0.09562 0.0587 .
# Residuals          53    301.61  5.6908         0.81348         
# Total              58    370.77                 1.00000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.Nycterophilia.dist.mat ~ CollectionLocality + BroadLocality, data=metadata.Nycterophilia.df, permutations=9999)


adonis(mm.Nycterophilia.dist.mat ~ BartWol, data=metadata.Nycterophilia.df, permutations=9999, method="euclidean", strata=metadata.Nycterophilia.df$BatFlySex)
# Blocks:  strata 
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)  
# BartWol    2     33.70  16.851  2.7996 0.0909 0.0331 *
#   Residuals 56    337.07   6.019         0.9091         
# Total     58    370.77                 1.0000         
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.Nycterophilia.dist.mat), metadata.Nycterophilia.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups     1  133.6 133.598 2.8759    999  0.072 .
# Residuals 57 2647.9  46.454                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.Nycterophilia.dist.mat), metadata.Nycterophilia.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     3  152.24  50.746 1.0404    999  0.354
# Residuals 55 2682.64  48.775

beta <- betadisper(dist(mm.Nycterophilia.dist.mat), metadata.Nycterophilia.df$BatFlySex)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq    F N.Perm Pr(>F)
# Groups     1   49.5  49.498 0.93    999  0.317
# Residuals 57 3033.7  53.223

beta <- betadisper(dist(mm.Nycterophilia.dist.mat), metadata.Nycterophilia.df$BatSex)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     1    4.77   4.766 0.0876    999  0.781
# Residuals 57 3102.37  54.428  

beta <- betadisper(dist(mm.Nycterophilia.dist.mat), metadata.Nycterophilia.df$AMCC)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups    34 1503.64  44.225 1.6974    999  0.169
# Residuals 24  625.29  26.054

beta <- betadisper(dist(mm.Nycterophilia.dist.mat), metadata.Nycterophilia.df$BartWol)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups     2  158.51  79.256 1.9039    999  0.143
# Residuals 56 2331.16  41.628 

#MM.BatCaveNycterophilia PERMANOVAs
adonis(mm.BatCaveNycterophilia.dist.mat ~ CollectionLocality, data=metadata.BatCaveNycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# CollectionLocality   4     458.0 114.490  3.3111 0.11593 0.0061 **
#   Residuals          101    3492.3  34.577         0.88407          
# Total              105    3950.3                 1.00000          
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatCaveNycterophilia.dist.mat ~ BroadLocality, data=metadata.BatCaveNycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# BroadLocality   1     354.5  354.55  10.255 0.08975  7e-04 ***
#   Residuals     104    3595.7   34.57         0.91025           
# Total         105    3950.3                 1.00000           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(mm.BatCaveNycterophilia.dist.mat ~ SampleSource, data=metadata.BatCaveNycterophilia.df, permutations=9999, method="euclidean")
# Permutation: free
# Number of permutations: 9999
# 
# Terms added sequentially (first to last)
# 
# Df SumsOfSqs MeanSqs F.Model   R2 Pr(>F)    
# SampleSource   2    3436.7 1718.37  344.66 0.87  1e-04 ***
#   Residuals    103     513.5    4.99         0.13           
# Total        105    3950.3                 1.00           
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatCaveNycterophilia.dist.mat), metadata.BatCaveNycterophilia.df$SampleSource)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)  
# Groups      2  292.4 146.211 2.9433    999  0.053 .
# Residuals 103 5116.6  49.676                       
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

beta <- betadisper(dist(mm.BatCaveNycterophilia.dist.mat), metadata.BatCaveNycterophilia.df$CollectionLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      4   8119  2029.7 1.4951    999  0.223
# Residuals 101 137116  1357.6

beta <- betadisper(dist(mm.BatCaveNycterophilia.dist.mat), metadata.BatCaveNycterophilia.df$BroadLocality)
permutest(beta)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 999
# 
# Response: Distances
# Df Sum Sq Mean Sq      F N.Perm Pr(>F)
# Groups      1   2664  2664.3 1.9154    999   0.18
# Residuals 104 144663  1391.0

#MM.NoEndo PERMANOVAs
adonis(mm.NoEndo.dist.mat ~ SampleSource, data=metadata.NoEndo.df, permutations=9999, method="euclidean")
beta <- betadisper(dist(mm.NoEndo.dist.mat), metadata.NoEndo.df$SampleSource)
permutest(beta)

adonis(mm.NoEndo.dist.mat ~ CollectionLocality, data=metadata.NoEndo.df, permutations=9999, method="euclidean")
beta <- betadisper(dist(mm.NoEndo.dist.mat), metadata.NoEndo.df$CollectionLocality)
permutest(beta)

adonis(mm.NoEndo.dist.mat ~ BroadLocality, data=metadata.NoEndo.df, permutations=9999, method="euclidean")
beta <- betadisper(dist(mm.NoEndo.dist.mat), metadata.NoEndo.df$BroadLocality)
permutest(beta)

################################Barplot of relative abundances###################################################
#https://joey711.github.io/phyloseq/plot_bar-examples.html
#https://github.com/joey711/phyloseq/issues/901

#These datasets are for creating bar plots of the average microbiome community composition for each Sample Type in each Collection Locality
#Add a new Sample Source column that distinguishes male/female bat flies and male/female bats
#Made this in excel, so just remake the MM object. 
metadata<- read.csv("phyloseq/metadata_phyloseq_filtered.csv", header = TRUE, row.names = 1)
metadata.df <- as.data.frame(metadata)
META<-sample_data(metadata)
MM<- phyloseq(OTU,TAX,META,phy_tree)

metadata.25 <- read.csv("phyloseq/metadata_phyloseq_filtered25.csv", header = TRUE, row.names = 1)
metadata.25.df <- as.data.frame(metadata.25)
META.25<-sample_data(metadata.25)
MM.25<- phyloseq(OTU.25,TAX.25,META.25,phy_tree.25)

#Split phyloseq object by CollectionLocality
MM_Gitana <- subset_samples(MM, CollectionLocality=="La Gitana, Sierra Cacachilas, Baja California Sur")
MM_Carmen <- subset_samples(MM, CollectionLocality=="Carmen, Loreto, Baja California Sur")
MM_Chivato <- subset_samples(MM, CollectionLocality=="Chivato, Sierra Cacachilas, Baja California Sur")
MM_Fabrica <- subset_samples(MM, CollectionLocality=="Cueva de la Fabrica, Coquimatlan, Colima, MX")
MM_Panchito <- subset_samples(MM, CollectionLocality=="Isla Panchito, near Chamela, Jalisco, MX")

MM.25_Gitana <- subset_samples(MM.25, CollectionLocality=="La Gitana, Sierra Cacachilas, Baja California Sur")
MM.25_Carmen <- subset_samples(MM.25, CollectionLocality=="Carmen, Loreto, Baja California Sur")
MM.25_Chivato <- subset_samples(MM.25, CollectionLocality=="Chivato, Sierra Cacachilas, Baja California Sur")
MM.25_Fabrica <- subset_samples(MM.25, CollectionLocality=="Cueva de la Fabrica, Coquimatlan, Colima, MX")
MM.25_Panchito <- subset_samples(MM.25, CollectionLocality=="Isla Panchito, near Chamela, Jalisco, MX")

#Within each Locality's phyloseq object, merge samples by the SampleSource column
MM_Gitana <- merge_samples(MM_Gitana, "SampleSource")
MM_Carmen <- merge_samples(MM_Carmen, "SampleSource")
MM_Chivato <- merge_samples(MM_Chivato, "SampleSource")
MM_Fabrica <- merge_samples(MM_Fabrica, "SampleSource")
MM_Panchito <- merge_samples(MM_Panchito, "SampleSource")

MM.25_Gitana <- merge_samples(MM.25_Gitana, "SampleSource")
MM.25_Carmen <- merge_samples(MM.25_Carmen, "SampleSource")
MM.25_Chivato <- merge_samples(MM.25_Chivato, "SampleSource")
MM.25_Fabrica <- merge_samples(MM.25_Fabrica, "SampleSource")
MM.25_Panchito <- merge_samples(MM.25_Panchito, "SampleSource")

#In the broader dataset, create a new phyloseq object with samples merged by SampleSourceSex,
#which differentiates between the sexes where applicable
MM_glom_SampSex <- merge_samples(MM, "SampleSourceSex")
MM.25_glom_SampSex <- merge_samples(MM.25, "SampleSourceSex")

#Replace the metadata, which is lost during merge and is important to plots
sample_data(MM_Gitana)$CollectionLocality <- "Gitana"
sample_data(MM_Gitana)$SampleSource <- row.names(sample_data(MM_Gitana))
sample_names(MM_Gitana) <- paste0(sample_names(MM_Gitana), sample_data(MM_Gitana)$CollectionLocality)

sample_data(MM_Carmen)$CollectionLocality <- "Carmen"
sample_data(MM_Carmen)$SampleSource <- row.names(sample_data(MM_Carmen))
sample_names(MM_Carmen) <- paste0(sample_names(MM_Carmen), sample_data(MM_Carmen)$CollectionLocality)

sample_data(MM_Chivato)$CollectionLocality <- "Chivato"
sample_data(MM_Chivato)$SampleSource <- row.names(sample_data(MM_Chivato))
sample_names(MM_Chivato) <- paste0(sample_names(MM_Chivato), sample_data(MM_Chivato)$CollectionLocality)

sample_data(MM_Fabrica)$CollectionLocality <- "Fabrica"
sample_data(MM_Fabrica)$SampleSource <- row.names(sample_data(MM_Fabrica))
sample_names(MM_Fabrica) <-paste0(sample_names(MM_Fabrica), sample_data(MM_Fabrica)$CollectionLocality)

sample_data(MM_Panchito)$CollectionLocality <- "Panchito"
sample_data(MM_Panchito)$SampleSource <- row.names(sample_data(MM_Panchito))
sample_names(MM_Panchito) <- paste0(sample_names(MM_Panchito), sample_data(MM_Panchito)$CollectionLocality)

sample_data(MM.25_Gitana)$CollectionLocality <- "Gitana"
sample_data(MM.25_Gitana)$SampleSource <- row.names(sample_data(MM.25_Gitana))
sample_names(MM.25_Gitana) <- paste0(sample_names(MM.25_Gitana), sample_data(MM.25_Gitana)$CollectionLocality)

sample_data(MM.25_Carmen)$CollectionLocality <- "Carmen"
sample_data(MM.25_Carmen)$SampleSource <- row.names(sample_data(MM.25_Carmen))
sample_names(MM.25_Carmen) <- paste0(sample_names(MM.25_Carmen), sample_data(MM.25_Carmen)$CollectionLocality)

sample_data(MM.25_Chivato)$CollectionLocality <- "Chivato"
sample_data(MM.25_Chivato)$SampleSource <- row.names(sample_data(MM.25_Chivato))
sample_names(MM.25_Chivato) <- paste0(sample_names(MM.25_Chivato), sample_data(MM.25_Chivato)$CollectionLocality)

sample_data(MM.25_Fabrica)$CollectionLocality <- "Fabrica"
sample_data(MM.25_Fabrica)$SampleSource <- row.names(sample_data(MM.25_Fabrica))
sample_names(MM.25_Fabrica) <-paste0(sample_names(MM.25_Fabrica), sample_data(MM.25_Fabrica)$CollectionLocality)

sample_data(MM.25_Panchito)$CollectionLocality <- "Panchito"
sample_data(MM.25_Panchito)$SampleSource <- row.names(sample_data(MM.25_Panchito))
sample_names(MM.25_Panchito) <- paste0(sample_names(MM.25_Panchito), sample_data(MM.25_Panchito)$CollectionLocality)

sample_data(MM_glom_SampSex)$SampleSourceSex <- row.names(sample_data(MM_glom_SampSex))

sample_data(MM.25_glom_SampSex)$SampleSourceSex <- row.names(sample_data(MM.25_glom_SampSex))

#Merge the phyloseq objects
MM_glom <- merge_phyloseq(MM_Gitana, MM_Chivato, MM_Carmen, MM_Panchito, MM_Fabrica)

MM.25_glom <- merge_phyloseq(MM.25_Gitana, MM.25_Chivato, MM.25_Carmen, MM.25_Panchito, MM.25_Fabrica)

#Convert aggregated phyloseq objects and the glom_SampSex objects into dataframes
MM_glom.df <-psmelt(MM_glom)
MM.25_glom.df <-psmelt(MM.25_glom)
MM_glom_SampSex.df <-psmelt(MM_glom_SampSex)
MM.25_glom_SampSex.df <-psmelt(MM.25_glom_SampSex)

#Convert taxonomic ranking to character
MM_glom.df$Genus <-as.character(MM_glom.df$Genus)
MM_glom.df$Family <-as.character(MM_glom.df$Family)
MM_glom.df$Class <- as.character(MM_glom.df$Class)
MM_glom.df$Phylum <- as.character(MM_glom.df$Phylum)

MM.25_glom.df$Genus <-as.character(MM.25_glom.df$Genus)
MM.25_glom.df$Family <-as.character(MM.25_glom.df$Family)
MM.25_glom.df$Class <- as.character(MM.25_glom.df$Class)
MM.25_glom.df$Phylum <- as.character(MM.25_glom.df$Phylum)

MM_glom_SampSex.df$Genus <-as.character(MM_glom_SampSex.df$Genus)
MM_glom_SampSex.df$Family <-as.character(MM_glom_SampSex.df$Family)
MM_glom_SampSex.df$Class <- as.character(MM_glom_SampSex.df$Class)
MM_glom_SampSex.df$Phylum <- as.character(MM_glom_SampSex.df$Phylum)

MM.25_glom_SampSex.df$Genus <-as.character(MM.25_glom_SampSex.df$Genus)
MM.25_glom_SampSex.df$Family <-as.character(MM.25_glom_SampSex.df$Family)
MM.25_glom_SampSex.df$Class <- as.character(MM.25_glom_SampSex.df$Class)
MM.25_glom_SampSex.df$Phylum <- as.character(MM.25_glom_SampSex.df$Phylum)

#Let's see how many Classes each sample type has within the combined phyloseq object
unique(MM_glom.df$Class)
unique(MM.25_glom.df$Class)
unique(MM_glom_SampSex.df$Class)
unique(MM.25_glom_SampSex.df$Class)

#The bat swabs have a lot of unique OTUs that are very low in abundance, however, these can be clustered by Order
#and increase the relative abundance of the more common Orders present
#To do this, we first create a Taxon column that contains the Order of the OTU for Cave and Bat samples,
#but the Genus of the OTU for BatFly Samples
#Second, we create a new column where we will sum the abundances for each unique Order or Genus within a 
#sample type and collection locality and remove duplicate rows
#Third, we calculate the relative abundances of Order and Genus within each unique Sample Type/Collection Locality
#combination.

#Create a new Taxon column that has the Bacterial Genus associated with bat flies, but the Bacterial Order for Bats and Cave
MM_glom.df$Taxon <- MM_glom.df$Genus
MM_glom.df$Taxon <-ifelse(MM_glom.df$SampleSource == "Cave", MM_glom.df$Class, MM_glom.df$Taxon)
MM_glom.df$Taxon <- ifelse(MM_glom.df$SampleSource == "Bat", MM_glom.df$Class, MM_glom.df$Taxon)

MM.25_glom.df$Taxon <- MM.25_glom.df$Genus
MM.25_glom.df$Taxon <- ifelse(MM.25_glom.df$SampleSource == "Cave", MM.25_glom.df$Class, MM.25_glom.df$Taxon)
MM.25_glom.df$Taxon <- ifelse(MM.25_glom.df$SampleSource == "Bat", MM.25_glom.df$Class, MM.25_glom.df$Taxon)

MM_glom_SampSex.df$Taxon <- MM_glom_SampSex.df$Genus
MM_glom_SampSex.df$Taxon <- ifelse(MM_glom_SampSex.df$SampleSourceSex == "Cave", MM_glom_SampSex.df$Class, MM_glom_SampSex.df$Taxon)
MM_glom_SampSex.df$Taxon <- ifelse(MM_glom_SampSex.df$SampleSourceSex == "Batfemale", MM_glom_SampSex.df$Class, MM_glom_SampSex.df$Taxon)
MM_glom_SampSex.df$Taxon <- ifelse(MM_glom_SampSex.df$SampleSourceSex == "Batmale", MM_glom_SampSex.df$Class, MM_glom_SampSex.df$Taxon)

MM.25_glom_SampSex.df$Taxon <- MM.25_glom_SampSex.df$Genus
MM.25_glom_SampSex.df$Taxon <- ifelse(MM.25_glom_SampSex.df$SampleSourceSex == "Cave", MM.25_glom_SampSex.df$Class, MM.25_glom_SampSex.df$Taxon)
MM.25_glom_SampSex.df$Taxon <- ifelse(MM.25_glom_SampSex.df$SampleSourceSex == "Batfemale", MM.25_glom_SampSex.df$Class, MM.25_glom_SampSex.df$Taxon)
MM.25_glom_SampSex.df$Taxon <- ifelse(MM.25_glom_SampSex.df$SampleSourceSex == "Batmale", MM.25_glom_SampSex.df$Class, MM.25_glom_SampSex.df$Taxon)

#Combine Orders and Genera that are the same within a SampleType and CollectionLocality combo and remove 
#the duplicate rows
MM_glom.df <- group_by(MM_glom.df, CollectionLocality, SampleSource, Taxon)
MM_glom_relabund <- summarize(MM_glom.df, count = sum(Abundance))

MM.25_glom.df <- group_by(MM.25_glom.df, CollectionLocality, SampleSource, Taxon)
MM.25_glom_relabund <- summarize(MM.25_glom.df, count = sum(Abundance))

MM_glom_SampSex.df <- group_by(MM_glom_SampSex.df, SampleSourceSex, Taxon)
MM_glom_SampSex_relabund <- summarize(MM_glom_SampSex.df, count = sum(Abundance))

MM.25_glom_SampSex.df <- group_by(MM.25_glom_SampSex.df, SampleSourceSex, Taxon)
MM.25_glom_SampSex_relabund <- summarize(MM.25_glom_SampSex.df, count = sum(Abundance))

#Remove rows that have a count of 0
MM_glom_relabund <- filter(MM_glom_relabund, count > 0)
MM.25_glom_relabund <- filter(MM.25_glom_relabund, count > 0)
MM_glom_SampSex_relabund <- filter(MM_glom_SampSex_relabund, count > 0)
MM.25_glom_SampSex_relabund <- filter(MM.25_glom_SampSex_relabund, count > 0)

#Convert counts to relative abundance within a SampleType-CollectionLocality Combo or by SampleSourceSex
MM_glom_relabund.df<- data.frame(MM_glom_relabund)
MM_glom_relabund.df<- group_by(MM_glom_relabund.df, CollectionLocality, SampleSource)
MM_glom_relabund.df <- mutate(MM_glom_relabund.df, RelAbund = count/sum(count))
MM_glom_relabund.df[1:20,1:5]

MM.25_glom_relabund.df<- data.frame(MM.25_glom_relabund)
MM.25_glom_relabund.df<- group_by(MM.25_glom_relabund.df, CollectionLocality, SampleSource)
MM.25_glom_relabund.df <- mutate(MM.25_glom_relabund.df, RelAbund = count/sum(count))
MM.25_glom_relabund.df[1:20,1:5]

MM_glom_SampSex_relabund.df<- data.frame(MM_glom_SampSex_relabund)
MM_glom_SampSex_relabund.df<- group_by(MM_glom_SampSex_relabund.df, SampleSourceSex)
MM_glom_SampSex_relabund.df <- mutate(MM_glom_SampSex_relabund.df, RelAbund = count/sum(count))
MM_glom_SampSex_relabund.df[1:20,1:4]

MM.25_glom_SampSex_relabund.df<- data.frame(MM.25_glom_SampSex_relabund)
MM.25_glom_SampSex_relabund.df<- group_by(MM.25_glom_SampSex_relabund.df, SampleSourceSex)
MM.25_glom_SampSex_relabund.df <- mutate(MM.25_glom_SampSex_relabund.df, RelAbund = count/sum(count))
MM.25_glom_SampSex_relabund.df[1:20,1:4]

#Aggegate all taxa with <1% abundance into an "Low Abundance" category
MM_glom_relabund.df$Taxon[MM_glom_relabund.df$RelAbund < 0.01] <- "LowAbundance"
MM.25_glom_relabund.df$Taxon[MM.25_glom_relabund.df$RelAbund < 0.01] <- "LowAbundance"
MM_glom_SampSex_relabund.df$Taxon[MM_glom_SampSex_relabund.df$RelAbund < 0.01] <- "LowAbundance"
MM.25_glom_SampSex_relabund.df$Taxon[MM.25_glom_SampSex_relabund.df$RelAbund < 0.01] <- "LowAbundance"

#Change NAs, uncultured, or uncultured bacterium in the Taxon column to "Unidentified"
MM_glom_relabund.df$Taxon[is.na(MM_glom_relabund.df$Taxon) | 
                            MM_glom_relabund.df$Taxon == "uncultured" |
                            MM_glom_relabund.df$Taxon == "uncultured bacterium"]<- "Unidentified"
MM.25_glom_relabund.df$Taxon[is.na(MM.25_glom_relabund.df$Taxon) | 
                               MM.25_glom_relabund.df$Taxon == "uncultured" |
                               MM.25_glom_relabund.df$Taxon == "uncultured bacterium"]<- "Unidentified"
MM_glom_SampSex_relabund.df$Taxon[is.na(MM_glom_SampSex_relabund.df$Taxon) | 
                            MM_glom_SampSex_relabund.df$Taxon == "uncultured" |
                            MM_glom_SampSex_relabund.df$Taxon == "uncultured bacterium"]<- "Unidentified"
MM.25_glom_SampSex_relabund.df$Taxon[is.na(MM.25_glom_SampSex_relabund.df$Taxon) | 
                               MM.25_glom_SampSex_relabund.df$Taxon == "uncultured" |
                               MM.25_glom_SampSex_relabund.df$Taxon == "uncultured bacterium"]<- "Unidentified"

#Get list of Ranks
ord<-MM_glom_relabund.df[with(MM_glom_relabund.df, order(RelAbund, decreasing=TRUE)),]
ord25<-MM.25_glom_relabund.df[with(MM.25_glom_relabund.df, order(RelAbund, decreasing=TRUE)),]
ord_SampSex<-MM_glom_SampSex_relabund.df[with(MM_glom_SampSex_relabund.df, order(RelAbund, decreasing=TRUE)),]
ord25_SampSex<-MM.25_glom_SampSex_relabund.df[with(MM.25_glom_SampSex_relabund.df, order(RelAbund, decreasing=TRUE)),]

unique(ord$Taxon)
unique(ord25$Taxon)
unique(ord_SampSex$Taxon)
unique(ord25_SampSex$Taxon)

#Write the order files to csv
write.csv(ord, "RelativeAbundancePlots/SampOrder.csv", row.names = FALSE)
write.csv(ord25, "RelativeAbundancePlots/25_SampOrder.csv", row.names = FALSE)
write.csv(ord_SampSex, "RelativeAbundancePlots/SampOrder_SampSex.csv", row.names = FALSE)
write.csv(ord25_SampSex, "RelativeAbundancePlots/25_SampOrder_SampSex.csv", row.names = FALSE)

#Change order of Taxon according to abundance
MM_glom_relabund.df$Taxon<- factor(MM_glom_relabund.df$Taxon, levels = c("endosymbionts3","Arsenophonus",
                                                                         "Bacilli","Gammaproteobacteria",
                                                                         "Bacteroidia","Actinobacteria",
                                                                         "Bartonella","Unidentified",
                                                                         "Alphaproteobacteria","Wolbachia",
                                                                         "Rhodothermia","Clostridia",
                                                                         "Deltaproteobacteria","Limnochordia",
                                                                         "S0134 terrestrial group",
                                                                         "Phycisphaerae","Sphingoaurantiacus",
                                                                         "Chloroflexia","Planctomycetacia",
                                                                         "Longimicrobia","Thermoleophilia",
                                                                         "Acidimicrobiia","Mollicutes",
                                                                         "Gemmatimonadetes","Microvirga",
                                                                         "Bosea","Ac37b","Verrucomicrobiae",
                                                                         "Subgroup 6","F0332","uncultured Bacillales bacterium",
                                                                         "PMMR1","Blastocatellia (Subgroup 4)","LowAbundance"))

MM.25_glom_relabund.df$Taxon<- factor(MM.25_glom_relabund.df$Taxon, levels = c("endosymbionts3","Bacilli","Arsenophonus",
                                                                               "Gammaproteobacteria","Bacteroidia",
                                                                               "Actinobacteria","Bartonella",
                                                                               "Alphaproteobacteria","Wolbachia",
                                                                               "Unidentified","Rhodothermia",
                                                                               "Clostridia","Deltaproteobacteria",
                                                                               "Limnochordia","Chloroflexia",
                                                                               "Phycisphaerae","Longimicrobia",
                                                                               "Sphingoaurantiacus","S0134 terrestrial group",
                                                                               "Acidimicrobiia","Mollicutes",
                                                                               "Thermoleophilia","Planctomycetacia",
                                                                               "Bosea","Subgroup 6","Microvirga",
                                                                               "Gemmatimonadetes","Ac37b","F0332",
                                                                               "bacterium LY17","uncultured Bacillales bacterium",
                                                                               "Oxyphotobacteria","PMMR1","LowAbundance"))

MM_glom_SampSex_relabund.df$Taxon<- factor(MM_glom_SampSex_relabund.df$Taxon, levels = c("endosymbionts3","Arsenophonus",
                                                                                         "Gammaproteobacteria","Bacilli",
                                                                                         "Actinobacteria","Bartonella",
                                                                                         "Wolbachia","Bacteroidia",
                                                                                         "Alphaproteobacteria",
                                                                                         "Unidentified","Clostridia",
                                                                                         "Rhodothermia","Deltaproteobacteria",
                                                                                         "Chloroflexia","Mollicutes",
                                                                                         "Thermoleophilia","Limnochordia",
                                                                                         "Planctomycetacia","S0134 terrestrial group",
                                                                                         "Ac37b","Phycisphaerae","LowAbundance"))

MM.25_glom_SampSex_relabund.df$Taxon<- factor(MM.25_glom_SampSex_relabund.df$Taxon, levels = c("endosymbionts3","Arsenophonus",
                                                                                               "Gammaproteobacteria","Bacilli",
                                                                                               "Actinobacteria","Bartonella",
                                                                                               "Wolbachia","Bacteroidia",
                                                                                               "Alphaproteobacteria","Clostridia",
                                                                                               "Rhodothermia","Unidentified",
                                                                                               "Deltaproteobacteria","Chloroflexia",
                                                                                               "Limnochordia","Mollicutes","Ac37b",
                                                                                               "Longimicrobia","Phycisphaerae",
                                                                                               "Thermoleophilia","LowAbundance"))

#Change order of CollectionLocality according to moving North to South
MM_glom_relabund.df$CollectionLocality<- factor(MM_glom_relabund.df$CollectionLocality, levels = c("Carmen", "Gitana", "Chivato", "Panchito", "Fabrica"))
MM.25_glom_relabund.df$CollectionLocality<- factor(MM.25_glom_relabund.df$CollectionLocality, levels = c("Carmen", "Gitana", "Chivato", "Panchito", "Fabrica"))

#Plot the average composition of the microbiome by Sample Type using Collection Locality as a facet
p_glom<- ggplot(data=MM_glom_relabund.df, aes(x=SampleSource, y=RelAbund, fill=Taxon))
p_glom + facet_wrap(~CollectionLocality) + geom_bar(stat="identity", position="fill") + 
  scale_fill_viridis_d(option= "magma",direction = -1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(limits=c("Cave", "Bat","Trichobius","Nycterophilia"))
ggsave("RelativeAbundancePlots/RelAbund_Sample_ColLoc.pdf", width=10, height=5, units="in")

p25_glom<- ggplot(data=MM.25_glom_relabund.df, aes(x=SampleSource, y=RelAbund, fill=Taxon))
p25_glom + facet_wrap(~CollectionLocality) + geom_bar(stat="identity", position="fill") + 
  scale_fill_viridis_d(option= "magma",direction = -1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_x_discrete(limits=c("Cave", "Bat","Trichobius","Nycterophilia"))
ggsave("RelativeAbundancePlots/25_RelAbund_Sample_ColLoc.pdf", width=10, height=5, units="in")

#Plot the average composition of the microbiome by collection locality, using Sample Type as a facet
#Change order of SampleSourceSex according to moving North to South
MM_glom_relabund.df$SampleSource<- factor(MM_glom_relabund.df$SampleSource, levels = c("Cave","Bat","Trichobius","Nycterophilia"))
MM.25_glom_relabund.df$SampleSource<- factor(MM.25_glom_relabund.df$SampleSource, levels = c("Cave","Bat","Trichobius","Nycterophilia"))

p_glom<- ggplot(data=MM_glom_relabund.df, aes(x=CollectionLocality, y=RelAbund, fill=Taxon))
p_glom + facet_wrap(~SampleSource) + geom_bar(stat="identity", position="fill") + 
  scale_fill_viridis_d(option= "magma",direction = -1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("RelativeAbundancePlots/RelAbund_ColLoc_Sample.pdf", width=10, height=5, units="in")

p25_glom<- ggplot(data=MM.25_glom_relabund.df, aes(x=CollectionLocality, y=RelAbund, fill=Taxon))
p25_glom + facet_wrap(~SampleSource) + geom_bar(stat="identity", position="fill") + 
  scale_fill_viridis_d(option= "magma",direction = -1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("RelativeAbundancePlots/25_RelAbund_ColLoc_Sample.pdf", width=10, height=5, units="in")

#Plot average microbiome compostion by Sample type including sex of bat and bat flies, but not collection locality
p_glom_SampSex<- ggplot(data=MM_glom_SampSex_relabund.df, aes(x=SampleSourceSex, y=RelAbund, fill=Taxon))
p_glom_SampSex + geom_bar(stat="identity", position="fill") + scale_fill_viridis_d(option= "magma",direction = -1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_x_discrete(limits=c("Cave", "Batfemale", "Batmale","Trichobiusfemale","Trichobiusmale","Nycterophiliafemale","Nycterophiliamale"))
ggsave("RelativeAbundancePlots/RelAbund_SampSex.pdf", width=8, height=6, units="in")

p25_glom_SampSex<- ggplot(data=MM.25_glom_SampSex_relabund.df, aes(x=SampleSourceSex, y=RelAbund, fill=Taxon))
p25_glom_SampSex + geom_bar(stat="identity", position="fill") + scale_fill_viridis_d(option= "magma",direction = -1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  scale_x_discrete(limits=c("Cave", "Batfemale", "Batmale","Trichobiusfemale","Trichobiusmale","Nycterophiliafemale","Nycterophiliamale"))
ggsave("RelativeAbundancePlots/25_RelAbund_SampSex.pdf", width=7, height=5, units="in")

################Random Forest Classification##################################################################################
#https://rpubs.com/phamdinhkhanh/389752
#http://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
#https://topepo.github.io/caret/index.html
#https://machinelearningmastery.com/how-to-estimate-model-accuracy-in-r-using-the-caret-package/

#We will run Random Forests on two datasets
#On the laptop:
#Create a dataset where the microbiome data is summarized as the PCoA coordinates
#Use Random Forests to attempt to classify the Sample Type, Collection Locality, whether a bat fly was collected
#from a specific bat (AMCC), and the sex of the bat fly (within datasets that look at the bat fly only)

#On Cuvier:
#Create a dataset where the microbiome members are each a column used to build the model
#Use Random Forests to perform the same classifications on this data

#LAPTOP RANDOM FORESTS
#We will use the PCoA axis loadings of philr-transformed data to summarize the microbiome (Abundance and OTU)
#So that the size of our data frames is smaller.

#Read in the dataset to be used for these random forests, or re-make it using the below code
mm.rf.df <- read.csv("RandomForest/MM_RF.csv", header=TRUE)
mm.25.rf.df <- read.csv("RandomForest/25_MM_RF.csv", header=TRUE)
mm.rf.Bat_Fly.df <- read.csv("RandomForest/MM_Bat_Fly_RF.csv", header=TRUE)
mm.Trichobius.rf.df<- read.csv("RandomForest/MM_Trichobius_RF.csv", header=TRUE)
mm.Nycterophilia.rf.df<- read.csv("RandomForest/MM_Nycterophilia_RF.csv", header=TRUE)
mm.Bat.rf.df <- read.csv("RandomForest/MM_Bat_RF.csv", header=TRUE)
mm.NoEndo.rf.df <- read.csv("RandomFores/MM_NoEndo_RF.csv", header=TRUE)

head(mm.pcoa)
head(mm.25.pcoa)
head(mm.BatFly.pcoa)
head(mm.Trichobius.pcoa)
head(mm.Nycterophilia.pcoa)
head(mm.NoEndo.pcoa)

#Add the complete data PCoA Axis 1 and Axis 2 to each sample in the sample metadata
head(metadata)
mm.rf.df<- merge(metadata, mm.pcoa$vectors[,c(1,2)], by="row.names")
colnames(mm.rf.df)[colnames(mm.rf.df) == "Row.names"] <- "Sample"
head(mm.rf.df)

mm.25.rf.df<- merge(metadata, mm.25.pcoa$vectors[,c(1,2)], by="row.names")
colnames(mm.25.rf.df)[colnames(mm.25.rf.df) == "Row.names"] <- "Sample"
head(mm.25.rf.df)

write.csv(mm.rf.df, file="RandomForest/MM_RF.csv", row.names=FALSE)
write.csv(mm.25.rf.df, file="RandomForest/25_MM_RF.csv", row.names=FALSE)

#Add the Bat + Fly PCoA Axis 1 and 2 to each sample in the Bat+Fly metadata
metadata_Bat_Fly <- rbind(metadata_Bat, metadata_BF)
mm.Bat_Fly.pcoa.vectors <- rbind(mm.BatFly.pcoa$vectors[,c(1,2)], mm.Bat.pcoa$vectors[,c(1,2)])
mm.rf.Bat_Fly.df <- merge(metadata_Bat_Fly, mm.Bat_Fly.pcoa.vectors, by="row.names")
colnames(mm.rf.Bat_Fly.df)[colnames(mm.rf.Bat_Fly.df) == "Row.names"] <- "Sample"

write.csv(mm.rf.Bat_Fly.df, file="RandomForest/MM_Bat_Fly_RF.csv", row.names=FALSE)

#Add the Trichobius PCoA data to the Trichobius metadata
metadata_Trich <- sample_data(MM.Trichobius.philr)
mm.Trichobius.rf.df <- merge(metadata_Trich, mm.Trichobius.pcoa$vectors[,c(1,2)], by="row.names")
colnames(mm.Trichobius.rf.df)[colnames(mm.Trichobius.rf.df) == "Row.names"] <- "Sample"
head(mm.Trichobius.rf.df)

write.csv(mm.Trichobius.rf.df, file="RandomForest/MM_Trichobius_RF.csv", row.names=FALSE)

#Add the Nycterophilia PCoA data to the Nycterophilia metadata
metadata_Nyct<- sample_data(MM.Nycterophilia.philr)

mm.Nycterophilia.rf.df <- merge(metadata_Nyct, mm.Nycterophilia.pcoa$vectors[,c(1,2)], by="row.names")
colnames(mm.Nycterophilia.rf.df)[colnames(mm.Nycterophilia.rf.df) == "Row.names"] <- "Sample"
head(mm.Nycterophilia.rf.df)

write.csv(mm.Nycterophilia.rf.df, file="RandomForest/MM_Nycterophilia_RF.csv", row.names=FALSE)

#Add the Bat PCoA data to the Bat metadata
mm.Bat.rf.df <- merge(metadata_Bat, mm.Bat.pcoa$vectors[,c(1,2)], by="row.names")
colnames(mm.Bat.rf.df)[colnames(mm.Bat.rf.df) == "Row.names"] <- "Sample"
head(mm.Bat.rf.df)

write.csv(mm.Bat.rf.df, file="RandomForest/MM_Bat_RF.csv", row.names=FALSE)

#Add the NoEndo PCoA data to the NoEndo metadata
metadata_NoEndo<- sample_data(MM.NoEndo.philr)

mm.NoEndo.rf.df <- merge(metadata_NoEndo, mm.NoEndo.pcoa$vectors[,c(1,2)], by="row.names")
colnames(mm.NoEndo.rf.df)[colnames(mm.NoEndo.rf.df) == "Row.names"] <- "Sample"
head(mm.NoEndo.rf.df)

write.csv(mm.NoEndo.rf.df, file="RandomForest/MM_NoEndo_RF.csv", row.names=FALSE)

#Classify by Sample Type
#First split the data into test data (20%) and training data (80%)
#Make sure that both data sets contain Cave Swabs, Bat Swabs (M and F), and Bat Flies of both species (M and F)
train_list <- createDataPartition(mm.rf.df$SampleSource, p=0.8, list=FALSE, times=1)

train_SampleType <- mm.rf.df[train_list,]
test_SampleType <- mm.rf.df[-train_list,]

train_list <- createDataPartition(mm.NoEndo.rf.df$SampleSource, p=0.8, list=FALSE, times=1)

train_SampleTypeNoEndo <- mm.NoEndo.rf.df[train_list,]
test_SampleTypeNoEndo <- mm.NoEndo.rf.df[-train_list,]

#Classify by Collection Locality
train_list <- createDataPartition(mm.rf.df$CollectionLocality, p=0.8, list=FALSE, times=1)

train_ColLoc <- mm.rf.df[train_list,]
test_ColLoc <- mm.rf.df[-train_list,]

train_list <- createDataPartition(mm.NoEndo.rf.df$CollectionLocality, p=0.8, list=FALSE, times=1)

train_ColLocNoEndo <- mm.NoEndo.rf.df[train_list,]
test_ColLocNoEndo <- mm.NoEndo.rf.df[-train_list,]

#Classify the Bat+Fly data by AMCC
#Remove AMCCs that are only represented by a single individual
ToRemove<-c("249595","249636","249642","249685","249687","249696","249715","249723","249726","249735","249746","249773","249785","255984","255986","255999","256000","256012","256013","256014","256015","256016","256048","256049","256052","256054","256061")
mm.rf.Bat_Fly.df<- mm.rf.Bat_Fly.df[!(mm.rf.Bat_Fly.df$AMCC %in% ToRemove),]
mm.rf.Bat_Fly.df<- mm.rf.Bat_Fly.df[!(is.na(mm.rf.Bat_Fly.df$AMCC)),]

train_list <- createDataPartition(factor(mm.rf.Bat_Fly.df$AMCC), p=0.7, list=FALSE, times=1)

train_AMCC <- mm.rf.Bat_Fly.df[train_list,]
test_AMCC <- mm.rf.Bat_Fly.df[-train_list,]

#Classify by BatFly Sex
#Remove NAs
mm.Trichobius.rf.df <- mm.Trichobius.rf.df[!(is.na(mm.Trichobius.rf.df$BatFlySex)),]

train_list <- createDataPartition(mm.Trichobius.rf.df$BatFlySex, p=0.8, list=FALSE, times=1)

train_FlySexTri <- mm.Trichobius.rf.df[train_list,]
test_FlySexTri <- mm.Trichobius.rf.df[-train_list,]

mm.Nycterophilia.rf.df <- mm.Nycterophilia.rf.df[!(is.na(mm.Nycterophilia.rf.df$BatFlySex)),]

train_list <- createDataPartition(mm.Nycterophilia.rf.df$BatFlySex, p=0.8, list=FALSE, times=1)

train_FlySexNyct <- mm.Nycterophilia.rf.df[train_list,]
test_FlySexNyct <- mm.Nycterophilia.rf.df[-train_list,]

#Classify by Wolbachia and Bartonella infection status
train_list <- createDataPartition(mm.Trichobius.rf.df$BartWol, p=0.8, list=FALSE, times=1)

train_FlyInfTri <- mm.Trichobius.rf.df[train_list,]
test_FlyInfTri <- mm.Trichobius.rf.df[-train_list,]

train_list <- createDataPartition(mm.Nycterophilia.rf.df$BartWol, p=0.8, list=FALSE, times=1)

train_FlyInfNyct <- mm.Nycterophilia.rf.df[train_list,]
test_FlyInfNyct <- mm.Nycterophilia.rf.df[-train_list,]

#Classify by Collection Locality
train_list <- createDataPartition(mm.Trichobius.rf.df$CollectionLocality, p=0.8, list=FALSE, times=1)

train_FlyColLocTri <- mm.Trichobius.rf.df[train_list,]
test_FlyColLocTri <- mm.Trichobius.rf.df[-train_list,]

train_list <- createDataPartition(mm.Nycterophilia.rf.df$CollectionLocality, p=0.8, list=FALSE, times=1)

train_FlyColLocNyct <- mm.Nycterophilia.rf.df[train_list,]
test_FlyColLocNyct <- mm.Nycterophilia.rf.df[-train_list,]

#Classify by Bat CollectionLocality
#Remove NAs
mm.Bat.rf.df <- mm.Bat.rf.df[!(is.na(mm.Bat.rf.df$CollectionLocality)),]

train_list <- createDataPartition(mm.Bat.rf.df$CollectionLocality, p=0.8, list=FALSE, times=1)

train_Bat <- mm.Bat.rf.df[train_list,]
test_Bat <- mm.Bat.rf.df[-train_list,]

#Classify by Bat Sex
mm.Bat.rf.df <- mm.Bat.rf.df[!(is.na(mm.Bat.rf.df$BatSex)),]

train_list <- createDataPartition(mm.Bat.rf.df$BatSex, p=0.8, list=FALSE, times=1)

train_BatSex <- mm.Bat.rf.df[train_list,]
test_BatSex <- mm.Bat.rf.df[-train_list,]

#Remove Columns that are not used
train_SampleType<- select(train_SampleType, SampleSource, Axis.1, Axis.2, CollectionLocality)
head(train_SampleType)

test_SampleType<- select(test_SampleType, SampleSource, Axis.1, Axis.2, CollectionLocality)
head(test_SampleType)

train_SampleTypeNoEndo<- select(train_SampleTypeNoEndo, SampleSource, Axis.1, Axis.2, CollectionLocality)
head(train_SampleTypeNoEndo)

test_SampleTypeNoEndo<- select(test_SampleTypeNoEndo, SampleSource, Axis.1, Axis.2, CollectionLocality)
head(test_SampleTypeNoEndo)

train_ColLoc<- select(train_ColLoc, CollectionLocality, Axis.1, Axis.2, SampleSource)
head(train_ColLoc)

test_ColLoc<- select(test_ColLoc, CollectionLocality, Axis.1, Axis.2, SampleSource)
head(test_ColLoc)

train_ColLocNoEndo<- select(train_ColLocNoEndo, CollectionLocality, Axis.1, Axis.2, SampleSource)
head(train_ColLocNoEndo)

test_ColLocNoEndo<- select(test_ColLocNoEndo, CollectionLocality, Axis.1, Axis.2, SampleSource)
head(test_ColLocNoEndo)

train_FlySexTri <- select(train_FlySexTri, BatFlySex, Axis.1, Axis.2, CollectionLocality)
head(train_FlySexTri)

test_FlySexTri <- select(test_FlySexTri, BatFlySex, Axis.1, Axis.2, CollectionLocality)
head(test_FlySexTri)

train_FlySexNyct <- select(train_FlySexNyct, BatFlySex, Axis.1, Axis.2, CollectionLocality)
head(train_FlySexNyct)

test_FlySexNyct <- select(test_FlySexNyct, BatFlySex, Axis.1, Axis.2, CollectionLocality)
head(test_FlySexNyct)

train_FlyInfTri <- select(train_FlyInfTri, BartWol, Axis.1, Axis.2, CollectionLocality, BatFlySex)
head(train_FlyInfTri)

test_FlyInfTri <- select(test_FlyInfTri, BartWol, Axis.1, Axis.2, CollectionLocality, BatFlySex)
head(test_FlyInfTri)

train_FlyInfNyct <- select(train_FlyInfNyct, BartWol, Axis.1, Axis.2, CollectionLocality, BatFlySex)
head(train_FlyInfNyct)

test_FlyInfNyct <- select(test_FlyInfNyct, BartWol, Axis.1, Axis.2, CollectionLocality, BatFlySex)
head(test_FlyInfNyct)

train_FlyColLocTri <- select(train_FlyColLocTri, CollectionLocality, Axis.1, Axis.2, BatFlySex, BartWol)
head(train_FlyColLocTri)

test_FlyColLocTri <- select(test_FlyColLocTri, CollectionLocality, Axis.1, Axis.2, BatFlySex, BartWol)
head(test_FlyColLocTri)

train_FlyColLocNyct <- select(train_FlyColLocNyct, CollectionLocality, Axis.1, Axis.2, BatFlySex, BartWol)
head(train_FlyColLocNyct)

test_FlyColLocNyct <- select(test_FlyColLocNyct, CollectionLocality, Axis.1, Axis.2, BatFlySex, BartWol)
head(test_FlyColLocNyct)

train_Bat <- select(train_Bat, CollectionLocality, BatSex, Axis.1, Axis.2)
head(train_Bat)

test_Bat <- select(test_Bat, CollectionLocality, BatSex, Axis.1, Axis.2)
head(test_Bat)

train_BatSex <- select(train_BatSex, BatSex, CollectionLocality, Axis.1, Axis.2)
head(train_Bat)

test_BatSex <- select(test_BatSex, BatSex, CollectionLocality, Axis.1, Axis.2)
head(test_Bat)

#Convert the factors (categorical variables) to binary or numerical values 
str(train_SampleType)
str(train_ColLoc)
str(train_AMCC)
train_AMCC$AMCC<-factor(train_AMCC$AMCC)
test_AMCC$AMCC<-factor(test_AMCC$AMCC)
str(train_FlySexTri)
table(train_FlySexTri$BatFlySex)
str(train_FlySexNyct)
str(train_FlyInfTri)
str(train_FlyInfNyct)
str(train_FlyColLocTri)
str(train_FlyColLocNyct)
table(train_FlySexNyct$BatFlySex)
str(train_Bat)
str(train_BatSex)
str(train_SampleTypeNoEndo)

#Set the resampling method for the random forst to 10-fold CV
#Only use up-samping for the SampleType dataset and the FlySexNyct dataset
#Only use Class Probs for the Fly Sex datasets, which have two classes
fit_control_up <- trainControl(method="repeatedcv", repeats=5, number=10, sampling="up")
fit_control <- trainControl(method="repeatedcv", repeats=5, number=10)
fit_control_2_up <- trainControl(method="repeatedcv", repeats=5, number=10, sampling="up", classProbs=TRUE, summaryFunction = twoClassSummary)
fit_control_2 <- trainControl(method="repeatedcv", repeats=5, number=10, classProbs=TRUE, summaryFunction = twoClassSummary)

#Fit a random forest model
set.seed(825)
rf_SampleType<- train(SampleSource ~ ., data = train_SampleType, method = "ranger", metric="Kappa", trControl = fit_control_up)
rf_SampleType

set.seed(825)
rf_SampleTypeNoEndo<- train(SampleSource ~ ., data = train_SampleTypeNoEndo, method = "ranger", metric="Kappa", trControl = fit_control_up)
rf_SampleTypeNoEndo

set.seed(825)
rf_ColLoc<- train(CollectionLocality ~ ., data=train_ColLoc, method="ranger", metric="Kappa", trControl = fit_control)
rf_ColLoc

set.seed(825)
rf_ColLocNoEndo<- train(CollectionLocality ~ ., data = train_ColLocNoEndo, method = "ranger", metric="Kappa", trControl = fit_control)
rf_ColLocNoEndo

set.seed(825)
rf_AMCC <- train(AMCC ~ ., data=train_AMCC, method="ranger", metric="Kappa", trControl = fit_control)
rf_AMCC

set.seed(825)
rf_FlySexTri <- train(BatFlySex ~ ., data=train_FlySexTri, method="ranger", metric ="ROC", trControl=fit_control_2_up)
rf_FlySexTri

set.seed(825)
rf_FlySexNyct <- train(BatFlySex ~ ., data=train_FlySexNyct, method="ranger", metric="ROC", trControl=fit_control_2_up)
rf_FlySexNyct

set.seed(825)
rf_FlyInfTri <- train(BartWol ~ ., data=train_FlyInfTri, method="ranger", metric ="Kappa", trControl=fit_control)
rf_FlyInfTri

set.seed(825)
rf_FlyInfNyct <- train(BartWol ~ ., data=train_FlyInfNyct, method="ranger", metric="Kappa", trControl=fit_control)
rf_FlyInfNyct

set.seed(825)
rf_FlyColLocTri <- train(CollectionLocality ~ ., data=train_FlyColLocTri, method="ranger", metric ="Kappa", trControl=fit_control)
rf_FlyColLocTri

set.seed(825)
rf_FlyColLocNyct <- train(CollectionLocality ~ ., data=train_FlyColLocNyct, method="ranger", metric="Kappa", trControl=fit_control_up)
rf_FlyColLocNyct

set.seed(825)
rf_Bat <- train(CollectionLocality ~ ., data=train_Bat, method="ranger", metric="Kappa", trControl=fit_control)
rf_Bat

set.seed(825)
rf_BatSex <- train(BatSex ~ ., data=train_BatSex, method="ranger", metric="ROC", trControl=fit_control_2_up)
rf_BatSex

#Predict the outcome using the trained model on the test data
rf_SampleType_Test<- predict(rf_SampleType, test_SampleType)
confusionMatrix(rf_SampleType_Test, as.factor(test_SampleType$SampleSource))

rf_SampleTypeNoEndo_Test<- predict(rf_SampleTypeNoEndo, test_SampleTypeNoEndo)
confusionMatrix(rf_SampleTypeNoEndo_Test, as.factor(test_SampleTypeNoEndo$SampleSource))

rf_ColLoc_Test<- predict(rf_ColLoc, test_ColLoc)
confusionMatrix(rf_ColLoc_Test, as.factor(test_ColLoc$CollectionLocality))

rf_ColLocNoEndo_Test<- predict(rf_ColLocNoEndo, test_ColLocNoEndo)
confusionMatrix(rf_ColLocNoEndo_Test, as.factor(test_ColLocNoEndo$CollectionLocality))

# rf_AMCC_Test<- predict(rf_AMCC, test_AMCC)
# confusionMatrix(rf_AMCC_Test, as.factor(test_AMCC$AMCC))

rf_FlySexTri_Test<- predict(rf_FlySexTri, test_FlySexTri)
confusionMatrix(rf_FlySexTri_Test, as.factor(test_FlySexTri$BatFlySex))

rf_FlySexNyct_Test<- predict(rf_FlySexNyct, test_FlySexNyct)
confusionMatrix(rf_FlySexNyct_Test, as.factor(test_FlySexNyct$BatFlySex))

rf_FlyInfTri_Test<- predict(rf_FlyInfTri, test_FlyInfTri)
confusionMatrix(rf_FlyInfTri_Test, as.factor(test_FlyInfTri$BartWol))

rf_FlyInfNyct_Test<- predict(rf_FlyInfNyct, test_FlyInfNyct)
confusionMatrix(rf_FlyInfNyct_Test, as.factor(test_FlyInfNyct$BartWol))

rf_FlyColLocTri_Test<- predict(rf_FlyColLocTri, test_FlyColLocTri)
confusionMatrix(rf_FlyColLocTri_Test, as.factor(test_FlyColLocTri$CollectionLocality))

rf_FlyColLocNyct_Test<- predict(rf_FlyColLocNyct, test_FlyColLocNyct)
confusionMatrix(rf_FlyColLocNyct_Test, as.factor(test_FlyColLocNyct$CollectionLocality))

rf_Bat_Test<- predict(rf_Bat, test_Bat)
confusionMatrix(rf_Bat_Test, as.factor(test_Bat$CollectionLocality))

rf_BatSex_Test<- predict(rf_BatSex, test_BatSex)
confusionMatrix(rf_BatSex_Test, as.factor(test_BatSex$BatSex))

