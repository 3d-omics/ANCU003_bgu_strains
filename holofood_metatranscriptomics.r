setwd("/Users/anttonalberdi/github/holofood_metatranscriptomics/")

library(tidyverse)
library(readr)
library(stringr)
library(data.table)
library(ggplot2)
library(hilldiv)

####
####
# LOAD AND PREPARE DATA
####
####

#Gene counts
gene_counts <- read_tsv("data/gene_counts.txt")
colnames(gene_counts) <- gsub(" Read Count","",colnames(gene_counts))
nsamples=ncol(gene_counts)-1

#Gene annotations
gene_annotations <- read_tsv("data/annotations_caecum.tsv") %>%
  rename("Contig" = 1, "MAG" = 2) %>%
  select(-scaffold, -gene_position, -start_position, -end_position, -strandedness)
gene_annotations$Genus <- unlist(lapply(strsplit(gene_annotations$bin_taxonomy, ";"),function(l) l[6]))
gene_annotations$Family <- unlist(lapply(strsplit(gene_annotations$bin_taxonomy, ";"),function(l) l[5]))
gene_annotations$Order <- unlist(lapply(strsplit(gene_annotations$bin_taxonomy, ";"),function(l) l[4]))
gene_annotations$Class <- unlist(lapply(strsplit(gene_annotations$bin_taxonomy, ";"),function(l) l[3]))
gene_annotations$Phylum <- unlist(lapply(strsplit(gene_annotations$bin_taxonomy, ";"),function(l) l[2]))

#Merge gene count and annotation tables
gene_counts_annot <- inner_join(gene_counts, gene_annotations, by="Contig")

#Save annotated counts table
write_delim(gene_counts_annot,"data/gene_counts_annot.tsv")

#MAG stats
MAG_stats <- read_tsv("data/mag_stats.tsv")

#MAG counts (and genome-size corrected counts)
MAG_counts <- read_csv("data/mag_counts.csv")
MAG_counts_weighed <- MAG_counts
MAG_stats$correction_factor=median(MAG_stats$Size)/MAG_stats$Size
MAG_counts_weighed[,-1]=round(sweep(MAG_counts_weighed[,-1], MARGIN=1, MAG_stats$correction_factor, `*`),0)

#MAG taxonomy
MAG_taxonomy <- read_tsv("data/mag_taxonomy.tsv")
MAG_taxonomy$Domain <- gsub("d__","",unlist(lapply(strsplit(MAG_taxonomy$classification, ";"),function(l) l[1])))
MAG_taxonomy$Phylum <- gsub("p__","",unlist(lapply(strsplit(MAG_taxonomy$classification, ";"),function(l) l[2])))
MAG_taxonomy$Class <- gsub("c__","",unlist(lapply(strsplit(MAG_taxonomy$classification, ";"),function(l) l[3])))
MAG_taxonomy$Order <- gsub("o__","",unlist(lapply(strsplit(MAG_taxonomy$classification, ";"),function(l) l[4])))
MAG_taxonomy$Family <- gsub("f__","",unlist(lapply(strsplit(MAG_taxonomy$classification, ";"),function(l) l[5])))
#Sort
MAG_taxonomy <- MAG_taxonomy[with(MAG_taxonomy, order(Domain,Phylum,Class,Order,Family)), ]

#ENAtoHF mapping file
ENAtoHF <- read_csv("data/ENAtoHF.csv")

#Kegg mapping file
kegg_steps <- read_tsv("data/module_step_form.tsv")
kegg_module_hierarchy <- read_tsv("data/module_hierarchy.tsv")

####
#Tidy up tables
####

#Create intersect

#Filter samples from MAG counts
intersect <- intersect(pull(ENAtoHF,ena_id),colnames(MAG_counts))
MAG_counts <- MAG_counts[,c("Genome",intersect)]

intersect <- intersect(pull(ENAtoHF,ena_id),colnames(MAG_counts_weighed))
MAG_counts_weighed <- MAG_counts_weighed[,c("Genome",intersect)]

#Rename samples from MAG counts
colnames(MAG_counts) <- with(ENAtoHF, holofood_id[match(colnames(MAG_counts), ena_id)])
colnames(MAG_counts)[1] <- "Genome"

colnames(MAG_counts_weighed) <- with(ENAtoHF, holofood_id[match(colnames(MAG_counts_weighed), ena_id)])
colnames(MAG_counts_weighed)[1] <- "Genome"

####
####
# GENERAL STATISTICS
####
####

####
#Compute aggregated values
####
gene_counts_mag <- aggregate(gene_counts_annot[2:(nsamples+1)],by=list(gene_counts_annot$MAG),FUN=sum)
gene_counts_genus <- aggregate(gene_counts_annot[2:(nsamples+1)],by=list(gene_counts_annot$Genus),FUN=sum)
gene_counts_family <- aggregate(gene_counts_annot[2:(nsamples+1)],by=list(gene_counts_annot$Family),FUN=sum)
gene_counts_order <- aggregate(gene_counts_annot[2:(nsamples+1)],by=list(gene_counts_annot$Order),FUN=sum)
gene_counts_phylum <- aggregate(gene_counts_annot[2:(nsamples+1)],by=list(gene_counts_annot$Phylum),FUN=sum)
gene_counts_kegg <- aggregate(gene_counts_annot[2:(nsamples+1)],by=list(gene_counts_annot$kegg_id),FUN=sum)
gene_counts_cazy <- aggregate(gene_counts_annot[2:(nsamples+1)],by=list(gene_counts_annot$cazy_hits),FUN=sum)

####
# General stats
####

#MAGs in catalogue
length(unique(gene_annotations$MAG))
#Genes in catalogue
nrow(gene_annotations)
#Average genes per MAG
nrow(gene_annotations)/length(unique(gene_annotations$MAG))
#Genes with expression data
nrow(gene_counts_annot[rowSums(gene_counts_annot[,c(2:71)]) > 0,])
nrow(gene_counts_annot[rowSums(gene_counts_annot[,c(2:71)]) > 0,])/nrow(gene_annotations)*100

####
# General expression plots
####

#General expression plot per phylum
gene_counts_phylum_rel <- gene_counts_phylum
gene_counts_phylum_rel[,-1] <- tss(gene_counts_phylum_rel[,-1])
gene_counts_phylum_melt <- melt(gene_counts_phylum_rel)
colnames(gene_counts_phylum_melt) <- c("Taxon","Sample","Counts")
ggplot(gene_counts_phylum_melt,aes(x=Sample,y=Counts,fill=Taxon)) +
  geom_bar(stat = "identity") +
  theme_grey(base_size=8)+
  theme(legend.position="bottom", strip.text.y = element_text(angle = 0),axis.text.x=element_blank())

####
# MAG abundance vs. gene abundance plot
####

#Subset and compute relative values
MAG_counts_sub_tss <- as.data.frame(MAG_counts_weighed[,colnames(MAG_counts_weighed) %in% gsub("Group.1","Genome",colnames(gene_counts_mag))])
MAG_counts_sub_tss[,-1] <- tss(MAG_counts_sub_tss[,-1])
gene_counts_mag_tss <- gene_counts_mag
gene_counts_mag_tss[,-1] <- tss(gene_counts_mag_tss[,-1])
colnames(gene_counts_mag_tss)[1] <- "Genome"

#Ensure same samples are considered
MAG_counts_sub_tss <- MAG_counts_sub_tss[, intersect(colnames(MAG_counts_sub_tss),colnames(gene_counts_mag_tss))]
gene_counts_mag_tss <- gene_counts_mag_tss[, intersect(colnames(MAG_counts_sub_tss),colnames(gene_counts_mag_tss))]

#Melt and merge
MetaG_melt <- melt(MAG_counts_sub_tss)
MetaG_melt <- MetaG_melt[order(MetaG_melt[,2],MetaG_melt[,1]),]
MetaT_melt <- melt(gene_counts_mag_tss)
MetaT_melt <- MetaT_melt[order(MetaT_melt[,2],MetaT_melt[,1]),]
MetaG_MetaT <- cbind(MetaG_melt,MetaT_melt[,3])
colnames(MetaG_MetaT) <- c("MAG","Sample","MetaG","MetaT")

#Append taxonomy
MetaG_MetaT <- merge(MetaG_MetaT,MAG_taxonomy[,c("user_genome","Phylum")],by.x="MAG",by.y="user_genome")

#Overall correlation
summary(lm(MetaG~MetaT, data=MetaG_MetaT))
cor.test(~MetaG+MetaT, data=MetaG_MetaT)

#Per phylum linear regressions
phyla <- unique(MetaG_MetaT$Phylum)
slopes <- c()
for(p in phyla){
MetaG_MetaT_sub <- MetaG_MetaT[MetaG_MetaT$Phylum == p,]
slope <- lm(MetaT~MetaG, data=MetaG_MetaT_sub)[1]$coefficients[2]
slopes <- c(slopes,slope)
}
names(slopes) <- phyla
slopes

#Plot (linear regressions plotted by Phylum)
ggplot(MetaG_MetaT,aes(x=MetaG,y=MetaT,color=Phylum)) +
  geom_point(alpha = 0.2) +
  xlim(0, 0.2) +
  ylim(0, 0.2) +
  geom_smooth(method='lm') +
  theme_minimal()

####
# Module data (overall)
####

#Merge count and KEGG KO data
gene_counts_ko <- merge(gene_counts,gene_annotations[,c("Contig","kegg_id")],by="Contig")
# Remove genes without annotation
gene_counts_ko <- gene_counts_ko[!is.na(gene_counts_ko$kegg_id),]
# Aggregate genes by KO
gene_counts_ko_agg <- aggregate(gene_counts_ko[,c(2:(ncol(gene_counts_ko)-1))],by=list(gene_counts_ko[,ncol(gene_counts_ko)]),FUN=sum)
colnames(gene_counts_ko_agg)[1] <- "ko"

#Merge count and KEGG Module data
gene_counts_module <- merge(gene_counts_ko_agg,kegg_steps[,c("ko","module")],by="ko")
gene_counts_module <- gene_counts_module[!is.na(gene_counts_module$module),]
gene_counts_module_agg <- aggregate(gene_counts_module[,c(2:(ncol(gene_counts_module)-1))],by=list(gene_counts_module[,ncol(gene_counts_module)]),FUN=sum)
colnames(gene_counts_module_agg)[1] <- "module"

rownames(gene_counts_module_agg) <- gene_counts_module_agg[,1]
gene_counts_module_agg <- gene_counts_module_agg[,-1]
heatmap(as.matrix(t(gene_counts_module_agg)))

####
# Module data (per MAG)
####

#Merge count and KEGG KO data
gene_counts_mag_ko <- merge(gene_counts,gene_annotations[,c("Contig","MAG","kegg_id")],by="Contig")
gene_counts_mag_ko <- gene_counts_mag_ko[!is.na(gene_counts_mag_ko$kegg_id),]
#Aggregate by MAG and KO (it takes a while)
gene_counts_mag_ko_agg <- aggregate(gene_counts_mag_ko[,c(2:(ncol(gene_counts_mag_ko)-2))],by=list(gene_counts_mag_ko[,"MAG"],gene_counts_mag_ko[,"kegg_id"]),FUN=sum)
#Get sample means
gene_counts_mag_ko_agg$mean <- rowSums(gene_counts_mag_ko_agg[,c(3:ncol(gene_counts_mag_ko_agg))])
#Aggregate by module
gene_counts_mag_ko_module <- merge(gene_counts_mag_ko_agg,kegg_steps[,c("ko","module")],by.x="Group.2",by.y="ko")
gene_counts_mag_ko_module_agg <- aggregate(gene_counts_mag_ko_module[,(ncol(gene_counts_mag_ko_module)-1)],by=list(gene_counts_mag_ko_module[,"Group.1"],gene_counts_mag_ko_module[,"module"]),FUN=sum)
colnames(gene_counts_mag_ko_module_agg) <- c("MAG","Module","Expression")
#Add category information
gene_counts_mag_ko_module_agg_cat <- merge(gene_counts_mag_ko_module_agg,kegg_module_hierarchy[,c(1:3)],by="Module")
gene_counts_mag_ko_module_agg_cat <- merge(gene_counts_mag_ko_module_agg_cat,MAG_taxonomy[,c("user_genome","Phylum","Order")],by.x="MAG",by.y="user_genome")

#Sort MAGs modules and categories
gene_counts_mag_ko_module_agg_cat$MAG <- as.factor(gene_counts_mag_ko_module_agg_cat$MAG)
gene_counts_mag_ko_module_agg_cat$MAG <- factor(gene_counts_mag_ko_module_agg_cat$MAG, levels = MAG_taxonomy$user_genome)
gene_counts_mag_ko_module_agg_cat$Module <- as.factor(gene_counts_mag_ko_module_agg_cat$Module)
gene_counts_mag_ko_module_agg_cat$Module <- factor(gene_counts_mag_ko_module_agg_cat$Module, levels = kegg_module_hierarchy$Module)
gene_counts_mag_ko_module_agg_cat$Category <- as.factor(gene_counts_mag_ko_module_agg_cat$Category)
gene_counts_mag_ko_module_agg_cat$Category <- factor(gene_counts_mag_ko_module_agg_cat$Category, levels = unique(kegg_module_hierarchy$Category))
gene_counts_mag_ko_module_agg_cat$Domain <- as.factor(gene_counts_mag_ko_module_agg_cat$Domain)
gene_counts_mag_ko_module_agg_cat$Domain <- factor(gene_counts_mag_ko_module_agg_cat$Domain, levels = unique(kegg_module_hierarchy$Domain))

#Modify to logarithmic scale
gene_counts_mag_ko_module_agg_cat_log <- gene_counts_mag_ko_module_agg_cat
gene_counts_mag_ko_module_agg_cat_log$Expression <- log(gene_counts_mag_ko_module_agg_cat_log$Expression)
gene_counts_mag_ko_module_agg_cat_log[gene_counts_mag_ko_module_agg_cat_log == "-Inf"] <- 0

#Plot (MAG vs Module)
ggplot(gene_counts_mag_ko_module_agg_cat_log, aes(Module, MAG, fill=Expression)) +
  geom_tile() +
  scale_fill_gradientn(colours=rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da","#f4f4f4")), na.value="#000000") +
  facet_grid(Phylum ~ Domain, switch = "both", scales = "free", space = "free") +
  theme(strip.text.y.left = element_text(angle = 0),panel.spacing = unit(0.1, "lines"),axis.text.x=element_blank(), axis.text.y=element_blank())

#Plot (Order vs Category)
ggplot(gene_counts_mag_ko_module_agg_cat_log, aes(Category, Order, fill=Expression)) +
    geom_tile() +
    scale_fill_gradientn(colours=rev(c("#d53e4f", "#f46d43", "#fdae61", "#fee08b", "#e6f598", "#abdda4", "#ddf1da","#f4f4f4")), na.value="#000000") +
    facet_grid(Phylum ~ Domain, switch = "both", scales = "free", space = "free") +
    theme(strip.text.y.left = element_text(angle = 0),panel.spacing = unit(0.1, "lines"),axis.text.x=element_blank(), axis.text.y=element_blank())

####
# Pathway analysis
####

library(pathview)

#Test
gene_counts_sub <- gene_counts[gene_counts$Contig %like% "ERR4836918_bin.11",c("Contig","CC17.07F1a")]
gene_annotations_sub <- gene_annotations[which(gene_annotations$MAG=="ERR4836918_bin.11"),]
gene_merged_sub <- merge(gene_counts_sub,gene_annotations_sub[,c("Contig","kegg_id")],by="Contig")
gene_merged_kegg <- gene_merged_sub[!is.na(gene_merged_sub$kegg_id),]
gene_merged_kegg_agg <- aggregate(gene_merged_kegg[,2],by=list(gene_merged_kegg[,3]),FUN=sum)
kegg_vector <- as.numeric(gene_merged_kegg_agg[,2])
names(kegg_vector) <- gene_merged_kegg_agg[,1]

#Percentage of KEGG annotated genes
nrow(gene_merged_kegg)/nrow(gene_merged_sub)*100

#Percentage of KEGG annotated reads
sum(gene_merged_kegg[,2])/sum(gene_merged_sub[,2])*100

#Plot
pathview(gene.data = kegg_vector/mean(kegg_vector), pathway.id = "00010", species="ko")

####
# Individual genome expression profile
####

library(DESeq2)
library(OmicCircos)
#https://bioconductor.org/packages/devel/bioc/vignettes/OmicCircos/inst/doc/OmicCircos_vignette.pdf

#Subset target genome
genome="ERR4836918_bin.11"
gene_annotations2 <- read_tsv("data/gene_annotations.tsv.gz") %>%
  rename("Contig" = 1, "MAG" = 2)
ERR4836918 <- gene_annotations2[gene_annotations2$MAG == genome,]
print(as_tibble(ERR4836918),n=50)
counts <- read_tsv("data/counts.tsv.gz")
ERR4836918_counts <- inner_join(counts, ERR4836918, by="Contig")
ERR4836918_counts[2:11] <- rlog(as.matrix(ERR4836918_counts[2:11]))

#Number of contigs
length(unique(ERR4836918$scaffold))
#Number of genes
length(unique(ERR4836918$Contig))

#Format genome (segments)
segments <- aggregate(ERR4836918$end_position,by=list(ERR4836918$scaffold),FUN=max)
segments$chromStart <- rep(0,nrow(segments))
colnames(segments)[1:2] <- c("chrom","chromEnd")
segments <- segments[,c("chrom","chromStart","chromEnd")]

#Log expression per sample
mapping <- cbind(ERR4836918_counts$scaffold,
ERR4836918_counts$start_position,
ERR4836918_counts$Contig,
ERR4836918_counts[2:11])
colnames(mapping)[1:3] <- c("chr","po","NAME")

#Mean expression
mapping2 <- as.data.frame(cbind(ERR4836918_counts$scaffold,
ERR4836918_counts$start_position,
ERR4836918_counts$Contig,
rowMeans(ERR4836918_counts[2:11])))
colnames(mapping2) <- c("chr","po","NAME","MeanExp")

segment_db <- segAnglePo(seg.dat=as.data.frame(ERR4836918[,c("scaffold","start_position","end_position","Contig","rank")]),seg=unique(ERR4836918$scaffold))

colors   <- rainbow(seg.num, alpha=0.5);

pdf("ERR4836918.pdf",width=8,height=8)
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")
#Plot contig references
circos(R=340, type="chr", cir=segment_db, col=colors, print.chr.lab=TRUE, W=4, scale=TRUE);
#Plot gene expression boxplots
circos(R=260, cir=segment_db, W=80, mapping=mapping, col.v=4, type="box", B=TRUE, col=colors[1], lwd=0.1, scale=TRUE);
#Plot average gene expression
circos(R=220, cir=segment_db, W=40, mapping=mapping2, col.v=4, type="b", B=TRUE, col=colors[1], lwd=0.1, scale=TRUE);
dev.off()
