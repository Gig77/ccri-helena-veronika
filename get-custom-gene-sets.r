library(biomaRt)
human = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl") # GRCh37, v75
mouse = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="feb2014.archive.ensembl.org", dataset="mmusculus_gene_ensembl") # GRCm38, v75

d <- read.delim("/mnt/projects/chrisi/results/deseq/zuber+strobl-expressed-vs-notexpressed.deseq2.chipseq-annotated.fuka-ross-boer.tsv", stringsAsFactors = F)

topN <- 3000
minFC <- 1
maxFDR <- 0.01

gmt <- data.frame(name=character(0), descr=character(0), genes=character(0), stringsAsFactors = F)

#############################################################################################################
# Chrisi inducible E/R overexpression
#############################################################################################################
chrisi <- d[!is.na(d$hgnc_symbol) & !is.na(d$padj) & d$padj <= maxFDR & abs(d$log2FoldChange) >= minFC, c("hgnc_symbol", "log2FoldChange", "padj")]
chrisi <- chrisi[order(chrisi$padj),]
chrisi <- chrisi[!duplicated(chrisi$hgnc_symbol),]
chrisi.up <- chrisi[chrisi$log2FoldChange >= 0,]
chrisi.dn <- chrisi[chrisi$log2FoldChange <= 0,]

gmt[nrow(gmt)+1,] <- c("PORTSMOUTH_2014_ER_OVEREXPRESSING_UP", "", paste(sort(chrisi.up$hgnc_symbol[1:min(nrow(chrisi.up), topN)]), collapse="\t"))
gmt[nrow(gmt)+1,] <- c("PORTSMOUTH_2014_ER_OVEREXPRESSING_DN", "", paste(sort(chrisi.dn$hgnc_symbol[1:min(nrow(chrisi.dn), topN)]), collapse="\t"))

#############################################################################################################
# Fuka et al. (2011) E/R knockdown
#############################################################################################################
fuka <- d[!is.na(d$hgnc_symbol) & !is.na(d$fuka.Padj) & d$fuka.Padj <= maxFDR  & abs(d$fuka.logFC) >= minFC, c("hgnc_symbol", "fuka.logFC", "fuka.Padj")]
fuka <- fuka[order(fuka$fuka.Padj),]
fuka <- fuka[!duplicated(fuka$hgnc_symbol),]
fuka.up <- fuka[fuka$fuka.logFC >= 0,]
fuka.dn <- fuka[fuka$fuka.logFC <= 0,]

gmt[nrow(gmt)+1,] <- c("FUKA_2011_ER_KNOCKDOWN_UP", "", paste(sort(fuka.up$hgnc_symbol[1:min(nrow(fuka.up), topN)]), collapse="\t"))
gmt[nrow(gmt)+1,] <- c("FUKA_2011_ER_KNOCKDOWN_DN", "", paste(sort(fuka.dn$hgnc_symbol[1:min(nrow(fuka.dn), topN)]), collapse="\t"))

#############################################################################################################
# Boer TEL-AML vs. non-TEL-AML
#############################################################################################################

boer1 <- d[!is.na(d$hgnc_symbol) & !is.na(d$boer.adjPval.TAvs.mean.noTall) & d$boer.adjPval.TAvs.mean.noTall <= maxFDR  & abs(d$boer.TAvs.mean.noTall) >= minFC, c("hgnc_symbol", "boer.TAvs.mean.noTall", "boer.adjPval.TAvs.mean.noTall")]
boer1 <- boer1[order(boer1$boer.adjPval.TAvs.mean.noTall),]
boer1 <- boer1[!duplicated(boer1$hgnc_symbol),]
boer1.up <- boer1[boer1$boer.TAvs.mean.noTall >= 0,]
boer1.dn <- boer1[boer1$boer.TAvs.mean.noTall <= 0,]

gmt[nrow(gmt)+1,] <- c("BOER_2009_ER_VS_NOT_TALL_UP", "", paste(sort(boer1.up$hgnc_symbol[1:min(nrow(boer1.up), topN)]), collapse="\t"))
gmt[nrow(gmt)+1,] <- c("BOER_2009_ER_VS_NOT_TALL_DN", "", paste(sort(boer1.dn$hgnc_symbol[1:min(nrow(boer1.dn), topN)]), collapse="\t"))

#############################################################################################################
# Boer TEL-AML vs. rest
#############################################################################################################

boer2 <- d[!is.na(d$hgnc_symbol) & !is.na(d$boer2009.adjP.TA_vs_rest) & d$boer2009.adjP.TA_vs_rest <= maxFDR  & abs(d$boer2009.TA_vs_rest) >= minFC, c("hgnc_symbol", "boer2009.TA_vs_rest", "boer2009.adjP.TA_vs_rest")]
boer2 <- boer2[order(boer2$boer2009.adjP.TA_vs_rest),]
boer2 <- boer2[!duplicated(boer2$hgnc_symbol),]
boer2.up <- boer2[boer2$boer2009.TA_vs_rest >= 0,]
boer2.dn <- boer2[boer2$boer2009.TA_vs_rest <= 0,]

gmt[nrow(gmt)+1,] <- c("BOER_2009_ER_VS_REST_UP", "", paste(sort(boer2.up$hgnc_symbol[1:min(nrow(boer2.up), topN)]), collapse="\t"))
gmt[nrow(gmt)+1,] <- c("BOER_2009_ER_VS_REST_DN", "", paste(sort(boer2.dn$hgnc_symbol[1:min(nrow(boer2.dn), topN)]), collapse="\t"))

#############################################################################################################
# Ross
#############################################################################################################

ross <- d[!is.na(d$hgnc_symbol) & !is.na(d$ross.adjPval.TAvs.mean_noTALL) & d$ross.adjPval.TAvs.mean_noTALL <= maxFDR  & abs(d$ross.TAvs.mean_noTALL) >= minFC, c("hgnc_symbol", "ross.TAvs.mean_noTALL", "ross.adjPval.TAvs.mean_noTALL")]
ross <- ross[order(ross$ross.adjPval.TAvs.mean_noTALL),]
ross <- ross[!duplicated(ross$hgnc_symbol),]
ross.up <- ross[ross$ross.TAvs.mean_noTALL >= 0,]
ross.dn <- ross[ross$ross.TAvs.mean_noTALL <= 0,]

gmt[nrow(gmt)+1,] <- c("ROSS_2003_ER_VS_NOT_TALL_UP", "", paste(sort(ross.up$hgnc_symbol[1:min(nrow(ross.up), topN)]), collapse="\t"))
gmt[nrow(gmt)+1,] <- c("ROSS_2003_ER_VS_NOT_TALL_DN", "", paste(sort(ross.dn$hgnc_symbol[1:min(nrow(ross.dn), topN)]), collapse="\t"))

#############################################################################################################
# Tijssen-Gottwald (2011)
# http://www.ncbi.nlm.nih.gov/pubmed/21571218
#############################################################################################################
tj <- read.csv("/mnt/projects/chrisi/results/chipseq/Tijssen_all.genes.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T)
tj.runx1 <- sort(unique(as.vector(do.call(c, tj))))
tj.runx1 <- tj.runx1[tj.runx1 != ""]

gmt[nrow(gmt)+1,] <- c("TIJSSEN_2011_BOUND_BY_RUNX1", "http://www.ncbi.nlm.nih.gov/pubmed/21571218", paste(tj.runx1, collapse="\t"))

#############################################################################################################
# RUNX1 ChIP-seq Niebuhr et al. (2013)
# http://www.ncbi.nlm.nih.gov/pubmed/23704093
#############################################################################################################
ni <-  read.csv("/mnt/projects/chrisi/results/chipseq/Niebuhr_TableS3_Runx1 Peaks Called in ProB-Cells.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T)
#ni[1:5,]
#par(mfrow=c(2,2)); hist(ni$dist_tss);  hist(ni$score); plot(ni$dist_tss,ni$score, pch=20, cex=0.8); plot(ni$dist_tss,ni$score, xlim=c(-50000, 50000), pch=20, cex=0.8)

cutoffDist <- c(-5000, 1000)
nif <- ni[ which(ni$dist_tss > cutoffDist[1] & ni$dist_tss < cutoffDist[2]), ]
#hist(nif$dist); hist(nif$score, br=100); length( unique( nif$nearest.gene))
#plot(nif$dist_tss,nif$score, pch=20, cex=0.8);

scoreCutoff <- 100
nifs <- nif[which(nif$score > scoreCutoff),]
#hist(nifs$dist); hist(nifs$score, br=100); length( unique( nifs$nearest.gene))

mausIDs <- unique(nifs$nearest.gene)

## get orthologs
humOrt <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol", values=mausIDs, mart = mouse,
                 attributesL = c("hgnc_symbol", "entrezgene"), martL = human)
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")#
if( length(which(is.na(humOrt$entrezgene)))>0 ) { humOrt <- humOrt[-which(is.na(humOrt$entrezgene)),] }
rm(mausIDs)

allnieRunx1hu <- sort(unique(humOrt$hgnc[humOrt$hgnc != ""]))
uniqOrt <- unique(humOrt[,1:2])

By <- by( uniqOrt, uniqOrt$mgi_symbol, function(x) { y <-  paste(sort(x$hgnc),collapse="|"); y }, simplify=F )
dfOrt <- as.data.frame(do.call(rbind, By))
dfOrt$V1 <- gsub("^\\|", "", dfOrt$V1)

By <- by( ni, ni$nearest.gene, function(x) { y <-  paste(sort(x$dist_tss),collapse="|"); y }, simplify=F )
dfni <- as.data.frame(do.call(rbind, By))
#head(dfni)

ni1Line <- merge(dfni, dfOrt, by="row.names", all=F)
colnames(ni1Line) <- c("mausSym", "Niedist","menschSym")
#head(ni1Line)
#length(unique(ni1Line$menschSym))

gmt[nrow(gmt)+1,] <- c("NIEBUHR_2013_BOUND_BY_RUNX1", "http://www.ncbi.nlm.nih.gov/pubmed/23704093", paste(unique(sort(ni1Line$menschSym)), collapse="\t"))

#############################################################################################################
# Wilson-Goettgens (2009)
# http://www.ncbi.nlm.nih.gov/pubmed/20887958
#############################################################################################################
wi <- read.csv("/mnt/projects/chrisi/results/chipseq/Wilson_Gottgens_ChIPseq.txt",  stringsAsFactors=F, sep="\t", header=T, fill=T)

# bound by RUNX1
wiRunx1 <- unique(wi$Runx1); wiRunx1 <- wiRunx1[ -which(wiRunx1 == "") ] 
humOrt <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol", values=wiRunx1, mart = mouse,
                 attributesL = c("hgnc_symbol", "entrezgene"), martL = human )
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[-which(is.na(humOrt$entrezgene)),]
wihumOrt <- humOrt; colnames(wihumOrt) <- c("Wi_mgi_symbol", "Wi_hgnc", "Wi_entrezgene")
wiRunx1hu <- sort( unique(humOrt$hgnc[humOrt$hgnc != ""]) ) 

# bound by "heptad" (Pimkin et al, 2014), including FLI1, GATA2, LYL1, TAL1 (SCL), ERG, RUNX1, and LMO2
wiHeptad <- wi$Scl[wi$Scl %in% wi$Fli.1 & wi$Scl %in% wi$Gata2 & wi$Scl %in% wi$Lyl1 & wi$Scl %in% wi$Erg & wi$Scl %in% wi$Runx1 & wi$Scl %in% wi$Lmo2]
humOrt <- getLDS(attributes = c("mgi_symbol"),
                 filters = "mgi_symbol", values=wiHeptad, mart = mouse,
                 attributesL = c("hgnc_symbol", "entrezgene"), martL = human )
colnames(humOrt) <- c("mgi_symbol", "hgnc", "entrezgene")
humOrt <- humOrt[-which(is.na(humOrt$entrezgene)),]
wihumOrt <- humOrt; colnames(wihumOrt) <- c("Wi_mgi_symbol", "Wi_hgnc", "Wi_entrezgene")
wiHeptadhu <- sort( unique(humOrt$hgnc[humOrt$hgnc != ""]) ) 

gmt[nrow(gmt)+1,] <- c("WILSON_2009_BOUND_BY_RUNX1", "http://www.ncbi.nlm.nih.gov/pubmed/20887958", paste(wiRunx1hu, collapse="\t"))
gmt[nrow(gmt)+1,] <- c("WILSON_2009_BOUND_BY_HEPTAD", "http://www.ncbi.nlm.nih.gov/pubmed/20887958", paste(wiHeptadhu, collapse="\t"))

write.table(gmt, file="/mnt/projects/helena_veronika/data/gsea_custom_gene_sets.gmt", row.names=F, col.names=F, quote=F, sep="\t")
