# DEG analysis
library(DESeq2)
library(pheatmap)
# General
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(reshape2)
# Enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
# Custom
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")

wideScreen()

pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/mcguigan_and_nila_ram/analysis'
setwd(pdir)

# Directories
rdsdir <- file.path("data", "deseq2_rds", "rds")

# Intermediate/data files
cnts_f <- file.path("results", "data", "counts.csv")
meta_f <- file.path("results", "data", "meta.csv")
deseq_f <- file.path("results", "data", "deseq.rds")
deseq_vst_f <- file.path("results", "data", "deseq_vst.rds")
deseq_res_f <- file.path("results", "data", "deseq_res.rds")
deseq_gsea_f <- file.path("results", "data", "deseq_gsea.rds")
deg_q_thresh <- 0.05
max_pc <- 5

# plotting params
sig_fills <- c('A366_only'='#2b8cbe', 'OICR_only'='#e34a33', 'Intersect'='#aea40c') 
sig_cols <- setNames(rep("black",3), names(sig_fills))
sig_sizes <- setNames(c(1,1,1.5), names(sig_fills))
grp_cols <- c('DMSO'='#2b8cbe', 'A366'='#a50f15', 'OICR'='#fb6a4a') #blue, red, red

# gsea params:
geneset_oi <- NULL
# geneset_oi <- c('HALLMARK_G2M_CHECKPOINT',
#                 'HALLMARK_E2F_TARGETS',
#                 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
#                 'REACTOME_CELL_CYCLE_MITOTIC',
#                 'REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS',
#                 'GOBP_MITOTIC_SISTER_CHROMATID_SEGREGATION', 
#                 'GOBP_DNA_REPLICATION',
#                 'GOCC_CHROMOSOME_CENTROMERIC_REGION',
#                 'GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT', 
#                 'GOMF_EXTRACELLULAR_MATRIX_STRUCTURAL_CONSTITUENT_CONFERRING_TENSILE_STRENGTH',
#                 'GOMF_CHROMATIN_BINDING')
# Create symbol<->entrezID mapping
genome_gse <- org.Hs.eg.db
txby <- keys(genome_gse, 'SYMBOL')
sym2entrez_ids <- mapIds(genome_gse, keys=txby, column='ENTREZID',
                         keytype='SYMBOL', multiVals="first")
entrez2sym_ids <- setNames(names(sym2entrez_ids), sym2entrez_ids)

###################
#### Functions ####
#
venn <- function(df, pcols, qthresh){
  sig <- apply(df[,pcols], 2, function(i){names(which(i<qthresh))})
  all_sig <- unlist(sig) %>% unique()
  split_sig <- list("A366_only"=setdiff(sig[[1]], sig[[2]]),
                    "Intersect"=intersect(sig[[1]], sig[[2]]),
                    "OICR_only"=setdiff(sig[[2]], sig[[1]]))
  
  df$sig <- 'nonSig'
  for(split_id in names(split_sig)){
    df[split_sig[[split_id]],]$sig <- split_id
  }
  return(list("res"=df, "sig"=split_sig))
}

formatStackedOverlap <- function(split_sig){
  sapply(split_sig, length) %>% 
    as.data.frame %>% 
    rename_with(., ~ "Count") %>%
    rownames_to_column(., "sig") %>%
    mutate(grp='CAF')
}

plotStackedOverlap <- function(overlap_cnt){
  ggplot(overlap_cnt, aes(x=Count, y=grp, fill=sig)) +
    geom_bar(stat='identity', position='stack', alpha=0.75) +
    scale_fill_manual(values=sig_fills) +
    theme_classic() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(), 
          legend.position = 'none') 
}

scatterplotIt <- function(df, xval, yval, sig_cols=NULL,
                          sig_fills=NULL, sig_sizes=NULL){
  ggplot(df, aes_string(x=xval, y=yval, fill='sig', color='sig')) +
    geom_point(pch=21, alpha=0.75) +
    scale_color_manual(values=c(sig_cols,
                                'nonSig'='grey')) +
    scale_fill_manual(values=c(sig_fills,
                               'nonSig'='grey')) +
    scale_size_manual(values=c(sig_sizes,
                               'nonSig'=0.5)) +
    theme_classic() +
    xlim(-2.5,6) + ylim(-2.5,6) +
    theme(legend.position = 'right',
          legend.margin =margin(0,0,0,0)) +
    guides(fill=guide_legend(title="q < 0.05"),
           color=guide_legend(title="q < 0.05"))
}
##################################################################
#### 1. Extract the Counts-matrix and the Metadata for DESeq2 ####
if(!file.exists(cnts_f)){
  dat <- lapply(list.files(rdsdir, pattern='RDS'), function(rds_f){
    load(file.path(rdsdir, rds_f))
    list('meta'=meta, "cnt"=data)
  })
  names(dat) <- list.files(rdsdir, pattern='RDS') %>% gsub(".RDS", "",  .)
  
  cnt <- do.call(cbind, lapply(dat, function(i) i$cnt)) %>%
    select(which(!duplicated(colnames(.)))) %>%
    as.data.frame
  
  meta <- do.call(rbind, lapply(dat, function(i) i$meta)) 
  dup_idx <- gsub("(_[0-9]).*", "\\1", rownames(meta)) %>%
    duplicated
  meta <- meta[!dup_idx,,drop=F]
  meta$sampletype <- gsub("_[0-9]*", "", rownames(meta))
  
  
  write.table(cnt, file=cnts_f, sep=",", 
              quote=F, row.names = T, col.names = T)
  write.table(meta, file=meta_f, sep=",", 
              quote=F, row.names = T, col.names = T)
}

################################
#### 2. Build DESeq2 object ####
if(!file.exists(deseq_f)){
  cnt <- read.csv(cnts_f, header=T, stringsAsFactors = F, check.names = F)
  meta <- read.csv(meta_f, header=T, stringsAsFactors = F, check.names = F)
  min_read_count <- 10
  
  # Combine all samples and create one DESeq object
  dds <- DESeqDataSetFromMatrix(countData=cnt[,rownames(meta)],
                                colData=meta,
                                design=as.formula('~sampletype'))
  keep <- rowSums(assay(dds)) >= min_read_count
  dds <- dds[keep,]
  
  # Set 'DMSO' as the base-level that Inhibitor1/2 will be compared against
  dds$sampletype <- relevel(dds$sampletype, ref = "DMSO")
  dds <- DESeq(dds)
  
  saveRDS(dds, file=deseq_f)
}

#####################
#### 3. Plot PCA ####
if(!exists("dds")) dds <- readRDS(file=deseq_f)
vsd <- vst(dds, blind=T)
# vsd <- rlog(dds, blind=T)

tcnts <- as.data.frame(t(assay(vsd)))
pca <- prcomp(tcnts, scale=F)
percent_var <- round(pca$sdev^2/sum(pca$sdev^2),2) %>%
  setNames(., paste0("PC", c(1:length(.))))

pca_x <- as.data.frame(cbind(pca$x[,paste0("PC", c(1:max_pc))],
                             "condition"=as.character(vsd$sampletype))) %>%
  rownames_to_column(., 'Name')
for(id in paste0("PC", c(1:max_pc))){
  pca_x[,id] <- as.numeric(as.character(pca_x[,id]))
}


pcs <- paste0("PC", c(1:max_pc))
pcs_l <- lapply(c(2:length(pcs)), function(i){c(pcs[[i-1]], pcs[[i]])})
ggps <- lapply(pcs_l, function(pc){
  ggplot(data=pca_x, aes_string(x=pc[1], y=pc[2], color="condition", label="Name")) +
    geom_point() + 
    xlab(paste0(pc[1], " (", percent_var[pc[1]], ")")) +
    ylab(paste0(pc[2], " (", percent_var[pc[2]], ")")) +
    theme_classic() +
    scale_color_manual(values=grp_cols) +
    geom_text_repel(aes(label=Name), size=3, max.overlaps=12,max.iter=1000000)  
})
pdf(file.path("results", "pca", "pca.pdf"), height = 4, width = 5)
ggps
dev.off()

saveRDS(vsd, file=deseq_vst_f)

######################################
#### 4. Generate DEGs and compare ####
if(!exists("dds")) dds <- readRDS(file=deseq_f)
if(!exists("vsd")) vsd <- readRDS(file=deseq_vst_f)

# Calculate differential genes
reslfc_inh1 <- lfcShrink(dds, coef=resultsNames(dds)[2], type="apeglm")
reslfc_inh2 <- lfcShrink(dds, coef=resultsNames(dds)[3], type="apeglm")
resl <- list("A366"=reslfc_inh1,
             "OICR"=reslfc_inh2)

# Merge differential genes into one dataframe
suffix <- sapply(strsplit(resultsNames(dds)[-1], split="_"), function(i) i[2]) %>%
  paste0(".", .)
.dfGene <- function(x) x %>% as.data.frame() %>% rownames_to_column(., "gene")
res <- full_join(.dfGene(reslfc_inh1), .dfGene(reslfc_inh2),
                 by='gene', suffix=suffix) %>% 
  column_to_rownames(., "gene")

# Venn the p-values, assign differential flag to the LFCresults
padj_col <- grep("padj", colnames(res), value=T)
res_tmp <- venn(df = res, pcols = padj_col, qthresh = deg_q_thresh)
res <- res_tmp$res
split_sig <- res_tmp$sig

# Setup for plotting
lfc_cols <- grep("log2Fold", colnames(res), value=T) %>%
  setNames(., gsub(".*\\.", "", .))

# stacked barplot in venn-diagram style showing comparison of significant genes
overlap_cnt <- formatStackedOverlap(split_sig)
gg_vennbar <- plotStackedOverlap(overlap_cnt)

# scatterplot showing comparison of LFC between significant genes
gg_vennpnt <- scatterplotIt(res, lfc_cols[1], lfc_cols[2], 
                            sig_cols, sig_fills, sig_sizes) +
  xlab(paste0("Log2 fold-change (", names(lfc_cols)[1], ")")) + 
  ylab(paste0("Log2 fold-change (", names(lfc_cols)[2], ")"))

pdf(file.path("results", "deg", "lfc_venn.pdf"), height = 6, width = 6)
plot_grid(gg_vennbar, gg_vennpnt, ncol=1, align='v', axis='lr', rel_heights=c(1,7))
dev.off()

resl <- c(resl, list("res"=res,
                     "all_sig"=all_sig_genes,
                     "split_sig"=split_sig))
saveRDS(resl, file=deseq_res_f)

################################################
#### 5.a) Generate GSEA and compare between ####
if(!exists("resl")) resl <- readRDS(deseq_res_f)

# Estimate the geneset enrichment per inhibitor
lfc_col <- 'log2FoldChange'
lfcs <- list("A366"=resl$res[,paste0(lfc_col, ".", "A366"), drop=F],
             "OICR"=resl$res[,paste0(lfc_col, ".", "OICR"), drop=F])
gseas <- lapply(lfcs, function(lfc_i){
  lfc_v <- setNames(lfc_i[,1],
                    sym2entrez_ids[rownames(lfc_i)])
  
  unlist(iterateMsigdb(species='Homo sapiens', fun=gseaFun, lfc_v=lfc_v),
         recursive=F)
})

# Merge GSEA dataframes between inhibitors by geneset, count overlap, and plot
gg_gseas <- lapply(names(gseas[[1]]), function(geneset){
  # setup
  gsea_inh1 <- gseas$A366[[geneset]]
  gsea_inh2 <- gseas$OICR[[geneset]]
  col_idx <- c(1,5,7,10,11)
  
  # merge gsea dataframes
  gsea_merge <- full_join(as.data.frame(gsea_inh1)[,col_idx],
                          as.data.frame(gsea_inh2)[,col_idx],
                          by='ID', suffix=c('.A366', '.OICR')) %>%
    column_to_rownames(., "ID")
  
  qcols <- grep("p.adjust", colnames(gsea_merge), value=T)
  gsea_tmp <- venn(gsea_merge, qcols, deg_q_thresh)
  gsea_merge <- gsea_tmp$res
  split_sig <- gsea_tmp$sig
  
  # stacked barplot in venn-diagram style showing comparison of significant genesets
  overlap_cnt <- formatStackedOverlap(split_sig)
  gg_vennbar <- plotStackedOverlap(overlap_cnt) + 
    ggtitle(geneset)

  # scatterplot showing comparison of NES between significant genesets
  nes_cols <- grep("NES", colnames(gsea_merge), value=T) %>%
    setNames(., gsub("^.*\\.", "",.))
  gg_vennpnt <- scatterplotIt(gsea_merge, nes_cols[1], nes_cols[2], 
                              sig_cols, sig_fills, sig_sizes) +
    xlab(paste0("NES (", names(nes_cols)[1], ")")) + 
    ylab(paste0("NES (", names(nes_cols)[2], ")")) +
    xlim(-3,3) + ylim(-3,3)
  
  gg_grid <- plot_grid(gg_vennbar, gg_vennpnt, ncol=1, align='v', axis='lr', rel_heights=c(1,6))
  return(list("gg"=gg_grid, "overlap"=split_sig, "gsea"=gsea_merge))
})
names(gg_gseas) <- names(gseas[[1]])

# Do the vizzy-vizzies
pdf(file=file.path("results", "gsea", "gsea_comp_scatter.pdf"), height = 12, width = 12)
plot_grid(plotlist=lapply(gg_gseas[-6], function(i) i$gg), ncol=2)
dev.off()

# Save the intermediate files
gsea_l <- list("gsea_by_inh"=gseas,
               "gsea_by_geneset"=lapply(gg_gseas,function(i) i$gsea),
               "sig_genesets"=lapply(gg_gseas,function(i) i$overlap))
saveRDS(gsea_l, file=deseq_gsea_f)


######################################################
#### 5.b) Barplot of select significant gene-sets ####
if(!exists("gsea_l")) gsea_l <- readRDS(deseq_gsea_f)
gsea_merge <- do.call(rbind, gsea_l$gsea_by_geneset[c(1:5)]) %>%
  as.data.frame 
geneset_ois <- sapply(gsea_l$sig_genesets, function(geneset_i) geneset_i$Intersect)

# label the GSEA results by their geneset and msigdb level
gsea_merge <- cbind(gsea_merge,
                    strsplit(rownames(gsea_merge), split="\\.") %>%
                      do.call(rbind,.) %>% 
                      as.data.frame %>%
                      rename_with(., ~ c("C1", "C2", "geneset")) %>%
                      mutate(C=paste0(C1, ".", C2)))

# validation check
if(is.null(geneset_oi)) warning(cat(paste0(
  "> No genesets specified, defaulting to top5 per geneset, 
  please list from the following list: \n",
  paste(paste0("\t- ", unlist(geneset_ois)), collapse="\n"))))
if(is.null(geneset_oi)) {
  geneset_oi <- gsea_merge %>% 
    mutate(pcum=abs(1-(p.adjust.A366 + p.adjust.OICR))) %>%
    filter((p.adjust.A366 < deg_q_thresh) | (p.adjust.OICR < deg_q_thresh)) %>%
    group_by(C) %>%
    top_n(., 5, pcum) %>% 
    select(geneset)
  geneset_oi <- geneset_oi$geneset
}
  
nes_cols <- grep("NES", colnames(gsea_merge), value=T) %>%
  setNames(., gsub("^.*\\.", "", .)) 
gsea_melt <- gsea_merge %>% 
  filter(geneset %in% geneset_oi) %>%
  select(nes_cols, C, geneset, sig) %>%
  mutate(geneset = gsub("_", " ", geneset) %>%
           gsub("HALLMARK |GOBP |REACTOME |GOCC |GOMF ", "", .) %>%
           str_to_sentence) %>%
  melt %>% 
  rename_with(., ~ c("msig", "Geneset", "Overlap", "Inhibitor", "NES")) 

gg_bar <- ggplot(gsea_melt, aes(y=Geneset, x=NES, color=Overlap, 
                                fill=Inhibitor, group=Inhibitor)) +
  geom_hline(aes(yintercept=Geneset), linetype='dashed', color='grey') + 
  geom_bar(stat='identity', position='dodge', size=0.7, width=0.8) +
  scale_fill_manual(values=c('A366'='grey1',
                             'OICR'='grey50')) +
  scale_color_manual(values=c(sig_fills,
                             'nonSig'='grey')) +
  facet_grid(msig ~ ., scales='free', space='free',switch = "y") +
  theme_classic() +
  xlim(-3,3) + ylim(-3,3) +
  geom_vline(xintercept=0, color='black') +
  theme(axis.text.y = element_text(size=6, hjust=0),
        strip.text.y.left = element_text(angle=90, size=8),
        legend.position = 'right',
        legend.margin =margin(0,0,0,0)) +
  guides(fill=guide_legend(title="q < 0.05")) +
  ylab("") + xlab("NES") +
  scale_y_discrete(labels =  scales::wrap_format(40))

pdf(file=file.path("results", "gsea", "gsea_bar.pdf"), height = 6, width = 6)
gg_bar
dev.off()
  
##########################################################
#### 5.c) DEG Heatmap of significant genes & genesets ####
if(!exists("resl")) resl <- readRDS(deseq_res_f)
if(!exists("gsea_l")) gsea_l <- readRDS(deseq_gsea_f)
if(!exists("vsd")) vsd <- readRDS(file=deseq_vst_f)

##############
#### TEST ####





dds_vst <- vst(sub_dds, blind=T)
pdf(file.path(PDIR, "results", "manual", paste0("deg-heatmap_", ctype, "-", cycle, ".pdf")),
    height=15)
rownames(sub_coldata) <- sub_coldata$Samples
df <- sub_coldata[,c('RECIST', 'irAE', 'AGE')]
sub_topn <- if(topn < nrow(reslfc_sig)) topn else nrow(reslfc_sig)
dds_assay <- t(apply(assay(dds_vst)[reslfc_sig[1:sub_topn,]$ensg, sub_coldata$Samples], 1, scale))
colnames(dds_assay) <- sub_coldata$Samples
rownames(dds_assay) <- reslfc_sig[1:sub_topn,]$symbol
if(any(is.na(rownames(dds_assay)))) {
  na_idx <- is.na(rownames(dds_assay))
  rownames(dds_assay)[na_idx] <- reslfc_sig[1:sub_topn,]$ensg[na_idx]
}

order_idx <- order(sub_coldata$RECIST)
print(pheatmap(dds_assay[,order_idx], 
               cluster_rows=TRUE, show_rownames=TRUE,
               cluster_cols=FALSE, annotation_col=df[order_idx,]))
dev.off()


cnts <- counts(dds,normalized=TRUE)
meta_df <- as.data.frame(colData(dds)[,c("lnstatus","treatment", "celltype")])
order_idx <- with(meta_df, order(celltype, lnstatus, treatment))
meta_df_tdln <- split(meta_df, meta_df$lnstatus)$TDLN
order_idx_tdln <- with(meta_df_tdln, order(celltype, lnstatus, treatmen

#########################################
#### 4. Calculate differential genes ####
# Calculate differential 
