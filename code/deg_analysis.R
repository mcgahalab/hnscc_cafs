# McGuigan RNAseq analysis of G9a/GLP inhibitor A366 and OICR on CAFs

# DEG analysis
library(DESeq2)
library(pheatmap)
library(ComplexHeatmap)
# General
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(reshape2)
library(data.table)
library(RColorBrewer)
# Enrichment
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
# Custom
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")

if(exists("wideScreen")) wideScreen()
pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/mcguigan_and_nila_ram/analysis'
setwd(pdir)


# Load in ChEA3 genesets
cheadir <- '/cluster/projects/mcgahalab/ref/chea3'
cheaf <- list.files(cheadir, pattern=".gmt$")
if(!file.exists(file.path(cheadir, "chea.rds"))){
  chea <- lapply(file.path(cheadir, cheaf), function(pathtochea){
    print(pathtochea)
    cdf <- read.table(pathtochea, sep="\t", header=F, stringsAsFactors = F, 
               check.names = F, fill=T) %>%
      mutate(V1=make.unique(V1)) %>%
      column_to_rownames(., "V1") %>% 
      apply(., 1, function(i) list(i[i!=''])) %>%
      unlist(., recursive = F)
    
    cdf[sapply(cdf, length) > 5]
  }) %>%
    setNames(., gsub(".gmt", "", cheaf))
  
  saveRDS(chea, file=file.path(cheadir, "chea.rds"))
} else {
  chea <- readRDS(file.path(cheadir, "chea.rds"))
}

# Directories
rdsdir <- file.path("data", "deseq2_rds", "rds")

# Intermediate/data files
cnts_f <- file.path("results", "data", "counts.csv")
meta_f <- file.path("results", "data", "meta.csv")
deseq_f <- file.path("results", "data", "deseq.rds")
deseq_vst_f <- file.path("results", "data", "deseq_vst.rds")
deseq_res_f <- file.path("results", "data", "deseq_res.rds")
deseq_gsea_f <- file.path("results", "data", "deseq_gsea.rds")
deseq_tfea_f <- file.path("results", "data", "deseq_tfea.rds")
deg_q_thresh <- 0.05
max_pc <- 5

# plotting params
sig_fills <- c('A366_only'='#2b8cbe', 'OICR_only'='#e34a33', 'Intersect'='#aea40c') 
sig_cols <- setNames(rep("black",3), names(sig_fills))
sig_sizes <- setNames(c(1,1,1.5), names(sig_fills))
grp_cols <- c('DMSO'='#2b8cbe', 'A366'='#a50f15', 'OICR'='#fb6a4a') #blue, red, red
# Generate 29 unique colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_sets <- c('Set1', 'Set2', 'Set3')
col_v = unlist(mapply(brewer.pal, qual_col_pals[col_sets,]$maxcolors, col_sets))
overlap_key <- c('nonSig'='-', 'OICR_only'='*', 'A366_only'='**', 'Intersect'='')

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

# create a GRCh38 txdb mapping

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
  overlap_cnt <- overlap_cnt %>%
    mutate("fraction"=round(Count/sum(Count),2))
  ggplot(overlap_cnt, aes(x=Count, y=grp, fill=sig, label=fraction)) +
    geom_bar(stat='identity', position='stack', alpha=0.75) +
    geom_text(size = 3, position = position_stack(vjust = 0.5)) +
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
  print(geneset)
  # setup
  gsea_inh1 <- gseas$A366[[geneset]]
  gsea_inh2 <- gseas$OICR[[geneset]]
  col_idx <- c(1,5,6,7,10,11)
  
  # merge gsea dataframes
  gsea_merge <- full_join(as.data.frame(gsea_inh1)[,col_idx],
                          as.data.frame(gsea_inh2)[,col_idx],
                          by='ID', suffix=c('.A366', '.OICR')) %>%
    column_to_rownames(., "ID")
  qcols <- grep("p.adjust", colnames(gsea_merge), value=T)
  # for(q in qcols){
  #   pcol <- gsub('p.adjust', 'pvalue', q)
  #   gsea_merge[,q] <- p.adjust(gsea_merge[,pcol], method='fdr')
  # }
  
  # Split into Intersect, A366_only, OICR_only
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
qcols <- grep("p.adjust", colnames(gsea_merge), value=T)
for(q in qcols){
  pcol <- gsub('p.adjust', 'pvalue', q)
  gsea_merge[,q] <- p.adjust(gsea_merge[,pcol], method='fdr')
}
# Split into Intersect, A366_only, OICR_only
gsea_tmp <- venn(gsea_merge, qcols, deg_q_thresh)
gsea_merge <- gsea_tmp$res
# geneset_ois <- sapply(gsea_l$sig_genesets, function(geneset_i) geneset_i$Intersect)

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
  please list from the following list: \n")))
gsea_merge <- gsea_merge %>% 
  mutate(fishersp=apply(gsea_merge[,c('pvalue.A366', 'pvalue.OICR')], 1, function(i)
    metap::sumlog(i)$p)) %>%
  mutate(p.adjust.fishers=p.adjust(fishersp, method='fdr'))
if(is.null(geneset_oi)) {
  geneset_oi <- gsea_merge %>% 
    mutate(pcum=abs(1-(p.adjust.A366 + p.adjust.OICR)),
           dir=((NES.A366 > 0) | (NES.OICR > 0))) %>%
    # filter((p.adjust.A366 < deg_q_thresh) | (p.adjust.OICR < deg_q_thresh)) %>%
    filter(p.adjust.fishers < deg_q_thresh) %>%
    group_by(C, dir) %>%
    top_n(., 5, p.adjust.fishers) %>% 
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
           stringr::str_to_sentence(.)) %>%
  melt %>% 
  rename_with(., ~ c("msig", "Geneset", "Overlap", "Inhibitor", "NES"))
gsea_melt <- gsea_melt %>%
  mutate(Geneset=factor(Geneset, 
                        levels=unique(gsea_melt[order(gsea_melt$NES),'Geneset'])),
         overlap=overlap_key[Overlap])

gg_bar <- ggplot(gsea_melt, aes(y=Geneset, x=NES, #color=Overlap,
                                fill=Inhibitor, group=Inhibitor,
                                label = overlap)) +
  geom_hline(aes(yintercept=Geneset), linetype='dashed', color='grey', size=0.25) + 
  geom_bar(stat='identity', position='dodge', size=0.7, width=0.8) +
  scale_fill_manual(values=c('A366'='grey1',
                             'OICR'='grey50')) +
  # scale_color_manual(values=c(sig_fills,
  #                            'nonSig'='grey')) +
  geom_text(x=2.8, hjust = 0, vjust=0.5) +
  facet_grid(msig ~ ., scales='free', space='free',switch = "y") +
  theme_classic() +
  xlim(-3,3) + ylim(-3,3) +
  geom_vline(xintercept=0, color='black') +
  theme(axis.text.y = element_text(size=6, hjust=1),
        strip.text.y.left = element_text(angle=90, size=8),
        legend.position = 'right',
        legend.margin =margin(0,0,0,0)) +
  guides(fill=guide_legend(title="q < 0.05")) +
  ylab("") + xlab("NES") +
  scale_y_discrete(labels =  scales::wrap_format(50))

pdf(file=file.path("results", "gsea", "gsea_bar.pdf"), height = 7, width = 6)
gg_bar
dev.off()
write.table(gsea_merge, file.path("results", "gsea", "gsea_merge.tsv"),
            sep="\t", col.names = T, row.names = F, quote = F)
gsea_merge[grep("centriole", gsea_merge$geneset, ignore.case = T), ] 

########################################################
#### 6. DEG Heatmap of significant genes & genesets ####
if(!exists("resl")) resl <- readRDS(deseq_res_f)
if(!exists("gsea_l")) gsea_l <- readRDS(deseq_gsea_f)
if(!exists("vsd")) vsd <- readRDS(file=deseq_vst_f)

## >1. Set up geneset dataframe
gsea_merge <- do.call(rbind, gsea_l$gsea_by_geneset[c(1:5)]) %>%
  as.data.frame 
geneset_ois <- sapply(gsea_l$sig_genesets, function(geneset_i) geneset_i$Intersect)
gsea_merge <- cbind(gsea_merge,
                    strsplit(rownames(gsea_merge), split="\\.") %>%
                      do.call(rbind,.) %>% 
                      as.data.frame %>%
                      rename_with(., ~ c("C1", "C2", "geneset")) %>%
                      mutate(C=paste0(C1, ".", C2)))

if(is.null(geneset_oi)) warning(cat(paste0(
  "> No genesets specified, defaulting to top5 per geneset, 
  please list from the following list: \n",
  paste(paste0("\t- ", unlist(geneset_ois)), collapse="\n"))))
if(is.null(geneset_oi)) {
  gsea_sig <- gsea_merge %>% 
    mutate(pcum=abs(1-(p.adjust.A366 + p.adjust.OICR)),
           C2=factor(C2, levels=c('GO:MF', 'Base', 'CP:REACTOME', 'GO:BP', 'GO:CC'))) %>%
    filter((p.adjust.A366 < deg_q_thresh) | (p.adjust.OICR < deg_q_thresh)) %>%
    arrange(., C2) %>%
    group_by(C) %>%
    top_n(., 5, pcum)
} else {
  gsea_sig <- gsea_merge[which(gsea_merge$geneset %in% geneset_oi),]
}

# Identify the genes in the top genesets or genesets of interest
geneset_genes <- apply(gsea_sig, 1, function(gsea_row_i){
  # Collapse core_enriched genes, split, and convert entrez->symbol
  entrez_i <- gsea_row_i[c('core_enrichment.OICR', 'core_enrichment.A366')] %>%
    paste(., collapse="/") %>%
    strsplit(., split="\\/")
  
  sym_i <- entrez2sym_ids[entrez_i[[1]]] %>% unique
}) %>% 
  setNames(., gsea_sig$geneset)
geneset_genes <- lapply(geneset_genes, data.frame) %>%
  rbindlist(., idcol='Geneset') %>%
  rename_with(., ~c("Geneset", "Gene")) %>%
  filter(!duplicated(Gene))

## >2. Set up DEGs
all_sig_genes <- resl$all_sig
split_sig_genes <- lapply(resl$split_sig, data.frame) %>%
  rbindlist(., idcol='split') %>%
  rename_with(., ~c("Overlap", "Gene"))

# Combine Overlap genes with Geneset Genes
geneset_split_genes <- full_join(geneset_genes, split_sig_genes, by='Gene') %>%
  data.frame %>% 
  filter(!is.na(Overlap)) %>%
  arrange(., Overlap, Geneset) %>%
  column_to_rownames(., "Gene")

## >3. Set up variance-scaled gene-expression matrix and scale it
tcnts <- as.data.frame(assay(vsd))[all_sig_genes,] 
z_cnts <- apply(tcnts, 1, scale) %>%
  as.data.frame %>%
  mutate(samples=colnames(tcnts)) %>%
  column_to_rownames(., 'samples') %>%
  t %>% as.data.frame
meta <- vsd@colData %>%
  as.data.frame %>% 
  mutate(sampletype=factor(sampletype, levels=c('DMSO', 'A366', 'OICR'))) %>%
  arrange(sampletype)

## >4. Do the Vizsters
unique_genesets <- unique(geneset_split_genes$Geneset)
unique_genesets <- unique_genesets[!is.na(unique_genesets)]
genes_geneset <- split(geneset_split_genes, f=is.na(geneset_split_genes$Geneset))
genes_geneset_overlap <- split(geneset_split_genes, f=geneset_split_genes$Overlap)
ann_colors <- list("sampletype"=grp_cols,
                   "Overlap"=sig_fills,
                   "Geneset"=setNames(col_v[1:length(unique_genesets)], unique_genesets))

pdf(file.path("results", "deg", "deg_heatmap_geneset.pdf"), width = 9, height = 9)
gene_meta_i <- genes_geneset$`FALSE`
ComplexHeatmap::pheatmap(z_cnts[match(rownames(gene_meta_i), rownames(z_cnts)),
                match(rownames(meta), colnames(z_cnts))], 
         cluster_rows=FALSE, show_rownames=TRUE, cluster_cols=FALSE, 
         annotation_col=meta[,1,drop=F], fontsize_row = 4,
         annotation_row=gene_meta_i, annotation_colors=(ann_colors))
dev.off()

pdf(file.path("results", "deg", "deg_heatmap_noGeneset.pdf"), width = 4.5, height = 9)
gene_meta_i <- genes_geneset$`TRUE`
ComplexHeatmap::pheatmap(z_cnts[match(rownames(gene_meta_i), rownames(z_cnts)),
                      match(rownames(meta), colnames(z_cnts))], 
               cluster_rows=TRUE, show_rownames=FALSE, cluster_cols=FALSE, 
               annotation_col=meta[,1,drop=F],
               annotation_row=gene_meta_i, annotation_colors=(ann_colors))
dev.off()

pdf(file.path("results", "deg", "deg_heatmap.intersect.pdf"), width = 4.5, height = 9)
ComplexHeatmap::pheatmap(z_cnts[match(rownames(genes_geneset_overlap$Intersect), rownames(z_cnts)),
                                match(rownames(meta), colnames(z_cnts))], 
                         cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE, 
                         annotation_col=meta[,1,drop=F],
                         annotation_colors=(ann_colors), fontsize_row=3)
dev.off()

####################################################
#### 7. Analyze enrichment of TF motifs in DEGs ####
if(!exists("resl")) resl <- readRDS(deseq_res_f)


# Estimate the geneset enrichment per inhibitor
lfc_col <- 'log2FoldChange'
lfcs <- list("A366"=resl$res[,paste0(lfc_col, ".", "A366"), drop=F],
             "OICR"=resl$res[,paste0(lfc_col, ".", "OICR"), drop=F])

# gseaFunction with an increased max GSSize
gseaFun <- function(msig_ds, lfc_v){
  gsea <- tryCatch({
    GSEA(sort(na.omit(lfc_v), decreasing = T), 
         TERM2GENE = msig_ds, pvalueCutoff = 1, maxGSSize=3000)
  }, error=function(e){NULL})
  return(gsea)
}

# Run GSEA on the TF databases from ChEA3
tfeas <- lapply(lfcs, function(lfc_i){
  lfc_v <- setNames(lfc_i[,1], rownames(lfc_i))
  
  tfea <- lapply(chea, function(chea_i){
    tryCatch({
      gsea_x <- iterateMsigdb(species='Homo sapiens', fun=gseaFun, lfc_v=lfc_v, 
                     msig_lvls=list("custom"=chea_i)) %>%
      unlist(., recursive=F)
      
      keep_idx <- sapply(gsea_x, nrow) >0
      lapply(gsea_x[keep_idx], as.data.frame) %>%
        do.call(rbind, .) %>% 
        as.data.frame %>%
        arrange(p.adjust) %>%
        mutate("TF"=names(chea_i)[keep_idx])
    }, error=function(e){NULL})
  })
  tfea_df <- rbindlist(tfea, idcol="database")[,c(1,2,5:11)] %>%
    arrange(p.adjust)
  return(tfea_df)
})
saveRDS(tfeas, file=file.path("results", "deg", "tfea.rds"))


# Ensure proper p.adjusted values
tfeas <- readRDS(file=file.path("results", "deg", "tfea.rds"))
tfeas <- lapply(tfeas, function(tfea){
  tfea %>% 
    group_by(database) %>%
    mutate(p.adjust = p.adjust(pvalue, method='fdr'))
})
col_sel <- c('database', 'ID', 'NES', 'p.adjust')

# Group together the two inhibitors and combine the p-values (and fdr adjust)
tfea_merge <- purrr::reduce(tfeas, full_join, by=c('database', 'ID'), 
                            suffix=paste0(".", names(tfeas))) %>%
  filter(!is.na(pvalue.A366) & !is.na(pvalue.OICR))
pcols <- grep("pvalue", colnames(tfea_merge), value=T)
padj_col <- grep("p.adjust", colnames(tfea_merge), value=T)

fishersp <- apply(tfea_merge, 1, function(i) metap::sumlog(as.numeric(i[pcols]))$p)
tfea_merge <- tfea_merge %>%
  as.data.frame %>%
  mutate(fishersp=fishersp,
         dir=((NES.A366 > 0) | (NES.OICR > 0))) %>%
  mutate(p.adjust.fishers=p.adjust(fishersp, method='fdr'))
rownames(tfea_merge) <- apply(tfea_merge[,c('database', 'ID')], 1, paste, collapse=".")

# Annotate interesect, A366_only, OICR_only
res_tmp <- venn(df = tfea_merge, pcols = padj_col, qthresh = deg_q_thresh)
tfea_merge <- res_tmp$res

# Extract the top significant TFs across both inhibitors
regulon_oi <- tfea_merge %>% 
  filter(p.adjust.fishers < deg_q_thresh) %>%
  mutate(C=gsub("_.*", "", ID)) %>%
  group_by(dir) %>%
  slice_min(., p.adjust.fishers, n=5) %>%
  select(ID)
regulon_oi <- regulon_oi$ID %>% 
  gsub("_.*", "", .) %>%
  unique

# Select and merge the top TFEA values, calculate SD around each
top_tfea_l <- sapply(regulon_oi, function(i) {
  i <- i %>% 
    gsub("_.*", "", .) %>%
    gsub("$", "(_.*)?$", .) %>%
    paste0("^", .)
  grep(i, x=tfea_merge$ID, value=T, perl=T)
})
tfea_merge_summ <- lapply(top_tfea_l, function(tfea_id){
  tfea_merge[which(tfea_merge$ID %in% tfea_id),,drop=F] %>%
    group_by(database) %>%
    summarise(NES.A366.mean=mean(NES.A366),
              NES.OICR.mean=mean(NES.OICR),
              dir=any(dir),
              p.adjust.fishers=p.adjust.fishers,
              sig=sig)
}) %>%
  rbindlist(., idcol="TF") %>%
  mutate(overlap=overlap_key[sig])

topn_tfs <- tfea_merge_summ %>% 
  group_by(TF) %>%
  summarise(NES=abs(sum(NES.A366.mean, NES.OICR.mean)),
            dir=any(dir)) %>%
  group_by(dir) %>%
  slice_max(., NES, n=5)
tfea_merge_summ <- tfea_merge_summ %>%
  filter(TF %in% topn_tfs$TF)

# Melt and plot
nes_cols <- grep("NES", colnames(tfea_merge_summ), value=T) %>%
  setNames(., gsub("^.*\\.", "", .)) 
melt_tfea <- tfea_merge_summ %>% 
  select(-p.adjust.fishers) %>%
  melt %>% 
  mutate(variable=gsub("^.*\\.(.*)\\..*", "\\1", variable), 
         database=gsub("_Coexpression", "", database) %>%
           gsub("_ChIP-seq", "", .)) %>%
  rename_with(., ~c("TF", "Database", "Overlap", "overlap", "Inhibitor", "NES"))


gg_bar <- ggplot(melt_tfea, aes(y=Database, x=NES, label = overlap, #color=Overlap,
                                fill=Inhibitor, group=Inhibitor)) +
  geom_hline(aes(yintercept=Database), linetype='dashed', color='grey', size=0.25) + 
  geom_bar(stat='identity', position='dodge', size=0.7, width=0.8) +
  geom_text(x=2.8, hjust = 0, vjust=0.5) +
  scale_fill_manual(values=c('A366'='grey1',
                             'OICR'='grey50')) +
  #scale_color_manual(values=c(sig_fills,
  #                            'nonSig'='grey')) +
  facet_grid(TF ~ ., scales='free', space='free',switch = "y") +
  theme_classic() +
  xlim(-3,3) + ylim(-3,3) +
  geom_vline(xintercept=0, color='black') +
  theme(axis.text.y = element_text(size=6, hjust=0),
        strip.text.y.left = element_text(angle=0, size=8),
        legend.position = 'right',
        legend.margin =margin(0,0,0,0)) +
  guides(fill=guide_legend(title="q < 0.05")) +
  ylab("") + xlab("NES") +
  scale_y_discrete(labels =  scales::wrap_format(40))

pdf(file=file.path("results", "tf", "tfea_bar.pdf"), height = 7, width = 6)
gg_bar
dev.off()

saveRDS(list("tfea"=tfeas, "tfea_merge"=tfea_merge), 
        file=deseq_tfea_f)

# tfea_merge[grep("SU", tfea_merge$ID, ignore.case = T),]

############################################
#### 8. Write out CSV for DEG/GSEA/TFEA ####
if(!exists("resl")) resl <- readRDS(deseq_res_f)
if(!exists("tfeal")) tfeal <- readRDS(deseq_tfea_f)
if(!exists("gseal")) gseal <- readRDS(deseq_gsea_f)

gsea_merge <- do.call(rbind, gseal$gsea_by_geneset[c(1:5)]) %>%
  as.data.frame %>%
  rownames_to_column(.,var="geneset")
gsea_merge <- gsea_merge[,-grep("core_enrichment", colnames(gsea_merge))]

tfea_merge <- tfeal$tfea_merge

deg_merge <- resl$res

dir.create(file.path("results", "tables"))
.writeTbl <- function(x,f){
  write.table(x, file=f, quote=F, sep=",", row.names = T, col.names = T)
}
.writeTbl(deg_merge, file.path("results", "tables", "deg.csv"))
.writeTbl(gsea_merge, file.path("results", "tables", "gsea.csv"))
.writeTbl(tfea_merge, file.path("results", "tables", "tfea.csv"))
#### Old JASPAR2020 Code ####
stop("DO NOT RUN")
# Get JASPAR2020 motifs for homo sapiens
motifs <- TFBSTools::getMatrixSet(JASPAR2020, list("species"="Homo sapiens",
                                                   "collection"=))
names(motifs) <- paste(names(motifs), TFBSTools::name(motifs),sep = "_")
# Get list of all promoters for all genes
promoters <- getPromoters(TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, by='gene', 
                          upstream = 1000, downstream = 1000) %>%
  keepStandardChromosomes(., pruning.mode='coarse')
txdb_genes <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb_genes$symbol <- entrez2sym_ids[txdb_genes$gene_id]
ov_idx <- findOverlaps(promoters, txdb_genes)

# Find max number of segments that overlap one reference bin
max_n <- max(table(queryHits(ov_idx)))
# Populate segmean matrix, iterating over multiple intervals per ref_bin
promoter_mat <- data.frame(matrix(nrow=length(promoters), ncol=max_n))

ov_idx_i <- ov_idx
idx <- 1
dup_status <- TRUE
while(dup_status){
  # Split based on bin with multiple overlaps 
  dup_status <- any(duplicated(queryHits(ov_idx_i)))
  spl_ov <- split(ov_idx_i, duplicated(queryHits(ov_idx_i)))
  if(dup_status) ov_idx_i <- spl_ov[['TRUE']]
  
  # Fill in segmean_matrix with all non-duplicated entries
  non_dup_ov <- spl_ov[['FALSE']]
  promoter_mat[queryHits(non_dup_ov),idx] <- elementMetadata(txdb_genes)[subjectHits(non_dup_ov),2]
  idx <- idx+1
}

# Reduce promoters into one value
promoter_mat <- as.matrix(promoter_mat)
elementMetadata(promoters) <- apply(promoter_mat, 1, function(i) paste(na.omit(i), collapse=",")) %>% 
  as.data.table

# match motifs to promoters
motif_ix <- matchMotifs(motifs, promoters,
                        genome = BSgenome.Hsapiens.UCSC.hg38)

# Create motif-specific genesets
motif_geneset <- apply(assay(motif_ix), 2, function(i) {
  genes <- elementMetadata(promoters[which(i)])[,1] %>%
    paste(., collapse=",") %>%
    strsplit(., split=",")
  return(genes[[1]])
})
