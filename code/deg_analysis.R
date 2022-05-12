library(DESeq2)
library(dplyr)
library(tidyr)
library(tibble)

wideScreen()

pdir <- '/cluster/projects/mcgahalab/data/mcgahalab/mcguigan_and_nila_ram/analysis'
rdsdir <- file.path("data", "deseq2_rds", "rds")


cnts_f <- file.path("results", "data", "counts.csv")
meta_f <- file.path("results", "data", "meta.csv")
deseq_f <- file.path("results", "data", "deseq.rds")
deseq_vst_f <- file.path("results", "data", "deseq_vst.rds")
max_pc <- 5
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
  
  # Set 'Control' as the base-level that Inhibitor1/2 will be compared against
  dds$sampletype <- relevel(dds$sampletype, ref = "Control")
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
  rownames_to_column(., 'Name') %>%
  mutate(Sample=gsub("_.*", "", Name))
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
    scale_color_manual(values=c(`Control`='#fc8d59', 
                                `Inhibitor1`='#43a2ca', `Inhibitor2`='#253494')) +
    geom_text_repel(aes(label=Name), size=3, max.overlaps=12,max.iter=1000000)  
})
pdf(file.path("results", "pca", "pca.pdf"), height = 4, width = 5)
ggps
dev.off()

#####################
#### 3. Plot PCA ####


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
overall_reslfc_ln <- lfcShrink(dds, coef=coef, type="apeglm")


# Group-wise pairings between samples
metal <- split(meta, meta$sampletype)
pairs <- combn(names(metal), 2)
ddsl <- apply(pairs, 2, function(i){
  meta_i <- rbind(metal[[i[1]]], metal[[i[2]]])
  dds_i <- DESeqDataSetFromMatrix(countData=cnt[,rownames(meta_i)],
                                  colData=meta_i,
                                  design=as.formula('~sampletype'))
  return(dds_i)
}) %>%
  setNames(., apply(pairs, 2, paste, collapse="_"))



