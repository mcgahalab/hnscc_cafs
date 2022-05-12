# Used to generate Figure 2-like UMAP plot of CAL-27 and CAL-33 cancer cell lines
# in comparison to the TCGA tumors from https://www.nature.com/articles/s41467-020-20294-x
git_repo <- '/Users/rquevedo/git/hnscc_cafs'
umap <- read.csv(file.path(git_repo, "ref", "celligner_suppdata.csv"), header=T) 
tissue_colors <- c(`other`='#999999',
                   `lung`='#1b9e77',
                   `upper_aerodigestive`='#e7298a',
                   `esophagus`='#d95f02',
                   `cervix`='#66a61e')

# Isolate lineages near the CAL-27/CAL-33 location on UMAP
umap$lineage_hnc = rep("other", nrow(umap))
hnc_idx <- which(umap$lineage %in% grep("other", names(tissue_colors), value = T, invert = T))
umap$lineage_hnc[hnc_idx] <-umap$lineage[hnc_idx]

# Isolate the CAL-27/33 cell lines in a separate dataframe
umap_cal <- umap[grep("(^CAL27)|(^CAL33)", umap$sampleID_CCLE_Name),] %>%
  as.data.frame() %>%
  mutate(sampleID=gsub("_.*", "", sampleID_CCLE_Name))


# Generate the umap plot
pdf(file.path(git_repo, "figures", "cal27_33_umap.pdf"), width = 6, height = 5)
ggplot(umap, aes(x=UMAP_1, y=UMAP_2, fill=lineage_hnc, size=type, color=type)) +
  geom_point(pch=21, alpha=0.5)  +
  scale_color_manual(values=c(`CL`='grey40', `tumor`='white')) +
  scale_size_manual(values=c(`CL`=0.6, `tumor`=0.5)) +
  theme_classic() + 
  scale_fill_manual(values=tissue_colors) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)),
         fill=guide_legend(ncol=1, byrow=T)) +
  theme(legend.position = 'right', 
        text=element_text(size=8),
        legend.margin =margin(0,0,0,0)) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  geom_point(data=umap_cal, aes(x=UMAP_1, y=UMAP_2), pch=21, fill='yellow', color='black', size=2, alpha=1) +
  geom_text_repel(data=umap_cal, aes(x=UMAP_1, y=UMAP_2, label=sampleID), color='black', size=4)
dev.off()
