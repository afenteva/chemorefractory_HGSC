library(readr)


deg = read_csv("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/cluster_IF-10%_gene_weights_median.csv")
head(deg)

deg$median = deg$median*(-1)
deg = as.data.frame(deg)

sel_sigs = deg[,c("median", "gene")]
sel_sigs = sel_sigs[order(sel_sigs$median, decreasing = T),]
head(sel_sigs)  



#BiocManager::install("clusterProfiler")
library(clusterProfiler)
organism = "org.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

library(msigdbr)
msigdbr_species()
co = msigdbr_collections()

gene_list = sel_sigs$median
names(gene_list) = sel_sigs$gene
head(gene_list)

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

em2 <- GSEA(gene_list, TERM2GENE = m_t2g)
head(em2)


library(ggplot2)
library(dplyr)
library(stringr)

## count the gene number
gene_count<- em2@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)

## merge with the original dataframe
dot_df<- left_join(em2@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)

dot_df$type = "upregulated"
dot_df$type[dot_df$NES < 0] = "downregulated"

dot_df = dot_df[dot_df$GeneRatio > 0.4 & dot_df$p.adjust < 0.05,]

# write.csv(dot_df, "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/gsea_plot_if10%.csv", row.names = F)


theme_lengend_bottom <- function(){ 
  
  theme_classic() +
    theme(
      legend.position = "bottom",
      #legend.title = element_blank(),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(-5,-0,0,-100),
      
      axis.text = element_text(colour = "black", size = 8),
      axis.title = element_text(colour = "black", size = 10),
      strip.text = element_text(colour = "black", size = 10),
      legend.text = element_text(colour = "black", size = 8),
      plot.title = element_text(size = 10)
      
    )}

# dot_df = read_csv("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/gsea_plot_cluster1.csv")
dot_df$ID = sub("HALLMARK_", "", dot_df$ID)
dot_df$Description = sub("HALLMARK_", "", dot_df$Description)

dot_df$Description = gsub("\\_", " ", dot_df$Description)
dot_df$Description = tolower(dot_df$Description)
dot_df = dot_df[order(dot_df$GeneRatio, decreasing = F),]
dot_df  
dot_df$Description = factor(dot_df$Description, levels = unique(dot_df$Description))

p <- ggplot(dot_df, aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  geom_vline(xintercept = 0.4, linetype = "dashed") +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +
  facet_grid(.~type) +
  ylab("HALLMARKS") +
  ggtitle("GSEA in cluster1\n(cisplatin less vs. more sensitive)") +
  theme_lengend_bottom() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 8, vjust = 0.75)) +
  guides(size = guide_legend(title.position="top", title.hjust = 0.5))

p
