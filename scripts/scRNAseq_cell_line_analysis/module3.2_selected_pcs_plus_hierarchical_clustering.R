
# 1) select the PCs that separate more vs less sensitive & discard PCs that separate cell cycle phases

files_pc = list.files(path = "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/pca_resistant_5%", pattern = "dfpc", full.names = T)
# files_pc = files_pc[grepl("setseed", files_pc)]
files_pc = files_pc[!grepl("exp025", files_pc)]
files_pc = files_pc[!grepl("exp028", files_pc)]

library(effsize)

res1 = dplyr::bind_rows(lapply(files_pc, function(f){
  df=read_csv(f)
  n = sub(".*resistant_.*%\\/.*dfpc\\_\\_(exp.*)\\_([A-Z0-9\\-]*)\\.csv", "\\1", f)
  d = sub(".*resistant_.*%\\/.*dfpc\\_\\_(exp.*)\\_([A-Z0-9\\-]*)\\.csv", "\\2", f)
  
  test_bypc = dplyr::bind_rows(lapply(paste0("PC", 1:10), function(p){
    red_df = df[,c(colnames(df)[1:6], p)]
    colnames(red_df)[7] = "PCx"
    
    c_general = cohen.d(red_df$PCx, f = factor(red_df$class_general, levels = c("treated", "control")))
    c_general = c_general$estimate
    
    r = red_df[red_df$class_general != "control",]
    c_resi = cohen.d(r$PCx, f = factor(r$class, levels = c("more_sensitive","less_sensitive")))
    c_resi = c_resi$estimate
    
    table(red_df$Phase)
    r1 = red_df[red_df$Phase %in% c("G1", "G2M"),]
    c_1 = cohen.d(r1$PCx, f = factor(r1$Phase, levels = c("G1","G2M")))
    c_1 = c_1$estimate
    
    r2 = red_df[red_df$Phase %in% c("G1", "S"),]
    c_2 = cohen.d(r2$PCx, f = factor(r2$Phase, levels = c("G1","S")))
    c_2 = c_2$estimate
    
    r3 = red_df[red_df$Phase %in% c("S", "G2M"),]
    c_3 = cohen.d(r3$PCx, f = factor(r3$Phase, levels = c("S","G2M")))
    c_3 = c_3$estimate
    
    # library(ggplot2)
    # ggplot(red_df, aes(x = class, y = PCx)) +
    #    geom_boxplot()
    # ggplot(red_df, aes(x = Phase, y = PCx)) +
    #   geom_boxplot()
    # ggplot(red_df, aes(x = class_general, y = PCx)) +
    #     geom_boxplot()
    
    return(data.frame(PC = p, c_general= c_general, c_resi = c_resi, 
                      c_G1_G2M = c_1, c_G1_S = c_2, c_S_G2M = c_3, stringsAsFactors = F))
  }))
  
  test_bypc$plate = n
  test_bypc$drug = d
  test_bypc$name= paste0(n, "__", d)
  return(test_bypc)
}))


selected = res1
head(selected)
selected$fin_name = paste0(selected$name, "__", selected$PC)

library(ggplot2)
ggplot(selected, aes(x = c_resi)) +
  geom_histogram(alpha = 0.5, col = "steelblue", fill = "steelblue", bins = 50) +
  geom_vline(xintercept = c(-0.5,0.5), linetype = "dashed", alpha = 0.5) +
  xlab("effect size sensitive vs resistant cells") +
  theme_classic()

summary(res1$c_general)
summary(res1$c_resi)

# select the PCs that separate more vs less sensitive 
head(res)
newres1 = res[,selected$fin_name[selected$c_resi>= 0.5]]
newres2 = res[,selected$fin_name[selected$c_resi<= -0.5]]
newres2 = -1*newres2
newres = cbind(newres1, newres2)
co = cor(newres)
dim(co)

library(ComplexHeatmap)
Heatmap(co, show_row_names = T, show_column_names = T, 
        row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6))

# discard PCs that separate cell cycle phases
head(newres)
dim(newres)
newres = newres[, ! colnames(newres) %in% c(selected$fin_name[abs(selected$c_G1_G2M) >= 0.5]) ]
newres = newres[, ! colnames(newres) %in% c(selected$fin_name[abs(selected$c_G1_S) >= 0.5]) ]
newres = newres[, ! colnames(newres) %in% c(selected$fin_name[abs(selected$c_S_G2M) >= 0.5]) ]


# 2) calculate similarity between selected PCs
nmf_intersect <- cor(newres)
dim(nmf_intersect)

saveRDS(nmf_intersect, "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/hc_pc_resistant_5%.rds")

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("white", "#d9d9d9", "yellow", "orange", "red", "maroon"))

Heatmap(nmf_intersect, col = col_fun, 
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8))

dim(nmf_intersect)

# based on the heatmap select PCs in the same cluster:
# IFN
sel_sigs1 = c("exp005__A__PC8", "exp003__A__PC9", "exp029__A__PC7", 
              "exp004__A__PC9") # 1k vs 10k

sel_sigs1 = c("exp005__A__PC8", "exp003__A__PC9", "exp029__A__PC7", 
              "exp004__A__PC9", "exp001R__A__PC7") # 1%

# EMT
sel_sigs1 = c("exp030__A__PC9", "exp014__A__PC9", "exp015__A__PC8",
              "exp012_26__A__PC8", "exp005__A__PC10")

sel_sigs1 = c("exp001R__A__PC9", "exp004__A__PC10", "exp003__A__PC10")


# save the median across the PCs in the same cluster
library(matrixStats)
sel_sigs1
sel_sigs = as.data.frame(newres[,sel_sigs1])

sel_sigs$median = rowMedians(as.matrix(sel_sigs))
sel_sigs$gene = rownames(sel_sigs)
sel_sigs = sel_sigs[,c("median", "gene")]
sel_sigs = sel_sigs[order(sel_sigs$median, decreasing = T),]
head(sel_sigs)  

write.csv(sel_sigs, "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/cluster_IF-10%_gene_weights_median.csv", row.names = F)



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

# m_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
#  dplyr::select(gs_name, gene_symbol)
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol)
head(m_t2g)

# em <- enricher(gene, TERM2GENE=m_t2g)
# head(em)

em2 <- GSEA(gene_list, TERM2GENE = m_t2g)
head(em2)
em2[1:10,1:5]


require(DOSE)
dotplot(em2, showCategory=10, split=".sign") + facet_grid(.~.sign) +
  theme(axis.text.y = element_text(size = 6))


y <- as.data.frame(em2)
head(y)
y = y[y$enrichmentScore > 0,]
y$Description = sub("HALLMARK\\_", "", y$Description)

colnames(y)
y = y[1:10,c("Description", "setSize", "enrichmentScore", "pvalue","p.adjust")]
rownames(y) = y$Description
y

sentSplit <- function(string, tolerance = 0.5, collapse = "\n") {
  d = paste(strwrap(strsplit(string, "_"), width = 40), collapse = collapse)
  d = sub("^c\\(", "", d)
  d = sub("\\)", "", d)
  d = gsub("\\\"", "", d)
  d = gsub(",", "",d)
  return(d)
}

for (i in 1:nrow(y)){
  y$Description[i] = sentSplit(as.character(y$Description[i]))
}
y$Description
y = y[1:10,]
y = y[order(y$enrichmentScore, decreasing = F),]
y$Description = factor(y$Description, levels = y$Description)

ggplot(y, aes(x = enrichmentScore, y = Description)) + 
  geom_point(aes(size = setSize, color = p.adjust)) +
  scale_colour_gradient(limits=c(0, 0.10), low="red") +
  ylab(NULL) +
  ggtitle("GSEA HALLMARK") +
  theme_classic() +
  theme(axis.text = element_text(size = 8))


####################################################################

# Daria decoupler

# install.packages("BiocManager")
# BiocManager::install("decoupleR")
# BiocManager::install("OmnipathR")

library(decoupleR)

net <- get_progeny(organism = 'human', top = 500)
net


library(matrixStats)
sel_sigs1
sel_sigs = as.data.frame(newres[,sel_sigs1])

sel_sigs$median = rowMedians(as.matrix(sel_sigs))
sel_sigs$gene = rownames(sel_sigs)
sel_sigs = sel_sigs[,c("median", "gene")]
sel_sigs = sel_sigs[order(sel_sigs$median, decreasing = T),]
head(sel_sigs)  

head(sel_sigs)
deg = sel_sigs
deg$median = -1*deg$median
deg$gene = NULL
colnames(deg)[1] = "t"
head(deg)


# Run mlm
contrast_acts <- run_mlm(mat=deg, net=net, .source='source', .target='target',
                         .mor='weight', minsize = 5)
contrast_acts

# write.csv(contrast_acts, "results/contrast_cluster_2.csv", row.names = F)
# write.csv(sel_sigs, "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/cluster1_gene_weights_median.csv", row.names = F)

# Plot
ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways")
