
library(Seurat)
library(readr)
library(dplyr)

# drug signature gene weights
# select 1000 top vs 10k bottom genes per signature 
sig1000 = read_csv("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/cohend_pergene/classify_top_1k_and_bottom_10k_genes.csv") 
sig1000 = sig1000[,c("gene", "type")]
head(sig1000)
table(sig1000$type)

separate_resistant_cells = function(expr, name){
  drug = "A"
  counts = as.data.frame(t(as.matrix(expr@assays$RNA$counts)))
  counts_top = counts[,colnames(counts) %in% sig1000$gene[sig1000$type == "top_1000"]]
  counts_bottom = counts[,colnames(counts) %in% sig1000$gene[sig1000$type == "bottom_1000"]]
  
  meta = expr@meta.data
  meta = meta[,c("cell", "Treatment", "Conc", "Time")]
  meta$sum_top = rowSums(counts_top)
  meta$sum_bottom = rowSums(counts_bottom)
  meta$ratio = (meta$sum_top)/(meta$sum_bottom)
  
  #write.csv(meta, "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/exp002R_cells_assignment_for_plot.csv", row.names = F)
  # library(ggplot2)
  # ggplot(meta[meta$Treatment %in% c("A", "Ctrl"),], 
  #        aes(x = Conc, y = log(ratio), col = Time)) +
  #    geom_boxplot() + 
  #   geom_point(alpha = 0.05, position = position_jitterdodge()) +
  #    facet_grid(cols = vars(Treatment), scales = "free", space = "free")
  
  meta = meta[meta$Treatment == drug & meta$Conc == "200%" & meta$Time == "24h",]
  meta$quintile = dplyr::ntile(meta$ratio, 4)
  meta$class = "rest"
  meta$class[meta$quintile == 4] = "more_sensitive"
  meta$class[meta$quintile == 1] = "less_sensitive"
  meta = meta[meta$class != "rest",]
  return(meta[,c("cell", "class")])
}

# calculate the resistant signatures  
# read scRNA-seq data 
list_plates = c("exp001R", "exp002R", "exp003", "exp004", "exp005", "exp012_26", 
                "exp014", "exp015", "exp025", "exp028", "exp029", "exp030")


lapply(list_plates, function(name){
  print(name)
  
  expr_ccle <- readRDS(paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module1_demultiplex.QC/seurat_obj/", name, ".rds"))
  cells = separate_resistant_cells(expr = expr_ccle, name = name)
  head(cells)
  
  # filter:select control & treated cells
  meta = expr_ccle@meta.data
  meta = meta[meta$Treatment == "Ctrl" | meta$Conc == "0%" | meta$cell %in% cells$cell,]
  meta = merge(meta[,c("Treatment", "Conc", "cell", "Phase")], cells, by = "cell", all = T)
  meta$class[is.na(meta$class)] = "control"
  meta$class_general = ifelse(meta$class == "control", "control", "treated")
  table(meta$class, meta$class_general)
  
  
  # PCA ********************************************
  cpm_mat = as.matrix(expr_ccle@assays$RNA$data)
  cpm_mat = cpm_mat[,meta$cell[meta$class_general == "control"]]
  CP100K_log <- log2((cpm_mat) + 1)
  CP100K_log <- CP100K_log - rowMeans(CP100K_log)
  CP100K_log[1:5,1:5]
  newmat = as.data.frame(t(CP100K_log))
  newmat$cell = rownames(newmat)
  newmat[1:5,1:5]
  
  newmat = merge(meta, newmat, by = "cell", all = F)
  newmat[1:5,1:10]
  table(newmat$class_general)  
  
  # PCA in control cells
  m1 = as.matrix(newmat[,7:ncol(newmat)])
  
  set.seed(100)
  pc = prcomp(m1, scale. = F, center = F)
  df_pc = cbind(newmat[,1:6], pc$x[,1:20])
  
  # Project treated data onto the PC coordinates system
  cpm_matB = as.matrix(expr_ccle@assays$RNA$data)
  cpm_matB = cpm_matB[,meta$cell[meta$class_general == "treated"]]
  CP100K_logB <- log2((cpm_matB) + 1)
  CP100K_logB <- CP100K_logB - rowMeans(CP100K_logB)
  CP100K_logB[1:5,1:5]
  newmatB = as.data.frame(t(CP100K_logB))
  newmatB$cell = rownames(newmatB)
  newmatB[1:5,1:5]
  
  newmatB = merge(meta, newmatB, by = "cell", all = F)
  newmatB[1:5,1:10]
  table(newmatB$class_general)
  mB1 = as.matrix(newmatB[,7:ncol(newmatB)])
  pro = predict(object = pc, mB1)
  df_pcB = cbind(newmatB[,1:6], pro[,1:20])
  
  # merge df_pc and projected data
  df_pc_comb = rbind(df_pc, df_pcB)
  
  # get gene PC weights (rotation matrix)
  w = as.data.frame(pc$rotation[,1:20])
  w$gene = rownames(w)
  
  experim = paste0(name, "_A.csv")
  write.csv(w, paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/pca_resistant_10k/setseed_rot__",experim), row.names = F)
  write.csv(df_pc_comb, paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/pca_resistant_10k/setseed_dfpc__",experim), row.names = F)
  
  
})



############################################################################
# plot resistant versus sensitive for 1 example

head(meta)

r = meta[meta$Treatment %in% c("A", "Ctrl"),]
r$class = "rest"
r$class[r$cell %in% meta$cell[meta$class == "resistant"]] = "resistant"
r$class[r$cell %in% meta$cell[meta$class == "sensitive"]] = "sensitive"

ggplot(r, aes(x = Conc, y = log(ratio))) +
  geom_point(aes(col = class, alpha = class), position = position_jitter()) +
  geom_boxplot(outlier.alpha = 0, alpha = 0.5) + 
  scale_alpha_manual(values = c(0.6,0.6, 0.6)) +
  scale_color_manual(values = c("tomato", "darkgray", "#5ab3d6")) +
  ylab("log ratio(counts_top100/counts_bottom100)\nCISPLATIN marker genes") +
  facet_grid(cols = vars(Treatment), scales = "free", space = "free") +
  theme_classic() +
  theme(legend.position = "top",
        legend.title = element_blank())
