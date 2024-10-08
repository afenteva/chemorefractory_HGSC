
library(Seurat)
library(parallel)


# read scRNA-seq data 

list_exp = c("exp001R", "exp002R", "exp003", "exp004", "exp005", "exp012_26", 
             "exp014", "exp015", "exp025", "exp028", "exp029", "exp030")



lapply(list_exp, function(name){
  
  print(name)
  drug = "A"
  expr_ccle <- readRDS(paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module1_demultiplex.QC/seurat_obj/", name, ".rds"))
  expr_ccle@assays$RNA$counts[1:5,1:5]
  expr_ccle@assays$RNA$scale.data[1:5,1:5]
  expr_ccle@assays$RNA$data[1:5,1:5]
  class(expr_ccle)
  
  # filter
  meta = expr_ccle@meta.data
  table(meta$Treatment)
  meta = meta[meta$Treatment %in% c(drug, "Ctrl"),]
  dim(meta)
  table(meta$Treatment)
  meta$group = "control"
  meta$group[meta$Treatment == "A" & meta$Conc %in% c("10%", "50%", "200%")] = "drugA"
  table(meta$group, meta$Conc)
  
  cpm_mat = as.matrix(expr_ccle@assays$RNA$scale.data)
  cpm_mat = cpm_mat[,meta$cell]
  dim(cpm_mat)
  cpm_mat[1:5,1:10]
  hist(as.matrix(cpm_mat))
  
  cpm_mat = cpm_mat[rowSums(cpm_mat) != 0,]
  
  mat = as.data.frame(t(cpm_mat))
  mat$cell = rownames(mat)
  mat = merge(meta[,c("cell", "group")], mat, by = "cell", all = F)
  mat[1:5,1:10]
  
  list_genes = colnames(mat[,3:ncol(mat)])
  df_fin = dplyr::bind_rows(lapply(list_genes, function(g){
    sub = mat[,c("group", g)]
    colnames(sub)[2] = "ge"
    cohendA = cohen.d(d = sub$ge, f = factor(sub$group, levels = c("drugA","control")))
    # ggplot(sub, aes(x = group, y = ge)) +
    #   geom_boxplot()
    return(data.frame(gene = g, cohend = cohendA$estimate, stringsAsFactors = F))
  }))
  
  df_fin$name = name
  
  # save output
  write.csv(df_fin, paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/cohend_pergene/cohend_", name, "_", drug, ".csv"), row.names = F)
  
  remove(expr_ccle)
  remove(meta)
  remove(cpm_mat)
  gc()
  
})# lapply

