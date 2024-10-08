# ******************************* module 1 *********************************************
# input: cell Ranger output matrix
# method: apply demultiplexing to keep the singlets, apply QC filtering and add metadata
# output: seurat object with counts and metadata
# **************************************************************************************

library(Matrix)
library(Seurat)

# packageurl <- "https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.1.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")


# 1) LOAD the UMI matrix (cell ranger output)
# give path to the folder with the matrix.mnt, barcodes and features matrix
data_dir = "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/cell_ranger_output/Rong_exp030_72_merged/filtered_feature_bc_matrix"
meta1 = readxl::read_excel("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/cell_ranger_output/cell_hashing_experiment_design_two_vehicles_exp009-016.xlsx")
name = "exp030"
pbmc <- Read10X(data_dir)

# Subset RNA and HTO counts by joint cell barcodes
pbmc.umis = pbmc$`Gene Expression`
pbmc.htos = pbmc$`Antibody Capture`

# Confirm that the HTO have the correct names
rownames(pbmc.htos)
pbmc.hto_row <- pbmc.htos[grep('row', rownames(pbmc.htos)), ]
pbmc.hto_col <- pbmc.htos[grep('column', rownames(pbmc.htos)), ]

# Setup Seurat object
pbmc.hashtag <- CreateSeuratObject(counts = pbmc.umis)



# 2) DEMULTIPLEX 
# Add HTO data as a new assay independent from RNA
pbmc.hashtag[["HTO_row"]] <- CreateAssayObject(counts = pbmc.hto_row)
pbmc.hashtag[["HTO_col"]] <- CreateAssayObject(counts = pbmc.hto_col)

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO_row", normalization.method = "CLR")
pbmc.hashtag <- NormalizeData(pbmc.hashtag, assay = "HTO_col", normalization.method = "CLR")

# demultiplex
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO_row", positive.quantile = 0.99, seed = 1001)
pbmc.hashtag <- HTODemux(pbmc.hashtag, assay = "HTO_col", positive.quantile = 0.99, seed = 1001)

# demultiplexing plots
p1 <- plot(HTOHeatmap(pbmc.hashtag, assay = "HTO_row"))
p2 <- plot(HTOHeatmap(pbmc.hashtag, assay = "HTO_col"))

# Keep the singlets
pbmc.hashtag$HTO_row_classification.global[1:10]
pbmc.hashtag$HTO_col_classification.global[1:10]
pbmc.singlet = subset(x = pbmc.hashtag, subset = HTO_row_classification.global == "Singlet" & 
                        HTO_col_classification.global == "Singlet")
total_cells = ncol(pbmc.hashtag)
total_singlets = ncol(pbmc.singlet)


# 3) QC ANALYSIS
as.data.frame(pbmc.singlet@assays$RNA$counts[1:10, 1:5])
head(pbmc.singlet@meta.data, 10)

# percentage of mitocondrial genes
pbmc.singlet <- PercentageFeatureSet(pbmc.singlet, "^MT-", col.name = "percent_mito")
# percentage of ribosomal genes
pbmc.singlet <- PercentageFeatureSet(pbmc.singlet, "^RP[SL]", col.name = "percent_ribo")
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
pbmc.singlet <- PercentageFeatureSet(pbmc.singlet, "^HB[^(P)]", col.name = "percent_hb")

# plot QC
# feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
# VlnPlot(pbmc.singlet, group.by = "HTO_row_classification", features = feats, pt.size = 0.1, ncol = 3, point_alpha = 0.1) +
#   NoLegend()
# FeatureScatter(pbmc.singlet, "nCount_RNA", "nFeature_RNA", 
#                group.by = "HTO_row_classification", pt.size = 0.5)

# calculate RNA assay QC metrics
# n_genes metrics
median_nFeat = summary(pbmc.singlet@meta.data$nFeature_RNA)[["Median"]]
mean_nFeat = summary(pbmc.singlet@meta.data$nFeature_RNA)[["Mean"]]
sd_nFeat = sd(pbmc.singlet@meta.data$nFeature_RNA)
mad_nFeat = mad(pbmc.singlet@meta.data$nFeature_RNA, constant = 1)

# n_counts metrics
median_nCount = summary(pbmc.singlet@meta.data$nCount_RNA)[["Median"]]
mean_nCount = summary(pbmc.singlet@meta.data$nCount_RNA)[["Mean"]]
sd_nCount = sd(pbmc.singlet@meta.data$nCount_RNA)
mad_nCount = mad(pbmc.singlet@meta.data$nCount_RNA, constant = 1)

# pct_mt genes metrics
median_pct_mt = summary(pbmc.singlet@meta.data$percent_mito)[["Median"]]
mean_pct_mt = summary(pbmc.singlet@meta.data$percent_mito)[["Mean"]]
mad_pct_mt = mad(pbmc.singlet@meta.data$percent_mito, constant = 1)

# filter out low quality cells
pbmc.singlet <- subset(
  pbmc.singlet, 
  subset = nFeature_RNA > median_nFeat-3*mad_nFeat &
    nFeature_RNA < mean_nFeat+3*sd_nFeat &
    nCount_RNA < mean_nCount+2*sd_nCount &
    nCount_RNA > median_nCount-2*mad_nCount &
    percent_mito < median_pct_mt+3*mad_pct_mt
)

dim(pbmc.singlet)
total_cells_afterQC = ncol(pbmc.singlet)

# plot filtered QC
# feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
# VlnPlot(pbmc.singlet, group.by = "HTO_row_classification", 
#         features = feats, pt.size = 0.1, ncol = 3) +  NoLegend()

# filter genes
dim(pbmc.singlet)

# Filter MALAT1
pbmc.singlet <- pbmc.singlet[!grepl("MALAT1", rownames(pbmc.singlet)), ]
# Filter Mitocondrial genes
pbmc.singlet <- pbmc.singlet[!grepl("^MT-", rownames(pbmc.singlet)), ]
# filter genes that are present in at least a certain amount of cells (0.1% of cells)
n = 0.001*ncol(pbmc.singlet)
selected_f <- rownames(pbmc.singlet)[Matrix::rowSums(pbmc.singlet) > n]
pbmc.singlet <- subset(pbmc.singlet, features = selected_f)
dim(pbmc.singlet)
total_num_genes = nrow(pbmc.singlet)


# 4) CELL CYCLE SCORE
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

pbmc.singlet <- NormalizeData(pbmc.singlet)
pbmc.singlet <- ScaleData(pbmc.singlet)
pbmc.singlet <- CellCycleScoring(pbmc.singlet, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
table(pbmc.singlet$Phase)


# 5) ADD METADATA (drug treatment information)
pbmc.singlet@meta.data[1:5,]
pbmc.singlet$rank = 1:ncol(pbmc.singlet)

meta = as.data.frame(pbmc.singlet@meta.data)
meta$name = paste0(meta$HTO_row_classification, "-", meta$HTO_col_classification)
meta$cell = rownames(meta)
meta1$name = paste0(meta1$row, "-", meta1$column)
table(unique(meta$name) %in% meta1$name)
meta = merge(meta, meta1, by.x = "name", by.y = "name", all = F)
meta = meta[order(meta$rank),]
head(meta)
table(colnames(pbmc.singlet) == meta$cell)
pbmc.singlet@meta.data = meta
pbmc.singlet@meta.data[1:5,]


# 6) SAVE DATA
# save seurat object
#pbmc.singlet$plate = name
pbmc.singlet
name
saveRDS(pbmc.singlet, file = paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module1_demultiplex.QC/seurat_obj/", name, ".rds"))
dim(sce.filtered)

# save qc matrix
df = data.frame(feat = c("mean_nCount", "median_nCount", "mad_nCount", "sd_nCount", 
                         "mean_nFeat", "median_nFeat", "mad_nFeat", "sd_nFeat",
                         "mean_pct_mt", "median_pct_mt", "mad_pct_mt", 
                         "total_cells", "total_singlets", "total_cells_afterQC", "total_num_genes"),
                value = c(mean_nCount, median_nCount, mad_nCount, sd_nCount, 
                          mean_nFeat, median_nFeat, mad_nFeat, sd_nFeat,
                          mean_pct_mt, median_pct_mt, mad_pct_mt, 
                          total_cells, total_singlets, total_cells_afterQC, total_num_genes), 
                stringsAsFactors = F)
head(df)
write.csv(df, paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module1_demultiplex.QC/qc_plots/summary_QC_", name, ".csv"), row.names = F)

# save figs
head(meta)
freq = as.data.frame(table(meta$row, meta$column))
freq = tidyr::spread(freq, Var2, Freq)
freq = freq[,c("Var1", paste0("column", 1:12))]
rownames(freq) = freq$Var1
mat = freq[,2:ncol(freq)]

library(ComplexHeatmap)
library(grid)
library(cowplot)
ht = grid.grabExpr(draw(Heatmap(mat, name = "# cells",cluster_columns = F, cluster_rows = F, row_names_side = "left",
                                cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                                  grid.text(mat[i, j], x, y)
                                })))

a = plot_grid(p1, p2, nrow = 1)
p = plot_grid(a,ht, nrow = 2, rel_heights = c(1,1.5))
save_plot(p, filename = paste0("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module1_demultiplex.QC/qc_plots/plot_", name, ".png"),
          base_height = 9, base_width = 9)
