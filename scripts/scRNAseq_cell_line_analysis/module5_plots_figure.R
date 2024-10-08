
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(readr)



#install.packages("extrafont")
library("extrafont")
font_import()
#dev.new(family = "Arial")

theme_lengend_bottom <- function(){ 
  
  theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.margin = margin(0,0,0,0),
      legend.box.margin = margin(-5,0,0,0),
      
      axis.text = element_text(colour = "black", size = 8),
      axis.title = element_text(colour = "black", size = 10),
      legend.text = element_text(colour = "black", size = 8),
      plot.title = element_text(size = 10)
      
    )}

names = data.frame(name = c("exp001R",     "exp002R", "exp003",    "exp004",   "exp005",   "exp012_26",        "exp014",        "exp015",   "exp025", "exp028", "exp029", "exp030"), 
                   newname = c("COV362_R1","Kuramochi_R1","COV362_R2","COV362_R3","COV362_R4", "Kuramochi_R5", "Kuramochi_R2", "Kuramochi_R3", "COV318", "GP5d", "COV362_R5", "Kuramochi_R4" ))
names



# fig D
meta = read_csv("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/exp001R_cells_assignment_for_plot.csv")
head(meta)

sub = meta[meta$Treatment == "A" & meta$Conc == "200%" & meta$Time == "24h",]
sub$quintile = dplyr::ntile(sub$ratio, 4)

meta$class = "rest"
meta$class[ meta$cell %in% sub$cell[sub$quintile == 4] ] = "more_sensitive"
meta$class[ meta$cell %in% sub$cell[sub$quintile == 1] ] = "less_sensitive"

meta$Conc = factor(meta$Conc, levels = c("Ctrl", "0%", "10%", "50%", "200%"))

df_text = data.frame(x = c("200%", "200%"), y = c(-0.8,-1.8), 
                     label = c("25% more\nsensitive", "25% less\nsensitive"),
                     col = c("more_sensitive", "less_sensitive"), 
                     vjust="inward",hjust="inward")#hjust = c(0.5, 0.5))

figd2 = ggplot(meta[meta$Treatment %in% c("A", "Ctrl") & meta$Time %in% c("24h", "Ctrl"),], 
               aes(x = Conc, y = log(ratio))) +
  geom_point(aes(col = class), size = 1, position = position_jitter(), alpha = 0.75,) +
  geom_boxplot(alpha = 0.5) + 
  geom_label(data = df_text, aes(x = x, y = y, label = label, 
                                 col = col, hjust = hjust), 
             lineheight = 0.75, show.legend = FALSE, size = 3) +
  scale_color_manual(values = c("#fc8d59", "#91cf60", "gray")) +
  ggtitle("COV362 R1") +
  ylab("log ratio (counts in cisplatin gene\n markers vs in background)") +
  xlab("cisplatin concentration") +
  theme_lengend_bottom() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.spacing.x = unit(0.05, 'cm'),
        legend.title = element_blank(),
        plot.margin = margin(5,5,5,5),
        legend.position = c(0.27,0.85))
figd2

# ******************************************************************************


df_pc1 = read_csv("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/pca_resistant_10k/dfpc__exp003_A.csv")

df1 = df_pc1
df1$label = "B"
df1$group = df1$class

df2 = df_pc1
df2$label = "A"
df2$group = df2$class_general

df3 = df_pc1
df3$label = "C"
df3$group = df3$Phase

df = dplyr::bind_rows(list(df1,df2,df3))
head(df)


fige1 = ggplot(df, aes(x = group, y = PC9, col = group)) +
  geom_point(alpha = 0.05, position = position_jitter(width = 0.25)) +
  geom_boxplot(alpha= 0.5)  +
  facet_grid(cols = vars(label), scales = "free", space = "free") +
  scale_color_manual(values = c("gray", "gray", "gray", "#92c5de", "#fc8d59", "#91cf60",  "tomato") ) +
  ylab("COV362_R2 PC9 activities") +
  ggtitle("(1) separate more vs.\nless sensitive cells") +
  theme_lengend_bottom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))


fige2 = ggplot(df, aes(x = group, y = PC2, col = group)) +
  geom_point(alpha = 0.05, position = position_jitter(width = 0.25)) +
  geom_boxplot(alpha= 0.5)  +
  facet_grid(cols = vars(label), scales = "free", space = "free") +
  scale_color_manual(values = c("gray", "gray", "gray", "#92c5de", "#fc8d59", "#91cf60",  "tomato") ) +
  ylab("COV362_R2 PC2 activities") +
  ggtitle("(2) separate phases\nof the cell cycle") +
  theme_lengend_bottom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
#fige2

fige3 = ggplot(df, aes(x = group, y = PC6, col = group)) +
  geom_point(alpha = 0.05, position = position_jitter(width = 0.25)) +
  geom_boxplot(alpha= 0.5)  +
  facet_grid(cols = vars(label), scales = "free", space = "free") +
  scale_color_manual(values = c("gray", "gray", "gray", "#92c5de", "#fc8d59", "#91cf60",  "tomato") ) +
  ylab("COV362_R2 PC6 activities") +
  ggtitle("(3) do not separate\nour categories") +
  theme_lengend_bottom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
#fige3


fige = plot_grid(fige1, fige2, fige3, nrow = 1)
fige
fig_empty = plot.new() 

down = plot_grid(fig_empty, fige, nrow = 1, rel_widths = c(1,1.5), labels = c("D", ""), 
                 label_size = 12)
down

# ******************************************************************************


nmf_intersect = readRDS("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/hc_pc_resistant10k.rds")

rownames(nmf_intersect) = sub("__A__", "--", rownames(nmf_intersect))
colnames(nmf_intersect) = sub("__A__", "--", colnames(nmf_intersect))

names
for(i in 1:12){
  sub = names[i,]
  rownames(nmf_intersect) = sub(sub$name, sub$newname, rownames(nmf_intersect))
  colnames(nmf_intersect) = sub(sub$name, sub$newname, colnames(nmf_intersect))
  
}

nmf_intersect = nmf_intersect[!grepl("COV318", rownames(nmf_intersect)),]
nmf_intersect = nmf_intersect[,!grepl("COV318", colnames(nmf_intersect))]

#col_fun = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("white", "#d9d9d9", "yellow", "orange", "red", "maroon"))
col_fun = colorRamp2(c(0, 0.2, 0.4, 0.6, 0.8, 1), c("white","#edf8e9","#bae4b3","#31a354", "#006d2c", "#00441b"))

df_class = data.frame(cell_line = sub("([A-Za-z0-9]*)\\_.*", "\\1", rownames(nmf_intersect)))
df_class

set.seed(2)

ht_f = grid.grabExpr(draw(
  Heatmap(nmf_intersect, col = col_fun, 
          column_title = "gene expression patterns of inter-individual\nvariability in non-treated cells agreement across replicates",
          column_title_gp = gpar(fontsize = 10),
          show_row_names = F,
          show_column_names = F,
          name = "Pearson correlation",
          
          row_names_gp = gpar(fontsize = 8),
          column_names_gp = gpar(fontsize = 8),
          column_names_rot = 45,
          
          heatmap_legend_param = list(
            legend_direction = "horizontal", 
            legend_width = unit(4, "cm"),
            title_gp = gpar(fontsize = 8),
            labels_gp = gpar(fontsize = 8))
          
  ) +
    Heatmap(df_class, show_column_names = F,
            name = "cell line",
            heatmap_legend_param = list(
              title_gp = gpar(fontsize = 8),
              labels_gp = gpar(fontsize = 8))
    )
  , 
  show_heatmap_legend = T, 
  heatmap_legend_side = "bottom") )#grid.grab

plot_grid(ht_f)


##################################################################################


contrast_acts1 = read_csv("results/contrast_cluster_ifn.csv")
contrast_acts1


fig_gsea1 = ggplot(contrast_acts1, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) +
  ylab("pathway activity") +
  ggtitle("less sensitive vs. more sensitive cells\ngene expression signature cluster1") +
  theme_lengend_bottom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

fig_gsea1  


contrast_acts2 = read_csv("results/contrast_cluster_emt.csv")
contrast_acts2


fig_gsea2 = ggplot(contrast_acts2, aes(x = reorder(source, score), y = score)) + 
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) +
  ylab("pathway activity") +
  ggtitle("less sensitive vs. more sensitive cells\ngene expression signature cluster2") +
  theme_lengend_bottom() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

fig_gsea2


fig_g = plot_grid(fig_gsea1, fig_gsea2, nrow = 2)
fig_g

##################################################################################




up = plot_grid(fig_empty, figd2, nrow = 1, labels = c("", "C"), label_size = 12,
               rel_widths = c(2,1))
up

down = plot_grid(fig_empty, fige, nrow = 1, rel_widths = c(1,1.65), labels = c("D", ""), 
                 label_size = 12)
down

fig_g = plot_grid(fig_gsea1, fig_empty, nrow = 2, labels = c("G", "H"), label_size = 12)
fig_g

set.seed(2)
fig_down2 = plot_grid(ht_f, fig_g, labels = c("F", ""), label_size = 12,nrow = 1)
fig_down2


plot_grid(up, down, fig_down2, nrow = 3,
          rel_heights = c(1,1,1.75))


pdf("../Fig_Rong/prueba_A4_arial.pdf", height = 11.69, width = 8.27,
    family = "ArialMT", 
    paper = "special", onefile = FALSE)

plot_grid(up, down, fig_down2, nrow = 3,
          rel_heights = c(1,1,1.75))

dev.off()

