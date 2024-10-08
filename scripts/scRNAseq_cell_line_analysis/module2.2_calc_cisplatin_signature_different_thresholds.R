


files = list.files("/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/cohend_pergene", "cohend", full.names = T)
files

fin = dplyr::bind_rows(lapply(files, function(f){
  df = read_csv(f)
  df = df[!is.na(df$cohend),]
  return(df)
}))

head(fin)
newfin = tidyr::spread(fin, name, cohend)
head(newfin)
na_sum = rowSums(is.na(newfin))
table(na_sum == 0)

newfin = newfin[na_sum == 0,]
head(newfin)

co = cor(newfin[,2:ncol(newfin)])
co = round(co,2)

library(circlize)
col_fun = colorRamp2(c(0.25, 0.5, 0.75, 1), c("#edf8e9","#bae4b3","#31a354", "#006d2c"))

library(ComplexHeatmap)
Heatmap(co, col = col_fun,
        cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
          grid.text(co[i, j], x, y)
        })



head(newfin)
newfin$mean = rowMeans(newfin[,2:ncol(newfin)])
head(newfin)
newfin = newfin[order(newfin$mean, decreasing = F),]
head(newfin)
newfin$type = "bottom"
# newfin$type[1:1000] = "bottom_1000"
newfin = newfin[order(newfin$mean, decreasing = T),]
head(newfin)
newfin$type[1:round(0.1*nrow(newfin))] = "top_10"
table(newfin$type)
# newfin = newfin[newfin$type!= "rest",]
write.csv(newfin, "/g/strcombio/fsupek_decider/msalvadores/KI-Rong/module_fig/cohend_pergene/classify_top_10%.csv", row.names = F)
