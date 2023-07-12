## Tengjia Jiang, 11.07.2023
## complexheatmap plot parameters
library(ggsci)
library(magrittr)
library(wesanderson)
library(ComplexHeatmap)
library(maftools)
library(scales)
library(dplyr)

workdir <- "./"
laml_clin <- readRDS(paste0(workdir,"revised_input_star/laml_clin.Rds"))

## maftools-----------------
## maftools plot mutation color 
vc_cols <- pal_npg("nrc")(8)
names(vc_cols) = c(
  'Frame_Shift_Del',
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Multi_Hit',
  'Frame_Shift_Ins',
  'In_Frame_Ins',
  'Splice_Site',
  'In_Frame_Del'
)
print(vc_cols)
show_col(vc_cols)

## maftools group color 
Location_col <- structure(c("#AB4B52","#1F78B4","#273046"), names=c("DT","PX","third"))
Location_col <- list(Location = Location_col)

Primary_loc_col <- structure(c("#899DA4","#FAD510"), names=c("Left","Right"))
Primary_loc_col <- list(Primary_loc = Primary_loc_col)

## ComplexHeatmap----------
## ComplexHeatmap oncoplot mutation color
col_own <- structure(c(pal_npg("nrc")(8),"#812321","#E5B0B0","#003987","#BACEE9"), 
                     names = c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation','Multi_Hit',
                               'Frame_Shift_Ins','In_Frame_Ins','Splice_Site','In_Frame_Del',
                               'Amplification', 'Gain', 'Deletion', 'Loss'))
print(col_own)

## ComplexHeatmap group color
res_col <- structure(c("#AB4B52","#1F78B4","#273046"), names=c("DT","PX","third"))
loc_col <- structure(c("#899DA4","#FAD510"), names=c("Left","Right"))
annotationColor_own <- list(Location = res_col, Primary_loc = loc_col)

clin_CH <- laml_clin 

ha_bl <- HeatmapAnnotation(
  df = clin_CH[,"Location", drop=F],
  col = list(Location = res_col,
             Site = loc_col)
)
alter_fun_own = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "#FFFFFF", col = "#CCCCCC"))
  },
  Frame_Shift_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Frame_Shift_Del"], col = NA))
  },
  Missense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Missense_Mutation"], col = NA))
  },
  Nonsense_Mutation = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Nonsense_Mutation"], col = NA))
  },
  Multi_Hit = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Multi_Hit"], col = NA))
  },
  Frame_Shift_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Frame_Shift_Ins"], col = NA))
  },
  In_Frame_Ins = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["In_Frame_Ins"], col = NA))
  },
  Splice_Site = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Splice_Site"], col = NA))
  },
  In_Frame_Del = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["In_Frame_Del"], col = NA))
  },
  Amplification = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Amplification"], col = NA))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Gain"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Deletion"], col = NA))
  },
  Loss = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_own["Loss"], col = NA))
  }
)
## legend
heatmap_legend_param_own = list(title = "Alternations",
                                at = c('Frame_Shift_Del', 'Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins', 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Amplification', 'Gain', 'Deletion', 'Loss'), 
                                labels = c('Frame_Shift_Del', 'Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins', 'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del','Amplification', 'Gain', 'Deletion', 'Loss'))

