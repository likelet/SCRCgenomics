## Tengjia Jiang, 11.07.2023
## CNV landscape of top 30 mutant genes

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(maftools)
library(dplyr)
library(tidyr)
library(magrittr)
library(reshape2)
library(ggthemr)
library(gg.gap)

workdir <- "./"
laml_clin <- readRDS(paste0(workdir,"revised_input_star/laml_clin.Rds"))
source(paste0(workdir,"revised_scripts_star/complexheatmap_parameters.R"))
mat_processed_cnv <- readRDS(paste0(workdir,"revised_input_star/mat_processed_cnv.Rds"))
mat_processed_cnv[mat_processed_cnv %in% c("nonCNV")] <- ""

altered_nums_cnv <- apply(mat_processed_cnv, 1, function(x){length(x)-sum(x=="")})
drivergene <- read.delim(paste0(workdir,"input/DriverDBv3_COAD.txt"),header = T,sep = "\t",stringsAsFactors = F)
select_cnvgene <- laml_CNV[laml_CNV$Gene %in% drivergene$gene, c("Gene","Hugo_Symbol","Pathway","OG_TSG")] %>% unique()
select_gene_freq <- as.data.frame(altered_nums_cnv)
select_gene_freq$Hugo_Symbol <- rownames(select_gene_freq)
select_cnvgene <- merge(select_cnvgene, select_gene_freq, by = "Hugo_Symbol",all.x = T)
select_cnvgene <- select_cnvgene %>% dplyr::filter(altered_nums_cnv/51 > 0.05)
select_cnvgene <- select_cnvgene[order(select_cnvgene$OG_TSG, -select_cnvgene$altered_nums_cnv,decreasing = F),]
mat_processed_cnv_driver <- mat_processed_cnv[select_cnvgene$Hugo_Symbol,]

## oncoplot------------
pct_num_cnv = rowSums(apply(mat_processed_cnv_driver, 1:2, function(x)any(x!="")))/ncol(mat_processed_cnv_driver)
pct_cnv = paste0(round(pct_num_cnv * 100, digits = 0), "%")
gp_fontsize <- 15

## mutation type fraction
CnvtypeBarplot <- function(margin_value,classified){
  cnv_type <- apply(mat_processed_cnv_driver, margin_value, function(x){as.data.frame(table(x))})
  cnv_type <- do.call(rbind,cnv_type)
  cnv_type$Hugo_Sample <- rownames(cnv_type)
  colnames(cnv_type)[1] <- "cnv_type"
  cnv_type$cnv_type <- as.character(cnv_type$cnv_type)
  cnv_type$Freq <- as.numeric(cnv_type$Freq)
  cnv_type <- cnv_type[cnv_type$cnv_type != "",]
  cnv_type$Hugo_Sample <- sub("[.]\\d$","",cnv_type$Hugo_Sample)
  cnv_type <- reshape2::dcast(Hugo_Sample ~ cnv_type,data = cnv_type,value.var = "Freq")
  rownames(cnv_type) <- cnv_type$Hugo_Sample
  cnv_type[is.na(cnv_type)] <- 0
  if(classified == 2){
    barplot_order <- c('Amplification', 'Deletion')
  }else if(classified == 4){
    barplot_order <- c('Amplification', 'Gain', 'Deletion', 'Loss')
  }
  if(margin_value == 1){
    cnv_type <- cnv_type[rownames(mat_processed_cnv_driver),-1]
  }else{
    cnv_type <- cnv_type[,-1]
    nucnv_sample <- setdiff(colnames(mat_processed_cnv_driver),rownames(cnv_type))
    noncnv_mat <- matrix(0,nrow=length(nucnv_sample), ncol=classified)
    colnames(noncnv_mat) <- barplot_order
    rownames(noncnv_mat) <- nucnv_sample
    cnv_type <- rbind(cnv_type, noncnv_mat) %>% as.matrix()
    cnv_type <- cnv_type[colnames(mat_processed_cnv_driver),]
  }
  cnv_type <- cnv_type[,barplot_order]
}
cnv_type_perSample <- CnvtypeBarplot(margin_value = 2, classified = 4)
cnv_type_perGene <- CnvtypeBarplot(margin_value = 1, classified = 4)

col_cnv <- structure(c("#812321","#E5B0B0","#003987","#BACEE9"), 
                     names=c('Amplification', 'Gain', 'Deletion', 'Loss'))

pdf(paste0(workdir,"revised_output_star/CNV_out/driver_gene_cnv_processed_FACET.pdf"),width = 15, height = 5)
Heatmap(mat_processed_cnv_driver, 
        column_split = laml_clin$reset_name, 
        cluster_rows = F, cluster_columns = F,
        col = col_own, 
        column_title = "Top 30 events of CNV",
        show_column_names = F,
        heatmap_legend_param = heatmap_legend_param_own,
        gap = unit(5, "mm"), rect_gp = gpar(col = "white"),
        na_col = "white",border = T,
        # add top annotation
        top_annotation = HeatmapAnnotation(barplot = anno_barplot(cnv_type_perSample,
                                                                  ylim = c(0,max(rowSums(cnv_type_perSample))),
                                                                  border = F, height = unit(4, "cm"),bar_width =0.8,
                                                                  # show barplot title
                                                                  show_annotation_name = F,
                                                                  gp = gpar(fill = col_cnv,col="white"),  axis = T)),
        # add clinical feature
        bottom_annotation = HeatmapAnnotation(df = laml_clin[,c("Location","Primary_loc")],
                                              col = annotationColor_own,
                                              show_annotation_name = T, show_legend = T),
        # add cnv variation frequency
        left_annotation = rowAnnotation(pct = anno_text(pct_cnv, just = "right", 
                                                        location = unit(1, "npc"), 
                                                        gp = gpar(fontsize = gp_fontsize), 
                                                        width = max_text_width(pct_cnv,  gp = gpar(fontsize = gp_fontsize)) 
                                                        + unit(1, "mm")), show_annotation_name = FALSE),
        # add cnv variation type bar plot
        right_annotation = rowAnnotation(barplot = anno_barplot(cnv_type_perGene, 
                                                                ylim = c(0, max(rowSums(cnv_type_perGene))), 
                                                                axis_param = list(side="top"),
                                                                border = F,width = unit(4, "cm"),bar_width =0.8,
                                                                show_annotation_name = F,
                                                                gp = gpar(fill = col_cnv,col="white"),  axis = T)))
dev.off()
