## Tengjia Jiang, 2023 July
## mutation landscape of top 30 mutant genes

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

workdir <- "/home/rstudio/storage-and-archive/multiple_primaryCRC/revised_process/multiPrimaryCRC/"
mat_processed_mut <- readRDS(paste0(workdir,"revised_input_star/mat_processed_mut.Rds"))
laml <- readRDS(paste0(workdir,"revised_input_star/laml.Rds"))
laml_clin <- readRDS(paste0(workdir,"revised_input_star/laml_clin.Rds"))

top_mutated_gene <- laml@gene.summary$Hugo_Symbol %>% head(30)
select_mutgene <- intersect(rownames(mat_processed_mut), top_mutated_gene)
mat_processed_mut_driver <- mat_processed_mut[select_mutgene,]

altered_nums_mut <- apply(mat_processed_mut_driver, 1, function(x){length(x)-sum(x=="")})
NUMBER_GENES<- 30
slice_indx_mut <- sort(altered_nums_mut, decreasing = TRUE)[1:NUMBER_GENES]
mat_processed_mut_driver <-  mat_processed_mut_driver[names(slice_indx_mut),]

## oncoplot-----------
## cauculate mut frequency
pct_num_mut = rowSums(apply(mat_processed_mut_driver, 1:2, function(x)any(x!="")))/ncol(mat_processed_mut_driver)
pct_mut = paste0(round(pct_num_mut * 100, digits = 0), "%")
gp_fontsize <- 15

## mutation type fraction
MuttypeBarplot <- function(margin_value){
  mut_type <- apply(mat_processed_mut_driver, margin_value, function(x){as.data.frame(table(x))})
  mut_type <- do.call(rbind,mut_type)
  mut_type$Gene <- rownames(mut_type)
  colnames(mut_type)[1] <- "Mut_type"
  mut_type$Mut_type <- as.character(mut_type$Mut_type)
  mut_type$Freq <- as.numeric(mut_type$Freq)
  mut_type <- mut_type[mut_type$Mut_type != "",]
  mut_type <- separate(mut_type,col = "Gene",into = "Hugo_Sample",sep = "[.]",remove = T)
  mut_type <- reshape2::dcast(Hugo_Sample ~ Mut_type,data = mut_type,value.var = "Freq")
  rownames(mut_type) <- mut_type$Hugo_Sample
  mut_type[is.na(mut_type)] <- 0
  barplot_order <- c('Frame_Shift_Del','Missense_Mutation','Nonsense_Mutation','Multi_Hit',
                     'Frame_Shift_Ins','In_Frame_Ins','Splice_Site','In_Frame_Del')
  if(margin_value == 1){
    mut_type <- mut_type[rownames(mat_processed_mut_driver),-1]
  }else{
    mut_type <- mut_type[,-1]
    numut_sample <- setdiff(colnames(mat_processed_mut_driver),rownames(mut_type))
    nonmut_mat <- matrix(0,nrow=length(numut_sample), ncol=8)
    colnames(nonmut_mat) <- barplot_order
    rownames(nonmut_mat) <- numut_sample
    mut_type <- rbind(mut_type, nonmut_mat) %>% as.matrix()
    mut_type <- mut_type[colnames(mat_processed_mut_driver),]
  }
  mut_type <- mut_type[,barplot_order]
  return(mut_type)
}
mut_type_perSample <- MuttypeBarplot(margin_value = 2)
mut_type_perGene <- MuttypeBarplot(margin_value = 1)


pdf(paste0(workdir,"revised_output_star/somatic_mut/top_gene_mutation_processed.pdf"),width = 20, height = 10)
Heatmap(mat_processed_mut_driver, 
        column_split = laml_clin$Sample, 
        cluster_rows = F, cluster_columns = F,
        col = col_own, 
        column_title = "Top 30 events of mutation",
        show_column_names = F,
        heatmap_legend_param = heatmap_legend_param_own,
        gap = unit(5, "mm"), rect_gp = gpar(col = "white"),
        na_col = "white",border = T,
        # add top annotation
        top_annotation = HeatmapAnnotation(barplot = anno_barplot(mut_type_perSample,
                                                                  ylim = c(0,max(rowSums(mut_type_perSample))),
                                                                  border = F, height = unit(4, "cm"),bar_width =0.8,
                                                                  # show barplot title
                                                                  show_annotation_name = F,
                                                                  gp = gpar(fill = col_own,col="white"),  axis = T)),
        # add clinical feature
        bottom_annotation = HeatmapAnnotation(df = laml_clin[,c("Location","Primary_loc")],
                                              col = annotationColor_own,
                                              show_annotation_name = T, show_legend = T),
        # add mutation frequency
        left_annotation = rowAnnotation(pct = anno_text(pct_mut, just = "right", 
                                                        location = unit(1, "npc"), 
                                                        gp = gpar(fontsize = gp_fontsize), 
                                                        width = max_text_width(pct_mut,  gp = gpar(fontsize = gp_fontsize)) 
                                                        + unit(1, "mm")), show_annotation_name = FALSE),
        # add mutation type bar plot
        right_annotation = rowAnnotation(barplot = anno_barplot(mut_type_perGene, 
                                                                ylim = c(0, max(rowSums(mut_type_perGene))), 
                                                                axis_param = list(side="top"),
                                                                border = F,width = unit(4, "cm"),bar_width =0.8,
                                                                show_annotation_name = F,
                                                                gp = gpar(fill = col_own,col="white"),  axis = T)))
dev.off()

## gene mutation share ----------
mut_location_prop <- laml@data %>% 
  dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  dplyr::filter(Hugo_Symbol %in% rownames(mat_processed_mut_driver)) %>%
  merge(., laml_clin[,c("Tumor_Sample_Barcode","Location","Sample")], by = "Tumor_Sample_Barcode", all.y = T) %>% 
  unique() %>%
  dplyr::filter(Location %in% c("DT","PX")) %>%
  dcast(Hugo_Symbol ~ Sample, value.var = "Location") %>%
  reshape2::melt(id.vars = "Hugo_Symbol") %>%
  dplyr::mutate(value = case_when(
    value == 0 ~ "Neither",
    value == 1 ~ "Either",
    value == 2 ~ "Both")) %>%
  dplyr::filter(!is.na(Hugo_Symbol))

LocShareMuttype <- function(gene){
  loc_type <- mut_location_prop %>% 
    dplyr::filter(Hugo_Symbol == gene) 
  loc_type <- table(loc_type$value) %>% as.data.frame() %>%
    dplyr::mutate(Gene = gene) %>%
    dplyr::rename(Mutant_type = Var1) %>%
    dplyr::select(Gene, everything())
  return(loc_type)
}
loc_muttype_res <- do.call(rbind, lapply(as.list(rownames(mat_processed_mut_driver)), LocShareMuttype))

pdf(paste0(workdir, "revised_output_star/somatic_mut/top_gene_mutation_share.pdf"), width = 3,height = 6)
ggthemr("fresh")
ggplot(loc_muttype_res, aes(x=Gene, y=Freq,  fill = Mutant_type)) + 
  geom_bar(stat = "identity",width = 0.6)+
  scale_x_discrete(limits = rev(rownames(mat_processed_mut_driver)))+
  scale_y_continuous(expand = c(0, 0),
                     position = "left",
                     breaks = seq(0, 55, 10)) +
  scale_fill_manual(
    name = NULL,
    values = c("#db6400","#0e918c","#595b83"),
    breaks = c("Both","Either","Neither")
  ) +
  coord_flip() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15))+
  labs(x = NULL, y = "Patients distribution", fill = NULL,
       title = "Multi-samples mutation distribution") 
ggthemr_reset()
dev.off()

## driver gene difference between left and right sides of colon-------
drivergene <- read.delim(paste0(workdir,"revised_input_star/DriverDBv3_COAD.txt"),header = T,sep = "\t",stringsAsFactors = F)
select_mutgene <- intersect(rownames(mat_processed_mut), drivergene$gene)
mat_processed_mut_driver <- mat_processed_mut[select_mutgene,]

altered_nums_mut <- apply(mat_processed_mut_driver, 1, function(x){length(x)-sum(x=="")})
NUMBER_GENES<- 30
slice_indx_mut <- sort(altered_nums_mut, decreasing = TRUE)[1:NUMBER_GENES]
mat_processed_mut_driver <-  mat_processed_mut_driver[names(slice_indx_mut),]


select_driver <- rownames(mat_processed_mut_driver)
left_maf <- subsetMaf(maf = laml,
                      clinQuery = "Primary_loc == 'Left' ")
right_maf <- subsetMaf(maf = laml,
                       clinQuery = "Primary_loc == 'Right' ")
left_mut <- left_maf@gene.summary %>%
  dplyr::filter(Hugo_Symbol %in% select_driver) %>%
  dplyr::mutate(Mut_rate = round(AlteredSamples/nrow(left_maf@variants.per.sample)*100, 1)) %>%
  dplyr::select(Hugo_Symbol,Mut_rate) %>%
  dplyr::mutate(group = "Left")
right_mut <- right_maf@gene.summary %>%
  dplyr::filter(Hugo_Symbol %in% select_driver) %>%
  dplyr::mutate(Mut_rate = round(AlteredSamples/nrow(right_maf@variants.per.sample)*100, 1)) %>%
  dplyr::select(Hugo_Symbol,Mut_rate) %>%
  dplyr::mutate(group = "Right")
left_right_mut <- rbind(left_mut, right_mut) %>%
  dplyr::mutate(label_prop = paste0(Mut_rate,"%")) %>% dplyr::arrange(desc(group),desc(Mut_rate))
mut_order <- left_right_mut %>% dplyr::filter(group == "Left") %>% dplyr::select(Hugo_Symbol)

pdf(paste0(workdir, "revised_output_star/somatic_mut/driver_left_right_compare.pdf"), width = 7,height = 3)
ggthemr("fresh")
ggplot(left_right_mut, aes(x=Hugo_Symbol, y=Mut_rate, fill=group)) +
  geom_bar(stat = "identity",colour="#393e46", width = 0.7, position = "dodge")+
  xlab("") +
  ylab("Mutation frequency")+
  scale_fill_manual(name = "",
                    values = c("#899DA4","#FAD510"),
                    breaks = c("Left","Right"))+
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 1, size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,80,10))+
  scale_x_discrete(limits = mut_order$Hugo_Symbol)
ggthemr_reset()
dev.off()

mafCompare(m1 = left_maf, m2 = right_maf, m1Name = 'Left', m2Name = 'Right', minMut = 5) %>%
  .$results %>%
  dplyr::filter(pval < 0.05) %>%
  dplyr::filter(Hugo_Symbol %in% select_driver)


## TMB distribution----------
all_mut <- laml@variant.classification.summary[,1:10] %>% 
  reshape2::melt(., id.vars = "Tumor_Sample_Barcode") %>%
  merge(., laml_clin, by = "Tumor_Sample_Barcode", all = T) %>%
  dplyr::arrange(Tumor_Sample_Barcode)

all_TMB <- aggregate(value~Tumor_Sample_Barcode, data = all_mut, FUN = sum) %>%
  dplyr::mutate(TMB = value/37) %>%
  merge(., laml_clin, by = "Tumor_Sample_Barcode", all = T) %>%
  dplyr::arrange(Tumor_Sample_Barcode) %>%
  dplyr::filter(Location != "third")

pdf(paste0(workdir, "revised_output_star/somatic_mut/all_TMB_barplot.pdf"), width = 9,height = 4)
ggthemr("fresh")
all_tmb_plot <- ggplot(all_TMB, aes(x=reset_name, y=TMB, fill=Location)) +
  geom_bar(stat = "identity", color = "#393e46")+
  scale_fill_manual(breaks = c("DT","PX"), values = c("#AB4B52","#1F78B4"))+
  xlab("Patients") +
  ylab("non-synonymous (mut/Mbp)")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))+
  scale_y_continuous(expand = c(0,0), breaks = seq(0,200,20))+
  scale_x_discrete(limits = levels(all_TMB$Sample), labels = unique(all_TMB$reset_name))
gg.gap(plot = all_tmb_plot, 
       segments = c(80,110), #截断区域的纵坐标
       ylim = c(0,140), #纵坐标的总长度
       tick_width = c(20,30), #上层和下层分别的纵坐标间隔值
       rel_heights = c(0.4,0.00001,0.1)) # 底层、截断层、上层在纵轴上所占的比例
ggthemr_reset()
dev.off()
