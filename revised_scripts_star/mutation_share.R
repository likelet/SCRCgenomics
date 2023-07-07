## Tengjia Jiang, 2023 July
## mutation share between paired-samples of one patient

library(magrittr)
library(dplyr)
library(plyr)
library(tidyr)
library(VennDiagram)
library(ggsci)
library(scales)
library(Rmisc)
library(gridExtra)
library(grid)
library(gg.gap)
library(ggplot2)
library(ggthemr)
library(RColorBrewer)
library(survival)
library(survminer)
library(ggrepel)

workdir <- "/home/rstudio/storage-and-archive/multiple_primaryCRC/revised_process/multiPrimaryCRC/"
laml <- readRDS(paste0(workdir,"revised_input_star/laml.Rds"))
laml_clin <- readRDS(paste0(workdir,"revised_input_star/laml_clin.Rds"))
all_cli <- readRDS(paste0(workdir,"revised_input_star/all_cli.Rds"))
ClonalityWithdraw_res <- readRDS(paste0(workdir,"revised_input_star/ClonalityWithdraw_res.Rds"))

## mutation share details venn plots------
mutation_share <- laml@data %>% as.data.frame() %>% 
  select(Hugo_Symbol,Start_Position,Reference_Allele,Tumor_Seq_Allele2,Tumor_Sample_Barcode) %>% 
  transmute(
    mutID = paste(Hugo_Symbol, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep = ":"),
    Tumor_Sample_Barcode = Tumor_Sample_Barcode
  ) %>% 
  separate(col = "Tumor_Sample_Barcode", into = "patient", sep = "_", remove = F) %>%
  merge(., unique(laml_clin[c("Sample","reset_name")]), by.x = "patient", by.y = "Sample", all = F) 

## venn plots
show_col(pal_startrek("uniform")(2))
pal_startrek("uniform")(2)

vennPlot <- function(patient){
  samples <- mutation_share[mutation_share$patient == patient, "Tumor_Sample_Barcode"] %>% unique() %>% sort()
  reset_name <- mutation_share[mutation_share$patient == patient, "reset_name"] %>% unique()
  DT <- mutation_share[mutation_share$patient == patient & 
                         mutation_share$Tumor_Sample_Barcode == samples[grep("DT",samples)], "mutID"]
  PX <- mutation_share[mutation_share$patient == patient & 
                         mutation_share$Tumor_Sample_Barcode == samples[grep("PX",samples)], "mutID"]
  other <- mutation_share[mutation_share$patient == patient & 
                            mutation_share$Tumor_Sample_Barcode == samples[grep("third",samples)], "mutID"]
  if(length(samples) == 2){
    venn.plot <- venn.diagram(
      list(DT= DT, PX = PX),
      filename = NULL, 
      lwd = 3,
      main = reset_name, # title
      main.pos = c(0.5,0.2),
      main.cex = 2,
      col = "transparent",
      fill = c("#AB4B52","#1F78B4"), 
      alpha = 0.6,
      label.col = "black",
      cex = 1.5,
      fontfamily = "serif",
      fontface = "bold",
      cat.col = c("#AB4B52","#1F78B4"), 
      cat.cex = 2,
      cat.fontfamily = "serif",
      cat.fontface = "bold",
      margin = 0.1,
      cat.dist = c(0.03, 0.03),
      cat.pos = c(3,9)
    )
  } else{
    venn.plot <- venn.diagram(
      list(DT= DT, PX = PX, other = other),
      filename = NULL, #设为空
      lwd = 3,
      main = reset_name, # title
      main.pos = c(0.5,0.2),
      main.cex = 2,
      col = "transparent",
      fill = c("#AB4B52","#1F78B4","#9896f1"), #填充颜色类别
      alpha = 0.6,
      label.col = "black",
      cex = 1.5,
      fontfamily = "serif",
      fontface = "bold",
      cat.col = c("#AB4B52","#1F78B4","#B697C5"), # 分类类别颜色
      cat.cex = 2,
      cat.fontfamily = "serif",
      cat.fontface = "bold",
      margin = 0.1,
      cat.dist = c(0.03, 0.03,0.03),
      cat.pos = c(3,9,6)
    )
  }
  venn.plot
}

for (patient in unique(laml_clin$Sample)) {
  vennPlot(patient = patient)
}


## combine venn plots
gs <- lapply(unique(mutation_share$patient), function(patient){
  grobTree(vennPlot(patient))
  
})

pdf(file=paste0(workdir, "revised_output_star/mutation_share/combind_VennDiagram.pdf"), width = 20,height = 25)
grid.arrange(grobs = gs, ncol = 6,
             top = "mutation shared by multi-samples")
dev.off()

## mutation share bar plot----
share_gene <- function(patient){
  samples <- mutation_share[mutation_share$patient == patient, "Tumor_Sample_Barcode"] %>% unique() %>% sort()
  DT <- mutation_share[mutation_share$patient == patient & 
                         mutation_share$Tumor_Sample_Barcode == samples[grep("DT",samples)], "mutID"]
  PX <- mutation_share[mutation_share$patient == patient & 
                         mutation_share$Tumor_Sample_Barcode == samples[grep("PX",samples)], "mutID"]
  share_mut <- intersect(DT, PX) %>% as.data.frame()
  if(length(intersect(DT, PX)) > 0){
    share_da <- data.frame(patient = patient, share_mut = share_mut)
  }
}

share_GENE <-do.call(rbind, lapply(as.list(unique(laml_clin$Sample)), share_gene))
colnames(share_GENE) <- c("patient","id")
gene_share <- share_GENE %>% separate(col = "id", into = "gene", remove = F)

drivergene <- read.delim(paste0(workdir,"revised_input_star/DriverDBv3_COAD.txt"),header = T,sep = "\t",stringsAsFactors = F)
gene_share <- gene_share %>% dplyr::mutate(isDriver = ifelse(gene %in% drivergene$gene, "Driver", "no"))

mutation_share_upset <- mutation_share %>% 
  separate(col = "Tumor_Sample_Barcode", into = c(NA,"Location"), remove = T)

HeapPercentBar <- function(pat){
  intersect_type <- mutation_share_upset %>% 
    dplyr::filter(patient == pat) %>% 
    xtabs(~Location+mutID, .) %>% 
    as.data.frame() %>%
    dplyr::filter(Location %in% c("DT","PX")) %>%
    reshape2::dcast(mutID ~ Location, value.var = "Freq") %>% 
    dplyr::mutate(Patient = pat)
  intersect_type <- intersect_type %>% 
    dplyr::mutate(instersect = case_when(
      DT == PX & DT == 1 ~ "PX_DT",
      DT == 1 & PX == 0  ~ "DT only",
      DT == 0 & PX == 1  ~ "PX only"))
  intersect_type_table <- table(intersect_type$instersect) %>% 
    as.data.frame() %>%
    dplyr::mutate(Patient = pat)
  return(intersect_type_table)
}


intersect_type_table <- do.call(rbind, lapply(as.list(unique(mutation_share_upset$patient)), HeapPercentBar))
colnames(intersect_type_table) <- c("intersect_type","Freq","Patient")
intersect_type_table$intersect_type <- factor(intersect_type_table$intersect_type, ordered = T, levels = c("DT only","PX only","PX_DT"))

itt <- intersect_type_table %>% 
  ddply("Patient", transform, percent_freq = round(Freq / sum(Freq) * 100,2)) %>%
  ddply("Patient", transform, percent_freq_label = round(cumsum(percent_freq),2)) %>%
  merge(., unique(laml_clin[c("Sample","reset_name")]), by.x = "Patient", by.y = "Sample", all = T)


pat_order <- itt %>% dplyr::select(reset_name,intersect_type,percent_freq) %>%
  dplyr::filter(!is.na(percent_freq)) %>%
  dplyr::mutate(intersect_type = case_when(intersect_type == "DT only" ~  "DT_only",
                                           intersect_type == "PX only" ~  "PX_only",
                                           intersect_type == "PX_DT" ~  "PX_DT")) %>%
  reshape2::dcast(.,reset_name~intersect_type,value.var = "percent_freq") %>%
  dplyr::arrange(desc(PX_DT),desc(DT_only))


pdf(paste0(workdir, "revised_output_star/mutation_share/share_proportion_barplot.pdf"), width = 7,height = 2)
ggthemr("fresh")
ggplot(itt, aes(x=reset_name, y=percent_freq, fill=intersect_type)) +
  geom_bar(stat = "identity")+
  xlab("Patients") +
  ylab("Intersect type proportion")+
  scale_fill_manual(name = "Intersect type",
                    values = c("#65ADC2","#233B43","#E84646"),
                    breaks = c("DT only","PX only","PX_DT"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 7),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(limits = pat_order$reset_name)
ggthemr_reset()
dev.off()

## clonal driver gene locus-------
drivergene <- read.delim(paste0(workdir,"revised_input_star/DriverDBv3_COAD.txt"),header = T,sep = "\t",stringsAsFactors = F)
clonality_patient <- ClonalityWithdraw_res %>% 
  merge(., laml_clin[c("BAM_name","reset_SamLocation","Location","reset_name")], by.x = "sample", by.y = "BAM_name", all.x = T) %>% tidyr::unite(., col = patient_gene, c("reset_name","Hugo_Symbol"), sep = "_", remove = F) %>%
  tidyr::unite(., col = patient_gene_id, c("patient_gene","Start_position"), sep = "_", remove = T) 

gene_share_patient <- gene_share %>% 
  merge(., laml_clin[c("Sample","reset_SamLocation","Location","reset_name")], by.x = "patient", by.y = "Sample", all.x = T) %>%
  tidyr::unite(., col = patient_gene, c("reset_name","gene"), sep = "_", remove = F) %>%
  tidyr::separate(., col = id, into = c(NA,"Start_position", NA,NA), remove = F, sep = ":") %>%
  tidyr::unite(., col = patient_gene_id, c("patient_gene","Start_position"), sep = "_", remove = T) 

geneshare_clonality <- clonality_patient %>%
  dplyr::filter(patient_gene_id %in% unique(gene_share_patient$patient_gene_id)) %>%
  dplyr::filter(Hugo_Symbol %in% drivergene$gene) 
table(geneshare_clonality$reset_name)

plotList <- list()
for (samp in c("P48","P11","P37","P2")) {
  geneshareColTb <- geneshare_clonality %>% 
    dplyr::filter(reset_name == samp)%>%
    dplyr::select(Location,patient_gene_id,cancer_cell_frac) %>%
    reshape2::dcast(patient_gene_id~Location, value.var = "cancer_cell_frac") %>%
    tidyr::separate(., col = patient_gene_id, into = c(NA,"Gene",NA), remove = T, sep = "_")
  
  rowshare <- nrow(geneshareColTb)
  print(samp)
  print(rowshare)
  p <- ggplot(geneshareColTb, aes(x=DT, y=PX))+
    geom_point(size=8, shape=21, colour="black", fill="#9D1A20")+
    geom_abline(intercept=0, slope=1, colour="black",linetype="dashed")+
    theme_classic()+
    xlab(samp)+
    geom_text_repel(aes(DT, PX, label = Gene))+
    labs(y="")+
    scale_y_continuous(limits = c(0,1), breaks = c(0,0.5,1))+
    scale_x_continuous(limits = c(0,1), breaks = c(0,0.5,1))+
    theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_blank(),
          plot.title = element_text(size=15,hjust = 0.5),
          axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15))+
    guides(color=F)
  plotList[[samp]] <- p
}

pdf(paste0(workdir,"revised_output_star/clonality/Clonality_compare_of_sharegene.pdf"),width = 5, height = 5)
multiplot(plotlist = plotList, cols = 2)
dev.off()

