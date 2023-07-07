## Tengjia Jiang, 2023 July
## sequencing data quality

library(ggplot2)
library(magrittr)
library(ggthemr)
library(reshape2)
library(dplyr)
library(stringr)
library(tidyr)
workdir <- "/home/rstudio/storage-and-archive/multiple_primaryCRC/revised_process/multiPrimaryCRC/"

## location data-------
loc_data <- read.delim(paste0(workdir,"revised_input_star/laml_clin.txt"),header = T,sep = "\t")

## clinical data-------
normal_clin <- read.delim(paste0(workdir,"revised_input_star/CA_blood_pair_location.list"),
                          header = T, sep = "\t" , stringsAsFactors = F)
normal_clin <- merge(normal_clin,loc_data[,c("SamLocation","reset_name")],by="SamLocation",all = T)

tumor_bam <- data.frame(bam=normal_clin$Tumor_Sample_Barcode, type="tumor",reset_name=normal_clin$reset_name)
normal_bam  <- data.frame(bam=normal_clin$Normal_Sample_Barcode, type="normal",reset_name=normal_clin$reset_name)

all_bam <- rbind(tumor_bam,normal_bam) %>% unique(.)
## mean coverage (sequencing depth)-------------
qualimap_depth <- read.delim(paste0(workdir,"revised_input_star/multiqc_qualimap_bamqc_genome_results.txt"),
                           sep = "\t", header = T) 
qualimap_depth$Sample <- gsub("_sort_dedup_realigned_recal","",x = qualimap_depth$Sample)
qualimap_depth <- merge(qualimap_depth,all_bam,by.x="Sample",by.y="bam",all.y=T) 
write.table(qualimap_depth, paste0(workdir,"revised_output_star/quality_assess/multiqc_qualimap_bamqc_genome_results_filter.txt"),
            col.names = T, row.names = F, sep = "\t",quote = F)

## mean sequencing depth for tumors
tumor_qualimap_depth <- merge(qualimap_depth,loc_data[,c("reset_SamLocation","BAM_name")],by.x="Sample",by.y="BAM_name")
tumor_qualimap_depth$order <- gsub("P","",tumor_qualimap_depth$reset_name)
tumor_qualimap_depth$order <- as.numeric(tumor_qualimap_depth$order)
tumor_qualimap_depth <- tidyr::separate(tumor_qualimap_depth,col = "reset_SamLocation",into = c(NA,"order2"),remove = F)
tumor_qualimap_depth <- tumor_qualimap_depth %>% dplyr::arrange(order,order2)
tumor_qualimap_depth$reset_SamLocation <- factor(tumor_qualimap_depth$reset_SamLocation, levels = tumor_qualimap_depth$reset_SamLocation)

pdf(paste0(workdir,"revised_output_star/quality_assess/sequencing_mean_coverage.pdf"),height = 4, width = 9)
ggplot(tumor_qualimap_depth, aes(x=reset_SamLocation , y=mean_coverage)) +
  geom_bar(stat = "identity",fill="#D97677",colour="#2D3A3F" )+
  ylab("Mean depth(X)")+
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        panel.grid=element_blank(),panel.background=element_rect( fill='transparent',color='white'), #网格线设置为不显示 
        panel.border=element_blank(), #边界线设置为空 
        axis.line = element_line(colour = "black"))+
  scale_y_continuous(expand = c(0,0))
dev.off()

## histogram plot of comparison between tumors and normals
pdf(paste0(workdir,"revised_output_star/quality_assess/sequencing_coverage_distribution.pdf"),height = 3, width = 3.5)
ggplot(qualimap_depth, aes(x = mean_coverage, color = type, fill=type))+
  geom_histogram(fill="white", alpha=0.5, position="identity")+
  scale_color_manual(breaks = c("tumor","normal"),values = c("#D37838","#4EA585"))+
  xlab("Mean coverage (X)")+
  ylab("N")+
  theme_classic()
dev.off()

## genome fraction coverage----------
genome_fraction <- read.delim(paste0(workdir, "revised_input_star/genome_frction_coverage.txt"),
                              header = F,sep = " ")
colnames(genome_fraction) <- c("Sample","Genome_fraction_30x")
genome_fraction <- genome_fraction %>% 
  dplyr::filter(Sample %in% all_bam$bam) %>% 
  dplyr::arrange(desc(Genome_fraction_30x))
## bar plot
pdf(paste0(workdir,"revised_output_star/quality_assess/genome_fraction_coverage.pdf"),height = 2.6, width = 5.5)
ggplot(genome_fraction, aes(x=reorder(Sample,-Genome_fraction_30x) , y=Genome_fraction_30x))+
  geom_bar(stat = "identity")+
  xlab("")+
  ylab("Fraction of reference (%)")+
  ggtitle("Genome fraction covered by at least 30 reads")+
  geom_hline(aes(yintercept=70),color="#990000",linetype="dashed")+
  theme_classic()+
  theme(axis.text.x = element_blank(),   #移除刻度线标签
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12))+
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,100)) #设置Y轴范围
dev.off()


## sequence quality--------------
sequence_quality <- read.delim(paste0(workdir, "revised_input_star/multiqc_general_stats.txt"),
                              header = T,sep = "\t")
 
sequence_quality <- tidyr::separate(sequence_quality,col = "Sample",into = c("Sample",NA),sep = "-") %>%
  dplyr::filter(Sample %in% all_bam$bam)


sequence_data <- merge(qualimap_depth,sequence_quality,by = "Sample") %>% 
  merge(.,genome_fraction,by="Sample") %>% 
  dplyr::select(-bam_file) %>%
  dplyr::arrange(factor(reset_name,levels = paste0("P",c(1:51)))) %>%
  dplyr::select(reset_name,Sample,type,everything()) 
write.table(sequence_data,
            paste0(workdir,"revised_output_star/quality_assess/sequence_data_quality.txt"),
            col.names = T, row.names = F, sep = "\t",quote = F)
