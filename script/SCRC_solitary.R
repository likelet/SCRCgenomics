## Tengjia Jiang, 11.07.2023
## SCRC & solitary comparison

library(magrittr)
library(data.table)
library(dplyr)
library(table1)
library(htmltools)
library(survival)
library(survminer)
library(ggsci)
library(Rmisc)
library(dplyr)
library(scales)
library(sampling)
library(ggplot2)
library(ggthemr)
library(gg.gap)
library(cowplot)

workdir <- "./"
laml <- readRDS(paste0(workdir,"input/laml.Rds"))
laml_clin <- readRDS(paste0(workdir,"input/laml_clin.Rds"))
all_cli <- readRDS(paste0(workdir,"input/all_cli.Rds"))
tcga_laml <- readRDS(paste0(workdir,"input/tcga_laml.Rds"))
followup_data <- readRDS(paste0(workdir,"input/followup_data.Rds"))
path_sum <- readRDS(paste0(workdir,"input/path_sum.Rds"))
oncokb_maf <- readRDS(paste0(workdir,"input/oncokb_maf.Rds"))
oncokb_maf_tcga <- readRDS(paste0(workdir,"input/oncokb_maf_tcga.Rds"))


## SCRC data
scrc_surv <- followup_data %>%
  dplyr::filter(reset_name %in% unique(laml_clin$reset_name)) %>%
  dplyr::select(reset_name,OS_status,OS) %>%
  dplyr::mutate(group="SCRC")

## TCGA data
Covariate_dcast_TCGA <- tcga_laml@data %>%
  dplyr::filter(Hugo_Symbol %in% drivergene$gene) %>%
  dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode) %>%
  reshape2::dcast(., Tumor_Sample_Barcode ~ Hugo_Symbol) %>%
  as.data.frame()

## survival analysis--------
## TMB
tcga_TMB <- tcga_laml@variants.per.sample %>% as.data.frame()

OS_data_TCGA <- tcga_laml@clinical.data %>% 
  dplyr::select(Tumor_Sample_Barcode,days_to_last_followup,vital_status) %>%
  dplyr::mutate(Osstatus = as.numeric(vital_status)) %>%
  dplyr::mutate(OS = as.numeric(days_to_last_followup)/30) 
tcga_os <- OS_data_TCGA %>%
  dplyr::select(Tumor_Sample_Barcode,Osstatus,OS) %>%
  dplyr::mutate(group = "TCGA") 
colnames(tcga_os) <- c("reset_name", "OS_status", "OS", "group")

SCRC_solitary_tcga <- rbind(scrc_surv, tcga_os) %>%
  dplyr::mutate(group = factor(group, levels = c("SCRC","TCGA")))

OS_diff_SCRC_soli_tcga <- survfit(Surv(OS, OS_status)~ group, data = SCRC_solitary_tcga)
survdiff(Surv(OS, OS_status)~ group, data = SCRC_solitary_tcga)


pdf(paste0(workdir,"output/SCRC_solitary/SCRC_solitary_survival_tcga.pdf"), width = 3, height = 4)
p3 <- ggsurvplot(OS_diff_SCRC_soli_tcga, data = SCRC_solitary_tcga,
                 fun = "pct",conf.int = F, pval = T,
                 tables.height = 0.3,
                 risk.table = T,
                 risk.table.y.text.col = T, 
                 risk.table.y.text = F, 
                 palette = c("#BC3C29FF","#20854EFF"),
                 legend.title = "",
                 legend.labs = c("SCRC", "TCGA"),
                 xlab = "Overall survival (months)")
print(p3,newpage = F)
dev.off()

## SCRC & solitary CRC mutation rate ------
pt.vs.rt_ss_seprate <- mafCompare(m1 = laml, m2 = tcga_laml, m1Name = 'SCRC', m2Name = 'TCGA', minMut = 5)
mut_freq <- pt.vs.rt_ss_seprate$results %>% 
  dplyr::filter(Hugo_Symbol %in% drivergene$gene) %>%
  dplyr::filter(adjPval < 0.25) %>%
  dplyr::filter(pval < 0.05) %>%
  dplyr::mutate(SCRC_mut = round(SCRC/103*100,digits = 0)) %>%
  dplyr::mutate(TCGA_mut = round(TCGA/489*100,digits = 0))
diff_mut_gene_seprate <- mut_freq %>%
  dplyr::select(Hugo_Symbol) 
path_order_seprate <- merge(path_sum, diff_mut_gene_seprate, by.x = "Gene", by.y = "Hugo_Symbol") %>% dplyr::arrange(Pathway)
mut_freq <- mut_freq %>% dplyr::filter(Hugo_Symbol %in% path_order_seprate$Gene)

pdf(paste0(workdir, "output/SCRC_solitary/SCRC_TCGA_mut_persample.pdf"), width = 4,height = 5)
coBarplot(m1 = laml, m2 = tcga_laml, 
          genes = path_order_seprate$Gene, 
          m1Name = "SCRC", m2Name = "TCGA",
          colors = vc_cols)
dev.off()

## TMB rainplot---------------
tcga_tmb <- tcga_laml@variants.per.sample %>% 
  dplyr::mutate(logTMB = log10(Variants/37+1)) %>%
  dplyr::mutate(TMB = Variants/37) %>%
  dplyr::mutate(group = "TCGA")
scrc_tmb <- laml@variants.per.sample %>% 
  dplyr::mutate(logTMB = log10(Variants/37+1)) %>%
  dplyr::mutate(TMB = Variants/37) %>%
  dplyr::mutate(group = "SCRC")
tcga_scrc_tmb <- rbind(tcga_tmb, scrc_tmb)


tcga_scrc_tmb_p <- ggplot(tcga_scrc_tmb,aes(x=group,y=logTMB, fill = group))+
  geom_flat_violin(aes(fill = group),
                   position = position_nudge(x = .05, y = 0),adjust =2, alpha = 1,trim = T,colour = NA)+
  geom_point(aes(x = group, y = logTMB,colour = group), 
             position = position_jitter(width = .05),size = 1, shape = 20)+
  geom_boxplot(aes(x = group, y = logTMB),
               position = position_nudge(x = -.1, y = 0), outlier.shape = NA, alpha = 1, width = .1, colour = "BLACK") +
  ylab('logTMB')+xlab('')+theme_cowplot()+guides(fill = "none", colour = "none") +
  scale_color_manual(breaks = c("SCRC","TCGA"), values = c("#CE364F","#86A875"))+
  scale_fill_manual(breaks = c("SCRC","TCGA"), values = c("#CE364F","#86A875"))+
  scale_y_continuous(breaks = c(seq(0,3,0.5)), expand = c(0,0))+
  scale_x_discrete(breaks = c("SCRC","TCGA"), labels = c("SCRC\n(N=103)","TCGA\n(N=489)"))
ggsave(paste0(workdir,"output/SCRC_solitary/TCGA_SCRC_tmb_raincloud.pdf"), width = 4, height = 4.5)

t.test(tcga_tmb$TMB, scrc_tmb$TMB, paired = F, alternative = "two.sided")

## pathway mutant rate-------
## SCRC data
path_oncokb <- merge(oncokb_maf, path_sum, by.x = "Hugo_Symbol", by.y = "Gene")

## tcga data
path_oncokb_tcga <- merge(oncokb_maf_tcga, path_sum, by.x = "Hugo_Symbol", by.y = "Gene")

## SCRC pathway mutation frequency
PathMut <- function(path){
  mut_count <- path_oncokb %>% dplyr::filter(Pathway == path) %>% dplyr::select(Tumor_Sample_Barcode) %>% unique() %>% nrow()
  prop <- round(mut_count/103*100, 1)
  pathmut <- data.frame(Path = path, Prop = prop, group = "SCRC", label_prop = paste0(prop,"%"), count = mut_count)
  return(pathmut)
}
pathmut_scrc <- do.call(rbind, lapply(as.list(unique(path_oncokb$Pathway)), PathMut)) %>% dplyr::mutate(Prop = Prop*(-1))

## TCGA pathway mutation frequency
PathMutCK <- function(path){
  mut_count <- path_oncokb_tcga %>% dplyr::filter(Pathway == path) %>% dplyr::select(Tumor_Sample_Barcode) %>% unique() %>% nrow()
  prop <- round(mut_count/489*100, 1)
  pathmut <- data.frame(Path = path, Prop = prop, group = "TCGA", label_prop = paste0(prop,"%"), count = mut_count)
  return(pathmut)
}
pathmut_tcga <- do.call(rbind, lapply(as.list(unique(path_oncokb_tcga$Pathway)), PathMutCK))

## barplot
pathmut_compare_tcga <- rbind(pathmut_scrc, pathmut_tcga) 
pathmut_order_tcga <- pathmut_compare_tcga %>% dplyr::filter(group == "SCRC") %>% dplyr::arrange(Prop)

pdf(paste0(workdir, "output/SCRC_solitary/SCRC_TCGA_pathway.pdf"), width = 4,height = 4)
ggplot(pathmut_compare_tcga, aes(Path, Prop, fill = group)) +
  geom_col(position = position_dodge(width = 0), 
           width = 1.2, size = 0.3, colour = 'black') +
  geom_text(aes(y=Prop, label=label_prop),
            vjust=1.5, colour="black")+
  coord_flip()+
  scale_fill_manual(breaks = c("SCRC","TCGA"), values = c("#CE364F","#86A875"))+
  scale_x_discrete(limits = rev(pathmut_order_tcga$Path))+
  labs(y = 'Percent of case',x = "") +
  geom_hline(yintercept = 0, size = 0.3) +
  guides(fill=F)+  
  scale_y_continuous(breaks = seq(-90, 90, 20), labels = as.character(abs(seq(-90, 90, 20))), limits = c(-90, 90)) +
  annotate('text', label = 'SCRC (N=103)', 1, -60) +
  annotate('text', label = 'TCGA (N=489)', 1, 60) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.title = element_blank())
dev.off()
tcga_laml@clinical.data %>% head()

## chi-square test
PathMutDiff <- function(path){
  print(path)
  scrc_path_mut <- pathmut_compare_tcga %>% dplyr::filter(Path == path & group == "SCRC") 
  scrc_path_mut <- scrc_path_mut$count
  tcga_path_mut <- pathmut_compare_tcga %>% dplyr::filter(Path == path & group == "TCGA") 
  tcga_path_mut <- tcga_path_mut$count
  
  P <- fisher.test(rbind(c(scrc_path_mut, 103-scrc_path_mut), c(tcga_path_mut, 489-tcga_path_mut))) %>% .$p.value
  table_bind <- data.frame(Path = path, P.value = P)
  return(table_bind)
}

do.call(rbind, lapply(as.list(unique(pathmut_compare_tcga$Path)), PathMutDiff)) %>% dplyr::filter(P.value < 0.05)


