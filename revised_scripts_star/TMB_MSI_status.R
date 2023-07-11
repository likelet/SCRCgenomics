## Tengjia Jiang, 11.07.2023
## TMB & MSI status

library(dplyr)
library(ggplot2)
library(magrittr)
library(RColorBrewer)
library(ggthemr)
library(VennDiagram)
library(ggsci)
library(survminer)
library(survival)
library(plyr)
library(tidyr)
library(survival)
library(survminer)
library(gtsummary)
library(table1)
library(htmltools)
library(finalfit)

workdir <- "./"
laml <- readRDS(paste0(workdir,"revised_input_star/laml.Rds"))
laml_clin <- readRDS(paste0(workdir,"revised_input_star/laml_clin.Rds"))
all_cli <- readRDS(paste0(workdir,"revised_input_star/all_cli.Rds"))
followup_data <- readRDS(paste0(workdir,"revised_input_star/followup_data.Rds"))
msi_score_res <- readRDS(paste0(workdir, "revised_input_star/msi_score_res.Rds"))

## TMB------------
tmb_compare <- laml_clin %>%
  dplyr::select(Tumor_Sample_Barcode,Location,reset_name) %>%
  merge(., laml@variants.per.sample, all.x = T) %>%
  dplyr::filter(Location %in% c("PX","DT")) %>%
  dplyr::mutate(Variants = ifelse(is.na(Variants),0,Variants)) %>%
  dplyr::mutate(logTMB = log10(Variants/37+1)) %>%
  dplyr::mutate(TMB = Variants/37)
P_test_tmb <- reshape2::dcast(tmb_compare, reset_name ~Location, value.var = "TMB")

tmb_reclass <- tmb_compare %>% 
  dplyr::mutate(TMB = ifelse(TMB > 10, "TMB_H", "TMB_L")) %>%
  reshape2::dcast(reset_name~Location, value.var = "TMB") %>%
  dplyr::mutate(CoTMB = case_when(DT == PX & DT == "TMB_H" ~ "TMB_H only",
                                  DT == PX & DT == "TMB_L" ~ "TMB_L only",
                                  DT != PX ~ "TMB_L/TMB_H")) 

## multi-samples TMB Status 
tmb_score_table <- table(tmb_reclass$CoTMB) %>% as.data.frame()
colnames(tmb_score_table) <- c("tmb_type","Num")

tmbLabel <- tmb_score_table$Num
tmbLabel <- paste0(tmbLabel, "(", round(tmb_score_table$Num / sum(tmb_score_table$Num ) * 100, 2), "%)") 

pdf(paste0(workdir, "revised_output_star/somatic_mut/tmb_score_summary_pie.pdf"), width = 5,height = 5)
ggthemr("fresh")
ggplot(tmb_score_table, aes(x = "",y = Num, fill = tmb_type)) +
  geom_bar(stat="identity",width=1) +
  coord_polar(theta="y",direction=1)+
  scale_fill_manual(values = c("#144d53","#307672","#e4eddb"),
                    breaks = c("TMB_H only","TMB_L/TMB_H","TMB_L only"))+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  labs(x="",y="",fill="TMB type")+
  ggtitle(label ="Multi-samples TMB Status",subtitle=NULL)+
  guides(fill = guide_legend(reverse = F)) +
  geom_text(aes(x=1.2,y=(sum(Num)-cumsum(Num)+Num/2),label=tmbLabel),
            size=3,color = "white",size = 16, face = "bold")
ggthemr_reset()
dev.off()

## MSS status-------
msi_score <- msi_score_res %>% 
  dplyr::mutate(MSIstatus = ifelse(MSIscore >= 20, "MSI_H", "MSS")) %>% 
  separate(.,col = "reset_SamLocation", into = c("Sample","Location"), sep = "_", remove = T) %>%
  dplyr::select(Sample, Location, MSIstatus) %>% 
  dcast(., Sample~Location, value.var = "MSIstatus") %>% 
  dplyr::select(Sample, DT, PX)

msi_score <- msi_score %>% 
  mutate(CoMSI = case_when(DT == PX & DT == "MSI_H" ~ "MSI_H only",
                           DT == PX & DT == "MSS" ~ "MSS only",
                           DT != PX & !is.na(PX) ~ "MSS/MSI_H",
                           is.na(PX) ~ "MSS only")) %>%
  mutate(MSI_2group = case_when(DT == "MSI_H" | PX == "MSI_H" ~ "MSI_H",
                                DT == "MSS" & PX == "MSS" ~ "MSS",
                                DT == "MSS" & is.na(PX) ~ "MSS"))

#### pie plot
msi_score_table <- table(msi_score$CoMSI) %>% as.data.frame()
colnames(msi_score_table) <- c("MSI_type","Num")

msiLabel <- msi_score_table$Num
msiLabel <- paste0(msiLabel, "(", round(msi_score_table$Num / sum(msi_score_table$Num ) * 100, 2), "%)") 
pdf(paste0(workdir, "revised_output_star/MSI_status/MSI_score_summary_pie.pdf"), width = 5,height = 5)
ggthemr("fresh")
ggplot(msi_score_table, aes(x = "",y = Num, fill = MSI_type)) +
  geom_bar(stat="identity",width=1) +
  coord_polar(theta="y",direction=1)+
  scale_fill_manual(values = c("#de95ba","#ffd9e8","#7f4a88"),
                    breaks = c("MSS/MSI_H","MSS only","MSI_H only"))+
  labs(x="",y="",fill="MSI type")+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  ggtitle(label ="Multi-samples MSI Status",subtitle=NULL)+
  guides(fill = guide_legend(reverse = T)) +
  geom_text(aes(x=1.2,y=(sum(Num)-cumsum(Num)+Num/2),label=msiLabel),
            size=3,color = "white",size = 16, face = "bold")
ggthemr_reset()
dev.off()

## MSS status & tumor samples location
msi_score_loc <- all_cli %>% select(reset_name,Tumor_sites) %>% merge(msi_score, .,by.x="Sample",by.y="reset_name")
msi_score_loc <- xtabs(~CoMSI+Tumor_sites, msi_score_loc) %>% as.data.frame()
msi_score_loc <- ddply(msi_score_loc, "Tumor_sites", transform,
                       percent_freq = round(Freq/sum(Freq)*100,2)) %>%
  ddply(., "Tumor_sites", transform, label_y_percent = round(cumsum(percent_freq),2)) %>%
  ddply(., "Tumor_sites", transform, label_y = cumsum(Freq)) %>%
  dplyr::mutate(myLabel = paste0(Freq," (",percent_freq,"%)"))

msi_score_loc$CoMSI <- factor(msi_score_loc$CoMSI,
                              levels = rev(c("MSI_H only","MSS only","MSS/MSI_H")))

#### stacked bar plot
pdf(paste0(workdir, "revised_output_star/MSI_status/MSI_score_location_bar_percentage.pdf"), width = 7,height = 4)
ggthemr("fresh")
ggplot(msi_score_loc, aes(x=Tumor_sites, y=percent_freq,  fill = CoMSI)) + 
  geom_bar(stat = "identity", width = 0.6)+
  geom_text(aes(y=label_y_percent, label=myLabel), vjust=1.5, colour="white")+
  scale_fill_manual(
    name = NULL,
    values = c("#de95ba","#ffd9e8","#7f4a88"),
    breaks = c("MSS/MSI_H","MSS only","MSI_H only")
  ) +
  scale_x_discrete(breaks = c("Right only", "Left only", "Right and left"), 
                   labels = c("Right only\n(N=9)", "Left only\n(N=14)", "Right and left\n(N=28)"))+
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15))+
  labs(x = NULL, y = "Patients distribution percent (%)", fill = NULL,
       title = "Multi-samples location distribution") +
  scale_y_continuous(expand = c(0,0))
ggthemr_reset()
dev.off()

msi_score_loc_p <- all_cli %>% select(reset_name,Tumor_sites) %>% merge(msi_score, .,by.x="Sample",by.y="reset_name")
xtabs(~CoMSI+Tumor_sites, msi_score_loc_p) %>% chisq.test()

## survival analysis
msi_cli <- merge(followup_data,msi_score,by.x="reset_name",by.y="Sample")
msi_cli$CoMSI <- factor(msi_cli$CoMSI, ordered = T, levels = c("MSS only","MSS/MSI_H","MSI_H only"))

OS_diff_msi_reclass <- survfit(Surv(OS, OS_status)~ CoMSI, data = msi_cli)
PFS_diff_msi_reclass <- survfit(Surv(PFS, PFS_status)~ CoMSI, data = msi_cli)
survdiff(Surv(OS, OS_status)~ CoMSI, data = msi_cli)
survdiff(Surv(PFS, PFS_status)~ CoMSI, data = msi_cli)

msi_splots_reclass <- list()
msi_splots_reclass[[1]] <- ggsurvplot(OS_diff_msi_reclass, data = msi_cli,
                                      fun = "pct",conf.int = F, pval = T,
                                      risk.table = T,
                                      palette = c("#ffd9e8","#de95ba","#7f4a88"),
                                      legend.labs = c("MSS only","MSS/MSI_H","MSI_H only"),
                                      legend.title = "",
                                      risk.table.y.text.col = T, 
                                      risk.table.y.text = F, 
                                      xlab = "Overall survival (months)")

msi_splots_reclass[[2]] <- ggsurvplot(PFS_diff_msi_reclass, data = msi_cli,
                                      fun = "pct",conf.int = F, pval = T,
                                      risk.table = T,
                                      palette = c("#ffd9e8","#de95ba","#7f4a88"),
                                      legend.labs = c("MSS only","MSS/MSI_H","MSI_H only"),
                                      legend.title = "",
                                      risk.table.y.text.col = T, 
                                      risk.table.y.text = F, 
                                      xlab = "Progression-free survival (months)")
msi_surres_reclass <- arrange_ggsurvplots(msi_splots_reclass, print = T, ncol = 2, nrow = 1, risk.table.height = 0.3)
ggsave(paste0(workdir,"revised_output_star/MSI_status/msi_reclass_survival_plot.pdf"),
       plot = msi_surres_reclass, width = 7, height = 5)

## MSS status & left/right location
msi_pri_loc <- msi_score_res %>% 
  merge(., laml_clin) %>% 
  dplyr::mutate(MSIgroup = ifelse(MSIscore >= 20, "MSI-H", "MSS")) %>% 
  dplyr::mutate(Primary_loc = as.character(Primary_loc)) %>%
  dplyr::mutate(Primary_loc = ifelse(is.na(Primary_loc), "Left", Primary_loc)) 

msi_pri_loc <- xtabs(~MSIgroup+Primary_loc, data = msi_pri_loc) %>% 
  as.data.frame() %>%
  ddply(., "Primary_loc", transform,
        percent_freq = round(Freq/sum(Freq)*100,2)) %>%
  ddply(., "Primary_loc", transform, label_y_percent = round(cumsum(percent_freq),2)) %>%
  dplyr::mutate(myLabel = paste0(Freq," (",percent_freq,"%)"))

msi_pri_loc$MSIgroup <- factor(msi_pri_loc$MSIgroup, ordered = T,levels = rev(c("MSI-H","MSS")))

pdf(paste0(workdir, "revised_output_star/MSI_status/MSI_location_bar_percentage.pdf"), width = 3.5,height = 3.5)
ggthemr("fresh")
ggplot(msi_pri_loc, aes(x=Primary_loc, y=percent_freq,  fill = MSIgroup)) + 
  geom_bar(stat = "identity", width = 0.6)+
  geom_text(aes(y=label_y_percent, label=myLabel), 
            vjust=1.5, colour="white")+
  scale_fill_manual(
    name = NULL,
    values = c("#7FCEB1","#4B3B56"),
    breaks = c("MSI-H","MSS")
  ) +
  theme(axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15))
ggthemr_reset()
dev.off()

xtabs(Freq~MSIgroup+Primary_loc, data = msi_pri_loc) %>% chisq.test()

