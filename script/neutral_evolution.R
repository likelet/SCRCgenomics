## Tengjia Jiang, 11.07.2023
## neutral evolution analysis

library(neutralitytestr)
library(ggthemr)
library(survival)
library(survminer)
library(ggsci)
library(Rmisc)
library(dplyr)
library(scales)
library(GenomicRanges)
library(magrittr)
library(tidyr)

workdir <- "./"
NeutralityR2_res <- readRDS(paste0(workdir,"input/NeutralityR2_res.Rds"))
laml_clin <- readRDS(paste0(workdir,"input/laml_clin.Rds"))
math_res_rmCNV <- readRDS(paste0(workdir,"input/math_res_rmCNV.Rds"))
followup_data <- readRDS(paste0(workdir,"input/followup_data.Rds"))

## pie plot-------
NeutralIdentical <- function(pat){
  Neutral_status <- NeutralityR2_res[NeutralityR2_res$reset_name == pat,"Neutral_status"] %>% unique()
  
  if(all(Neutral_status == "Yes") == TRUE){
    NeutralStatus <- data.frame(Patient = pat, Neutral_status = "Neutral only")
  } else if(all(Neutral_status == "No") == TRUE){
    NeutralStatus <- data.frame(Patient = pat, Neutral_status = "non_Neutral only")
  } else{
    NeutralStatus <- data.frame(Patient = pat, Neutral_status = "Different")
  }
}

NeutralStatus <- do.call(rbind, lapply(as.list(unique(NeutralityR2_res$reset_name)), NeutralIdentical))
NeutralStatus <- table(NeutralStatus$Neutral_status) %>% as.data.frame()
colnames(NeutralStatus) <- c("NeutralStatus","Num")


myLabel <- NeutralStatus$Num
myLabel <- paste0(myLabel, "(", round(NeutralStatus$Num / sum(NeutralStatus$Num ) * 100, 2), "%)") 

pdf(paste0(workdir, "output/evolution/neutralStatus.pdf"), width = 5,height = 5)
ggthemr("fresh")
ggplot(NeutralStatus, aes(x = "",y = Num, fill = NeutralStatus)) +
  geom_bar(stat="identity",width=1) +
  coord_polar(theta="y",direction=1)+
  scale_fill_manual(values = c("#BC3C29FF","#E18727FF","#20854EFF"),
                    breaks = c("non_Neutral only","Different","Neutral only"))+
  labs(x="",y="",fill="NeutralStatus")+
  ggtitle(label ="multi-samples neutral evolution Status",subtitle=NULL)+
  guides(fill = guide_legend(reverse = T)) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank())+
  geom_text(aes(x=1.2,y=(sum(Num)-cumsum(Num)+Num/2),label=myLabel),
            size=3,color = "white",size = 16, face = "bold")
ggthemr_reset()
dev.off()

## box plot---------
evol_MATH <- merge(NeutralityR2_res, math_res_rmCNV, by = "reset_SamLocation")

pdf(paste0(workdir,"output/evolution/neutral_MATH_compare.pdf"),width = 5, height = 7)
ggplot(evol_MATH, aes(x=Neutral_status,y=MATH,color=Neutral_status))+
  geom_boxplot(outlier.size=0, size=0.9, width=0.6,fill="white")+
  geom_point(size=5)+
  scale_color_manual(breaks = c("No","Yes"), values = c("#CE364F","#86A875"))+
  ylab("MATH")+
  scale_x_discrete(breaks = c("No","Yes"), labels = c("non_Neutral","Neutral"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 25),
        plot.title = element_text(size=28,hjust = 0.5),
        axis.text.x = element_text(size=25),
        axis.text.y = element_text(size=25))+
  guides(color=F)
dev.off()

neutral_math <- evol_MATH %>% dplyr::filter(Neutral_status == "Yes") %>% dplyr::select(MATH)
nonneutral_math <- evol_MATH %>% dplyr::filter(Neutral_status == "No") %>% dplyr::select(MATH)

wilcox.test(neutral_math$MATH, nonneutral_math$MATH, paired = F, alternative = "greater")

## survival analysis------
PaNeutralStatus <- do.call(rbind, lapply(as.list(unique(NeutralityR2_res$reset_name)), NeutralIdentical)) %>%
  mutate(Neutral_status_reclas =  ifelse(Neutral_status %in% c("Different","Neutral only"), 
                                         "other", "non_Neutral only"))

surv_CRC <- unique(followup_data) %>%
  merge(., PaNeutralStatus,by.x="reset_name",by.y="Patient")

OS_diff <- survfit(Surv(OS, OS_status)~ Neutral_status, data = surv_CRC)
PFS_diff <- survfit(Surv(PFS, PFS_status)~ Neutral_status, data = surv_CRC)
OS_diff_reclass <- survfit(Surv(OS, OS_status)~ Neutral_status_reclas, data = surv_CRC)
PFS_diff_reclass <- survfit(Surv(PFS, PFS_status)~ Neutral_status_reclas, data = surv_CRC)

# 3 groups
splots <- list()
show_col(pal_nejm("default")(7))
splots[[1]] <- ggsurvplot(OS_diff, data = surv_CRC,
                          fun = "pct",conf.int = F, pval = T,
                          risk.table = T,
                          palette = c("#E18727FF","#20854EFF","#BC3C29FF"),
                          legend.labs = c("Different","Neutral only","non_Neutral only"),
                          risk.table.y.text.col = T, 
                          risk.table.y.text = F, 
                          xlab = "OS (months)")


splots[[2]] <- ggsurvplot(PFS_diff, data = surv_CRC,
                          fun = "pct",conf.int = F, pval = T,
                          risk.table = T,
                          palette = c("#E18727FF","#20854EFF","#BC3C29FF"),
                          legend.labs = c("Different","Neutral only","non_Neutral only"),
                          risk.table.y.text.col = T, 
                          risk.table.y.text = F, 
                          xlab = "PFS (months)")
surres <- arrange_ggsurvplots(splots, print = T, ncol = 2, nrow = 1, risk.table.height = 0.3)
ggsave(paste0(workdir,"output/evolution/survival_3_neutral_pattern.pdf"),
       plot = surres, width = 7, height = 5)

# 2 groups
splots_reclass <- list()
splots_reclass[[1]] <- ggsurvplot(OS_diff_reclass, data = surv_CRC,
                                  fun = "pct",conf.int = F, pval = T,
                                  risk.table = T,
                                  palette = c("#BC3C29FF","#7876B1FF"),
                                  legend.labs = c("Non-Neutral only","other"),
                                  risk.table.y.text.col = T, 
                                  risk.table.y.text = F, 
                                  xlab = "OS (months)")


splots_reclass[[2]] <- ggsurvplot(PFS_diff_reclass, data = surv_CRC,
                                  fun = "pct",conf.int = F, pval = T,
                                  risk.table = T,
                                  palette = c("#BC3C29FF","#7876B1FF"),
                                  legend.labs = c("Non-Neutral only","other"),
                                  risk.table.y.text.col = T, 
                                  risk.table.y.text = F, 
                                  xlab = "PFS (months)")
surres_reclass <- arrange_ggsurvplots(splots_reclass, print = T, ncol = 2, nrow = 1, risk.table.height = 0.3)
ggsave(paste0(workdir,"output/evolution/survival_2_neutral_pattern.pdf"),
       plot = surres_reclass, width = 7, height = 5)


