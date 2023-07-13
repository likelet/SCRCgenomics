## Tengjia Jiang, 11.07.2023
## mutation signature

library(magrittr)
library(dplyr)
library(tidyr)
library(NMF)
library(maftools)
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library(pheatmap)
library(BSgenome)
library(plyr)
library(wesanderson)
library(MutationalPatterns)


workdir <- "./"
laml_clin <- readRDS(paste0(workdir,"input/laml_clin.Rds"))
mut_mat <- readRDS(paste0(workdir,"input/NMF_rank_all.Rds"))
estimate <- readRDS(paste0(workdir,"input/estimate_all.Rds"))


pdf(paste0(workdir,"output/signature/NMF_rank_all.pdf"))
plot(estimate)
dev.off()

## extract signature
k.best <- 3
nmf_res <- extract_signatures(mut_mat, rank = k.best, nrun = 10)
colnames(nmf_res$signatures) <- paste0("Signature_",LETTERS[1:k.best])
rownames(nmf_res$contribution) <- paste0("Signature_",LETTERS[1:k.best])

## visualize signature
pdf(paste0(workdir,"output/signature/signature_all.pdf"), width = 5, height = 2.5)
plot_96_profile(nmf_res$signatures)
dev.off()


## similarity with COSMIC signature
pattern.og30.cosm <- compareSignatures(nmfRes = nmf_res, sig_db = "legacy")


## signature barplot-----
sig_prop <- nmf_res$contribution %>% 
  as.data.frame() %>% 
  dplyr::mutate(Signature = rownames(.)) %>%
  reshape2::melt(., id.vars = "Signature") %>%
  ddply("variable", transform, percent_freq = round(value / sum(value) * 100,2)) %>%
  ddply("variable", transform, percent_freq_label = round(cumsum(percent_freq),2)) %>%
  merge(., unique(laml_clin[c("reset_SamLocation","Primary_loc","reset_name")]), by.x = "variable", by.y = "reset_SamLocation", all = T) 
sig_prop$reset_name <- sub(pattern = "P",replacement = "",x = sig_prop$reset_name)
sig_prop$reset_name <- as.integer(sig_prop$reset_name)
sig_prop <- sig_prop %>%
  dplyr::mutate(Primary_loc = as.character(Primary_loc)) %>%
  dplyr::mutate(Primary_loc = ifelse(is.na(Primary_loc), "Multi", Primary_loc)) %>%
  dplyr::filter(!is.na(Signature)) %>%
  dplyr::arrange(reset_name)

sig_prop$variable <- factor(sig_prop$variable, levels = unique(sig_prop$variable))


pdf(paste0(workdir, "output/signature/signature_bar_percentage_left_right_all.pdf"), width = 10,height = 3)
ggplot(sig_prop, aes(x=variable, y=percent_freq,  fill = Signature)) + 
  geom_bar(stat = "identity", width = 1, color = "#393e46")+
  scale_fill_manual(
    name = NULL,
    values = c("#82CADE","#FE9E83","#DEDBB9"),
    breaks = c("Signature_A","Signature_B","Signature_C")
  ) +
  scale_y_continuous(expand = c(0,0))+
  theme(axis.text.x = element_text(size = 7, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15))+
  labs(x = NULL, y = "Fraction of signature (%)", fill = NULL) 
dev.off()


