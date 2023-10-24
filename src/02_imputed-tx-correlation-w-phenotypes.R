################################################################################
#        checking correlations between imputed-tx and SLI phenotypes           #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggh4x);library(ggpubr);library(ggrepel)
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SLI"
setwd(project.dir)
################################################################################
factors <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors_resid.csv")
sample.mappings <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/pheno_age_corrected.csv") %>%
  select(id, sample = core_id) %>%
  filter(sample %in% factors$sample)
################################################################################
# tissue.ls <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
# # tissue.ls <- c("Excitatory", "Astrocytes", "Endothelial", "Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb")
# # i=8
# tissues.output <- list()
# registerDoMC(cores = 6)
# foreach(i=1:length(tissue.ls)) %dopar% {
#   tissue <- tissue.ls[i]
#   tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
#   tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% sample.mappings$id,] %>% 
#     as.data.frame() %>%
#     rownames_to_column("id")
#   rm(tissue.tx.r); gc()
#   tx.fa <- inner_join(factors, inner_join(sample.mappings,tissue.tx))
#   registerDoMC(cores = 4)
#   
#   # build the model as gene_expression ~ all factors
#   # tissue.output <- foreach(j=1:ncol(tx.fa), .combine = rbind) %dopar% {
#   #   gene <- colnames(tx.fa)[j]
#   #   if (gene %in% colnames(tissue.tx)[-1]) {
#   #     model <- glm(y ~ Factor1 + Factor2 + Factor3 + Factor4 + Factor5 + Factor6 + Factor7, 
#   #                 data = tx.fa %>% 
#   #                   select(sample, id, starts_with("Factor"), y = gene) %>%
#   #                   mutate(y = scale(y, scale = T, center = T)[,1]))
#   #     df <- summary(model)$coefficients %>%
#   #       as.data.frame() %>%
#   #       rownames_to_column("variable") %>%
#   #       mutate(gene = gene)
#   #   }else {
#   #     return(NULL)
#   #   }
#   # }
#   # write_csv(tissue.output, paste0("data/derivatives/factors-predicting-gene/", tissue, ".csv"))
#   # factor.output <- foreach(j=1:7, .combine = rbind) %dopar% {
#   #   factor <- paste0("Factor", j)
#   #   model <- glm(y ~ ., 
#   #                data = tx.fa %>% 
#   #                  select(factor, colnames(tissue.tx)[-1]) %>%
#   #                  mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]) %>%
#   #                  rename(y = 1))
#   #   df <- summary(model)$coefficients %>%
#   #     as.data.frame() %>%
#   #     rownames_to_column("variable") %>%
#   #     mutate(factor = factor)
#   # }
#   # write_csv(factor.output, paste0("data/derivatives/genes-predicting-factor/", tissue, ".csv"))
#   
#   # build the model as Factor_score ~ all genes expression
#   factor.output <- corr.func(tx.fa %>% 
#                                select(starts_with("Factor")) %>%
#                                mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]), 
#                              tx.fa %>% 
#                                select(colnames(tissue.tx)[-1]) %>%
#                                mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]))
#   write_csv(factor.output, paste0("data/derivatives/gene-by-factor-correlations/", tissue, ".csv"))
#   
# }
# 
# ################################################################################
# # factors predicting genes
# files <- list.files("data/derivatives/factors-predicting-gene")
# tx.fa.corr <- foreach(i=1:length(files), .combine = rbind) %dopar% {
#   f <- read_csv(paste0("data/derivatives/factors-predicting-gene/", files[i])) %>% 
#     mutate(tissue = sub(".csv", "", files[i]))
#   return(f)
# }
# tx.fa.corr.filt <- tx.fa.corr %>%
#   filter(grepl("Brain", tissue)) %>%
#   group_by(tissue) %>%
#   mutate(FDR = p.adjust(`Pr(>|t|)`))
# ################################################################################
# # gene by factor correlations
# files <- list.files("data/derivatives/gene-by-factor-correlations")
# ge.fa.corr <- foreach(i=1:length(files), .combine = rbind) %dopar% {
#   f <- read_csv(paste0("data/derivatives/gene-by-factor-correlations/", files[i])) %>% 
#     mutate(tissue = sub(".csv", "", files[i]))
#   return(f)
# }
# ge.fa.corr.filt <- ge.fa.corr %>%
#   filter(grepl("Brain", tissue), pval < 0.05) %>%
#   group_by(tissue) %>%
#   mutate(FDR = p.adjust(pval))
################################################################################
################################################################################
################################################################################
################################################################################
# get correlations between imputed tx and raw phenotypes
# independent section
phenotypes <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/pheno_imputed.csv")
tissue.ls <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
# tissue.ls <- c("Excitatory", "Astrocytes", "Endothelial", "Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb")
registerDoMC(cores=4)
foreach(i=1:length(tissue.ls)) %dopar% {
  tissue <- tissue.ls[i]
  tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
  tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% phenotypes$id,] %>% 
    as.data.frame() %>%
    rownames_to_column("id")
  rm(tissue.tx.r); gc()
  tx.phe <- inner_join(phenotypes, tissue.tx)
  tx.phe.corr <- corr.table(tx.phe %>% 
                             select(colnames(phenotypes)[-c(1:4)]) %>%
                             mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]),
                           tx.phe %>% 
                             select(colnames(tissue.tx)[-1]) %>%
                             mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1])) %>%
    filter(V1 %in% colnames(phenotypes), !V2 %in% colnames(phenotypes)) %>%
    mutate(FDR = p.adjust(pval, method = "fdr"))
  # write_csv(tx.phe.corr, paste0("data/derivatives/gene-by-phenotype-correlations/", tissue, ".csv"))
  # pdssave(tx.phe.corr, file = paste0("data/derivatives/gene-by-phenotype-correlations/", tissue, ".rds"))
  pdssave(tx.phe.corr %>% filter(pval < 0.05), 
          file = paste0("data/derivatives/gene-by-phenotype-correlations/", tissue, "-sig.rds"))
}

files <- list.files("data/derivatives/gene-by-phenotype-correlations", pattern = "Brain")
tx.phe.corr <- foreach(i=1:length(files), .combine = rbind) %dopar% {
  f <- pdsload(fname = paste0("data/derivatives/gene-by-phenotype-correlations/", files[i])) %>% 
    mutate(tissue = sub("-sig.rds.pxz", "", files[i]))
  return(f)
}
write_rds(tx.phe.corr, "data/derivatives/gene-by-phenotype-correlations/combined-brain-data.rds")
tx.phe.corr <- read_rds("data/derivatives/gene-by-phenotype-correlations/combined-brain-data.rds")
gc()
################################################################################
################################################################################
################################################################################
# random forest predicting factos scores from imputed gene expression
library(randomForest)
factors <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors_resid.csv")
sample.mappings <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/pheno_age_corrected.csv") %>%
  select(id, sample = core_id) %>%
  filter(sample %in% factors$sample)
#####
tissue.ls <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
# i=8
registerDoMC(cores = 4)
rf.importance <- foreach(i=7:16, .combine = rbind) %dopar% {
  tissue <- tissue.ls[i]
  tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
  tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% sample.mappings$id,] %>% 
    as.data.frame() %>%
    rownames_to_column("id")
  rm(tissue.tx.r); gc()
  ge.names <- data.frame(o = colnames(tissue.tx)[-1]) %>% 
    mutate(n = sub("-", "_", o))
  colnames(tissue.tx)[-1] <- ge.names$n
  tx.fa <- inner_join(factors, inner_join(sample.mappings,tissue.tx)) %>%
    select(starts_with("Factor"), colnames(tissue.tx)[-1])
  # Train the Random Forest model
  rf.model <- randomForest(as.formula(paste0(paste(paste0("Factor", 1:7), collapse = " + "), " ~ .")),
                           data = tx.fa, 
                           ntree = 100)
  # Get feature importance
  importance <- rf.model$importance %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    mutate(tissue = tissue)
  # plots
  # number of trees against MSE
  p1 <- data.frame(t = 1:100, MSE = rf.model$mse) %>%
    ggplot(aes(x=t, y=MSE)) +
    geom_line() +
    labs(x = "trees", title = tissue)
  p2 <- data.frame(o = rf.model$y, predicted = rf.model$predicted) %>%
    ggplot(aes(x = o, y = predicted)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "glm") +
    stat_cor()+
    labs(x = paste(paste0("Factor", 1:7), collapse = " + "))
  svg(width = 5.2, height = 6, filename = paste0("figs/RF/", tissue, ".svg"))
  print(patchwork::wrap_plots(p1,p2, ncol = 1))
  dev.off()
  rm(rf.model);rm(tissue.tx);gc()
  return(importance)
}
rf.importance %>%
  ggplot(aes(x = IncNodePurity)) +
  geom_histogram(bins = 50) +
  facet_wrap(~tissue, scales = "free", ncol = 2)
top5 <- rf.importance %>%
  group_by(tissue) %>%
  filter(IncNodePurity > as.numeric(quantile(IncNodePurity, probs = 0.95)))
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# check how many gene of Tanner's langugae list are there with significant pearson corr
tk.genes <- read.csv("/wdata/common/SLI_WGS/gene-lists/Language-Literatue-Review-Genes.csv")
right_join(tx.phe.corr %>% rename(gene = V2),
           tk.genes %>% select(gene =1),
           relationship = "many-to-many") %>%
  # tx.phe.corr %>% rename(gene = V2) %>%
  distinct(V1, gene, tissue, .keep_all = T) %>%
  group_by(tissue, V1) %>%
  drop_na() %>%
  dplyr::summarise(count = n()) %>%
  mutate(grade = sub("_.*", "", V1),
         phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V1)),
         tissue = sub("Brain_", "", tissue)) %>%
  ggplot(aes(x = tissue, y=grade, fill = count, label = ifelse(count>5,count,""))) +
  # geom_bar(width = 0.5) +
  geom_tile() +
  geom_text(size = 3) +
  # geom_bar(stat = "identity") +
  facet_grid2(rows = vars(phenotype),
              scales = "free", space = "free")+
  labs(x="", y="") +
  my.guides + scale_fill_gradient2(low = "black", high = "#800000") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.text.x.top = element_text(angle = 90)
        )
length(unique(tk.genes$Gene.Region))

inner_join(tx.phe.corr %>% rename(gene = V2),
           tk.genes %>% select(gene =1),
           relationship = "many-to-many") %>%
  distinct(V1, gene, tissue, .keep_all = T) %>%
  drop_na() %>%
  mutate(grade = sub("_.*", "", V1),
         phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V1)),
         tissue = sub("Brain_", "", tissue)) %>%
  ggplot(aes(x = tissue, y = grade))+
  geom_tile(aes(fill = r)) +
  facet_grid2(rows = vars(phenotype),axes = "all" ,remove_labels = "all",
              cols = vars(gene),
              scales = "free", space = "free")+
  my.guides + scale_fill_gradient2(low = "black", high = "#800000") +
  labs(x = "", y="", title = "correlation between imputed tx for TK language genes and phenotypes per tissue") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.text.x.top = element_text(angle = 90),
        axis.text.x.bottom = element_text(size=7)
  )
################################################################################
# check available TK genes in important genes by RF
p2 <- top5 %>%
  mutate(TK = ifelse(gene %in% tk.genes$Gene.Region, 1, 0),
         tissue = sub("Brain_", "", tissue)) %>%
  group_by(tissue) %>%
  dplyr::summarise(count = n(), TK = sum(TK)) %>%
  ggplot(aes(y=tissue, x = TK)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "number of genes in TK language list", y="") +
  theme(axis.text.y.left = element_blank())
p1 <- top5 %>%
  mutate(TK = ifelse(gene %in% tk.genes$Gene.Region, 1, 0),
         tissue = sub("Brain_", "", tissue)) %>%
  group_by(tissue) %>%
  dplyr::summarise(count = n(), TK = sum(TK)) %>%
  ggplot(aes(y=tissue, x = count)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(x = "number of top 5% important genes", y="", 
       title = "Top 5% of important genes per tissue in RF model")
d <- top5 %>%
  mutate(TK = ifelse(gene %in% tk.genes$Gene.Region, 1, 0)%>%as.factor(),
         tissue = sub("Brain_", "", tissue))
p3 <- d %>%
  ggplot(aes(y = tissue, x = IncNodePurity)) +
  geom_boxplot(width = 0.5, outlier.size = 0) +
  geom_point(alpha = 0.3, size = 0.4, data = d %>% filter(TK==1), color = "red") +
  # geom_text(size = 2)+
  # scale_color_manual(values = c("black","red")) +
  geom_text_repel(aes(label = ifelse(TK==1, gene, "")), 
                  box.padding = 1, force = 0.25, max.overlaps = Inf, size = 2) +
  labs(y="") +
  theme(axis.text.y.left = element_blank())
patchwork::wrap_plots(p1+theme(axis.title.x.bottom = element_text(size = 7)),
                      p2+theme(axis.title.x.bottom = element_text(size = 7)),
                      p3, widths = c(1,1,3))

################################################################################


################################################################################

################################################################################
