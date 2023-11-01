################################################################################
#        checking correlations between imputed-tx and SLI phenotypes           #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggh4x);library(ggpubr);library(ggrepel);library(RGCCA)
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
# files <- list.files("data/derivatives/gene-by-factor-correlations", pattern = "Brain")
# ge.fa.corr <- foreach(i=1:length(files), .combine = rbind) %dopar% {
#   f <- read_csv(paste0("data/derivatives/gene-by-factor-correlations/", files[i])) %>%
#     mutate(tissue = sub(".csv", "", files[i])) %>%
#     mutate(FDR = p.adjust(pval, method = "fdr")) %>%
#     filter(pval<0.05)
#   return(f)
# }
# write_rds(ge.fa.corr, "data/derivatives/gene-by-factor-correlations/combined-sig-brain-data.rds")
ge.fa.corr <- read_rds("data/derivatives/gene-by-factor-correlations/combined-sig-brain-data.rds")
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
#####
tx.phe.corr <- read_rds("data/derivatives/gene-by-phenotype-correlations/combined-brain-data.rds")
gc()
################################################################################
################################################################################
################################################################################
# # random forest predicting factors scores from imputed gene expression
# library(randomForest)
# factors <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/factors/pheno_factors_resid.csv")
# sample.mappings <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/pheno_age_corrected.csv") %>%
#   select(id, sample = core_id) %>%
#   filter(sample %in% factors$sample)
# #####
# tissue.ls <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
# # i=8
# registerDoMC(cores = 4)
# # build a RF to predict factors from gene expression of all genes per tissue
# # extract feature importance per RF and return it
# rf.importance <- foreach(i=7:16, .combine = rbind) %dopar% {
#   tissue <- tissue.ls[i]
#   tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
#   tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% sample.mappings$id,] %>% 
#     as.data.frame() %>%
#     rownames_to_column("id")
#   rm(tissue.tx.r); gc()
#   ge.names <- data.frame(o = colnames(tissue.tx)[-1]) %>% 
#     mutate(n = sub("-", "_", o))
#   colnames(tissue.tx)[-1] <- ge.names$n
#   tx.fa <- inner_join(factors, inner_join(sample.mappings,tissue.tx)) %>%
#     select(starts_with("Factor"), colnames(tissue.tx)[-1])
#   # Train the Random Forest model
#   rf.model <- randomForest(as.formula(paste0(paste(paste0("Factor", 1:7), collapse = " + "), " ~ .")),
#                            data = tx.fa, 
#                            ntree = 100)
#   # Get feature importance
#   importance <- rf.model$importance %>%
#     as.data.frame() %>%
#     rownames_to_column("gene") %>%
#     mutate(tissue = tissue)
#   # plots
#   # number of trees against MSE
#   p1 <- data.frame(t = 1:100, MSE = rf.model$mse) %>%
#     ggplot(aes(x=t, y=MSE)) +
#     geom_line() +
#     labs(x = "trees", title = tissue)
#   # actual VS. predicted values
#   p2 <- data.frame(o = rf.model$y, predicted = rf.model$predicted) %>%
#     ggplot(aes(x = o, y = predicted)) +
#     geom_point(alpha = 0.6) +
#     geom_smooth(method = "glm") +
#     stat_cor()+
#     labs(x = paste(paste0("Factor", 1:7), collapse = " + "))
#   svg(width = 5.2, height = 6, filename = paste0("figs/RF/", tissue, ".svg"))
#   print(patchwork::wrap_plots(p1,p2, ncol = 1))
#   dev.off()
#   rm(rf.model);rm(tissue.tx);gc()
#   return(importance)
# }
# # histogram of feature importance per tissue for genes
# rf.importance %>%
#   ggplot(aes(x = IncNodePurity)) +
#   geom_histogram(bins = 50) +
#   facet_wrap(~tissue, scales = "free", ncol = 2)
# # keep top 5% of genes per tissue based on their feature importance (95% quantile)
# top5 <- rf.importance %>%
#   group_by(tissue) %>%
#   filter(IncNodePurity > as.numeric(quantile(IncNodePurity, probs = 0.95)))
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# check how many genes of TK list are there with sign gene phenotype correlations
tk.genes <- read.csv("/wdata/common/SLI_WGS/gene-lists/Language-Literatue-Review-Genes.csv")
# plot how many genes of TK list are there in your gene-phenotype sig correlation
right_join(tx.phe.corr %>% rename(gene = V2),
           tk.genes %>% select(gene =1),
           relationship = "many-to-many") %>%
  distinct(V1, gene, tissue, .keep_all = T) %>%
  group_by(tissue, V1) %>%
  drop_na() %>%
  dplyr::summarise(count = n()) %>%
  mutate(grade = sub("_.*", "", V1),
         phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V1)),
         tissue = sub("Brain_", "", tissue)) %>%
  ggplot(aes(x = tissue, y=grade, fill = count, label = ifelse(count>5,count,""))) +
  geom_tile() +
  geom_text(size = 3) +
  facet_grid2(rows = vars(phenotype),
              scales = "free", space = "free")+
  labs(x="", y="") +
  my.guides + scale_fill_gradient2(low = "black", high = "#800000") +
  theme(strip.text.y.right = element_text(angle = 0),
        strip.text.x.top = element_text(angle = 90)
        )

# plot heatmap for correlations found between tx and phenotype for language genes in all brain tissues
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
# # check available TK genes in important genes by RF
# # barplot for how many language genes found in top 5% important genes in RF models per tissue
# p2 <- top5 %>%
#   mutate(TK = ifelse(gene %in% tk.genes$Gene.Region, 1, 0),
#          tissue = sub("Brain_", "", tissue)) %>%
#   group_by(tissue) %>%
#   dplyr::summarise(count = n(), TK = sum(TK)) %>%
#   ggplot(aes(y=tissue, x = TK)) +
#   geom_bar(stat = "identity", width = 0.5) +
#   labs(x = "number of genes in TK language list", y="") +
#   theme(axis.text.y.left = element_blank())
# # barplot for how many genes are the top 5% in RF models per tissue
# p1 <- top5 %>%
#   mutate(TK = ifelse(gene %in% tk.genes$Gene.Region, 1, 0),
#          tissue = sub("Brain_", "", tissue)) %>%
#   group_by(tissue) %>%
#   dplyr::summarise(count = n(), TK = sum(TK)) %>%
#   ggplot(aes(y=tissue, x = count)) +
#   geom_bar(stat = "identity", width = 0.5) +
#   labs(x = "number of top 5% important genes", y="", 
#        title = "Top 5% of important genes per tissue in RF model")
# # plot the importance distribution in a boxplot for RF models per tissue
# # add labels if TK language genes are in these top 5% genes
# d <- top5 %>%
#   mutate(TK = ifelse(gene %in% tk.genes$Gene.Region, 1, 0)%>%as.factor(),
#          tissue = sub("Brain_", "", tissue))
# p3 <- d %>%
#   ggplot(aes(y = tissue, x = IncNodePurity)) +
#   geom_boxplot(width = 0.5, outlier.size = 0) +
#   geom_point(alpha = 0.3, size = 0.4, data = d %>% filter(TK==1), color = "red") +
#   geom_text_repel(aes(label = ifelse(TK==1, gene, "")), 
#                   box.padding = 1, force = 0.25, max.overlaps = Inf, size = 2) +
#   labs(y="") +
#   theme(axis.text.y.left = element_blank())
# patchwork::wrap_plots(p1+theme(axis.title.x.bottom = element_text(size = 7)),
#                       p2+theme(axis.title.x.bottom = element_text(size = 7)),
#                       p3, widths = c(1,1,3))

################################################################################
# check available TK genes in gene factors sig correlations
# plot how many genes of TK list are there in your gene-factor sig correlation
right_join(ge.fa.corr %>% rename(gene = V2),
           tk.genes %>% select(gene =1),
           relationship = "many-to-many") %>%
  distinct(V1, gene, tissue, .keep_all = T) %>%
  group_by(tissue, V1) %>%
  drop_na() %>%
  dplyr::summarise(count = n()) %>%
  mutate(phenotype = V1,
         tissue = sub("Brain_", "", tissue)) %>%
  ggplot(aes(x = tissue, y=phenotype, fill = count, label = ifelse(count>5,count,""))) +
  geom_tile() +
  geom_text(size = 3) +
  labs(x="", y="") +
  my.guides + scale_fill_gradient2(low = "black", high = "#800000")

# plot heatmap for correlations found between tx and factors for language genes in all brain tissues
inner_join(ge.fa.corr %>% rename(gene = V2),
           tk.genes %>% select(gene =1),
           relationship = "many-to-many") %>%
  distinct(V1, gene, tissue, .keep_all = T) %>%
  drop_na() %>%
  mutate(phenotype = V1,
         tissue = sub("Brain_", "", tissue)) %>%
  ggplot(aes(x = tissue, y = phenotype))+
  geom_tile(aes(fill = r)) +
  facet_grid2(axes = "all" ,remove_labels = "all",
              cols = vars(gene),
              scales = "free", space = "free")+
  my.guides + scale_fill_gradient2(low = "black", high = "#800000") +
  labs(x = "", y="", title = "correlation between imputed tx for TK language genes and factors per tissue") +
  theme(strip.text.x.top = element_text(angle = 90),
        axis.text.x.bottom = element_text(size=7))


################################################################################
# do a CCA using phenotypes, PGS, and imputed tx for language genes
phenotypes <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/pheno_imputed.csv")
sli.pgs.raw <- read_csv("/sdata/devGenes/merged_arrays/polygenic_scores/LDPred2_gathered_pgs_pc_corrected_long.csv") %>%
  mutate(pgs_name_b = pgs_name) %>%
  # mutate(pgs_name = str_replace_all(pattern = " ", replacement = "_", string = tolower(pgs_clean_name)))
  mutate(pgs_name = pgs_clean_name)
  
sli.pgs <- sli.pgs.raw %>%
  select(id = IID, pgs_name, pgs_pc_corrected) %>%
  pivot_wider(names_from = "pgs_name", values_from = "pgs_pc_corrected")
tk.genes <- read.csv("/wdata/common/SLI_WGS/gene-lists/Language-Literatue-Review-Genes.csv")
tissue.ls <- c("Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia")
registerDoMC(cores=4)
foreach(j=1:length(tissue.ls)) %dopar% {
  tissue <- tissue.ls[j]
  tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
  tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% phenotypes$id,] %>% 
    as.data.frame() %>%
    rownames_to_column("id") %>%
    select(id, any_of(tk.genes$Gene.Region))
  rm(tissue.tx.r); gc()
  
  # correct the imputed-tx for 20 genetic PCs
  tissue.tx.pc <- inner_join(sli.pgs.raw %>% select(id = IID, starts_with("pc")), 
                                  tissue.tx)
  tissue.tx.pc.corr <- data.frame(id = tissue.tx.pc$id, 
                                  lapply(tissue.tx.pc %>% select(colnames(tissue.tx)[-1]), 
                                         function(x) {
                                           model <- glm(as.formula(paste("y ~", paste(paste0("pc", 1:20), collapse = " + "))), 
                                                        data = data.frame(y = x, tissue.tx.pc%>%select(starts_with("pc"))))
                                           return(residuals(model))
                                           }))
  rm(list = c("tissue.tx.pc"));gc()
  tx.phe.pgs <- inner_join(inner_join(phenotypes, sli.pgs), tissue.tx.pc.corr) %>%
    mutate_at(.vars = -c(1:4), .funs = function(x) scale(x)[,1])
  
  # quick scatterplot check for correlations with genes
  # tx.phe %>%
  #   ggplot(aes(x=scale(HLCS), y = scale(g0_sentence_sentence_imitation))) +
  #   geom_point() +
  #   geom_smooth(method = "glm") +
  #   stat_cor()
  
  # tx.phe.cca.cv <- rgcca_cv(blocks = list(tx.phe.pgs%>%select(colnames(phenotypes)[-c(1:4)]),
  #                                         tx.phe.pgs%>%select(any_of(tk.genes$Gene.Region))), 
  #                           ncomp = 10, scale = F, verbose = T, method = "sgcca", 
  #                           response = 1, validation = "kfold", k = 10, n_cores = 4)
  # p1 <- plot(tx.phe.cca.cv) # cross-validation plot
  
  ################################################
  # CCA for 2 blocks ~ imputed-tx and phenotypes #
  ################################################
  tx.phe.cca <- rgcca(blocks = list(phenotypes = tx.phe.pgs%>%
                                      select(colnames(phenotypes)[-c(1:4)]),
                                    imputed_tx = tx.phe.pgs%>%
                                      select(any_of(tk.genes$Gene.Region))), 
                      ncomp = 10, verbose = T, method = "sgcca", 
                      # sparsity = as.numeric(tx.phe.cca.cv$best_params), 
                      response = 1)
  # needed diagnostic plots
  p1<- plot(tx.phe.cca, type = "ave", cex = 1) # average variance explained for components by block
  # component loadings
  loadings.plots <- list()
  for (i in 1:5) {
    loadings.plots[[i]] <- plot(tx.phe.cca, type = "weight", comp = i, display_order = FALSE, cex = 0.7)
  }
  p2<- patchwork::wrap_plots(loadings.plots[[1]],loadings.plots[[2]], loadings.plots[[3]],loadings.plots[[4]], 
                             loadings.plots[[5]], ncol = 2)
  # look at correlation between components 
  t1 <- data.frame(id= tx.phe.pgs$id,
                   ph_cc = tx.phe.cca$Y$phenotypes,
                   tx_cc = tx.phe.cca$Y$imputed_tx)
  p3 <- t1 %>%
    pivot_longer(cols = starts_with("tx_"), names_to = "tx_cc", values_to = "tx_cc_v") %>%
    pivot_longer(cols = starts_with("ph_"), names_to = "ph_cc", values_to = "ph_cc_v") %>%
    mutate(tx_cc = as.numeric(sub("tx_cc.comp", "", tx_cc)),
           ph_cc = as.numeric(sub("ph_cc.comp", "", ph_cc))) %>%
    ggplot(aes(x=tx_cc_v, y=ph_cc_v)) +
    geom_point(alpha=0.1, size=0.005) +
    geom_smooth(method = "glm")+
    facet_grid2(rows = vars(ph_cc), cols = vars(tx_cc), scales = "free") +
    labs(x="imputed-tx components", y = "phenotypes components") +
    theme(strip.text.y.right = element_text(angle = 0))
  p4 <- corr.table(t1%>%select(starts_with("tx")), t1%>%select(starts_with("ph"))) %>%
    filter(grepl("tx", V1), grepl("ph", V2)) %>%
    ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
    geom_tile()+
    redblu.col.gradient+my.guides+
    geom_text(size=2) +
    null_labs
  # visualize component loadings
  l1 <- data.frame(cc = as.factor(1:10),
                   t(tx.phe.cca$a$phenotypes),
                   t(tx.phe.cca$a$imputed_tx)) %>%
    pivot_longer(cols = colnames(phenotypes)[-c(1:4)], names_to = "ph_var", values_to = "ph_var_lo") %>%
    mutate(grade = sub("_.*", "", ph_var),
           phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", ph_var))) %>%
    pivot_longer(cols = any_of(tk.genes$Gene.Region), names_to = "tx_var", values_to = "tx_var_lo")
  p5_1 <- l1 %>%
    ggplot(aes(x=cc, y=tx_var, fill = tx_var_lo)) +
    geom_tile()+
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    null_labs
  p5_2 <- l1 %>%
    ggplot(aes(x=cc, y=grade, fill = ph_var_lo)) +
    geom_tile()+
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
    null_labs +
    theme(strip.text.y.right = element_text(angle = 0)) +
    null_labs
  p5 <- patchwork::wrap_plots(p5_1, p5_2, nrow = 1)
  # save plots in pdf
  # pdf(file = paste0("figs/CCA/tx-phenotypes/", tissue, ".pdf"), width = 13.5, height = 12.5)
  # print(p1);print(p2);print(p3);print(p4);print(p5)
  # dev.off()
  ggsave(p1, filename = paste0("figs/CCA/tx-phenotypes/p1-", tissue, ".png"), 
         bg = "white", width = 9, height = 5, units = "in")
  ggsave(p3, filename = paste0("figs/CCA/tx-phenotypes/p3-", tissue, ".png"), 
         bg = "white", width = 9, height = 9, units = "in")
  ggsave(p4, filename = paste0("figs/CCA/tx-phenotypes/p4-", tissue, ".png"), 
         bg = "white", width = 6, height = 5, units = "in")
  ggsave(p5, filename = paste0("figs/CCA/tx-phenotypes/p5-", tissue, ".png"), 
         bg = "white", width = 10, height = 10, units = "in")
  ######################################################
  # CCA for 3 blocks ~ imputed-tx, PGS, and phenotypes #
  ######################################################
  tx.phe.pgs.cca <- rgcca(blocks = list(phenotypes = tx.phe.pgs %>%
                                          select(colnames(phenotypes)[-c(1:4)]),
                                        imputed_tx = tx.phe.pgs %>%
                                          select(any_of(tk.genes$Gene.Region)),
                                        PGS = tx.phe.pgs %>%
                                          select(colnames(sli.pgs)[-1])), 
                          ncomp = 10, verbose = T, method = "sgcca", 
                          response = 1)
  # needed diagnostic plots
  # p1<- plot(tx.phe.pgs.cca, type = "ave", cex = 1) # average variance explained for components by block
  # # component loadings
  # loadings.plots <- list()
  # for (i in 1:5) {
  #   loadings.plots[[i]] <- plot(tx.phe.pgs.cca, type = "weight", comp = i, display_order = FALSE, cex = 0.7)
  # }
  # p2<- patchwork::wrap_plots(loadings.plots[[1]],loadings.plots[[2]], loadings.plots[[3]],loadings.plots[[4]], 
  #                            loadings.plots[[5]], ncol = 2)
  # look at correlation between components 
  t2 <- data.frame(id= tx.phe.pgs$id,
                   ph_cc = tx.phe.pgs.cca$Y$phenotypes,
                   tx_cc = tx.phe.pgs.cca$Y$imputed_tx,
                   pgs_cc = tx.phe.pgs.cca$Y$PGS)
  # p3_12 <- t2 %>%
  #   pivot_longer(cols = starts_with("tx_"), names_to = "tx_cc", values_to = "tx_cc_v") %>%
  #   pivot_longer(cols = starts_with("ph_"), names_to = "ph_cc", values_to = "ph_cc_v") %>%
  #   mutate(tx_cc = as.numeric(sub("tx_cc.comp", "", tx_cc)),
  #          ph_cc = as.numeric(sub("ph_cc.comp", "", ph_cc))) %>%
  #   ggplot(aes(x=tx_cc_v, y=ph_cc_v)) +
  #   geom_point(alpha=0.1, size=0.005) +
  #   geom_smooth(method = "glm")+
  #   facet_grid2(rows = vars(ph_cc), cols = vars(tx_cc), scales = "free") +
  #   labs(x="imputed-tx components", y = "phenotypes components") +
  #   theme(strip.text.y.right = element_text(angle = 0))
  # p3_13 <- t2 %>%
  #   pivot_longer(cols = starts_with("tx_"), names_to = "tx_cc", values_to = "tx_cc_v") %>%
  #   pivot_longer(cols = starts_with("pgs_"), names_to = "pgs_cc", values_to = "pgs_cc_v") %>%
  #   mutate(tx_cc = as.numeric(sub("tx_cc.comp", "", tx_cc)),
  #          pgs_cc = as.numeric(sub("pgs_cc.comp", "", pgs_cc))) %>%
  #   ggplot(aes(x=tx_cc_v, y=pgs_cc_v)) +
  #   geom_point(alpha=0.1, size=0.005) +
  #   geom_smooth(method = "glm")+
  #   facet_grid2(rows = vars(pgs_cc), cols = vars(tx_cc), scales = "free") +
  #   labs(x="imputed-tx components", y = "PGS components") +
  #   theme(strip.text.y.right = element_text(angle = 0))
  # p3_23 <- t2 %>%
  #   pivot_longer(cols = starts_with("pgs_"), names_to = "pgs_cc", values_to = "pgs_cc_v") %>%
  #   pivot_longer(cols = starts_with("ph_"), names_to = "ph_cc", values_to = "ph_cc_v") %>%
  #   mutate(pgs_cc = as.numeric(sub("pgs_cc.comp", "", pgs_cc)),
  #          ph_cc = as.numeric(sub("ph_cc.comp", "", ph_cc))) %>%
  #   ggplot(aes(x=pgs_cc_v, y=ph_cc_v)) +
  #   geom_point(alpha=0.1, size=0.005) +
  #   geom_smooth(method = "glm")+
  #   facet_grid2(rows = vars(ph_cc), cols = vars(pgs_cc), scales = "free") +
  #   labs(x="PGS components", y = "phenotypes components") +
  #   theme(strip.text.y.right = element_text(angle = 0))
  # p4_12 <- corr.table(t2%>%select(starts_with("tx")), t2%>%select(starts_with("ph"))) %>%
  #   filter(grepl("tx", V1), grepl("ph", V2)) %>%
  #   ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   redblu.col.gradient+my.guides+
  #   geom_text(size=2) +
  #   null_labs
  # p4_13 <- corr.table(t2%>%select(starts_with("tx")), t2%>%select(starts_with("pgs"))) %>%
  #   filter(grepl("tx", V1), grepl("pgs", V2)) %>%
  #   ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   redblu.col.gradient+my.guides+
  #   geom_text(size=2) +
  #   null_labs
  # p4_23 <- corr.table(t2%>%select(starts_with("pgs")), t2%>%select(starts_with("ph"))) %>%
  #   filter(grepl("pgs", V1), grepl("ph", V2)) %>%
  #   ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   redblu.col.gradient+my.guides+
  #   geom_text(size=2) +
  #   null_labs
  p4 <- patchwork::wrap_plots(p4_12, p4_13, p4_23)
  # visualize component loadings
  l2 <- cbind(cc = as.factor(1:10),
                   t(tx.phe.pgs.cca$a$phenotypes),
                   t(tx.phe.pgs.cca$a$imputed_tx),
                   t(tx.phe.pgs.cca$a$PGS)) %>%
    as.data.frame() %>%
    mutate(cc = as.factor(cc)) %>%
    pivot_longer(cols = colnames(phenotypes)[-c(1:4)], names_to = "ph_var", values_to = "ph_var_lo") %>%
    mutate(grade = sub("_.*", "", ph_var),
           phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", ph_var))) %>%
    pivot_longer(cols = any_of(tk.genes$Gene.Region), names_to = "tx_var", values_to = "tx_var_lo") %>%
    pivot_longer(cols = colnames(sli.pgs)[-1], names_to = "pgs_var", values_to = "pgs_var_lo")
    
  p5_1 <- l2 %>%
    ggplot(aes(x=cc, y=tx_var, fill = tx_var_lo)) +
    geom_tile()+
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    null_labs
  p5_2 <- l2 %>%
    ggplot(aes(x=cc, y=grade, fill = ph_var_lo)) +
    geom_tile()+
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
    null_labs +
    theme(strip.text.y.right = element_text(angle = 0))
  p5_3 <- l2 %>%
    ggplot(aes(x=cc, y=pgs_var, fill = pgs_var_lo)) +
    geom_tile()+
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    null_labs
  p5 <- patchwork::wrap_plots(p5_1, p5_2, p5_3, nrow = 1)
  # save plots in pdf
  # pdf(file = paste0("figs/CCA/tx-phenotypes-pgs/", tissue, ".pdf"), width = 13.5, height = 12.5)
  # print(p1);print(p2);print(p3_12);print(p3_13);print(p3_23);print(p4_12);print(p4_13);print(p4_23);print(p5)
  # dev.off()
  # ggsave(p1, filename = paste0("figs/CCA/tx-phenotypes-pgs/p1-", tissue, ".png"), 
  #        bg = "white", width = 9, height = 5, units = "in")
  # ggsave(p3_12, filename = paste0("figs/CCA/tx-phenotypes-pgs/p3_12-", tissue, ".png"), 
  #        bg = "white", width = 9, height = 9, units = "in")
  # ggsave(p3_13, filename = paste0("figs/CCA/tx-phenotypes-pgs/p3_13-", tissue, ".png"), 
  #        bg = "white", width = 9, height = 9, units = "in")
  # ggsave(p3_23, filename = paste0("figs/CCA/tx-phenotypes-pgs/p3_23-", tissue, ".png"), 
  #        bg = "white", width = 9, height = 9, units = "in")
  # ggsave(p4, filename = paste0("figs/CCA/tx-phenotypes-pgs/p4-", tissue, ".png"), 
  #        bg = "white", width = 10.2, height = 5, units = "in")
  ggsave(p5, filename = paste0("figs/CCA/tx-phenotypes-pgs/p5-", tissue, ".png"), 
         bg = "white", width = 11, height = 11, units = "in")
  #################################################
  # saving data ~ CCA loadings and sample weights #
  #################################################
  write_csv(t1, file = paste0("data/derivatives/CCA/tx-phenotypes/sample-cc-", tissue, ".csv"))
  write_csv(t2, file = paste0("data/derivatives/CCA/tx-phenotypes-pgs/sample-cc-", tissue, ".csv"))
  write_csv(l1, file = paste0("data/derivatives/CCA/tx-phenotypes/cc-loadings-", tissue, ".csv"))
  write_csv(l2, file = paste0("data/derivatives/CCA/tx-phenotypes-pgs/cc-loadings-", tissue, ".csv"))
  ####################
  # delete and clean #
  ####################
  rm(list = c("p1","p2","p3","p4","p5","loadings.plots", "tx.phe.cca", "tx.phe.pgs.cca","tissue.tx",
              "tx.phe.pgs", "t1", "t2", "l1", "l2", "p3_12", "p3_13", "p3_23", "p4_12", "p4_13", "p4_23",
              "p5_1", "p5_2", "p5_3", "tissue.tx.pc.corr"))
  gc()
}


################################################################################




