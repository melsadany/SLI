################################################################################
#        checking correlations between imputed-tx and SLI phenotypes           #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
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
tissue.ls <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
# tissue.ls <- c("Excitatory", "Astrocytes", "Endothelial", "Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb")
# i=8
tissues.output <- list()
registerDoMC(cores = 6)
foreach(i=1:length(tissue.ls)) %dopar% {
  tissue <- tissue.ls[i]
  tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
  tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% sample.mappings$id,] %>% 
    as.data.frame() %>%
    rownames_to_column("id")
  rm(tissue.tx.r); gc()
  tx.fa <- inner_join(factors, inner_join(sample.mappings,tissue.tx))
  registerDoMC(cores = 4)
  # tissue.output <- foreach(j=1:ncol(tx.fa), .combine = rbind) %dopar% {
  #   gene <- colnames(tx.fa)[j]
  #   if (gene %in% colnames(tissue.tx)[-1]) {
  #     model <- glm(y ~ Factor1 + Factor2 + Factor3 + Factor4 + Factor5 + Factor6 + Factor7, 
  #                 data = tx.fa %>% 
  #                   select(sample, id, starts_with("Factor"), y = gene) %>%
  #                   mutate(y = scale(y, scale = T, center = T)[,1]))
  #     df <- summary(model)$coefficients %>%
  #       as.data.frame() %>%
  #       rownames_to_column("variable") %>%
  #       mutate(gene = gene)
  #   }else {
  #     return(NULL)
  #   }
  # }
  # write_csv(tissue.output, paste0("data/derivatives/factors-predicting-gene/", tissue, ".csv"))
  # factor.output <- foreach(j=1:7, .combine = rbind) %dopar% {
  #   factor <- paste0("Factor", j)
  #   model <- glm(y ~ ., 
  #                data = tx.fa %>% 
  #                  select(factor, colnames(tissue.tx)[-1]) %>%
  #                  mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]) %>%
  #                  rename(y = 1))
  #   df <- summary(model)$coefficients %>%
  #     as.data.frame() %>%
  #     rownames_to_column("variable") %>%
  #     mutate(factor = factor)
  # }
  # write_csv(factor.output, paste0("data/derivatives/genes-predicting-factor/", tissue, ".csv"))
  factor.output <- corr.func(tx.fa %>% 
                               select(starts_with("Factor")) %>%
                               mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]), 
                             tx.fa %>% 
                               select(colnames(tissue.tx)[-1]) %>%
                               mutate_all(.funs = function(x) scale(x, scale = T, center = T)[,1]))
  write_csv(factor.output, paste0("data/derivatives/gene-by-factor-correlations/", tissue, ".csv"))
}

################################################################################
# factors predicting genes
files <- list.files("data/derivatives/factors-predicting-gene")
tx.fa.corr <- foreach(i=1:length(files), .combine = rbind) %dopar% {
  f <- read_csv(paste0("data/derivatives/factors-predicting-gene/", files[i])) %>% 
    mutate(tissue = sub(".csv", "", files[i]))
  return(f)
}
tx.fa.corr.filt <- tx.fa.corr %>%
  filter(grepl("Brain", tissue)) %>%
  group_by(tissue) %>%
  mutate(FDR = p.adjust(`Pr(>|t|)`))
################################################################################
# gene by factor correlations
files <- list.files("data/derivatives/gene-by-factor-correlations")
ge.fa.corr <- foreach(i=1:length(files), .combine = rbind) %dopar% {
  f <- read_csv(paste0("data/derivatives/gene-by-factor-correlations/", files[i])) %>% 
    mutate(tissue = sub(".csv", "", files[i]))
  return(f)
}
ge.fa.corr.filt <- ge.fa.corr %>%
  filter(grepl("Brain", tissue)) %>%
  group_by(tissue) %>%
  mutate(FDR = p.adjust(pval))
################################################################################



################################################################################

