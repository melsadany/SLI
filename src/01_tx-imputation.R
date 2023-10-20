################################################################################
#                   imputing tissue transcriptome for SLI data                 #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SLI"
setwd(project.dir)
################################################################################
source("/Dedicated/jmichaelson-wdata/msmuhammad/pipelines/tx-imputation/01_tx-imputation-functions.R")
################################################################################
genotypes.path <- "/Dedicated/jmichaelson-sdata/devGenes/merged_arrays/data/gsa2020_psycharray_merged"
# subset genotypes to keep the ones with known weight or eQTL effect
subset_genotypes(genotypes_path_base_name = genotypes.path, 
                 # tissue_type = "all",
                 tissue_type = c("Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"),
                 project_dir = project.dir, 
                 verbose = T, 
                 celltype = T,
                 celltype_type = c("Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb"),
                 tissue = T)
################################################################################
# impute transcriptome and save it
impute.tx(project_dir = project.dir, 
          # tissue_type = "all",
          tissue_type = c("Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"),
          celltype = T, 
          tissue = T, 
          celltype_type = c("Inhibitory", "Microglia", "Oligodendrocytes", "OPCs", "Pericytes", "pb"),
          verbose = T)

################################################################################



################################################################################