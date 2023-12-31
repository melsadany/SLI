################################################################################
#                    doing CCA for imputed-tx, PGS, and phenotypes             #
################################################################################
rm(list = ls())
gc()
source("/Dedicated/jmichaelson-wdata/msmuhammad/msmuhammad-source.R")
library(ggh4x);library(ggpubr);library(ggrepel);library(RGCCA);library(ggdendro);library(grid)
################################################################################
################################################################################
project.dir <- "/Dedicated/jmichaelson-wdata/msmuhammad/projects/SLI"
setwd(project.dir)
################################################################################
# do a CCA using phenotypes, PGS, and imputed tx for language genes
phenotypes <- read_csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/public/phenotype/pheno_imputed.csv")
sli.pgs.raw <- read_csv("/Dedicated/jmichaelson-sdata/devGenes/merged_arrays/polygenic_scores/LDPred2_gathered_pgs_pc_corrected_long.csv") %>%
  mutate(pgs_name_b = pgs_name) %>%
  mutate(pgs_name = pgs_clean_name)

sli.pgs <- sli.pgs.raw %>%
  select(id = IID, pgs_name, pgs_pc_corrected) %>%
  pivot_wider(names_from = "pgs_name", values_from = "pgs_pc_corrected")

NDD.genes <- readxl::read_xlsx("/Dedicated/jmichaelson-wdata/lcasten/tools/ref_data/simons/spark_2022_nature_genetics/41588_2022_1148_MOESM4_ESM.xlsx", 
                               sheet = 2, skip = 1)
tk.genes <- read.csv("/Dedicated/jmichaelson-wdata/common/SLI_WGS/gene-lists/Language-Literatue-Review-Genes.csv")
all.genes <- full_join(NDD.genes%>%select(gene=HGNC),
                       tk.genes%>%select(gene=Gene.Region))

tissue.ls <- c("Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia")
registerDoMC(cores=4)
foreach::foreach(j=1:length(tissue.ls)) %dopar% {
  tissue <- tissue.ls[j]
  # cross validation does not work for tissues 3,4,9 in the 3 blocks CCA
  tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
  tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% phenotypes$id,] %>% 
    as.data.frame() %>%
    rownames_to_column("id") %>%
    select(id, any_of(all.genes$gene))
  rm(tissue.tx.r); gc()
  
  # correct the imputed-tx for 20 genetic PCs
  tissue.tx.pc <- inner_join(sli.pgs.raw %>% 
                               distinct(IID, .keep_all = T) %>% 
                               select(id = IID, starts_with("pc")), 
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
    mutate_at(.vars = -c(1:4), .funs = function(x) scale(x)[,1]) %>%
    select(-core_id)
  ################################################
  # CCA for 2 blocks ~ imputed-tx and phenotypes #
  ################################################
  # tx.phe.cca.cv <- rgcca_cv(blocks = list(phenotypes = tx.phe.pgs%>%
  #                                           select(colnames(phenotypes)[-c(1:4)]),
  #                                         imputed_tx = tx.phe.pgs%>%
  #                                           select(any_of(all.genes$gene))), 
  #                           ncomp = 10, verbose = T, method = "sgcca",
  #                           response = 1, validation = "kfold", k = 10, n_cores = 4)
  # tx.phe.cca <- rgcca(tx.phe.cca.cv)
  # ####
  # write_rds(tx.phe.cca, file = paste0("data/derivatives/CCA/tx-phenotypes/CCA-", tissue, ".rds"))
  ################################################
  # CCA for 2 blocks ~ pgs and phenotypes #
  ################################################
  # pgs.phe.cca.cv <- rgcca_cv(blocks = list(phenotypes = tx.phe.pgs%>%
  #                                            select(colnames(phenotypes)[-c(1:4)]),
  #                                          pgs = tx.phe.pgs%>%
  #                                            select(colnames(sli.pgs)[-1])), 
  #                           ncomp = 10, verbose = T, method = "sgcca",
  #                           response = 1, validation = "kfold", k = 10, n_cores = 4)
  # pgs.phe.cca <- rgcca(pgs.phe.cca.cv)
  pgs.phe.cca <- rgcca(blocks = list(phenotypes = tx.phe.pgs%>%
                                       select(colnames(phenotypes)[-c(1:4)]),
                                     pgs = tx.phe.pgs%>%
                                       select(colnames(sli.pgs)[-1])), 
                       ncomp = 10, verbose = T, method = "sgcca",
                       response = 1)
  ####
  write_rds(pgs.phe.cca, file = paste0("data/derivatives/CCA/pgs-phenotypes/CCA-", tissue, ".rds"))
  ######################################################
  # CCA for 3 blocks ~ imputed-tx, PGS, and phenotypes #
  ######################################################
  # tx.phe.pgs.cca.cv <- rgcca_cv(blocks = list(phenotypes = tx.phe.pgs %>%
  #                                               select(colnames(phenotypes)[-c(1:4)]),
  #                                             imputed_tx = tx.phe.pgs %>%
  #                                               select(any_of(all.genes$gene)),
  #                                             PGS = tx.phe.pgs %>%
  #                                               select(colnames(sli.pgs)[-1])), 
  #                               ncomp = 10, verbose = T, method = "sgcca", 
  #                               response = 1, validation = "kfold", k = 10, n_cores = 1)
  # tx.phe.pgs.cca <- rgcca(tx.phe.pgs.cca.cv)
  # ####
  # write_rds(tx.phe.pgs.cca, file = paste0("data/derivatives/CCA/tx-phenotypes-pgs/CCA-", tissue, ".rds"))
}

#############
# CCA plots #
#############

foreach(j=1:length(tissue.ls)) %dopar% {
  tissue <- tissue.ls[j]
  # ################################################
  # # CCA for 2 blocks ~ imputed-tx and phenotypes #
  # ################################################
  # tx.phe.cca <- read_rds(paste0("data/derivatives/CCA/tx-phenotypes/CCA-", tissue, ".rds"))
  # tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
  # tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% phenotypes$id,] %>% 
  #   as.data.frame() %>%
  #   rownames_to_column("id") %>%
  #   select(id, any_of(all.genes$gene))
  # rm(tissue.tx.r); gc()
  # 
  # # correct the imputed-tx for 20 genetic PCs
  # tissue.tx.pc <- inner_join(sli.pgs.raw %>% 
  #                              distinct(IID, .keep_all = T) %>% 
  #                              select(id = IID, starts_with("pc")), 
  #                            tissue.tx)
  # tissue.tx.pc.corr <- data.frame(id = tissue.tx.pc$id, 
  #                                 lapply(tissue.tx.pc %>% select(colnames(tissue.tx)[-1]), 
  #                                        function(x) {
  #                                          model <- glm(as.formula(paste("y ~", paste(paste0("pc", 1:20), collapse = " + "))), 
  #                                                       data = data.frame(y = x, tissue.tx.pc%>%select(starts_with("pc"))))
  #                                          return(residuals(model))
  #                                        }))
  # rm(list = c("tissue.tx.pc"));gc()
  # tx.phe.pgs <- inner_join(inner_join(phenotypes, sli.pgs), tissue.tx.pc.corr) %>%
  #   mutate_at(.vars = -c(1:4), .funs = function(x) scale(x)[,1])
  # 
  # ####
  # # needed diagnostic plots
  # ####
  # # average variance explained for components by block
  # p1<- plot(tx.phe.cca, 
  #           type = "ave", 
  #           cex = 1) 
  # # look at correlation between components 
  # t1 <- data.frame(id= tx.phe.pgs$id,
  #                  ph_cc = tx.phe.cca$Y$phenotypes,
  #                  tx_cc = tx.phe.cca$Y$imputed_tx)
  # ####
  # ## correlations scatterplot
  # p3 <- t1 %>%
  #   pivot_longer(cols = starts_with("tx_"), names_to = "tx_cc", values_to = "tx_cc_v") %>%
  #   pivot_longer(cols = starts_with("ph_"), names_to = "ph_cc", values_to = "ph_cc_v") %>%
  #   mutate(tx_cc = factor(sub("tx_cc.comp", "", tx_cc), levels = c(1:10)),
  #          ph_cc = factor(sub("ph_cc.comp", "", ph_cc), levels = c(1:10))) %>%
  #   ggplot(aes(x=tx_cc_v, y=ph_cc_v)) +
  #   geom_point(alpha=0.1, size=0.005) +
  #   geom_smooth(method = "glm")+
  #   facet_grid2(rows = vars(ph_cc), cols = vars(tx_cc), scales = "free") +
  #   labs(x="imputed-tx components", y = "phenotypes components") +
  #   theme(strip.text.y.right = element_text(angle = 0))
  # ####
  # ## correlations heatmap
  # p4 <- corr.table(t1%>%select(starts_with("tx")), t1%>%select(starts_with("ph"))) %>%
  #   filter(grepl("tx", V1), grepl("ph", V2)) %>%
  #   mutate(V1 = factor(sub("tx_cc.comp", "", V1), levels = c(1:10)),
  #          V2 = factor(sub("ph_cc.comp", "", V2), levels = c(1:10))) %>%
  #   ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   redblu.col.gradient+my.guides+
  #   geom_text(size=2) +
  #   labs(x="imputed-tx components",
  #        y="phenotypes componenets")
  # ####
  # # visualize component loadings from the CCA opbject
  # l1 <- data.frame(cc = as.factor(1:10),
  #                  t(tx.phe.cca$a$phenotypes),
  #                  t(tx.phe.cca$a$imputed_tx)) %>%
  #   pivot_longer(cols = colnames(phenotypes)[-c(1:4)], names_to = "ph_var", values_to = "ph_var_lo") %>%
  #   mutate(grade = sub("_.*", "", ph_var),
  #          phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", ph_var))) %>%
  #   pivot_longer(cols = any_of(all.genes$gene), names_to = "tx_var", values_to = "tx_var_lo")
  # p5_1 <- l1 %>%
  #   distinct(cc, tx_var_lo, tx_var) %>%
  #   ggplot(aes(x=cc, y=tx_var, fill = tx_var_lo)) +
  #   geom_tile()+
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   null_labs+labs(title = "loadings from CCA") +
  #   theme(axis.text.y.left = element_text(size = 2))
  # p5_2 <- l1 %>%
  #   distinct(cc, grade, ph_var_lo, phenotype) %>%
  #   ggplot(aes(x=cc, y=grade, fill = ph_var_lo)) +
  #   geom_tile()+
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
  #   null_labs +
  #   theme(strip.text.y.right = element_text(angle = 0)) +
  #   null_labs 
  # p5 <- patchwork::wrap_plots(p5_1, p5_2, nrow = 1)
  # 
  # ####
  # # visualize sample component weights correlation with the raw block data
  # # components correlation with genes
  # b1 <- corr.table(tx.phe.pgs%>%
  #                    select(any_of(all.genes$gene)),
  #                  tx.phe.cca$Y$imputed_tx) %>%
  #   filter(V1 %in% all.genes$gene,
  #          V2 %in% colnames(tx.phe.cca$Y$imputed_tx)) %>%
  #   mutate(V2 = factor(sub("comp", "", V2), levels = as.character(c(1:10))))
  # # components correlation with phenotypes
  # b2 <- corr.table(tx.phe.pgs%>%
  #                    select(colnames(phenotypes)[-c(1:4)]),
  #                  tx.phe.cca$Y$phenotypes) %>%
  #   filter(V1 %in% colnames(phenotypes)[-c(1:4)],
  #          V2 %in% colnames(tx.phe.cca$Y$phenotypes)) %>%
  #   mutate(V2 = factor(sub("comp", "", V2), levels = as.character(c(1:10)))) %>%
  #   mutate(grade = sub("_.*", "", V1),
  #          phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V1))) 
  # # cluster components correlation with genes
  # b1.d <- as.dendrogram(hclust(d = dist(x = b1 %>% pivot_wider(names_from = "V2", values_from = "r", id_cols = "V1")%>%column_to_rownames("V1"))))
  # b1.d2 <- ggdendrogram(b1.d, rotate = T) +
  #   theme(axis.text.y.left = element_text(size=2))
  # order1 <- order.dendrogram(b1.d)
  # 
  # p6_11 <- b1 %>%
  #   mutate(V1 = factor(V1, levels = unique(V1)[order1])) %>%
  #   ggplot(aes(x=V2, y=V1, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   geom_text(size = 1) +
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   null_labs +labs(title = "correlation between CCA values and samples' raw data") +
  #   theme(axis.text.y.left = element_text(size = 2), 
  #         plot.title = element_text(size = 10))
  # # p6_12 <- print(b1.d2, 
  # #                vp = viewport(height = 1))
  # p6_1 <- cowplot::plot_grid(p6_11,
  #                            b1.d2,
  #                            rel_widths = c(2,0.5))
  # 
  # # grid.newpage()
  # # pp <- print(p6_11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1))
  # # pp2 <- print(b1.d2, 
  # #       vp = viewport(x = 0.90, y = 0.51, height =1.035, width = 0.2))
  # 
  # p6_2 <- b2 %>%
  #   ggplot(aes(x=V2, y=grade, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   geom_text(size = 3) +
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
  #   theme(strip.text.y.right = element_text(angle = 0)) +
  #   null_labs
  # p6 <- cowplot::plot_grid(p6_1, p6_2, nrow = 1)
  # ####
  # # save plots in pdf
  # # pdf(file = paste0("figs/CCA/tx-phenotypes/", tissue, ".pdf"), width = 13.5, height = 12.5)
  # # print(p1);print(p2);print(p3);print(p4);print(p5)
  # # dev.off()
  # ggsave(p1, filename = paste0("figs/CCA/tx-phenotypes/p1-", tissue, ".png"), 
  #        bg = "white", width = 9, height = 5, units = "in")
  # ggsave(p3, filename = paste0("figs/CCA/tx-phenotypes/p3-", tissue, ".png"), 
  #        bg = "white", width = 9, height = 9, units = "in")
  # ggsave(p4, filename = paste0("figs/CCA/tx-phenotypes/p4-", tissue, ".png"), 
  #        bg = "white", width = 6, height = 5, units = "in")
  # ggsave(p5, filename = paste0("figs/CCA/tx-phenotypes/p5-", tissue, ".png"), 
  #        bg = "white", width = 10, height = 10, units = "in")
  # ggsave(p6, filename = paste0("figs/CCA/tx-phenotypes/p6-", tissue, ".png"), 
  #        bg = "white", width = 12, height = 17, units = "in")
  ################################################
  # CCA for 2 blocks ~ pgs and phenotypes #
  ################################################
  pgs.phe.cca <- read_rds(paste0("data/derivatives/CCA/pgs-phenotypes/CCA-", tissue, ".rds"))
  tissue.tx.r <- pdsload(paste0("data/derivatives/imputed-tx/", tissue, ".rds.pxz"))
  tissue.tx <- tissue.tx.r[rownames(tissue.tx.r) %in% phenotypes$id,] %>% 
    as.data.frame() %>%
    rownames_to_column("id") %>%
    select(id, any_of(all.genes$gene))
  rm(tissue.tx.r); gc()
  
  # correct the imputed-tx for 20 genetic PCs
  tissue.tx.pc <- inner_join(sli.pgs.raw %>% 
                               distinct(IID, .keep_all = T) %>% 
                               select(id = IID, starts_with("pc")), 
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
  
  ####
  # needed diagnostic plots
  ####
  # average variance explained for components by block
  p1<- plot(pgs.phe.cca, 
            type = "ave", 
            cex = 1) 
  # look at correlation between components 
  t1 <- data.frame(id= tx.phe.pgs$id,
                   ph_cc = pgs.phe.cca$Y$phenotypes,
                   pgs_cc = pgs.phe.cca$Y$pgs)
  ####
  ## correlations scatterplot
  p3 <- t1 %>%
    pivot_longer(cols = starts_with("pgs_"), names_to = "pgs_cc", values_to = "pgs_cc_v") %>%
    pivot_longer(cols = starts_with("ph_"), names_to = "ph_cc", values_to = "ph_cc_v") %>%
    mutate(tx_cc = factor(sub("pgs_cc.comp", "", pgs_cc), levels = c(1:10)),
           ph_cc = factor(sub("ph_cc.comp", "", ph_cc), levels = c(1:10))) %>%
    ggplot(aes(x=pgs_cc_v, y=ph_cc_v)) +
    geom_point(alpha=0.1, size=0.005) +
    geom_smooth(method = "glm")+
    facet_grid2(rows = vars(ph_cc), cols = vars(pgs_cc), scales = "free") +
    labs(x="pgs components", y = "phenotypes components") +
    theme(strip.text.y.right = element_text(angle = 0))
  ####
  ## correlations heatmap
  p4 <- corr.table(t1%>%select(starts_with("pgs")), t1%>%select(starts_with("ph"))) %>%
    filter(grepl("pgs", V1), grepl("ph", V2)) %>%
    mutate(V1 = factor(sub("pgs_cc.comp", "", V1), levels = c(1:10)),
           V2 = factor(sub("ph_cc.comp", "", V2), levels = c(1:10))) %>%
    ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
    geom_tile()+
    redblu.col.gradient+my.guides+
    geom_text(size=2) +
    labs(x="pgs components",
         y="phenotypes componenets")
  ####
  # visualize component loadings from the CCA opbject
  l1 <- data.frame(cc = as.factor(1:10),
                   t(pgs.phe.cca$a$phenotypes),
                   t(pgs.phe.cca$a$pgs)) %>%
    pivot_longer(cols = colnames(phenotypes)[-c(1:4)], names_to = "ph_var", values_to = "ph_var_lo") %>%
    mutate(grade = sub("_.*", "", ph_var),
           phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", ph_var))) %>%
    pivot_longer(cols = any_of(colnames(sli.pgs)[-1]), names_to = "pgs_var", values_to = "pgs_var_lo")
  p5_1 <- l1 %>%
    distinct(cc, pgs_var_lo, pgs_var) %>%
    ggplot(aes(x=cc, y=pgs_var, fill = pgs_var_lo)) +
    geom_tile()+
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    null_labs+labs(title = "loadings from CCA")
  p5_2 <- l1 %>%
    distinct(cc, grade, ph_var_lo, phenotype) %>%
    ggplot(aes(x=cc, y=grade, fill = ph_var_lo)) +
    geom_tile()+
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
    null_labs +
    theme(strip.text.y.right = element_text(angle = 0)) +
    null_labs 
  p5 <- patchwork::wrap_plots(p5_1, p5_2, nrow = 1)
  
  ####
  # visualize sample component weights correlation with the raw block data
  # components correlation with genes
  b1 <- corr.table(tx.phe.pgs%>%
                     select(colnames(sli.pgs)[-1]),
                   pgs.phe.cca$Y$pgs) %>%
    filter(V1 %in% colnames(sli.pgs)[-1],
           V2 %in% colnames(pgs.phe.cca$Y$pgs)) %>%
    mutate(V2 = factor(sub("comp", "", V2), levels = as.character(c(1:10))))
  # components correlation with phenotypes
  b2 <- corr.table(tx.phe.pgs%>%
                     select(colnames(phenotypes)[-c(1:4)]),
                   pgs.phe.cca$Y$phenotypes) %>%
    filter(V1 %in% colnames(phenotypes)[-c(1:4)],
           V2 %in% colnames(pgs.phe.cca$Y$phenotypes)) %>%
    mutate(V2 = factor(sub("comp", "", V2), levels = as.character(c(1:10)))) %>%
    mutate(grade = sub("_.*", "", V1),
           phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V1))) 
  # cluster components correlation with genes
  b1.d <- as.dendrogram(hclust(d = dist(x = b1 %>% pivot_wider(names_from = "V2", values_from = "r", id_cols = "V1")%>%column_to_rownames("V1"))))
  b1.d2 <- ggdendrogram(b1.d, rotate = T) +
    theme(axis.text.y.left = element_text(size=2))
  order1 <- order.dendrogram(b1.d)
  
  p6_11 <- b1 %>%
    mutate(V1 = factor(V1, levels = unique(V1)[order1])) %>%
    ggplot(aes(x=V2, y=V1, fill = r, label = ifelse(pval<0.05, "*", ""))) +
    geom_tile()+
    geom_text(size = 1) +
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    null_labs +labs(title = "correlation between CCA values and samples' raw data") +
    theme(plot.title = element_text(size = 10))
  # p6_12 <- print(b1.d2, 
  #                vp = viewport(height = 1))
  p6_1 <- cowplot::plot_grid(p6_11,
                             b1.d2,
                             rel_widths = c(2,0.5))
  
  # grid.newpage()
  # pp <- print(p6_11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1))
  # pp2 <- print(b1.d2, 
  #       vp = viewport(x = 0.90, y = 0.51, height =1.035, width = 0.2))
  
  p6_2 <- b2 %>%
    ggplot(aes(x=V2, y=grade, fill = r, label = ifelse(pval<0.05, "*", ""))) +
    geom_tile()+
    geom_text(size = 3) +
    my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
    facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
    theme(strip.text.y.right = element_text(angle = 0)) +
    null_labs
  p6 <- cowplot::plot_grid(p6_1, p6_2, nrow = 1, rel_widths = c(2,1.2))
  ####
  # save plots in pdf
  # pdf(file = paste0("figs/CCA/pgs-phenotypes/", tissue, ".pdf"), width = 13.5, height = 12.5)
  # print(p1);print(p2);print(p3);print(p4);print(p5)
  # dev.off()
  ggsave(p1, filename = paste0("figs/CCA/pgs-phenotypes/p1-", tissue, ".png"), 
         bg = "white", width = 9, height = 5, units = "in")
  ggsave(p3, filename = paste0("figs/CCA/pgs-phenotypes/p3-", tissue, ".png"), 
         bg = "white", width = 9, height = 9, units = "in")
  ggsave(p4, filename = paste0("figs/CCA/pgs-phenotypes/p4-", tissue, ".png"), 
         bg = "white", width = 6, height = 5, units = "in")
  ggsave(p5, filename = paste0("figs/CCA/pgs-phenotypes/p5-", tissue, ".png"), 
         bg = "white", width = 10, height = 10, units = "in")
  ggsave(p6, filename = paste0("figs/CCA/pgs-phenotypes/p6-", tissue, ".png"), 
         bg = "white", width = 12, height = 17, units = "in")
  # ######################################################
  # # CCA for 3 blocks ~ imputed-tx, PGS, and phenotypes #
  # ######################################################
  # tx.phe.pgs.cca <- read_rds(paste0("data/derivatives/CCA/tx-phenotypes-pgs/CCA-", tissue, ".rds"))
  # # needed diagnostic plots
  # ####
  # # average variance explained for components by block
  # p1<- plot(tx.phe.pgs.cca, 
  #           type = "ave", 
  #           cex = 1) 
  # # look at correlation between components 
  # t2 <- data.frame(id= tx.phe.pgs$id,
  #                  ph_cc = tx.phe.pgs.cca$Y$phenotypes,
  #                  tx_cc = tx.phe.pgs.cca$Y$imputed_tx,
  #                  pgs_cc = tx.phe.pgs.cca$Y$PGS)
  # ####
  # ## correlation scatterplot
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
  # ####
  # # correlation heatmap
  # p4_12 <- corr.table(t2%>%select(starts_with("tx")), t2%>%select(starts_with("ph"))) %>%
  #   filter(grepl("tx", V1), grepl("ph", V2)) %>%
  #   mutate(V1 = factor(sub("tx_cc.comp", "", V1), levels = c(1:10)),
  #          V2 = factor(sub("ph_cc.comp", "", V2), levels = c(1:10))) %>%
  #   ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   redblu.col.gradient+my.guides+
  #   geom_text(size=2) +
  #   labs(x = "imputed_tx componenets",
  #        y = "phenotypes componenets")
  # p4_13 <- corr.table(t2%>%select(starts_with("tx")), t2%>%select(starts_with("pgs"))) %>%
  #   filter(grepl("tx", V1), grepl("pgs", V2)) %>%
  #   mutate(V1 = factor(sub("tx_cc.comp", "", V1), levels = c(1:10)),
  #          V2 = factor(sub("pgs_cc.comp", "", V2), levels = c(1:10))) %>%
  #   ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   redblu.col.gradient+my.guides+
  #   geom_text(size=2) +
  #   labs(x = "imputed_tx componenets",
  #        y = "PGS componenets")
  # p4_23 <- corr.table(t2%>%select(starts_with("pgs")), t2%>%select(starts_with("ph"))) %>%
  #   filter(grepl("pgs", V1), grepl("ph", V2)) %>%
  #   mutate(V1 = factor(sub("pgs_cc.comp", "", V1), levels = c(1:10)),
  #          V2 = factor(sub("ph_cc.comp", "", V2), levels = c(1:10))) %>%
  #   ggplot(aes(x=V1, y=V2, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   redblu.col.gradient+my.guides+
  #   geom_text(size=2) +
  #   labs(x = "PGS componenets",
  #        y = "phenotypes componenets")
  # p4 <- patchwork::wrap_plots(p4_12, p4_13, p4_23)
  # ####
  # # visualize component loadings
  # l2 <- cbind(cc = c(1:10),
  #             t(tx.phe.pgs.cca$a$phenotypes),
  #             t(tx.phe.pgs.cca$a$imputed_tx),
  #             t(tx.phe.pgs.cca$a$PGS)) %>%
  #   as.data.frame() %>%
  #   mutate(cc = factor(cc, levels = c(1:10))) %>%
  #   pivot_longer(cols = colnames(phenotypes)[-c(1:4)], names_to = "ph_var", values_to = "ph_var_lo") %>%
  #   mutate(grade = sub("_.*", "", ph_var),
  #          phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", ph_var))) %>%
  #   pivot_longer(cols = any_of(all.genes$gene), names_to = "tx_var", values_to = "tx_var_lo") %>%
  #   pivot_longer(cols = colnames(sli.pgs)[-1], names_to = "pgs_var", values_to = "pgs_var_lo")
  # 
  # p5_1 <- l2 %>%
  #   distinct(cc, tx_var, tx_var_lo) %>%
  #   ggplot(aes(x=cc, y=tx_var, fill = tx_var_lo)) +
  #   geom_tile()+
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   null_labs +
  #   theme(axis.text.y.left = element_text(size = 2))
  # p5_2 <- l2 %>%
  #   distinct(cc, grade, ph_var_lo, phenotype) %>%
  #   ggplot(aes(x=cc, y=grade, fill = ph_var_lo)) +
  #   geom_tile()+
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
  #   null_labs +
  #   theme(strip.text.y.right = element_text(angle = 0))
  # p5_3 <- l2 %>%
  #   distinct(cc, pgs_var_lo, pgs_var) %>%
  #   ggplot(aes(x=cc, y=pgs_var, fill = pgs_var_lo)) +
  #   geom_tile()+
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   null_labs
  # p5 <- patchwork::wrap_plots(p5_1, p5_2, p5_3, nrow = 1)
  # ####
  # # visualize sample component weights correlation with the raw block data
  # # components correlation with genes
  # b1 <- corr.table(tx.phe.pgs%>%
  #                    select(any_of(all.genes$gene)),
  #                  tx.phe.pgs.cca$Y$imputed_tx) %>%
  #   filter(V1 %in% all.genes$gene,
  #          V2 %in% colnames(tx.phe.pgs.cca$Y$imputed_tx)) %>%
  #   mutate(V2 = factor(sub("comp", "", V2), levels = as.character(c(1:10))))
  # # components correlation with phenotypes
  # b2 <- corr.table(tx.phe.pgs%>%
  #                    select(colnames(phenotypes)[-c(1:4)]),
  #                  tx.phe.pgs.cca$Y$phenotypes) %>%
  #   filter(V1 %in% colnames(phenotypes)[-c(1:4)],
  #          V2 %in% colnames(tx.phe.pgs.cca$Y$phenotypes)) %>%
  #   mutate(V2 = factor(sub("comp", "", V2), levels = as.character(c(1:10)))) %>%
  #   mutate(grade = sub("_.*", "", V1),
  #          phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V1))) 
  # b3 <- corr.table(tx.phe.pgs%>%
  #                    select(colnames(sli.pgs)[-1]),
  #                  tx.phe.pgs.cca$Y$PGS) %>%
  #   filter(V1 %in% colnames(sli.pgs)[-1],
  #          V2 %in% colnames(tx.phe.pgs.cca$Y$PGS)) %>%
  #   mutate(V2 = factor(sub("comp", "", V2), levels = as.character(c(1:10))))
  # # cluster components correlation with genes
  # b1.d <- as.dendrogram(hclust(d = dist(x = b1 %>% pivot_wider(names_from = "V2", values_from = "r", id_cols = "V1")%>%column_to_rownames("V1"))))
  # b1.d2 <- ggdendrogram(b1.d, rotate = T) +
  #   theme(axis.text.y.left = element_text(size=2))
  # order1 <- order.dendrogram(b1.d)
  # 
  # p6_11 <- b1 %>%
  #   mutate(V1 = factor(V1, levels = unique(V1)[order1])) %>%
  #   ggplot(aes(x=V2, y=V1, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   geom_text(size = 1) +
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   null_labs +labs(title = "correlation between CCA values and samples' raw data") +
  #   theme(axis.text.y.left = element_text(size = 2), 
  #         plot.title = element_text(size = 10))
  # # p6_12 <- print(b1.d2, 
  # #                vp = viewport(height = 1))
  # p6_1 <- cowplot::plot_grid(p6_11,
  #                            b1.d2,
  #                            rel_widths = c(2,0.5))
  # 
  # # grid.newpage()
  # # pp <- print(p6_11, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1))
  # # pp2 <- print(b1.d2, 
  # #       vp = viewport(x = 0.90, y = 0.51, height =1.035, width = 0.2))
  # 
  # p6_2 <- b2 %>%
  #   ggplot(aes(x=V2, y=grade, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   geom_text(size = 3) +
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   facet_grid2(rows = vars(phenotype), scales = "free", space = "free") +
  #   theme(strip.text.y.right = element_text(angle = 0)) +
  #   null_labs
  # 
  # # cluster components correlation with pgs
  # b3.d <- as.dendrogram(hclust(d = dist(x = b3 %>% pivot_wider(names_from = "V2", values_from = "r", id_cols = "V1")%>%column_to_rownames("V1"))))
  # b3.d2 <- ggdendrogram(b3.d, rotate = T) +
  #   theme(axis.text.y.left = element_text(size=2))
  # order3 <- order.dendrogram(b3.d)
  # 
  # 
  # p6_3 <- b3 %>%
  #   mutate(V1 = factor(V1, levels = unique(V1)[order3])) %>%
  #   ggplot(aes(x=V2, y=V1, fill = r, label = ifelse(pval<0.05, "*", ""))) +
  #   geom_tile()+
  #   geom_text(size = 3) +
  #   my.guides+scale_fill_gradient2(low = redblack.col[2], high = redblack.col[1], name = "")+
  #   theme(strip.text.y.right = element_text(angle = 0)) +
  #   null_labs
  # 
  # p6 <- cowplot::plot_grid(p6_1, p6_2, p6_3, nrow = 1)
  # 
  # ####
  # 
  # # save plots in pdf
  # # pdf(file = paste0("figs/CCA/tx-phenotypes-pgs/", tissue, ".pdf"), width = 13.5, height = 12.5)
  # # print(p1);print(p2);print(p3_12);print(p3_13);print(p3_23);print(p4_12);print(p4_13);print(p4_23);print(p5)
  # # dev.off()
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
  # ggsave(p5, filename = paste0("figs/CCA/tx-phenotypes-pgs/p5-", tissue, ".png"), 
  #        bg = "white", width = 11, height = 11, units = "in")
  # ggsave(p6, filename = paste0("figs/CCA/tx-phenotypes-pgs/p6-", tissue, ".png"), 
  #        bg = "white", width = 15, height = 17, units = "in")
  #################################################
  # saving data ~ CCA loadings and sample weights #
  #################################################
  # write_csv(t1, file = paste0("data/derivatives/CCA/tx-phenotypes/sample-cc-", tissue, ".csv"))
  # write_csv(t2, file = paste0("data/derivatives/CCA/tx-phenotypes-pgs/sample-cc-", tissue, ".csv"))
  # write_csv(l1, file = paste0("data/derivatives/CCA/tx-phenotypes/cc-loadings-", tissue, ".csv"))
  # write_csv(l2, file = paste0("data/derivatives/CCA/tx-phenotypes-pgs/cc-loadings-", tissue, ".csv"))
  ####################
  # delete and clean #
  ####################
  rm(list = c("p1","p3","p4","p5", "tx.phe.cca", "tx.phe.pgs.cca","tissue.tx",
              "tx.phe.pgs", "t1", "t2", "l1", "l2", "p3_12", "p3_13", "p3_23", "p4_12", "p4_13", "p4_23",
              "p5_1", "p5_2", "p5_3", "tissue.tx.pc.corr"))
  gc()
}
################################################################################
################################################################################
# plot correlation between the phenotypes raw scores themselves
# assume you read some random tissue tx and got the tx.phe.pgs
corr.table(tx.phe.pgs %>% select(colnames(phenotypes)[-c(1:4)]),
           y= data.frame(idk = matrix(nrow = nrow(tx.phe.pgs), 0, ncol = 1))) %>%
  filter(V1 != "idk", V2 != "idk") %>%
  mutate(grade = sub("_.*", "", V1),
         phenotype = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V1)),
         grade2 = sub("_.*", "", V2),
         phenotype2 = factor(sub("^[^_]+_([^_]+)_(\\w+)$", "\\2", V2))) %>%
  ggplot(aes(x = grade, y=grade2, fill = r, label = ifelse(pval < 0.05, "*", ""))) +
  geom_tile() +
  geom_text(size=2) +
  facet_grid2(rows = vars(phenotype2), cols = vars(phenotype), scales = "free", space = "free") +
  redblack.col.gradient + null_labs + my.guides +
  theme(strip.text.x.top = element_text(angle = 90),
        strip.text.y.right = element_text(angle = 0))
ggsave(filename = paste0("figs/phenotypes-corr.png"), 
              bg = "white", width = 8, height = 8, units = "in")
################################################################################


################################################################################


