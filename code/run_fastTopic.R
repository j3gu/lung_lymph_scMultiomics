#!/software/R-4.2.0-el7-x86_64/bin/R
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
set.seed(10)
workdir<-"~/projects/lung_lymph_scMultiomics/output/fastTopics"

# load data
counts <-
  readRDS(paste0(workdir, "/aggregated_lymph_scRNA_counts.RDS"))

# fit poission nmf with pre-defined k topics
k_topics <-12
main_numiter <-150

fit <-
  fit_poisson_nmf(
    counts,
    k = k_topics,
    init.method = "random",
    method = "em",
    numiter = 20,
    control = list(nc = 8)
  )

fit <-
  fit_poisson_nmf(
    counts,
    fit0 = fit,
    method = "scd",
    numiter = main_numiter,
    control = list(
      numiter = 4,
      nc = 8,
      extrapolate = TRUE
    )
  )

saveRDS(fit,
        sprintf("%s/fit_model_updates%d_k%d.RDS", workdir, k_topics, main_numiter))
