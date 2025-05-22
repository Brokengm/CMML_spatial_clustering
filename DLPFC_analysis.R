library(BayesSpace)
library(ggplot2)
library(patchwork)
library(mclust)

dlpfc <- readVisium('C:/Users/dyq37/Desktop/CMML_2/DLPFC')
# Preprocess the data
set.seed(101)
dlpfc <- scuttle::logNormCounts(dlpfc)
dec <- scran::modelGeneVar(dlpfc)
top <- scran::getTopHVGs(dec, n = 2000)

set.seed(102)
dlpfc <- scater::runPCA(dlpfc, subset_row=top)

## Add BayesSpace metadata
dlpfc <- spatialPreprocess(dlpfc, platform="Visium", skip.PCA=TRUE)

q <- 7  # Number of clusters
d <- 15  # Number of PCs

## Run BayesSpace clustering
set.seed(104)
dlpfc <- spatialCluster(dlpfc, q=q, d=d, platform='Visium',
                        nrep=50000, gamma=3, save.chain=FALSE)

## We recoded the cluster labels to match the expected brain layers (different from the tutorial)
labels <- dplyr::recode(dlpfc$spatial.cluster, `1` = 4, `2` = 1, `3` = 7, `4` = 5, `5` = 6, `6` = 2, `7` = 3)

## View results
clusterPlot(dlpfc, label=labels, palette=NULL, size=0.05) +
  scale_fill_viridis_d(option = "D", labels = 1:7) +
  labs(title="BayesSpace")

#collect different clustering results
clusters <- list()
clusters[["BayesSpace (t-distributed error)"]] <- dlpfc$spatial.cluster

## Louvain
set.seed(100)
g.jaccard = scran::buildSNNGraph(dlpfc, use.dimred="PCA", type="jaccard")
clusters[["Louvain"]] <- igraph::cluster_louvain(g.jaccard)$membership
# get name of different methods
methods <- names(clusters)

# to store different plots
plot_list <- vector("list", length(methods))
names(plot_list) <- methods

for (meth in methods) {
  if (meth == "Louvain") {
    colData(dlpfc)$method_label <- clusters[[meth]]
  # recode to match the expected brain layers,different from Bayesspace's
    colData(dlpfc)$method_label <- dplyr::recode(
    colData(dlpfc)$method_label, `1` = 1, `2` = 4, `3` = 2, `4` = 5, `5` = 3, `6` = 7, `7` = 6,`8` = 8)
    p <- clusterPlot(
         dlpfc,label = "method_label",npalette = NULL,nsize = 0.05) +
       scale_fill_viridis_d(option="D", labels=1:8) +
       labs(title = meth)

    plot_list[[meth]] <- p
  } else {
    colData(dlpfc)$method_label <- clusters[[meth]]
    colData(dlpfc)$method_label <- dplyr::recode(
      colData(dlpfc)$method_label, `1` = 4, `2` = 1, `3` = 7, `4` = 5, `5` = 6, `6` = 2, `7` = 3)
    p <- clusterPlot(
          dlpfc,label = "method_label",npalette = NULL,nsize = 0.05) +
        scale_fill_viridis_d(option="D", labels=1:7) +
        labs(title = meth)

    plot_list[[meth]] <- p
  }
}


png("Louvain.png", width = 800, height = 800)
plot_list[["Louvain"]]
dev.off()

png("bayesspace.png", width = 800, height = 800)
plot_list[["BayesSpace (t-distributed error)"]]
dev.off()

## get the ground truth (mannual annotation)

#Loading objects:sce
load('Human_DLPFC_Visium_processedData_sce_scran_spatialLIBD.Rdata')
sce <- sce[, sce$sample_name == "151673"]
colData(sce)$manual_ann <- colData(sce)$layer_guess_reordered

# rename coordinate，BayesSpace::clusterPlot()need these
colData(sce)$array_row           <- colData(sce)$row
colData(sce)$array_col           <- colData(sce)$col
colData(sce)$pxl_row_in_fullres  <- colData(sce)$imagerow
colData(sce)$pxl_col_in_fullres  <- colData(sce)$imagecol

#plot the ground truth
p1 <- clusterPlot(
  sce,
  label   = "manual_ann",
  palette = NULL,
  size    = 0.05
) +
  scale_fill_viridis_d(option = "D") +
  labs(
    title = "manual_annotation (ground truth)",
    fill  = "Layer"
  )

png("manual_annotation.png", width = 800, height = 800)
p1
dev.off()


#benchmark and visualization
# here I found that the order of the spots in the dlpfc object is not the same as in the sce object
# reorder the dlpfc object to match the order of the spots in the sce object
truth <- colData(sce)$manual_ann
names(truth) <- colnames(sce)

pred_bs   <- clusters[["BayesSpace (t-distributed error)"]]
names(pred_bs) <- colnames(dlpfc)

pred_lv   <- clusters[["Louvain"]]
names(pred_lv) <- colnames(dlpfc)


# find the common spots name between the two objects
common <- intersect(names(truth), names(pred_bs))
#reorder
truth_common <- truth[common]
pred_bs_common <- pred_bs[common]
pred_lv_common <- pred_lv[common]

levels_truth <- c(paste0("Layer", 1:6), "WM")
# convert to sce's label (layer1-6, WM) into numeric
truth_num <- as.integer(factor(truth_common, levels = levels_truth))

# read in results from stSME
stSME <- read.csv("louvain_labels.csv", row.names = 1, stringsAsFactors = FALSE)
pred_stsme <- stSME$louvain
names(pred_stsme) <- rownames(stSME)

# become numeric
pred_stsme_num <- as.integer(pred_stsme)

preds <- list(
  BayesSpace = pred_bs_common,
  Louvain    = pred_lv_common,
  stSME      = pred_stsme_num
)


# drop any NA’s
ok <- !is.na(truth_num)
truth_num <- truth_num[ok]
preds <- lapply(preds, `[`, ok)

ari <- sapply(preds, function(x) adjustedRandIndex(x, truth_num))

# Put into a data.frame and print
ari_df <- data.frame(
  method = names(ari),
  ARI    = unname(ari)
)
print(ari_df)

# visualize as a barplot
png("ARI_plot.png", width = 500, height = 500)
ggplot(ari_df, aes(x = reorder(method, ARI), y = ARI, fill = method)) +
  geom_col(show.legend = FALSE,width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(x = NULL, y = "Adjusted Rand Index (ARI)",
       title = "ARI vs. Manual Layer Annotation") +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text        = element_text(size = 16),
    axis.title.y     = element_text(size = 16)
  )
dev.off()



#Moran's I
library(spdep)
# extract the coordinates of the spots
coords <- cbind(
  colData(sce)$row,
  colData(sce)$col
)

coords <- coords[ok, ]


labels_list <- preds
# construct a neighbor matrix
# using dnearneigh to find neighbors of spots
nb <- dnearneigh(
  x         = coords,
  d1        = 0,
  d2        = 1.5,
  longlat   = FALSE
)
# convert to standard listw
lw <- nb2listw(nb, style="W", zero.policy=TRUE)

# compute Moran’s I
moran_res <- lapply(labels_list, function(lab) {
  lab_num <- as.integer(factor(lab))
  # compute general Moran’s I
  test <- moran.test(lab_num, lw, zero.policy=TRUE)
  c(
    `Moran’s I`    = unname(test$estimate["Moran I statistic"]),
    `Expectation`  = unname(test$estimate["Expectation"]),
    `Variance`     = unname(test$estimate["Variance"])
  )
})

# combine results into a data frame
moran_df <- do.call(rbind, moran_res)
print(moran_df)

# 1. convert into data.frame
moran_df <- as.data.frame(moran_df)           
moran_df$method <- rownames(moran_df)         
rownames(moran_df) <- NULL                   

plot_df <- moran_df[, c("method", "Morans I")]


str(plot_df)

# draw barplot
png("Moran's I_plot.png", width = 500, height = 500)
ggplot(plot_df, aes(x = reorder(method, `Morans I`), y = `Morans I`, fill = method)) +
  geom_col(show.legend = FALSE, width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    x     = NULL,
    y     = "Moran's I",
    title = "Spatial Autocorrelation of Cluster(Moran's I)"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text        = element_text(size = 16),
    axis.title.y     = element_text(size = 16)
  )
dev.off()


#NMI
library(aricode)
nmi_vals <- sapply(preds, function(pred) {
  NMI(truth_num, pred)
})

nmi_df <- data.frame(
  method = names(nmi_vals),
  NMI    = unname(nmi_vals)
)
print(nmi_df)

png("NMI_plot.png", width = 500, height = 500)
ggplot(nmi_df, aes(y = NMI,x = reorder(method, NMI), fill = method)) +
  geom_col(show.legend = FALSE, width = 0.5) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  labs(
    x     = NULL,
    y     = "Normalized Mutual Information (NMI)",
    title = "NMI vs. Manual Annotation"
  ) +
  theme_minimal() +
  theme(
    plot.title       = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text        = element_text(size = 16),
    axis.title.y     = element_text(size = 16)
  )

dev.off()
