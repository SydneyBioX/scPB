library(ELVAR)
library(Seurat)
library(igraph)
data_dir = "/albona/nobackup2/biostat/datasets/singlecell/scPB/simulation/data/sce_sim_test.rds"

sce_sim <- readRDS(data_dir) 
counts(sce_sim) <- logcounts(sce_sim)
sc_dataset <- as.Seurat(sce_sim)
names(sc_dataset@assays) <- "RNA"
DefaultAssay(sc_dataset) <- "RNA"


sc_dataset <- FindVariableFeatures(sc_dataset, selection.method = "vst", verbose=FALSE);
sc_dataset <- ScaleData(sc_dataset, verbose=FALSE)
sc_dataset <- RunPCA(sc_dataset, verbose=FALSE)
sc_dataset <- FindNeighbors(sc_dataset, dims = 1:8, k.param=20, verbose = FALSE);

adj_m <- as.matrix(sc_dataset@graphs$RNA_nn)
diag(adj_m) <- 0
gr <- graph.adjacency(adj_m, mode="undirected")
vertexN.v <- names(V(gr))
vertex_attr(gr, name = "Condition") <- sc_dataset@meta.data$Condition
comp <- components(gr)
print(is.connected(gr))

eva <- Eva_partitions(gr, resolution = 1, threshold = 0.0000001, alpha = 0.8, Vattr.name= "Condition")
table(eva[[1]], sc_dataset$scClassify_tumour_prediction_coarse)
table(eva[[1]], sc_dataset$Condition)
eva[[1]]
sce_sim$eva <- eva[[1]]
scater::plotUMAP(sce_sim, colour_by = "eva")
sigcl <- ProcessEVA(eva, sc_dataset, attrName = "Condition");
sigcl$pvalEnr
da <- DoDA(sc_dataset, sigcl, varDA = "Condition", varREP = "Sample")

print(da$stat)
print(da$df) # this is the matrix used to fit the negative binomial
