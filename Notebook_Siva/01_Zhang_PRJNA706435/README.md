Wild-type O. sativa (ZH11): Is a Japonica rice variety

Data Download:
PRJNA706435 (scRNA)

```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/039/SRR13853439/SRR13853439_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/039/SRR13853439/SRR13853439_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/040/SRR13853440/SRR13853440_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR138/040/SRR13853440/SRR13853440_2.fastq.gz
```

* Abstract:
Root development relies on the establishment of meristematic tissues that give rise to distinct cell types that differentiate across defined temporal and spatial gradients. Dissection of the developmental trajectories and the transcriptional networks that underlie them is therefore crucial for better understanding of the function of the root apical meristem in both dicots and monocots. Here we performed a single-cell RNA (scRNA) sequencing and chromatin accessibility survey of rice radicles. The temporal profiling of individual root tip cells enabled the reconstruction of continuous developmental trajectories of epidermal cells and ground tissues, and elucidation of the regulatory networks underlying cell fate determination of these cell lineages. We further identified characteristic processes, transcriptome profiles, and marker genes for these cell types and used them to reveal conserved and divergent root developmental pathways between dicots and monocots. Finally, we demonstrate the potential of the platform for functional genetic studies by using spatiotemporal modelling to identify a novel rice root meristematic mutant from a cell-specific gene cohort


#### Cellranger:
WD: `/work/gif4/Siva/2025_Reuben_Peters_scRNA_exploration/01_Cellranger/Zhang_2021`

```
#!/bin/bash
#SBATCH -N 1
#SBATCH -n 40
#SBATCH -t 1-00:00:00 
#SBATCH -p nova # partition nova
#SBATCH -J cellranger_count
#SBATCH -o cr_count.o%j # std output 
#SBATCH -e cr_count.e%j # std error


cellranger count \
        --id=SRR13853439 \
        --transcriptome=Nipponbare \
        --fastqs=$PWD \
        --sample=SRR13853439 \
        --expect-cells=5000 \
        --create-bam true \
        --localcores 24


cellranger count \
        --id=SRR13853440 \
        --transcriptome=Nipponbare \
        --fastqs=$PWD \
        --sample=SRR13853440 \
        --expect-cells=5000 \
        --create-bam true \
        --localcores 24

```

#### Job Details and Efficiency:
```bash
Job ID: 6077217
Cluster: nova
User/Group: csiva/domain users
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 40
CPU Utilized: 2-04:51:33
CPU Efficiency: 28.32% of 7-18:40:40 core-walltime
Job Wall-clock time: 04:40:01
Memory Utilized: 148.42 GB
```

I copied the results folder for both SRR13853439 and SRR13853440, excluding the BAM file to my laptop to do futher ananalysis:


WD: `/Users/Lab/ISU_GIF_Job/2025_Reuben_Peters/Zhang_2021`


`more Zhang_scRNAseq.R`

```bash
getwd()
"/Users/Lab/ISU_GIF_Job/2025_Reuben_Peters/Zhang_2021"
install.packages("Seurat")

library(ggplot2)
library(dplyr)
library(Seurat)
library(Matrix)

S439 <- "SRR13853439/raw_feature_bc_matrix/"
S439_obj <- Read10X(data.dir = S439)
S439_obj <- CreateSeuratObject(counts = S439_obj, 
                                 project = "Zhang_SC", 
                                 min.cells = 3, 
                                 min.features = 200)
S439_obj
# An object of class Seurat 
# 29563 features across 81495 samples within 1 assay 
# Active assay: RNA (29563 features, 0 variable features)
# 1 layer present: counts

## QC:
# Add metadata: number of genes and percentage of mitochondrial reads
S439_obj[["percent.mt"]] <- PercentageFeatureSet(S439_obj, pattern = "^MT-")
# Visualize QC metrics
VlnPlot(S439_obj, features = c("nFeature_RNA", 
                               "nCount_RNA", "percent.mt"), ncol = 3)
# Set QC thresholds (adjust based on dataset)
S439_obj <- subset(S439_obj, 
                   subset = nFeature_RNA > 200 
                   & nFeature_RNA < 6000 
                   & percent.mt < 10)
S439_obj
# An object of class Seurat 
# 29563 features across 78747 samples within 1 assay 
# Active assay: RNA (29563 features, 0 variable features)
# 1 layer present: counts
# Normalization and Scaling
S439_obj <- NormalizeData(S439_obj)
S439_obj <- FindVariableFeatures(S439_obj, 
                                 selection.method = "vst", 
                                 nfeatures = 2000)

S439_obj
# An object of class Seurat 
# 29563 features across 78747 samples within 1 assay 
# Active assay: RNA (29563 features, 2000 variable features)
# 2 layers present: counts, data

# Dimensionality Reduction

S439_obj <- ScaleData(S439_obj)
S439_obj <- RunPCA(S439_obj)
ElbowPlot(S439_obj)

S439_obj <- RunUMAP(S439_obj, dims = 1:10)
names(S439_obj@reductions)
DimPlot(S439_obj, reduction = "pca")
DimPlot(S439_obj, reduction = "umap")
 # Clustering and Marker
S439_obj <- FindNeighbors(S439_obj, dims = 1:10)
S439_obj <- FindClusters(S439_obj, resolution = 0.5)

# Identify marker genes for each cluster
markers <- FindAllMarkers(S439_obj, 
                          only.pos = TRUE, 
                          min.pct = 0.25, logfc.threshold = 0.25)
head(markers)

# Compare a cluster vs all others
cluster0_markers <- FindMarkers(S439_obj, 
                               ident.1 = 0, min.pct = 0.25)
head(cluster0_markers)

cluster1_markers <- FindMarkers(S439_obj, 
                                ident.1 = 1, min.pct = 0.25)
head(cluster1_markers)

cluster2_markers <- FindMarkers(S439_obj, 
                                ident.1 = 1, min.pct = 0.25)
head(cluster2_markers)
### Checking a specific gene
table(FetchData(S439_obj, vars = "LOC-Os04g44790") > 0)

VlnPlot(S113_obj, features = "LOC-Os04g44790", group.by = "seurat_clusters")
FeaturePlot(S439_obj, features = "LOC-Os04g44790")
```

