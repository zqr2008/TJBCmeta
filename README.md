TJBCmeta
Reproducible analysis code for the Tianjin Birth Cohort (TJBC) maternal–infant metagenomic study. The repository contains scripts used for metadata processing, microbiome diversity analysis, differential-abundance modelling, strain-level transmission inference, CAZy/KEGG functional analysis, benchmarking, and supporting tables/figures.
> **Important:** This repository contains analysis scripts and workflow templates. The individual-level cohort data, large intermediate objects, and raw metagenomic files are not bundled because they may contain sensitive human-subject information. To reproduce the full analyses, users must prepare the required input files locally and update the hard-coded project paths in the scripts.
---
Repository structure
```text
TJBCmeta/
├── README.md
├── trajectoryDefine.R
├── strainphlan/
│   └── StrainPhlAn.sh
├── PanPhlAn/
│   ├── step1.sh
│   ├── step2.sh
│   ├── singularity_envirment_step1.sh
│   └── singularity_envirment_step2.sh
└── publish_Script/
    ├── 00_documentation/
    │   ├── README_scripts__original.md
    │   └── update_publish_metadata.py
    ├── 01_metadata_qc_baseline_Metadata_QC_baseline/
    ├── 02_diversity_permanova_Diversity_PERMANOVA/
    ├── 03_differential_abundance_Differential_abundance_trajectory/
    ├── 04_transmission_and_strain_Transmission_strain_persistence/
    ├── 05_functional_cazy_kegg_Functional_CAZy_KEGG/
    ├── 06_benchmarking_methods_Benchmarking_methods/
    └── 07_supporting_tables_data_Supporting_tables_data/
```
The `publish_Script/` directory is the main publication-oriented code collection. Script filenames usually indicate the figure or supplementary table they generate, for example:
```text
new_maslin2__no_figure.R
draw_separeatevisit__Figure3C_Figure4B.R
instrain_transmission__Figure4A_Supplementary Table 11.R
cazy__Figure5A_Figure5B.R
somaticgrowth_with_GT_GH__Figure6A.R
```
Most published scripts begin with a `PUBLISH_ANNOTATION_START` block describing:
category,
figure/table linkage,
purpose,
main input files,
expected outputs,
static warnings.
---
System requirements
The original analysis environment was:
```text
CentOS 7.9.2009
R 4.2.2
```
Recommended practical environment:
Linux workstation or HPC node.
Conda/Mamba for command-line bioinformatics tools.
R 4.2.x for compatibility with the original scripts.
Singularity/Apptainer if running the provided HPC container wrapper scripts.
Sufficient storage for metagenomic data. Full cohort-level analyses require large storage and memory; figure-level R scripts usually require the processed `.rds`, `.rda`, `.tsv`, or `.xlsx` objects used by the manuscript.
---
Core external software
The current README lists the following command-line tools and R packages. The table below expands them with purpose and official links.
Software	Version used	Purpose in this project	Official link
fastp	0.23	FASTQ quality control and trimming	https://github.com/OpenGene/fastp
Bowtie2	2.4.4	Read alignment	https://github.com/BenLangmead/bowtie2
MetaPhlAn	4.0.6	Taxonomic profiling; marker database `mpa_vOct22_CHOCOPhlAnSGB_202212`	https://github.com/biobakery/MetaPhlAn
StrainPhlAn	MetaPhlAn 4 suite	Strain-level marker phylogeny and strain-sharing inference	https://github.com/biobakery/MetaPhlAn/wiki/Strain-Sharing-Inference
PanPhlAn	not specified	Pangenome profiling for selected species	https://github.com/biobakery/PanPhlAn
MEGAHIT	1.2.9	Metagenomic assembly	https://github.com/voutcn/megahit
MaxBin2	2.2.7	Genome binning	https://sourceforge.net/projects/maxbin2/
MetaBAT2	2.18	Genome binning	https://bitbucket.org/berkeleylab/metabat/src/master/
CONCOCT	1.1.0	Genome binning	https://github.com/BinPro/CONCOCT
DAS Tool	1.1.7	Consensus bin selection/refinement	https://github.com/cmks/DAS_Tool
dRep	3.6.2	Genome dereplication	https://github.com/MrOlm/drep
CheckM	1.2.4	MAG quality assessment	https://github.com/Ecogenomics/CheckM
QUAST	5.3.0	Assembly/bin quality metrics	https://github.com/ablab/quast
GTDB-Tk	2.5.2	Taxonomic classification of genomes	https://github.com/Ecogenomics/GTDBTk
dbCAN3	4.1.4	CAZyme annotation	https://github.com/linnabrown/dbcan
inStrain	1.10.0	Strain-level population profiling	https://github.com/MrOlm/inStrain
MaAsLin2	1.12.0	Multivariable differential-abundance modelling	https://github.com/biobakery/Maaslin2
vegan	2.7.1	Diversity analysis and PERMANOVA	https://cran.r-project.org/package=vegan
ape	5.8-1	Phylogenetic analysis	https://cran.r-project.org/package=ape
nnet	7.3-18	Multinomial/statistical modelling	https://cran.r-project.org/package=nnet
R stats	4.2.2	Base statistical tests/models	https://www.r-project.org/
---
R package dependencies
The R scripts use a broad set of packages. The main packages observed across the repository are:
```r
c(
  "tidyverse", "dplyr", "tidyr", "tibble", "readr", "readxl", "haven",
  "data.table", "ggplot2", "ggpubr", "ggrepel", "ggh4x", "pheatmap",
  "RColorBrewer", "viridis", "scales", "patchwork", "ComplexHeatmap",
  "circlize", "MetBrewer", "wesanderson", "hrbrthemes",
  "vegan", "phyloseq", "microbiome", "microViz", "ape", "picante",
  "zCompositions", "compositions", "cluster", "umap", "uwot",
  "Maaslin2", "rstatix", "lme4", "lmerTest", "nlme", "lcmm",
  "broom", "gtsummary", "tableone", "labelled", "sjmisc", "memisc",
  "randomForest", "ranger", "pROC", "caret", "smotefamily",
  "openxlsx", "writexl", "officer", "flextable", "Hmisc",
  "tidytext", "igraph", "nnet"
)
```
Not every package is required for every analysis. For example, MaAsLin2 differential-abundance scripts mainly need `Maaslin2`, `microViz`, `phyloseq`, `microbiome`, and `tidyverse`, whereas figure/table scripts may additionally need plotting and reporting packages.
---
Installation
1. Clone the repository
```bash
git clone https://github.com/zqr2008/TJBCmeta.git
cd TJBCmeta
```
2. Create a command-line metagenomics environment
Using `mamba` is recommended because the environment is large.
```bash
mamba create -n tjbcmeta-tools -c conda-forge -c bioconda \
  fastp=0.23.* \
  bowtie2=2.4.4 \
  metaphlan=4.0.6 \
  megahit=1.2.9 \
  maxbin2=2.2.7 \
  metabat2=2.18 \
  concoct=1.1.0 \
  dastool=1.1.7 \
  drep=3.6.2 \
  checkm-genome=1.2.4 \
  quast=5.3.0 \
  gtdbtk=2.5.2 \
  instrain=1.10.0 \
  samtools \
  python
```
Activate it with:
```bash
mamba activate tjbcmeta-tools
```
Some exact historical versions may be unavailable on current channels. If exact version resolution fails, create separate environments for conflicting tools, or use container images to freeze the original software stack.
3. Create an R environment
```bash
mamba create -n tjbcmeta-r -c conda-forge -c bioconda r-base=4.2.2 r-essentials
mamba activate tjbcmeta-r
```
Then install R packages:
```r
install.packages(c(
  "tidyverse", "data.table", "ggpubr", "ggrepel", "ggh4x", "pheatmap",
  "RColorBrewer", "viridis", "scales", "patchwork", "MetBrewer",
  "wesanderson", "hrbrthemes", "vegan", "ape", "picante",
  "zCompositions", "compositions", "cluster", "umap", "uwot",
  "rstatix", "lme4", "lmerTest", "nlme", "lcmm", "broom",
  "gtsummary", "tableone", "labelled", "sjmisc", "memisc",
  "randomForest", "ranger", "pROC", "caret", "smotefamily",
  "openxlsx", "writexl", "officer", "flextable", "Hmisc",
  "tidytext", "igraph", "nnet", "readxl", "haven"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "Maaslin2", "phyloseq", "microbiome", "ComplexHeatmap", "circlize"
))
```
4. Verify installation
```bash
Rscript -e 'pkgs <- c("tidyverse","vegan","phyloseq","microViz","Maaslin2","ggpubr","rstatix","lme4","lcmm"); sapply(pkgs, requireNamespace, quietly=TRUE)'
```
Expected output: a named logical vector with `TRUE` for installed packages.
For command-line tools:
```bash
fastp --version
bowtie2 --version
metaphlan --version
megahit --version
checkm --version
gtdbtk --version
inStrain --version
```
---
Input data organization
The original scripts use absolute local paths such as:
```text
C:/Users/zqr20/Documents/tjmeta/BIG_revision/
C:/Users/zqr20/OneDrive - BGI Hong Kong Tech Co., Limited/文档/tjmeta/BIG_revision/
```
Before running the scripts on another machine, edit these paths to your local project directory. A recommended local layout is:
```text
project/
├── 00_raw_data/
│   ├── fastq/
│   └── metadata/
├── 01_reference/
│   ├── metaphlan/
│   ├── gtdbtk/
│   └── dbcan/
├── 02_processed_objects/
│   ├── phyloseq_objects/
│   ├── transmission/
│   ├── cazy/
│   └── kegg/
├── 03_intermediate/
├── 04_analysis_results/
│   ├── differential_abundance_results/
│   ├── diversity_results/
│   ├── transmission_results/
│   ├── functional_results/
│   └── figures/
└── 05_tables/
```


How to run
A. Growth-trajectory modelling
The top-level script `trajectoryDefine.R` defines child growth trajectories using group-based trajectory modelling with `lcmm::hlme`.
```bash
Rscript trajectoryDefine.R
```
Main inputs:
SPSS `.sav` file containing child z-score data.
Required variables include `IDchild`, `agemos`, and `zwhz`.
Expected outputs:
fitted latent-class mixed models for quadratic and cubic models,
model summary table including class proportions,
average posterior probability summary,
selected child trajectory class labels.
Before running, update the input path in the script.
B. MaAsLin2 differential-abundance analysis
Example script:
```bash
Rscript publish_Script/03_differential_abundance_Differential_abundance_trajectory/new_maslin2__no_figure.R
```
Main inputs:
`02_processed_objects/phyloseq_objects/newdata_phylo_NMrevision_0305.rds`
Main expected outputs:
MaAsLin2 output directories under `04_analysis_results/differential_abundance_results/`,
`all_results.tsv`,
`significant_results.tsv`,
model diagnostics and residual files generated by MaAsLin2,
downstream heatmap-ready significant taxa tables.
C. Diversity and PERMANOVA analyses
Example scripts:
```bash
Rscript publish_Script/02_diversity_permanova_Diversity_PERMANOVA/revision_everydistance__FigureS3_Supplementary\ Table\ 4_Supplementary\ Table\ 5.R

Rscript publish_Script/02_diversity_permanova_Diversity_PERMANOVA/permonva_revise_05_mother_permanova_bubbles__Figure2C_FigureS4A_Supplementary\ Table\ 6.R
```
Expected outputs:
alpha/beta diversity summaries,
PERMANOVA tables,
bubble plots or figure-ready data for Figure 2 and Supplementary Figures/Tables.
D. StrainPhlAn strain-sharing workflow
The script `strainphlan/StrainPhlAn.sh` generates per-species shell scripts for StrainPhlAn analysis.
Prepare:
```bash
cd strainphlan
mkdir -p shell db_marker Accurateoutput Fastoutput Sample database
```
Required files:
```text
samples.marker.list
database/mpa_vOct22_CHOCOPhlAnSGB_202212.pkl
Sample/*.pkl
1w8.metadata.forR
1w8.metadata.forTrans
threshold_calc.R
strain_transmission.py
```
Run script generation:
```bash
bash StrainPhlAn.sh
```
Expected outputs:
```text
shell/<species>.accurate.sh
shell/<species>.fast.sh
```
Then run the generated scripts, for example:
```bash
bash shell/t__Bacteroides_thetaiotaomicron.accurate.sh
bash shell/t__Bacteroides_thetaiotaomicron.fast.sh
```
Expected StrainPhlAn outputs:
```text
db_marker/<species>/<species>.fna
Accurateoutput/<species>/RAxML_bestTree.<species>.StrainPhlAn4.tre
Accurateoutput/<species>/<species>_nGD.tsv
Accurateoutput/<species>/threshold.txt
Fastoutput/<species>/RAxML_bestTree.<species>.StrainPhlAn4.tre
Fastoutput/<species>/<species>_nGD.tsv
Fastoutput/<species>/threshold.txt
```
The generated workflow extracts marker genes, builds strain-level phylogenies, calculates pairwise normalized genetic distances, estimates species-specific thresholds, and detects strain-sharing events.
E. PanPhlAn pangenome workflow
The `PanPhlAn/` scripts are HPC-oriented wrappers originally configured for `Bacteroides_thetaiotaomicron`.
Main scripts:
```text
PanPhlAn/step1.sh
PanPhlAn/step2.sh
PanPhlAn/singularity_envirment_step1.sh
PanPhlAn/singularity_envirment_step2.sh
```
Before running, edit:
working directory,
sample list file,
FASTQ location,
PanPhlAn database/index path,
pangenome TSV path,
temporary directory,
output directory,
Singularity image path if using the container wrappers.
Expected outputs:
```text
shell/<species>_<sample>_step1.sh
shell/<species>_<sample>_sing.sh
Result/<species>_<sample>.tsv
```
F. Functional CAZy/KEGG analyses
Example scripts:
```bash
Rscript publish_Script/05_functional_cazy_kegg_Functional_CAZy_KEGG/cazy__Figure5A_Figure5B.R
Rscript publish_Script/05_functional_cazy_kegg_Functional_CAZy_KEGG/KEGG_organize__Figure5C.R
Rscript publish_Script/05_functional_cazy_kegg_Functional_CAZy_KEGG/Bifidobacterium_pseudocatenulatum_PCA__Figure5D.R
Rscript publish_Script/05_functional_cazy_kegg_Functional_CAZy_KEGG/somaticgrowth_with_GT_GH__Figure6A.R
```
Expected outputs:
CAZy heatmaps,
KEGG summary plots,
PCA plots for CAZyme profiles,
GH/GT associations with abundance or infant growth,
Supplementary Tables 14–17 depending on the script.
G. Benchmarking and supporting tables
Benchmarking scripts are in:
```text
publish_Script/06_benchmarking_methods_Benchmarking_methods/
```
Supporting-table scripts are in:
```text
publish_Script/07_supporting_tables_data_Supporting_tables_data/
```
Expected outputs:
benchmark summary tables,
prevalence calculations,
supplementary table-ready data.
---
