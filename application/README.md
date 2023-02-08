# File description

### GLMM analysis of *Arabidopsis halleri*

-   Ahal_FieldSurvey.R\
    Poisson GLMM of the flower number in wild *Arabidopsis halleri* in a natural population. Data are available via Dryad repository (<https://doi.org/10.5061/dryad.53k2d>).

### GLMM analysis of *Ischnura elegans*

-   damselfly_cage_data.R\
    Poisson GLMM of the egg number in *Ischnura elegans* in semi-field cages. Data are available upon request for Y. Takahashi.

### Branch number GWAS in *Arabidopsis thaliana*

-   pheno_branchNum.csv\
    Accession list and phenotype data for "BranchNo2019GWAS.R".CSV subset of Table S3.

-   reshapeSNP.R\
    R script to prepare genotype data for for "BranchNo2019GWAS.R".

-   subsetSNP.py\
    Python script to extract genotype data from .hdf available at AraGWAS Catalog (<https://aragwas.1001genomes.org>).

-   BranchNo2019GWAS.R\
    Field GWAS of the branch number using 199 accessions of *A. thaliana*.

-   BranchNo2019GWASfigure.R\
    Manhattan and QQ plots for the output from "BranchNo2019GWAS.R".

**Note that the others are input or intermediate files given as symbolic links**
