# Regression model of frequency-dependent selection  
Source code of Sato et al. "Detecting frequency-dependent selection using a genetic marker regression of fitness components".  

## Directory

### /application
R source code or input files for pilot tests using real data.  

- pheno_branchNum.csv  
CSV file including an accession list and branch number data of <i>Arabidopsis thaliana</i>.  


### /simulation
Source code for SLiM and GWAS simulations.  

SLiM codes (xxxx.slim) for  
- negative frequency-dependent selection (NFDS);  
- positive frequency-dependent selection (PFDS);  
- overdominance (OD);  
- and spatiotemporally varying selection (STVS).  

GWAS simulation using the SLiM outputs.  

/output  
Folder to store SLiM and GWAS output files.  

