# Regression model of frequency-dependent selection  

## Overview  
Source code and original data generated by Sato et al. "Detecting frequency-dependent selection using a marker-based regression of fitness components".  

## Directory

### /simulation
A folder including SLiM and GWAS simulations.  

SLiM codes (xxxx.slim) for  
- negative frequency-dependent selection (NFDS);  
- positive frequency-dependent selection (PFDS);  
- overdominance (OD);  
- and Spatio-temporally varying selection (STVS).  

GWAS simulation using the SLiM outputs.  

/output  
Folder to store SLiM and GWAS output files.  

### /application
A folder including R source codes and input files for pilot tests using <i>Arabidopsis</i>.  

- pheno_branchNum.csv  
CSV file including an accession list and branch number data of <i>Arabidopsis thaliana</i>.   
