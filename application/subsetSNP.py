####################################
# Subsetting the full imputed SNPs #
####################################

# (1) load libraries
import h5py
import pandas as pd
import numpy as np
import csv

# (2) load data
# data downloaded from AraGWAS catalog (https://aragwas.1001genomes.org/#/download-center)
input_file = "all_chromosomes_binary.hdf5"
h5file = h5py.File(input_file,"r")

snps = h5file["snps"].value
accs = h5file["accessions"].value
pos = h5file["positions"].value

acID = []
for i in accs:
    acID.append(int(i))

acID = np.array(acID)
gwasID = pd.read_csv("gwasIDlist.csv") # load gwasIDs

acc_list = []
for i in gwasID["GWASid"]:
    place = np.where(acID == i) # search a focal accession
    acc_list.append(int(place[0]))

# (2) export the subset data
sub_snps = snps[:,acc_list]
del snps
sub_snps = pd.DataFrame(sub_snps)
sub_snps.to_csv("sub_snps.csv") # export SNPs

pos = pd.DataFrame(pos)
pos.to_csv("positions.csv") # export positions
