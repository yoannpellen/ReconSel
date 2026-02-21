from scipy.stats import hypergeom
import argparse
import os
import sys
from scipy import stats
import math

# Description for script whith '-h'
parser = argparse.ArgumentParser(description = "Script to run the hypergeometric test after dN/dS pipeline")
optional_args = parser.add_argument_group(title = "Arguments:")

# Define arguments with help message
optional_args.add_argument('-map', dest = 'map', help = "Mapping file with statistics.")
optional_args.add_argument('-clade', dest = 'clade', default = "all", help = "Clade to run (see first column). Default 'all'. To run only the 'colony size' for example, put '_cs'")
optional_args.add_argument('-sf', dest = 'subfamily', default = "all", help = "Subfamily to run. Default 'all'. To run only the '9-Exons', put '9E'")
optional_args.add_argument('-branch', dest = 'branches', default = "all", help = "Species branch to focus on. Default 'all'.")
optional_args.add_argument('-pval', dest = 'pval', default = 0.05, help = "pvalue from Godon to consider positive selection. Default '0.05'.")

# Group all arguments in a list
args = parser.parse_args()

map_file = args.map
clade = args.clade
subfamily = args.subfamily
branches = args.branches
pvalue = float(args.pval)

with open("hypergeometrict_test.tsv", 'w') as output:
    output.write("Clade\tSpecies branch\tNormalized sum of the length of all gene tree branches\tTotal number of positive gene tree branches\tNormalized sum of the length of the gene tree branches mapped on the species branch of interest\tNumber of positive gene tree branches mapped on the species branch of interest\tpval\tqval\tlog2 ratio\tfold change\n")

pval_list = []
gene_list = []
sum_clade_gene_length = 0
count_gene_branches = 0
species_list = []
species_length = 0
count_species = 0
output_lines = []

with open(map_file, 'r') as statistics:
    next(statistics) # skip header
    for line in statistics:
        if f"{line.split(',')[1]}-{line.split(',')[2]}" not in gene_list:
            gene_list.append(f"{line.split(',')[1]}-{line.split(',')[2]}")
            sum_clade_gene_length += float(line.split(',')[3])
            count_gene_branches += 1
        if line.split(',')[4] not in species_list:
            species_list.append(f"{line.split(',')[4]}")
            species_length += float(line.split(',')[5])
            count_species += 1

if branches == "all":
    branches = ",".join(species_list)

average_sp_length = species_length/count_species
average_sum_clade_gene_length = sum_clade_gene_length/count_gene_branches

for branch in branches.split(','):
    M=sum_clade_gene_length # population size (sum of the length of all gene branches in the clade)
    n=0 # number of successes in the population (branches under positive selection in all gene branches)
    N=0 # sample size (sum of the length of all gene branches mapped on the species branch of interest)
    X=0 # number of drawn successes (branches under positive selection on the species branch of interest)
    
    with open(map_file, 'r') as statistics:
        next(statistics)
        unique_gene_branches = []
        for line in statistics:
            if (clade in line.split(',')[0] or clade == "all") and (subfamily == line.split(',')[1] or subfamily == "all"): # Clade and subfamily of interest
                # To use either the p-value or q-value as threshold to consider positive selection by Godon
                # pval => float(line.split(',')[9])
                # qval => float(line.split(',')[10][:-1])
                if pvalue >= float(line.split(',')[9]) and f"{line.split(',')[1]}-{line.split(',')[2]}" not in unique_gene_branches: # Any unique positive gene branch
                    n += 1
                    unique_gene_branches.append(f"{line.split(',')[1]}-{line.split(',')[2]}")
                if branch == line.split(',')[4]: # Species branch of interest
                    N += float(line.split(',')[3])
                    if pvalue >= float(line.split(',')[9]): # Positive gene branch on species branch of interest
                        X += 1
    
    N = round(N/(average_sum_clade_gene_length))
    M = round(M/average_sum_clade_gene_length)
    
    hg_pval = hypergeom.pmf(X, M, n, N)
    pval_list.append(hg_pval)
    
    output_lines.append(clade + "\t" + branch + "\t" + str(M) + "\t" + str(n) + "\t" + str(N) + "\t" + str(X) + "\t" + str(hg_pval))
    
qvalues = stats.false_discovery_control(pval_list)

for line, qval in zip(output_lines, qvalues):
    log2ratio = math.log2((int(line.split('\t')[5])/int(line.split('\t')[4]))/(int(line.split('\t')[3])/int(line.split('\t')[2])))
    fold_change = 2**float(log2ratio)
    
    with open("hypergeometrict_test.tsv", 'a') as output:
        output.write(f"{line}\t{qval}\t{log2ratio}\t{fold_change}\n")
    
