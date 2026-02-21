#!/bin/bash

# Bash script to extract the values for the null and alternative models from godon's results
# Creates a csv file for all subfamilies, with as header orthogp / branch / lnull / lalt
# Runs the R script '~/scripts/branchSite.lrt.fdr.r' to do statistics and FDR correction
# FDR correction is set at 0.05, within the R script
# Outputs a tree file with the branches under positive selection marked with a '+' (to view in FigTree -> 'Branch Labels > Display > label')

# Define source directory for all scripts based on the path for the current script
scripts_path="$(dirname "$(readlink -f "$0")")"

clade=$(basename "${PWD}")

for folder in $(ls -d ${PWD}/Subfamily_*); do 
	subfamily=$(basename "${folder}" | cut -c 11-)
	if [ ! -f "${folder}/godon/Godon_OK" ]; then
		continue
	fi
	cd ${folder}/godon/
	# Extracting statistics from the Godon output and building a CSV file
	python ${scripts_path}/godon_correction.py -f extract -godon godonBS.${subfamily}.out -csv ${clade}.godonBS.${subfamily}.out.csv
	
	if [ $? -ne 0 ]; then 
		echo "ERROR getting the statistics for ${subfamily}" >&2 
		exit 1
	fi
	cd ${folder}/../
done

# Merge all subfamilies into one file for FDR correction
echo "Subfamily,branch,lnull,lalt" > ${PWD}/${clade}.godonBS.all.csv
for file in ${PWD}/Subfamily_*/godon/${clade}.godonBS.*.out.csv; do
	tail -n +2 ${file} >> ${PWD}/${clade}.godonBS.all.csv
done

# Running R script to get the p-values and q-values
module load R/4.3.3 gcc/13.2.0
R --slave --args ${PWD}/${clade}.godonBS.all.csv ${PWD}/${clade}.godonBS.all.fdr.csv < ${scripts_path}/branchSite.lrt.fdr.r
if [ ${?} -ne 0 ]; then 
	echo "ERROR: Problem running R" >&2
	exit 1
fi


# Naming every internal branch of the subfamilies gene trees to use in Generax
for folder in $(ls -d ${PWD}/Subfamily_*); do 

	if [ ! -f "${folder}/godon/Godon_OK" ]; then
		continue
	fi
	
	subfamily=$(basename "${folder}" | cut -c 11-)
	
	python ${scripts_path}/godon_correction.py -f tree_label -godon ${folder}/godon/godonBS.${subfamily}.out -it ${folder}/raxml/RAxML_bestTree.${subfamily} -ot ${folder}/${clade}.godonBS.${subfamily}.named.tree
	if [ $? -ne 0 ]; then 
		echo "ERROR making tree for ${subfamily}" >&2 
		exit 1
	fi
done
