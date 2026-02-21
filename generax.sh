#!/bin/bash

prep() {
# Define source directory for all scripts based on the path for the current script
scripts_path="$(dirname "$(readlink -f "$0")")"
clade="$(basename "${PWD}")"
rm -rf generax_${clade}
mkdir generax_${clade}
cp ${species_tree} generax_${clade}
echo "[FAMILIES]" > generax_${clade}/families_${clade}.txt

for tree in $(ls ${PWD}/Subfamily_*/*.named.tree); do
	subfamily="$(cut -d'.' -f3 <<<"$(basename "${tree}")")"
	
	cp ${tree} generax_${clade}/generax.${subfamily}.tree
	cp Subfamily_${subfamily}/guidance2/MSA.PRANK.${subfamily}.nucl.NNN0.*.trimmed.aln generax_${clade}/generax.${subfamily}.fasta
	
	python ${scripts_path}/generax.py -f mapping -it generax_${clade}/generax.${subfamily}.tree -map generax_${clade}/mapping.${subfamily}.link
	
	# Create family file
	echo "- family_${subfamily}" >> generax_${clade}/families_${clade}.txt
	echo "starting_gene_tree = generax.${subfamily}.tree" >> generax_${clade}/families_${clade}.txt
	echo "alignment = generax.${subfamily}.fasta" >> generax_${clade}/families_${clade}.txt
	echo "mapping = mapping.${subfamily}.link" >> generax_${clade}/families_${clade}.txt
	echo "subst_model = LG+I" >> generax_${clade}/families_${clade}.txt
done
}

run() {
sbatch -n32 -N1 -J ${clade}.Generax -p preempt1d --mem 0 --wrap \
"cd generax_${clade}
rm -rf Generax *xml MSA.PRANK* slurm*
generax -f families_${clade}.txt -s ${species_tree} --prune-species-tree -r UndatedDL --unrooted-gene-tree --reconcile
cp GeneRax/reconciliations/*.xml ./"
}

while getopts f:t option
do 
    case "${option}"
        in
        f)function=${OPTARG};;
        t)species_tree=${OPTARG};;
    esac
done

case ${function} in

	prep)
		if [ -z ${species_tree} ]; then
			echo "Missing the species tree to be able to run the script, give its path with '-t'."
			exit 1
		elif [ ! -f ${species_tree} ]; then
			echo "File ${species_tree} does not exist."
			exit 1
		fi
		prep;;
	
	run)
		if [ -z ${species_tree} ]; then
			echo "Missing the species tree to be able to run the script, give its path with '-t'."
			exit 1
		elif [ ! -f ${species_tree} ]; then
			echo "File ${species_tree} does not exist."
			exit 1
		fi
		run;;
	
	*)
		echo "${function} does not exist."
		exit 1;;
esac