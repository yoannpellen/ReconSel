#!/bin/bash

# Runs Godon for all orthogroups for which RAxML finished successfully
# See https://github.com/idavydov/godon-tutorial for more details on parameters
# Uses 32 CPUs
# Uses M0 model to estimate branch length
# 4 categories for codon rate variation
# Test all branches

clade=$(basename "${PWD}")
for folder in $(ls -d ${PWD}/Subfamily_*); do
	subfamily=$(basename "${folder}" | cut -c 11-)
	jobname="${clade}_${subfamily}.Godon"
	
	# Skip subfamilies for which this job is already running
	if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "${jobname}"; then 
		continue
	fi
	
	if [ -f "${folder}/raxml/RAxML_bestTree.${subfamily}" ] && [ ! -f "${folder}/godon/Godon_OK" ]; then
		cd ${folder}
		rm -rf godon && mkdir godon
		sbatch -n32 -N1 -p preempt1d -J ${jobname} --mem=0 --wrap \
		"~/softwares/godon test BS \
		--procs 32 \
		--out ${folder}/godon/godonBS.${subfamily}.out \
		--m0-tree \
		--ncat-codon-rate 4 \
		--all-branches \
		${folder}/guidance2/MSA.PRANK.${subfamily}.nucl.NNN0.*.trimmed.aln \
		${folder}/raxml/RAxML_bestTree.${subfamily}
		
		if grep -q 'Running time: ' ${folder}/godon/godonBS.${subfamily}.out; then
			touch ${folder}/godon/Godon_OK
		fi"
		cd ..
	fi
done
