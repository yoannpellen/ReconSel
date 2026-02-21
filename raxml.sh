#!/bin/bash

# Run RAxML for all subfamilies
# GTRGAMMA model
# Random starting seed '123456'
# 100 bootstraps
# 32 CPUs

clade=$(basename ${PWD})
rm -f RAxML_fail.log
for folder in $(ls -d ${PWD}/Subfamily_*); do
	subfamily=$(basename "${folder}" | cut -c 11-)
	jobname="${clade}.${subfamily}.RAxML"
	
	# Skip subfamilies where this job is already running
	if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "${jobname}"; then 
		continue
	fi
	
	if [ -f "${folder}/guidance2/G_TRIM" ] && [ ! -f "${folder}/guidance2/BAD_ALN" ] && [ ! -f "${folder}/raxml/RAxML_bestTree.${subfamily}" ]; then
		cd ${folder}
		sbatch -n32 -N1 -p preempt1d -J ${jobname} --wrap \
		"module load raxml/1.2.2
		rm -rf raxml && mkdir raxml
		
		raxmlHPC-PTHREADS \
		-s ${folder}/guidance2/MSA.PRANK.${subfamily}.nucl.NNN0.*.trimmed.aln \
		-n ${subfamily} \
		-m GTRGAMMA \
		-p 123456 \
		-N 100 \
		-f d \
		-T 32 \
		-w ${folder}/raxml
		
		if [ -f '${folder}/raxml/RAxML_bestTree.${subfamily}' ]; then
			rm -f raxml/*RUN*
		else
			echo -e '${subfamily}\tFAILED' >> ../RAxML_fail.log
		fi"
		cd ..
	fi
done
