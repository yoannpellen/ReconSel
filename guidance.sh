#!/bin/bash

# Guidance2 script for ReconSel
# Each function successfully executed will create a file saying so, the next function can only be run if this file exists

aligning() {
	# Run Guidance2 for each subfamily separately
	# codon alignment using PRANK
	# Set to use 32 CPUs
	clade=$(basename "${PWD}")
	rm -f guidance_fail.txt
	for folder in $(ls -d ${PWD}/Subfamily_*); do
		sf=$(basename "${folder}" | cut -c 11-)
		jobname="${clade}.${subfamily}.Guidance_alignment"
		
		# Skip subfamilies where this job is already running
		if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "${jobname}"; then 
			continue
		fi
		
		if [ ! -f "${folder}/guidance2/ENDS_OK" ]; then
			rm -rf ${folder}/guidance2_err/ ${folder}/guidance2/ ${folder}/slurm*
			cd ${folder}
			sbatch -n32 -N1 -p preempt1d -J ${jobname} --wrap \
			"echo '${subfamily} Guidance2 with PRANK'
			module load prank
			
			perl ${guidance2_path}/www/Guidance/guidance.pl \
			--seqfile ${folder}/${subfamily}.nucl.trimmed.fna \
			--seqType codon \
			--outDir ${folder}/guidance2 \
			--program GUIDANCE2 \
			--msaProgram PRANK \
			--proc_num 32
			
			if [ ! -f '${folder}/guidance2/ENDS_OK' ]; then
				echo -e '${subfamily}\tFAILED' >> ../guidance_fail.log
			fi"
			cd ..
		fi
	done
}

masking() {
	# Starting from masking at 0.8, loops with -0.1 steps to mask less than 20% of the alignment
	subfamily=$(basename "${PWD}" | cut -c 11-)
	rm -f ${PWD}/guidance2/MSA.PRANK.${subfamily}*NNN* ${PWD}/guidance2/G_CORREC*
	cp ${PWD}/guidance2/MSA.PRANK.aln.With_Names ${PWD}/guidance2/MSA.PRANK.${subfamily}.nucl.aln
	for ((i=8; i>=0; i--)); do
		if [[ "${i}" -eq 0 ]]; then
			echo "The alignment of ${subfamily} has a bad score, no threshold masks less than 20%." > ${PWD}/guidance2/BAD_ALN
			break
		fi

		mask="0.${i}"
		perl ${guidance2_path}/www/Guidance/maskLowScoreResidues.pl \
		${PWD}/guidance2/MSA.PRANK.${subfamily}.nucl.aln \
		${PWD}/guidance2/MSA.PRANK.Guidance2_res_pair_res.scr \
		${PWD}/guidance2/MSA.PRANK.${subfamily}.nucl.NNN${mask}.aln \
		${mask} \
		nuc

		python ${scripts_path}/guidance_correction.py -f translate -nucl ${PWD}/guidance2/MSA.PRANK.${subfamily}.nucl.NNN${mask}.aln -prot ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN${mask}.aln
		masked=''
		masked=$(python ${scripts_path}/guidance_correction.py -f masking -prot ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN${mask}.aln)
		if [[ "${masked}" -lt 20 ]]; then
			echo "Final masking at ${mask} for subfamily ${subfamily}"
			break
		fi
		rm -f ${PWD}/guidance2/MSA.PRANK.${subfamily}.*.NNN*.aln
	done
	
	if [ ! -f "${PWD}/guidance2/BAD_ALN" ]; then
		# Checking for identical or fully masked sequences after translated to amino acids
		check_seq=''
		check_seq=$(python ${scripts_path}/guidance_correction.py -f seq_check -sf ${subfamily} -prot ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN${mask}.aln -log ${PWD}/G_checkup.${subfamily}.txt)
		echo "${check_seq}" >> ${PWD}/G_checkup.${subfamily}.txt
		touch ${PWD}/guidance2/G_CORREC
	fi
}

trimming() {
	# Delete sequences than can't be used after masking with Guidance2
	# Either two of them became identical, or one is completely masked
	clade=$(basename "${PWD}")
	for folder in $(ls -d ${PWD}/Subfamily_*); do
		subfamily=$(basename "${folder}" | cut -c 11-)
		jobname="${clade}.${subfamily}.Guidance_trimming"
		
		# Skip subfamilies where this job is already running
		if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "${jobname}"; then 
			continue
		fi
		
		if [ -f "${folder}/guidance2/G_CORREC" ] && [ ! -f "${folder}/guidance2/G_TRIM" ]; then
			cd ${folder}
			sbatch -n1 -N1 -p preempt1d,preempt7d -J ${jobname} --get-user-env --wrap \
			"python ${scripts_path}/guidance_correction.py -f trim -log G_checkup.${subfamily}.txt \
			&& touch ${folder}/guidance2/G_TRIM"
			cd ..
		fi
	done
}

cleaning() {
	# Delete all useless files after Guidance2 succeeds, basically keeping the alignments and score files
	# Can be runned even if all jobs aren't finished
	clade=$(basename "${PWD}")
	for folder in $(ls -d ${PWD}/Subfamily_*); do
		subfamily=$(basename "${folder}" | cut -c 11-)
		jobname="${clade}.${subfamily}.Guidance_cleaning"
		if [ -f "${folder}/guidance2/ENDS_OK" ]; then
			cd ${folder}/guidance2
			rm -rf BP COS.std *JalView* Seqs* *.tar.gz *html Sampled* *Without* *semphy.tree SampledOPVals.log ../slurm*.out slurm*.out
			cd ../..
		fi
	done
}

failing() {
	# Define source directory for all scripts based on the path for the current script
	guidance2_path=$1
	subfamily=$2
	scripts_path=$3 # It wouldn't automatically detect it the same way as other functions, so I give it manually

	echo "Extracting Guidance score from successful runs."
	${guidance2_path}/programs/msa_set_score/msa_set_score \
	${PWD}/guidance2/MSA.PRANK.aln \
	${PWD}/guidance2/MSA.PRANK.Guidance2 \
	-d ${PWD}/guidance2/BP/GUIDANCE2_MSA/ > ${PWD}/guidance2/MSA.PRANK.Guidance2.msa_set_score.std

	echo "Renaming sequences based on sf_${subfamily}/guidance2/Seq.Codes."
	python ${scripts_path}/guidance_correction.py -f rename -sf ${subfamily}

	echo "Masking unreliably aligned residues"
	for ((i=8; i>=1; i--)); do
		mask="0.${i}"
		perl ${guidance2_path}/www/Guidance/maskLowScoreResidues.pl \
		${PWD}/guidance2/MSA.PRANK.${subfamily}.fullname.aln \
		${PWD}/guidance2/MSA.PRANK.Guidance2_res_pair_res.scr \
		${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN.aln \
		${mask} \
		aa && touch ${PWD}/guidance2/ENDS_OK
		
		masked=""
		masked=$(python ${scripts_path}/guidance_correction.py -f masking -sf ${subfamily} -prot ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN.aln)
		if [[ "${masked}" -lt 20 ]]; then
			echo "Final masking at ${mask} for subfamily ${subfamily}."
			cp ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN.aln ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN${mask}.aln
			cp ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN.aln ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN${mask}.trimmed.aln \
			&& touch ${PWD}/guidance2/G_CORREC
			break
		fi
	done

	# Sometimes prot sequences end up being identical after masking by Guidance and it can cause error further in the analysis
	# This is to check if there is such case
	echo "Checking for identical sequences after masking and removing them."
	check=""
	check=$(python ${scripts_path}/guidance_correction.py -f seq_check -sf ${subfamily} -prot ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN${mask}.aln)
	if [[ ${check} != "" ]]; then
		python ${scripts_path}/guidance_correction.py -f delete -log G_checkup.${subfamily}.txt
	fi

	echo "Creating codon alignment based on prot alignment"
	python ${scripts_path}/guidance_correction.py -f reverse -prot ${PWD}/guidance2/MSA.PRANK.${subfamily}.prot.NNN${mask}.trimmed.aln \
	&& touch ${PWD}/guidance2/G_TRIM
	
	cd guidance2
	rm -rf ../guidance2_err BP COS.std *JalView* Seqs* *.tar.gz *html Sampled* *Without* *semphy.tree SampledOPVals.log ../slurm*.out
	cd ..
	echo "Done correcting failed Guidance run for subfamily ${subfamily}."
}

# Define source directory for all scripts based on the path for the current script
# Path to the main directory of Guidance2, should be called 'guidance.v2.02'
scripts_path="$(dirname "$(readlink -f "$0")")"

while getopts f:p: option
do 
    case "${option}"
        in
        f)function=${OPTARG};;
        p)guidance2_path=${OPTARG};;
    esac
done
# Remove the trailing "/" in the guidance path if it's there to be able to accept any form of path
guidance2_path="${guidance2_path%/}"

# Check which function to run based on the first argument when calling the script
case ${function} in

	align)
		aligning;;

	mask)
		# Easier to separate the loop and the call of the function through sbatch
		# The quote system in sbatch is giving me headaches
		clade=$(basename "${PWD}")
		for folder in $(ls -d ${PWD}/Subfamily_*); do
			subfamily=$(basename "${folder}" | cut -c 11-)
			jobname="${clade}.${subfamily}.Guidance_masking"
			
			# Skip subfamilies where this job is already running
			if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "${jobname}"; then 
				continue
			fi
			
			if [ -f "${folder}/guidance2/ENDS_OK" ] && [ ! -f "${folder}/guidance2/G_CORREC" ]; then
				cd ${folder}
				sbatch -n1 -N1 -p preempt1d -J ${jobname} --wrap \
				"bash ${scripts_path}/guidance.sh -f masking -p ${guidance2_path}"
				cd ..
			fi
		done;;

	masking)
		masking;;

	trim)
		trimming;;

	clean)
		cleaning;;

	fail)
		clade=$(basename "${PWD}")
		rm -f guidance_fail.txt
		scripts_path="$(dirname "$(readlink -f "$0")")"
		subfamily=$3
		cd sf_${subfamily}
		jobname="${clade}.${subfamily}.Guidance_fail"
		
		# Skip running subfamilies
		if squeue -u $USER -o "%.12i %.9P %.30j %.8u" | grep -q "${clade}.${subfamily}."; then 
			exit
		fi
		
		failing ${guidance2_path} ${subfamily} ${scripts_path};;
	
	*)
		echo "${function} does not exist."
		exit 1;;
	
esac
