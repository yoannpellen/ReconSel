# Script to create a working directory for each subfamily for all selected species
# Run as 'python PATH/TO/prepare.py' from the clade directory created before

import argparse
import os
import sys
import glob
from Bio import SeqIO

parser = argparse.ArgumentParser(description = "Script to prepare the working directory")
optional_args = parser.add_argument_group(title = "Arguments:")
optional_args.add_argument("-dir", dest = "species_directory", help = "Path to the directory containing all the species fasta files.")
optional_args.add_argument("-species", nargs="+", dest = "list_species", help = "List of species.")

args = parser.parse_args()

# List of species to select
# Each element of the list matches the beginning of the file name associated with that species
species_list = args.list_species

# Path for the directory containing all the sequences for the species to run
species_dir = args.species_directory
root = os.getcwd()
clade = os.path.basename(root)

prepare_logs = f"{clade}.prep.log"
with open(prepare_logs,'w') as log:
    log.write(f"List of {len(species_list)} species used in {clade}: {str(species_list)}\n")

for species in species_list:
    # Check if every species has a fasta file before going further
    species_file = any(glob.iglob(os.path.join(species_dir, f"{species}*")))
    if not species_file_:
        print(f"\nNo fasta file found for {species}, check if the corresponding fasta file exists in {species_dir}.\n")
        sys.exit()
    
    for file in os.listdir(species_dir):
        species_fasta = os.fsdecode(file)
        if species_fasta.startswith(species):
            
            with open(prepare_logs, 'a') as log:
                log.write(f"\nTaking sequences from {species_fasta}, species {species_list.index(species)}/{len(species_list)}\n")
            
            with open(f"{species_dir}{species_fasta}", 'r') as speciesfasta:
                for sequence in SeqIO.parse(speciesfasta, "fasta"):
                    stop = False
                    subfamily = sequence.id.split('__')[1]
                    subfamily_dir = os.getcwd() + "/Subfamily_" + subfamily + "/"
                    os.makedirs(os.path.dirname(subfamily_dir), exist_ok=True)
                    subfamily_fasta_out = subfamily_dir + subfamily + ".nucl.fna"
            
                    # Checking if the nucleotide sequence is divisible by 3 (to be able to use codon models), removing the sequence from the dataset if it's not
                    if len(sequence.seq)%3 != 0:
                        with open(prepare_logs, 'a') as log:
                            log.write(f"    {str(sequence.id)} discarded because not divisible by 3\n")
                        stop = True
                        continue
                        
                    # Checking for a stop codon at the end of the sequence and removing it if it's found
                    nucl_seq = str(sequence.seq)
                    nucl_seq = nucl_seq.upper()
                    if nucl_seq.endswith(("TAA", "TAG", "TGA")):
                        sequence.seq = sequence.seq[:-3]
                    nucl_seq = str(sequence.seq)
                    nucl_seq = nucl_seq.upper()
                    
                    # Checking for internal stop codon, and removing the sequence from the dataset if it has one
                    for i in range(0,len(nucl_seq),3):
                        codon = nucl_seq[i:i+3]
                        if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                            with open(prepare_logs, 'a') as log:
                                log.write(f"    {str(sequence.id)} discarded because of internal stop codon\n")
                            stop = True
                            break
                            
                    if stop == False:
                        with open(subfamily_fasta_out, 'a') as subfamily_out:
                            SeqIO.write(sequence, subfamily_out, 'fasta')

# Checking for identical sequences within each subfamily and removing them
# Identical sequences are problematic for Godon, causing the following fatal error "panic: division of zero by zero or infinity by infinity"
# Part 1: list all identical sequences in the log file
with open(prepare_logs, 'a') as log:
    log.write("\nListing identical sequences removed in each subfamily.\n")
                        
for path, dirs, files in os.walk(root):
    for file in files:
        if file.endswith(".nucl.fna"):
            subfamily = file.split('.')[0]
            subfamily_fasta = "Subfamily_" + subfamily + "/" + file
            # Opening the subfamily file twice, and compare all sequences to each others
            with open(subfamily_fasta, 'r') as subfamily_fasta1:
                id_list=[]
                for subfamily_seq1 in SeqIO.parse(subfamily_fasta1, "fasta"):
                    with open(subfamily_fasta, "r") as subfamily_fasta2:
                        for subfamily_seq2 in SeqIO.parse(subfamily_fasta2, "fasta"):
                            if subfamily_seq1.seq == subfamily_seq2.seq:
                                if subfamily_seq1.id != subfamily_seq2.id and subfamily_seq1.id not in id_list and subfamily_seq2.id not in id_list:
                                    id_list.append(subfamily_seq1.id)
                if id_list:
                    with open(prepare_logs, 'a') as log:
                        log.write(f"    {subfamily} -> identical sequences: " + str(id_list)[1:-1] + "\n")

# Part 2: reading the log file and creating a trimmed version of the subfamily fasta, without identical sequences
# This file will be the one used as a start for the rest of the pipeline
with open(prepare_logs, 'r') as log:
    logdata = log.read()
for path, dirs, files in os.walk(root):
    for file in files:
        if file.endswith(".nucl.fna"):
            subfamily = file.split('.')[0]
            reference_fasta = "Subfamily_" + subfamily + "/" + subfamily + ".nucl.fna"
            trimmed_fasta = "Subfamily_" + subfamily + "/" + subfamily + ".nucl.trimmed.fna"
            open(trimmed_fasta, 'w').close()
            with open(reference_fasta, 'r') as reference:
                for sequence in SeqIO.parse(reference, 'fasta'):
                    if sequence.id not in logdata:
                        with open(trimmed_fasta, 'a') as trimmed:
                            SeqIO.write(sequence, trimmed, 'fasta')

print("\nRun [egrep -c '>' Subfamily_*/*.nucl.trimmed.fna] for subfamilies count.")
print("Subfamilies with less than 4 sequences need to be merged with the closest one in the phylogeny to be able to run Guidance2.")
print("To merge subfamilies, refer to your phylogeny and put the small subfamilies with the closest one, using <merge_subfamilies.sh Subfamily1 Subfamily2 [Subfamily3...]>.\n")
