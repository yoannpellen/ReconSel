# Utility python script for the <guidance_multi.sh> bash script
# Performs various operations helpful for the positive selection pipeline

import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

# Description for script whith '-h'
parser = argparse.ArgumentParser(description = "Different functions needed to deal with Guidance2")
optional_args = parser.add_argument_group(title = "Arguments (depending on the function called)")

# Define arguments with help message
optional_args.add_argument('-f', dest = 'function', help = "Function to call (rename, masking, seq_check, delete, translate, reverse, fail)")
optional_args.add_argument('-sf', dest = 'subfamily', help = "Subfamily ID")
optional_args.add_argument('-nucl', dest = 'nucl_aln', help = "Nucleic acid alignment")
optional_args.add_argument('-prot', dest = 'prot_aln', help = "Amino acid alignment")
optional_args.add_argument('-log', dest = 'logf', help = "Sequence check up file")

# Group all arguments in a list
args = parser.parse_args()

def rename(subfamily):
    out_name = "guidance2/MSA.PRANK." + subfamily + ".fullname.aln"
    with open("guidance2/MSA.PRANK.aln", "r") as msa:
        msa_data = msa.read()
    with open("guidance2/Seqs.Codes", "r") as seq_code:
        for line in seq_code:
            OR = line.split("\t")[0]
            OR = ">"+OR
            code = line.split("\t")[1]
            code = ">"+code
            msa_data = msa_data.replace(code, OR+"\n")
    with open(out_name, "w") as msa_out:
        msa_out.write(msa_data)

def delete(logf):
    with open(logf, 'r') as log:
        for line in log:
            if line .startswith("Subfamily"):
                sf = line.split(":")[0][10:]
                mask = line.split(" ")[4][:-1]
                sf_file = "guidance2/MSA.PRANK." + subfamily + ".prot.NNN.aln"
                sf_file2 = "guidance2/MSA.PRANK." + subfamily + ".prot.NNN" + mask + ".trimmed.aln"
            elif line.startswith("   ->"):
                OR1 = line.split(" ")[4] + "_"
                OR2 = line.split(" ")[6] + "_"
                with open(sf_file, 'r') as sf_fasta:
                    for record in SeqIO.parse(sf_fasta, 'fasta'):
                        if record.id.startswith(OR1) or record.id.startswith(OR2):
                            os.system(rf"sed -i '/{record.id}/d' {sf_file2}")
                            os.system(rf"sed -i '/{record.seq}/d' {sf_file2}")

def reverse(prot_aln):
    # Builds a nucleotide alignment based on an amino acid alignment and the reference nucleotide fasta
    nucl_aln = prot_aln.replace(".prot.", ".nucl.")
    open(nucl_aln, "w").close()
    with open(prot_aln, 'r') as prot_file:
        for prot_record in SeqIO.parse(prot_file, 'fasta'):
            with open("guidance2/Seqs.Orig_DNA.fas.FIXED", 'r') as nucl_source:
                for nucl_record in SeqIO.parse(nucl_source, 'fasta'):
                    if prot_record.id == nucl_record.id:
                        seq = str(nucl_record.seq)
                        nuclcodon = [seq[i:i+3] for i in range(0,len(seq),3)]
                        newcodon = []
                        protflag = 0
                        for prot in prot_record.seq:
                            if prot == "-":
                                newcodon.append('---')
                            elif prot == "X":
                                newcodon.append('NNN')
                                protflag += 1
                            else:
                                newcodon.append(nuclcodon[protflag])
                                protflag += 1
                        nucl_record.seq = ''.join(newcodon)
                        with open(nucl_aln, 'a') as nucl:
                            nucl.write(">" + nucl_record.id + "\n")
                            nucl.write(str(nucl_record.seq) + "\n")

def masking(aln):
    # Checks the % of masked residues in the alignment
    percent_mask = 100
    masked_prot = 0
    total_prot = 0
    with open(aln, "r") as msa:
        for sequence in SeqIO.parse(msa, "fasta"):
            masked_prot += sequence.seq.count("X")
            total_prot += len(sequence.seq)
    percent_mask = masked_prot/total_prot
    print(int(percent_mask*100))

def translate(nucl, prot):
    # Translate nucleotide sequences into amino acids
    with open(nucl, 'r') as nucl_fasta:
        for sequence in SeqIO.parse(nucl_fasta, 'fasta'):
            seq_prot = Seq(str(sequence.seq))
            seq_prot = str(seq_prot.translate())
            with open(prot, 'a') as prot_fasta:
                prot_fasta.write(">" + sequence.id + "\n")
                prot_fasta.write(seq_prot + "\n")

def seq_check(aln, subfamily):
    # Checks sequences in the alignment right after masking to detect identical sequences
    # Identical sequences IDs are stored in a list
    # If "seq_1" and "seq_5" are identical, I enter "seq_1seq_5" and "seq_5seq_1" in the list
    # The print statements are added to the file "G_checkup.{subfamily}.txt"
    id_list=[]
    masked_prot = 0
    gap_prot = 0
    total_prot = 0
    mask = aln.split("NNN")[1][:3]
    with open(aln, 'r') as msa:
        for msa_sequence in SeqIO.parse(msa, "fasta"):
            if all(prot in "-X" for prot in str(msa_sequence.seq)):
                print(rf"   -> all masked: {msa_sequence.id} ")
                continue
            with open(aln, "r") as ref_msa:
                for ref_sequence in SeqIO.parse(ref_msa, "fasta"):
                    if msa_sequence.seq == ref_sequence.seq:
                        merge_name = rf"{msa_sequence.id}{ref_sequence.id}"
                        merge_name_reverse = rf"{ref_sequence.id}{msa_sequence.id}"
                        if msa_sequence.id != ref_sequence.id and merge_name not in id_list and merge_name_reverse not in id_list:
                            print(rf"   -> same seq: {msa_sequence.id} and {ref_sequence.id} ")
                            id_list.extend([merge_name,merge_name_reverse])
            masked_prot += msa_sequence.seq.count("X")
            gap_prot += msa_sequence.seq.count("-")
            total_prot += len(msa_sequence.seq)
    percent_mask = masked_prot/total_prot
    percent_gap = gap_prot/total_prot
    outlog = rf"G_checkup.{subfamily}.txt"
    open(outlog, 'w').close()
    with open(outlog, 'a') as log:
        log.write(rf"{subfamily}: masking at {mask}")
        log.write("\n")
        log.write(rf"   {percent_mask:.0%} of the alignment is masked")
        log.write("\n")
        log.write(rf"   {percent_gap:.0%} of gaps in the alignment")
        log.write("\n")
    if len(id_list) == 0:
        print("")

def trim(logf):
    # Creates a trimmed msa without the identical sequences listed in the log file created by the 'seq_check' function
    with open(logf) as log:
        first_line = log.readline()
    sf = first_line.split(":")[0]
    mask = first_line.split(" ")[3][:-1]
    with open(logf, 'r') as log:
        logdata = log.read()
    sf_in_prot = "guidance2/MSA.PRANK." + sf + ".prot.NNN" + mask + ".aln"
    sf_in_cds = "guidance2/MSA.PRANK." + sf + ".nucl.NNN" + mask + ".aln"
    sf_out_prot = "guidance2/MSA.PRANK." + sf + ".prot.NNN" + mask + ".trimmed.aln"
    open(sf_out_prot, 'w').close()
    sf_out_cds = "guidance2/MSA.PRANK." + sf + ".nucl.NNN" + mask + ".trimmed.aln"
    open(sf_out_cds, 'w').close()
    with open(sf_in_prot, 'r') as in_prot:
        for sequence in SeqIO.parse(in_prot, 'fasta'):
            if sequence.id not in logdata:
                with open(sf_out_prot, 'a') as out_prot:
                    SeqIO.write(sequence, out_prot, 'fasta')
    with open(sf_in_cds, 'r') as in_cds:
        for sequence in SeqIO.parse(in_cds, 'fasta'):
            if sequence.id not in logdata:
                with open(sf_out_cds, 'a') as out_cds:
                    SeqIO.write(sequence, out_cds, 'fasta')


if __name__ == "__main__":
    
    if args.function == "rename":
        rename(args.sf)
    elif args.function == "delete":
        delete(args.logf)
    elif args.function == "reverse":
        reverse(args.prot_aln)
    elif args.function == "masking":
        masking(args.prot_aln)
    elif args.function == "seq_check":
        seq_check(args.prot_aln, args.sf)
    elif args.function == "trim":
        trim(args.logf)
    elif args.function == "translate":
        translate(args.nucl_aln, args.prot_aln)
    else:
        print("ERROR: function '" + args.function + "' doesn't exist")
