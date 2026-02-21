# Utility python script for the <godon_stats.sh> bash script
# Mark branches under positive selection on the tree file
# Threshold to consider positive selection defined in qval_threshold

import re
import os
import sys
import argparse
from ete3 import Tree

# Description for script whith '-h'
parser = argparse.ArgumentParser(description = "Different functions needed to deal with godon_stats.sh")
optional_args = parser.add_argument_group(title = "Arguments (depending on the function called)")

# Define arguments with help message
optional_args.add_argument('-f', dest = 'function', help = "Function to call (tree_label,)")
optional_args.add_argument('-csv', dest = 'csv', help = "CSV file")
optional_args.add_argument('-godon', dest = 'godon', help = "Godon's output")
optional_args.add_argument('-it', dest = 'input_tree', help = "RAxML tree")
optional_args.add_argument('-ot', dest = 'output_tree', help = "Output tree")

# Group all arguments in a list
args = parser.parse_args()

def extract(godon, csv):
    with open(csv, 'w') as csv_out:
        csv_out.write(f"Subfamily,branch,lnull,lalt\n")
    subfamily = godon.split('.')[-2]
    with open(godon, 'r') as godon_output:
        for line in godon_output:
            if line.startswith("Testing branch "):
                branch_name = line.split()[2].split('/')[0]
            if line.startswith("Foreground branch: "):
                if line.split("#1")[0][-1] != ')':
                    # The next line was working, but still giving an error "SyntaxWarning: invalid escape sequence '\('", so I replace all the '(' by '>' in the tree just to get rid of the error message
                    # branch = re.split("\(|,", line.split("#1")[0])[-1]
                    pattern = line.split("#1")[0].replace("(", ">")
                    branch_name = re.split(">|,", pattern)[-1]
            if line.startswith("lnL0="):
                lnull = line.split(',')[0].split('=')[1]
                lalt = line.split(',')[1].split('=')[1]
                
                with open(csv, 'a') as csv_out:
                    csv_out.write(f"{subfamily},{branch_name},{lnull},{lalt}\n")

def tree_label(godon, input_tree, output_tree):
    ps = 0
    subfamily = str(godon).split('.')[1]
    dico_branch_num = {}

    with open(godon, 'r') as godon_file:
        for line_godon in godon_file:
            if line_godon.startswith("Testing branch "):
                branch_num = line_godon.split(" ")[2]
                branch_num = branch_num.split("/")[0]
                gtree = next(godon_file)[19:]
                godon_tree = Tree(gtree, format=1)
                
                count = 1
                for branch in godon_tree.traverse():
                    if "#1" in branch.name and branch.is_leaf():
                        dico_branch_num[count] = branch.name
                    elif "#1" in branch.name:
                        dico_branch_num[count] = branch_num
                    count += 1

    with open(input_tree, 'r') as intree:
        in_tree = intree.read()
    out_tree = Tree(in_tree, format=1)
    count = 1
    for branch in out_tree.traverse():
        if not branch.is_leaf():
            if count in dico_branch_num:
                branch.name = subfamily + "_gbranch" + dico_branch_num[count]
        count += 1
    out_tree.write(format=1, outfile=output_tree)


if __name__ == "__main__":
    
    if args.function == "extract":
        extract(args.godon, args.csv)
    elif args.function == "tree_label":
        tree_label(args.godon, args.input_tree, args.output_tree)
    else:
        print("ERROR: function '" + args.function + "' doesn't exist")