# Utility python script for the <generax.sh> bash script
# Various functions to prepare files for generax.sh

from ete3 import Tree
import argparse

# Description for script whith '-h'
parser = argparse.ArgumentParser(description = "Different functions needed to deal with Guidance2")
optional_args = parser.add_argument_group(title = "Arguments (depending on the function called)")

# Define arguments with help message
optional_args.add_argument('-f', dest = 'function', help = "Function to call")
optional_args.add_argument('-it', dest = 'intree', help = "Input tree")
optional_args.add_argument('-ot', dest = 'outtree', help = "Output tree")
optional_args.add_argument('-map', dest = 'mapping', help = "Name of mapping file")

# Group all arguments in a list
args = parser.parse_args()

def treelabel(intree, outtree):
    # Label every internal branch of the species tree with its number
    with open(intree, 'r') as in_tree:
        tree = in_tree.read()
    out_tree = Tree(tree, format=1)
    branch_num = 1
    for branch in out_tree.traverse():
        if not branch.is_leaf():
            branch.name = "spbranch"+str(branch_num)
        branch_num += 1
    out_tree.write(format=1, outfile=outtree)

def mapping(intree, mapping):
    # Create .link file for generax using only the gene tree file, getting the species name from the leaves name
    with open(intree, 'r') as intree:
        intree = intree.read()
    tree = Tree(intree, format=1)
    for leaf in tree.iter_leaves():
        species_name = leaf.name.split('__')[2]
        with open(mapping, 'a') as mapfile:
            mapfile.write(species_name+':'+leaf.name+'\n')

if __name__ == "__main__":
    if args.function == "mapping":
        mapping(args.intree, args.mapping)
    elif args.function == "treelabel":
        treelabel(args.intree, args.outtree)
    else:
        print("ERROR: function '" + str(args.function) + "' doesn't exist")
        