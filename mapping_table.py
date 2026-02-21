import os
import sys
from ete3 import Tree

# Define source directory for all scripts based on the path for the current script
scripts_path = os.path.dirname(__file__)

with open(f"{scripts_path}/320_tree_newick.reference_reconciliation.nwk", 'r') as speciestree:
    tree = speciestree.read()
species_tree = Tree(tree, format=1)
dic_species_branch_length = {}
for species_branch in species_tree.traverse():
    if "spbranch" in species_branch.name:
        species = species_branch.name.split('spbranch')[1]
    else:
        species = species_branch.name
    dic_species_branch_length[species] = str(species_branch.dist)

clade = os.path.basename(os.getcwd())
generax_dir = os.getcwd() + "/generax_" + clade + "/"
mapping = generax_dir + clade + ".csv"
godon_stats = clade + ".godonBS.all.fdr.csv"
mapping_stats = generax_dir + clade + ".mapping_table.csv"

# Making a first mapping file, without statistics
with open(mapping, 'w') as mapf:
    mapf.write("Clade,Subfamily,Gene branch,Gene branch length,Species branch,Species branch length\n")
for file_xml in os.listdir(generax_dir):
    xml = os.fsdecode(file_xml)
    if xml.endswith(".xml"):
        recGeneTree = False
        subfamily = xml.split('_')[1]
        xml = generax_dir + xml
        with open(xml, 'r') as xmlfile:
            for line_xml in xmlfile:
                if "<recGeneTree>" in line_xml:
                    recGeneTree = True
                if recGeneTree:
                    if '<name>' in line_xml and '_' in line_xml:
                        species_branch = next(xmlfile)
                        species_branch = next(xmlfile)
                        if 'species_0' in species_branch: # Skipping the root of the species tree
                            continue
                        if 'speciesLocation=' in species_branch:
                            gene_branch_xml = line_xml.split('<name>')[1]
                            gene_branch_xml = gene_branch_xml.split('<')[0]
                            subfamily_tree = generax_dir+rf"generax.{subfamily}.tree"
                            
                            with open(subfamily_tree, 'r') as genetree:
                                tree = genetree.read()
                            gene_tree = Tree(tree, format=1)
                            for gene_branch in gene_tree.traverse():
                                if gene_branch_xml == gene_branch.name:
                                    gene_branch_dist = str(gene_branch.dist)
                            
                            # Keeping only the number from the internal branches names
                            if "gbranch" in gene_branch_xml:
                                gene_branch_xml = gene_branch_xml.split('_gbranch')[1]
                            species_branch = species_branch.split('"')[1]
                            if 'spbranch' in species_branch:
                                species_branch = species_branch.split('spbranch')[1]
                            
                            with open(mapping, 'a') as mapf:
                                towrite = clade+","+subfamily+","+gene_branch_xml+","+gene_branch_dist+","+species_branch+","+dic_species_branch_length[species_branch]+"\n"
                                mapf.write(towrite)

# Building the mapping file with all statistics
with open(mapping_stats, 'w') as mapstatsf:
    mapstatsf.write("Clade,Subfamily,Gene branch,Gene branch length,Species branch,Species branch length,lnull,lalt,deltal2,pval,qval\n")

with open(mapping, 'r') as mapf:
    for line_map in mapf:
        clade = line_map.split(',')[0]
        subfamily_map = line_map.split(',')[1]
        gene_branch_map = line_map.split(',')[2]
        gene_branch_length = line_map.split(',')[3]
        species_branch = line_map.split(',')[4]
        species_branch_length = line_map.split(',')[5][:-1]
        
        with open(godon_stats, 'r') as godonstats:
            for line_godonstats in godonstats:
                subfamily_stats = line_godonstats.split(',')[0]
                gene_branch_stats = line_godonstats.split(',')[1]
                lnull = line_godonstats.split(',')[2]
                lalt = line_godonstats.split(',')[3]
                deltal2 = line_godonstats.split(',')[4]
                pval = line_godonstats.split(',')[5]
                qval = line_godonstats.split(',')[6][:-1]
                
                if subfamily_map == subfamily_stats and gene_branch_map == gene_branch_stats:
                    with open(mapping_stats, 'a') as map_stats:
                        map_stats.write(clade+","+subfamily_map+","+gene_branch_map+","+gene_branch_length+","+species_branch+","+species_branch_length+","+lnull+","+lalt+","+deltal2+","+pval+","+qval+"\n")

os.remove(mapping)