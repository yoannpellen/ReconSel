# **ReconSel**

ReconSel is a pipeline designed to test for enrichment in positive selection in a species phylogeny, using multi-copy gene families. This work was published in Pellen and Privman (2025)[^1].
It allows to run dN/dS analysis on multi-copy gene families (e.g. ants odorant receptors, cytochrome P450, orthogroups, etc... referred hereafter as subfamilies), and report the results on the associated species tree using gene tree species tree reconciliation.

---

## :hammer_and_wrench: General Setup

### Input data

> [!NOTE]
> ReconSel was originally developed for Pellen *et al.* (2026)[^2], to work on ants odorant receptors. Therefore, it was written to fit the naming convention we had in this project, and will require to first make sure your data have what the scripts will look for.
> Many Bash script are using the Slurm Workload Manager to submit CPU intensive tasks. Make sure the default parameters fit your system capacities if you use it too, or adapt the scripts by removing the *sbatch* command if you don't.

To be able to run ReconSel, you need only two things as input: 
1. A rooted species tree in newick format.
2. A fasta files for every species you want to use, with all the sequences for every multi-copy gene family for this species.

ReconSel uses three different informations from the input data: the sequence ID, subfamily ID and species ID.
 - The leaves of the species tree are the species ID
 - Each species must have its own *.fasta* file, with nucleotide sequences. The name of the fasta file must start with the species ID. It must contain all the sequences, for every subfamily of the species.
 - Each sequence in the fasta files must be formatted as such: *sequenceID__subfamilyID__speciesID*. The double underscore will be used as a separator by many scripts to get the necessary info, so all your IDs can include any character except "__". E.g. if you have a sequence called *CfloOR120*, in subfamily *9E*, for the species *Camponotus_floridanus*, the sequence ID will be *CfloOR120__9E__Camponotus_floridanus* in the file *Camponotus_floridanus\*.fasta*.

### :computer: Softwares
Make sure you have the following softwares installed and running:
   - [Guidance2](https://guidance.tau.ac.il/source)
   - [PRANK](https://ariloytynoja.github.io/prank-msa/)
   - [RAxML](https://github.com/amkozlov/raxml-ng)
   - [Godon](https://bitbucket.org/Davydov/godon/src/master/)
   - [R](https://www.r-project.org/)
   - [GeneRax](https://github.com/BenoitMorel/GeneRax)

### :snake: Python dependencies
The following Python modules need to be installed:
   - [Biopython](https://biopython.org/wiki/Download)
   - [ete3](https://etetoolkit.org/download/)
   - [SciPy](https://scipy.org/install/)

### :file_cabinet: Directories structure

1. Keep all the scripts within the same directory.
2. Create a main working directory for each group of species you selected. This directory will be referred to as *Clade*.
3. Run all scripts from this main directory.
4. Some scripts can require minor edits before use:
   - List of species to select
   - Thresholds for statistics
   - Settings for Slurm

---

## :card_index_dividers: 1. Data Preparation

This script sets the necessary folders and files to running the pipeline. It takes as parameter the path to your species fasta files, and the list of species you selected.

```bash
python Path/to/prepare.py -dir <Path/to/species/fasta/folder> -species <species_ID_1> <species_ID_2> [<species_ID_3> ...]
```

This script creates a separate folder for each subfamily represented in your clade, and automatically removes:
- Sequences with internal stop codons
- Sequences which length is not divisible by 3
- Final stop codons

> [!WARNING]
> *Guidance2* needs at least 4 sequences per alignment to be able to work. 
> To check how many sequences you have in each subfamily, run:
> ```bash
> egrep -c ">" Subfamily_*/*.nucl.trimmed.fna
> ```
>
>In case you have subfamilies with less than 4 sequences, and if it's possible in your dataset, merge them with the closest one using:
> ```bash
> bash merge_subfamilies.sh <subfamily1> <subfamily2> [subfamily3 ...]
> ```

---

## :aquarius: 2. Multiple Sequence Alignment

This first step will generate a multiple sequence alignment for each subfamily, using the software *PRANK* as aligner, and *Guidance2* will calculate an alignment score, before masking unreliably aligned residues.

Every step must be launched individually, after the previous one finished. A file is created at the end of a successful run, to be used as a flag for the next step to be able to start.

1. `align`
   Runs Guidance2 with PRANK to generate codon alignments.
   Creates `ENDS_OK` flag.

   ```bash
   bash Path/to/guidance.sh -f align -p <Path/to/guidance.v2.02>
   ```

2. `mask`
   Masks low-confidence residues based on residue scores from Guidance2.
   As we don't want to loose too much of the alignment after masking unreliably aligned residues, we start with a threshold of 0.8 for the Guidance score, and if more than 20% of the alignment is masked, we then loop with -0.1 steps until less than 20% of the alignment is masked.
   Aborts if the percentage of masked residues exceeds 20% at threshold 0.1.
   Creates `G_CORREC` flag.

   ```bash
   bash Path/to/guidance.sh -f mask -p <Path/to/guidance.v2.02>
   ```

   >[!TIP]
   > If the masking aborts, it means that the alignment has a low score, most likely because there are too many differences between the sequences. You can try to subdivide the subfamily further based on groups of more similar sequences you can identify.

3. `trim` 
   Removes sequences that might have become identical after masking.
   - Creates `G_TRIM` flag

   ```bash
   bash Path/to/guidance.sh -f trim
   ```


4. `clean`  
   Cleans unnecessary files generated by *Guidance2*, keeping only the alignments and score files.

   ```bash
   bash Path/to/guidance.sh -f clean
   ```

---

## :deciduous_tree: 3. Phylogenetic Gene Trees Construction

Builds a phylogenetic tree for each subfamily, using the masked codon alignment generated by *Guidance*.
The default parameters in the script are:
- `-f d`: new rapid hill­climbing algorithm
- `-m GTRGAMMA`: nucleotide model
- `-p 123456`: random seed for parsimony inferences
- `-N 100`: bootstrap
- `-T 32`: number of CPUs allocated

Creates an output tree called `RAxML_bestTree*` wich will be used as a success flag.
Run as:
```bash
bash Path/to/raxml.sh
```

---

## 🔍 4. dN/dS inferences

Runs dN/dS inferences with the branch-site model, using as input the masked codon alignment from *Guidance2* and the associated gene tree from *RAxML*.
The default parameters in the script are:
- `test BS`: testing for positive selection using the branch-site model
- `--m0-tree`: use M0 model to estimate branch length
- `--ncat-codon-rate 4`: number of categories for the codon rate variation
- `--all-branches`: test all branches of the gene tree
- `-procs 32`: 32 CPUs allocated

Run as:
```bash
bash Path/to/godon.sh
```

---

## :bar_chart: 5. Godon Statistics

Run as:
```bash
bash Path/to/godon_stats.sh
```

---

## :deciduous_tree::evergreen_tree: 6. Gene Tree – Species Tree Reconciliation

### Step 1: Label Internal Branches of the Species Tree
```bash
python PATH/TO/generax.py -f treelabel -it <species_tree> -ot <new_species_tree>
```

### Step 2: Prepare Input for Generax 
Run as:
```bash
bash Path/to/generax.sh -f prep -t <new_species_tree>
```

### Step 3: Run Generax
Run as:
```bash
bash Path/to/generax.sh -f run -t <new_species_tree>
```

### Step 4: Build Mapping Table
Run as:
```bash
python Path/to/mapping_table.py
```

---

## :chart_with_upwards_trend: 7. Enrichment Test

Run as:
```bash
python Path/to/hypergeometric_test.py -map <mapping_file> []
```
---

## :book: References

[^1]: Pellen Y, Privman E. 2025. Testing for positive selection in multi-copy gene families using a reconciliation approach. :2025.08.18.670524. Available from: https://www.biorxiv.org/content/10.1101/2025.08.18.670524v1

[^2]: Pellen Y, Vizueta J, Beck E, Leibig J, Schrader L, Privman E. 2025. Adaptive evolution of odorant receptors is associated with elaborations of social organization in ants. :2025.08.20.670513. Available from: https://www.biorxiv.org/content/10.1101/2025.08.20.670513v2

