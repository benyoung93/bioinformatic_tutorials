# Species Tree Inference from Single Copy Orthologs

In this tutorial we’ll build a species tree from single-copy orthologs. These genes are ideal because they avoid problems with duplicates and give us a consistent set of markers across species.

We’ll use two complementary approaches:

Concatenation – combine all alignments into one big dataset and infer a tree with maximum likelihood.

Coalescent methods (e.g. ASTRAL) – first build gene trees, then summarize them into a species tree that accounts for gene tree conflict.

Finally, we’ll look at quartet support values (q1, q2, q3), which tell us how much the gene trees agree with each branch of the species tree versus supporting alternative histories. This gives a clearer picture of which parts of the tree are strongly supported and where there’s conflict.

## Step 0. General Setup 

All my written programs are in the bin directory of this tutorial page.  
Additionally, I have provided the conda environment yaml files for programs I used as I know my written programs work with these. 
Finally, below is the general directroy layout (identified using `tree -d`) for this analysis. You are welcome to use a similar makeup, or rename, whatever floats your ⛵️.  


## Step 1. Generating Cleaned and Aligned Single Copy Orthologs  

This step here is broadly the same as Step 1 - 4 in the Phylogenomic tutorial. As such, I do not go into detail here but you can refer to other tutorial if you want more information. 

```
cd /scratch/alpine/beyo2625/species_tree_tut
bin/rename.py \
proteomes \
modified_proteomes \
analysis_lists
```

Now things are cleaned and looking good we run `proteinortho`.  

```
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=xxxxx
#SBATCH --partition=xxxxx
#SBATCH --account=xxxxx
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --job-name=protortho
#SBATCH --error=/scratch/alpine/beyo2625/species_tree_tut/prot_ortho/protortho.err
#SBATCH --output=/scratch/alpine/beyo2625/species_tree_tut/prot_ortho/protortho.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxxxx
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=2

module purge 
eval "$(conda shell.bash hook)"
conda activate proteinorthoold_env

cd /scratch/alpine/beyo2625/species_tree_tut/prot_ortho

proteinortho6.pl \
-project=spectree \
-cpus=20 \
--temp=/scratch/alpine/beyo2625/temp_space \
/scratch/alpine/beyo2625/species_tree_tut/modified_proteomes/*.fasta
```

Next, we use `summarise_proteinortho.py` to identify the different sets of single copy orthologs (SCOs) from the `proteinortho` run. 

```{bash summarising proteinortho results}
cd /scratch/alpine/beyo2625/species_tree_tut
bin/summarise_proteinortho.py \
--extension .fasta \
--output_dir prot_ortho \
prot_ortho/spectree.proteinortho.tsv \
modified_proteomes
```
```
✅ Found 80 input FASTA files
100% (80 species): 69 SCOs written to prot_ortho/SCO_all.txt
95% (≥76 species): 335 SCOs written to prot_ortho/SCO_great95.txt
90% (≥72 species): 377 SCOs written to prot_ortho/SCO_great90.txt
85% (≥68 species): 388 SCOs written to prot_ortho/SCO_great85.txt
80% (≥64 species): 394 SCOs written to prot_ortho/SCO_great80.txt
75% (≥60 species): 398 SCOs written to prot_ortho/SCO_great75.txt
70% (≥56 species): 403 SCOs written to prot_ortho/SCO_great70.txt
65% (≥52 species): 420 SCOs written to prot_ortho/SCO_great65.txt
60% (≥48 species): 448 SCOs written to prot_ortho/SCO_great60.txt
55% (≥44 species): 476 SCOs written to prot_ortho/SCO_great55.txt
50% (≥40 species): 527 SCOs written to prot_ortho/SCO_great50.txt
```

Then we *grab* the SCOs from the results



## Step 2. Generating Concatenated Alignment and Partition File 


## Step 3. Identyfying the Best Model for each Gene Tree and Generating Gene Trees


## Step 4. Adding Best Models to Partition File


## Step 5. Collapsing Poorly Supported Branches


## Step 6. Building the Species Tree


## Step 7. Quartet Support Analysis


## Step 8. Visualising Tree with Quartets
