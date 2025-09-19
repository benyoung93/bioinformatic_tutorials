# Phylogenomic Tutorial  

So this phylogenomic tutorial will be using proteomes. Namely, you have *n* input proteomes, and you want to identify single copy orthologs (SCOs) and generate a phyogenomic tree from them. In most cases, you should get a nice set of SCOs between your species (100-300). This can not be the case however. As such, this tutorial also shows you a case in if you want to identify SCOs in >x% of species (e.g. >95% of species) as this can greatly bump up you number of SCOs, and programs like `RAxML` still work well with gaps present. Obviously you do not want massive stretches of gaps, but when you run RAxML it will tell you the percentage of gaps, and anything <15% is pretty good. Okay let us proceed.  

## Step 1: Proteome Cleaning and Preperation

So when downloading proteomes they can have  
* different endings (e.g. .fa, .fasta, .faa)
* different numbers of underscores/characters in the main body
* different naming conventions of the proteins within

This can cause alot of problems for the programs we are using, so it makes sense to standardise and clean everything that is going into the pipeline.  

I wrote a python script called `rename.py` which automates all of this for you. This python program renames all your files with  
* a standard prefix
* allows you to choose the suffix (i.e. .fasta)
* generates a mapping file of the new names to the old proteome names (important for when visualising your tree)

```
./rename.py --help
usage: rename.py [-h] [-p PREFIX]
                 input_dir output_dir_fastas output_dir_summary

Standardize FASTA filenames and headers for phylogenomics projects, generate
species list, separate summary files, and show progress.

positional arguments:
  input_dir             Directory containing input FASTA files (only FASTAs
                        should be present)
  output_dir_fastas     Directory where cleaned/renamed FASTA files will be
                        saved
  output_dir_summary    Directory where summary files (mapping, cleanup,
                        species list) will be saved

optional arguments:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Prefix for renamed files and headers (default:
                        'sample')
```

If doing the following  
```
./rename.py \
--prefix fungi \
protein_fastas \
modified_fastas \
analysis_lists
```

This will make a directory called `modified_fastas` with the renamed files. It will also make an analysis directory to have all the mapping files and cleanup files produced.  

## Step 2: Identyfying Single Copy Orthologs

For this step we use `proteinortho` - https://gitlab.com/paulklemm_PHD/proteinortho  

I have included in the `conda_envs` directory in this tutorial the yaml files that I know all my programs work with the ouputs of things from programs.  

For this tutorial I am running the pipeline on the following. 
* SCO in all species
* SCO in >95% of species  

In the tutorial, I will only give the commands for the SCO all, but will provide the ouputs of both to demonstrate programs. You will therfore see one command, but duplicated outputs. I have included my `slurm` scripts for this, but this will depend on the scheduling system you use, as well as account and partitions etc etc.  

```
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=normal
#SBATCH --partition=amilan
#SBATCH --account=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --job-name=protortho
#SBATCH --error=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/prot_ortho/protortho.err
#SBATCH --output=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/prot_ortho/protortho.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXX
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=2

module purge 
eval "$(conda shell.bash hook)"
conda activate proteinortho_env

cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/prot_ortho

proteinortho6.pl \
-project=auto_test_tree \
-cpus=20 \
--temp=/scratch/alpine/beyo2625/temp_space \
../modified_fastas/*.fasta
```

As you can see, we use the `modified_fasta` directory that was output in the `rename.py` script with all the cleaned proteomes.  

Following the completion of `proteinortho`, we then use my written script `sumarise_proteinortho.py` to  
* identify SCO in all and percentages from 50% up to 95% in increments of 5
* give the desired extension specified in the `rename.py` script for the modified proteomes. 

```
usage: sumarise_proteinortho.py [-h] [-e EXTENSION] [-o OUTPUT_DIR] tsv_file fasta_dir

Process ProteinOrtho output to extract SCOs at different percentage thresholds.

positional arguments:
  tsv_file              ProteinOrtho .tsv file
  fasta_dir             Directory with input FASTA files

options:
  -h, --help            show this help message and exit
  -e EXTENSION, --extension EXTENSION
                        FASTA file extension (default: .faa)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Directory for output SCO files
```
```
cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline
./process_proteinortho.py \
-e fasta \
-o prot_ortho \
prot_ortho/auto_test_tree.proteinortho.tsv \
modified_fastas
```

This writes all the results to the `prot_ortho` directory created when running `proteinortho`, and specifies the modified fasta files from `rename.py` as the input.  

This will generate a set of files that you need in the next step to make directories of all the SCO for downstream steps.  
