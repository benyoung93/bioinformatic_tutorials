# Phylogenomic Tutorial  

So this phylogenomic tutorial will be using proteomes. Namely, you have *n* input proteomes, and you want to identify single copy orthologs (SCOs) and generate a phyogenomic tree from them. In most cases, you should get a nice set of SCOs between your species (100-300). This can not be the case however. As such, this tutorial also shows you a case in if you want to identify SCOs in >x% of species (e.g. >95% of species) as this can greatly bump up you number of SCOs, and programs like `RAxML` still work well with gaps present. Obviously you do not want massive stretches of gaps, but when you run RAxML it will tell you the percentage of gaps, and anything <15% is pretty good. Okay let us proceed.  

## Step 1: Proteome Cleaning and Preperation

So when downloading proteomes they can have  
* different endings (e.g. .fa, .fasta, .faa)
* different numbers of underscores/characters in the main body
* different naming conventions of the proteins within
* Be gunzipped, tar archived etc etc etc. 

This can cause alot of problems for the programs we are using, so it makes sense to standardise and clean everything that is going into the pipeline.  

I wrote a python script called `rename.py` which automates all of this for you. This python program renames all your files with  
* a standard prefix
* allows you to choose the suffix (i.e. .fasta)
* generates a mapping file of the new names to the old proteome names (important for when visualising your tree)

The one caveat here is that the input directory **MUST** only contain the proteomes you want to process, and also not have anything that is tar archived, gzipped etc. So prior to running the `./rename.py`, make sure everyting is unzipped, unarchived etc. 

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
```{bash}
./rename.py \
--prefix fungi \
protein_fastas \
modified_fastas \
analysis_lists
```

This will make a directory called `modified_fastas` with the renamed proteomes from the directory `protein_fastas`. The renamed proteomes will be `fungi01, fungi01 .... funginn`. It will also make an analysis directory to have all the mapping files and cleanup files produced.  

## Step 2: Identyfying and Selecting Single Copy Orthologs

For this step we use `proteinortho` - https://gitlab.com/paulklemm_PHD/proteinortho  

I have included in the `conda_envs` directory in this tutorial the yaml files that I know all my programs work with the generated outputs.   

For this tutorial I am running the pipeline on the following. 
* SCO in all species
* SCO in >95% of species  

I will only give the commands for the SCO all, but will provide the ouputs fromthe >95% to demonstrate programs and how the outputs differ. You will therfore see one command, but duplicated outputs. I have included my `slurm` scripts for this, but you will need to use your own HPC cluster scheduling tool if running on different supercomputers. 

```{bash}
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

Following the completion of `proteinortho`, we then use my written script `sumarise_proteinortho.py` to identify the SCO in all proteomes, and SCO percentages from 50% up to 95% in increments of 5. 

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
-e fasta \ ## the extension of the fasta files in the modified fasta directory
-o prot_ortho \
prot_ortho/auto_test_tree.proteinortho.tsv \
modified_fastas
```

This writes all the results to the `prot_ortho` directory created when running `proteinortho`, and specifies the modified fasta files from `rename.py` as the input.  

This will generate athe SCO in all and different percentages that you need in the next step. I recommend inspecting these to see what they look like, as this can affect your choice in the donwstream analysis. For here, we are going with SCO all and SCO > 95% species.  

At this point I like to generate all the downstream working directories. So whatever your working firectory is, i like to `cd` and do the following. This just means we do not need to worry about it when running the pipeline.  

```
cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline ## going to the project directory
mkdir sco_all sco_g95 \ ## making the SCO directories
sco_align_all sco_align_g95 \ ## making directories for the aligned SCOs
sco_align_all/loop_err_out sco_align_g95/loop_err_out \ ## directories in the aligned for the loop scripts for alignment
sco_clean_all sco_clean_g95 \ ## making diresctories for the cleaned SCOs
sco_comp_all sco_comp_g95 \ ## making the directories for completed SCOs with all species
raxml_all raxml_g95 ## making the raxml files. 
```

Now we want to select SCO and generate fasta files with each SCO in it. We will use the `proteinortho_grab_proteins.pl` in the `proteinortho` environmnet to do this.  

```
#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --qos=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --account=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --job-name=grabprots_90perc
#SBATCH --error=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/loop_err_out/grabprots.err
#SBATCH --output=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/loop_err_out/grabprots.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2

## activating the conda environment
module purge 
eval "$(conda shell.bash hook)"
conda activate proteinortho_env

cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline

## getting the SCOs in all species
proteinortho_grab_proteins.pl \
-exact \
-cpus=10 \
-tofiles=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_all \ ## path to the created directory to save SCO fastas
prot_ortho/SCO_all.txt \ ## generated SCO_all.txt file specyfying which SCOs we want
modified_fastas/*.fasta ## directory with the renamed proteomes

## getting the 95% SCOs
proteinortho_grab_proteins.pl \
-exact \
-cpus=10 \
-tofiles=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_g95 \
prot_ortho/SCO_great95.txt \
modified_fastas/*.fasta
```

If you then look in the `SCO_all` or `SCO_great95` they will have fasta files present, and the number of fasta files will equal rows - 1 in the .txt file specified.  

## Step 2: Aligning SCOs

The next step is to align all the sequences in each SCO with themselves. This uses `muscle` and is straightforward. The script below allows you to submit 1 file to the supercomputer, and then it submits a single job for each SCO that needs to be aligned.  

```
#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --qos=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --job-name=muscle
#SBATCH --error=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_align_all/loop_err_out/muscle.err
#SBATCH --output=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_align_all/loop_err_out/muscle.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX

# ---- USER CONFIG ----
## Edit these so they match your paths for things
BASE_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline
FASTA_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_all
ALIGN_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_align_all
LOG_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_align_all/loop_err_out
# ---------------------

cd "$FASTA_DIR" || exit 1
PALMATA=$(ls *.fasta | sed 's/\.fasta$//')

echo "Files going through muscle5 alignment:"
echo "$PALMATA"

for PALPAL in $PALMATA; do
    JOB_SCRIPT="$LOG_DIR/${PALPAL}_muscle.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --account=ucb423_asc2
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --job-name=${PALPAL}_muscle
#SBATCH --error=$LOG_DIR/${PALPAL}_muscle.err
#SBATCH --output=$LOG_DIR/${PALPAL}_muscle.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=beyo2625@colorado.edu

module purge
eval "\$(conda shell.bash hook)"
conda activate muscle_env

cd $BASE_DIR
muscle \\
    -super5 sco_all/${PALPAL}.fasta \\
    -output sco_align_all/${PALPAL}.align
EOF

    sbatch "$JOB_SCRIPT"
done
```

This will then create the same number of fasta files as in the `SCO_all` directory, but they will be aligned. To run this for the `SCO_great95` you would just substitute paths in to the respective files for processing.  

## Step 3. Trimming Samples

We want to clean the alignments we generated, and we will use `trimal` for this. This is a very quick and easy step so I will not go into this to much. Again, it will generate the same number of fasta files as in the `SCO_all` directory, and the `SCO_align_all` directory.  

```
mamba activate trimal_env

## for all SCOs
cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_align_all
PALMATA=$(ls *.align | sed 's/\.align$//') # making variable of all the files in the aligned directory

cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline

for PALPAL in $PALMATA ## loop script to iterate through all the files and generate cleaned ones
do
trimal \
-in sco_align_all/"$PALPAL".align \
-out sco_clean_all/"$PALPAL".trimmed \
-gappyout
done
```
