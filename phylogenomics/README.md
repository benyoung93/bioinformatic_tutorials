# Phylogenomic Tutorial  

So this phylogenomic tutorial will be using proteomes. Namely, you have *n* input proteomes, and you want to identify single copy orthologs (SCOs) and generate a phyogenomic tree from them. In most cases, you should get a nice set of SCOs between your species (100-300). This can not be the case however. As such, this tutorial also shows you a case in if you want to identify SCOs in >x% of species (e.g. >95% of species) as this can greatly bump up you number of SCOs, and programs like `RAxML` still work well with gaps present. Obviously you do not want massive stretches of gaps, but when you run RAxML it will tell you the percentage of gaps, and anything <15% is pretty good. Okay let us proceed.  

Quick side note, for the written programs in the `bin` directory, all you need to do is copy and paste these into terminal when in `nano`. You can then save them as a name, make them executable (`chmod +x [programname.py]`) and use them (`./[programname.py]`. The `./` is just saying that the executable is in the current working ddirectory, you can also hard path (`path/to/executable/[programname.py]` if you want. These programs also work without needing to load any environments or have specific tools present which is nice.

PUT IN STUFF ABOUT CONDA ENVS

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

This will make a directory called `modified_fastas` with the renamed proteomes from the directory `protein_fastas`. The renamed proteomes will be `fungi01, fungi01 .... funginn`. It will also make an analysis directory to have all the mapping files and cleanup files produced. One file produced is a mapping between the new name and the old proteome name. This is very useful for adding to a metadata file which you can then import into your tree visualisation program and specify your tip labels.  

## Step 2: Identyfying and Selecting Single Copy Orthologs

For this step we use `proteinortho` - https://gitlab.com/paulklemm_PHD/proteinortho  

I have included in the `conda_envs` directory in this tutorial the yaml files that I know all my programs work with the generated outputs.   

For this tutorial I am running the pipeline on the following. 
* SCO in all species
* SCO in >95% of species  

I will only give the commands for the SCO all, but will provide the ouputs fromthe >95% to demonstrate programs and how the outputs differ. You will therfore see one command, but duplicated outputs. I have included my `slurm` scripts for this, but you will need to use your own HPC cluster scheduling tool if running on different supercomputers. 

For CU Boulder, this link goes to the RC Documentation for the flags - https://curc.readthedocs.io/en/latest/clusters/alpine/alpine-hardware.html

```{bash}
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --account=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=40G
#SBATCH --job-name=protortho
#SBATCH --error=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/prot_ortho/protortho.err
#SBATCH --output=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/prot_ortho/protortho.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX
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
./sumarise_proteinortho.py \
-e fasta \ ## the extension of the fasta files in the modified fasta directory
-o prot_ortho \
prot_ortho/auto_test_tree.proteinortho.tsv \
modified_fastas
```

And the output will print to the terminal.  

```
✅ Found 69 input FASTA files
100% (69 species): 9 SCOs written to prot_ortho/SCO_all.txt
95% (≥66 species): 365 SCOs written to prot_ortho/SCO_great95.txt
90% (≥63 species): 837 SCOs written to prot_ortho/SCO_great90.txt
85% (≥59 species): 1337 SCOs written to prot_ortho/SCO_great85.txt
80% (≥56 species): 1625 SCOs written to prot_ortho/SCO_great80.txt
75% (≥52 species): 1987 SCOs written to prot_ortho/SCO_great75.txt
70% (≥49 species): 2171 SCOs written to prot_ortho/SCO_great70.txt
65% (≥45 species): 2410 SCOs written to prot_ortho/SCO_great65.txt
60% (≥42 species): 2541 SCOs written to prot_ortho/SCO_great60.txt
55% (≥38 species): 2691 SCOs written to prot_ortho/SCO_great55.txt
50% (≥35 species): 2780 SCOs written to prot_ortho/SCO_great50.txt
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

And then if we use the `tree` command we can see out directory layout.  

```
cd 
tree -d ## only print directories
.
├── analysis_lists
├── conda_envs
├── loop_err_out
├── modified_fastas
├── protein_fastas
├── prot_ortho
├── raxml_all
├── raxml_g95
├── sco_align_all
│   └── loop_err_out
├── sco_align_g95
│   └── loop_err_out
├── sco_all
├── sco_clean_all
├── sco_clean_g95
├── sco_comp_all
├── sco_comp_g95
└── sco_g95
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
#SBATCH --job-name=grabprots_all
#SBATCH --error=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_all/grabprots.err
#SBATCH --output=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_all/grabprots.out
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
```

If you then look in the `SCO_all` or `SCO_great95` they will have fasta files present, and the number of fasta files will equal rows - 1 in the .txt file specified.  

## Step 3: Aligning SCOs

The next step is to align all the sequences in each SCO with themselves. This uses `muscle` and is straightforward. The script below allows you to submit 1 job to the supercomputer, and then it submits a single job for each SCO that needs to be aligned.  

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
## Edit these so they match your paths for things wooooooooo 🫠
BASE_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline
FASTA_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_all
ALIGN_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_align_all
LOG_DIR=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/sco_align_all/loop_err_out
# ---------------------

## This is making a variable with all your sample names (i.e. fungi01, fungi02, fungi03 .... etc etc). 
cd "$FASTA_DIR" || exit 1
PALMATA=$(ls *.fasta | sed 's/\.fasta$//')

echo "Files going through muscle5 alignment:"
echo "$PALMATA"

for PALPAL in $PALMATA; do
    JOB_SCRIPT="$LOG_DIR/${PALPAL}_muscle.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --account=xxxxx
#SBATCH --partition=xxxxx
#SBATCH --qos=xxxxx
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --job-name=${PALPAL}_muscle
#SBATCH --error=$LOG_DIR/${PALPAL}_muscle.err
#SBATCH --output=$LOG_DIR/${PALPAL}_muscle.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxxxx

module purge
eval "\$(conda shell.bash hook)"
conda activate muscle_env

cd $BASE_DIR
muscle \\
    -super5 $FASTA_DIR/${PALPAL}.fasta \\
    -output $ALIGN_DIR/${PALPAL}.align
EOF

    sbatch "$JOB_SCRIPT"
done
```

This will then create the same number of fasta files as in the `SCO_all` directory, but they will be aligned. To run this for the `SCO_great95` you would just substitute paths in to the respective files for processing.  

## Step 4. Trimming Samples

We want to clean the alignments we generated, and we will use `trimal` for this. This is a very quick and easy step so I will not go into this to much. So easy that you do not even have to submit a job (and and run on the login node hehehehhe). Again, it will generate the same number of fasta files as in the `SCO_all` directory, and the `SCO_align_all` directory.  

```
mamba activate trimal_env ## activate the trimal environment

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

## Step 5. Completing the Missing Species in SCOs

If you are running true SCO (i.e. SCO_all as we are doing in this tutorial) you do not need to run this step. Saying that, how I have written the python programs allow this to be included and for ease of use you can include it.  

So for SCO that do not have all species present (i.e. you are using >95% SCO and thus there will be SCOs with < the total number of species included) you need to coplete them or when you concatenate it will create different sequence lengths and cause super whacky trees generated using `RAxML`.  

I wrote the python script `complete_missing_species.py` that will do this for you, and print the ouput saying for each SCO what species are added. For each SCO, it identifies the length of the sequence (which, after trimming will be the same for all the sequences in the file) and then adds ---- padding sequences for the missing species (as they have non alignment).  

The `species_file` for this step we generated when we ran `rename.py`. It is in the generated `analysis_list` directory that was made wehn you ran `rename.py` and it should be called `species_list.fa`.  

```
./complete_missing_species.py --help
usage: complete_missing_species.py [-h] --species_file SPECIES_FILE
                                   --input_dir INPUT_DIR
                                   [--input_ext INPUT_EXT] --output_dir
                                   OUTPUT_DIR [--output_ext OUTPUT_EXT]
                                   [--line_length LINE_LENGTH]
                                   [--missing_report MISSING_REPORT]

Fill missing species in orthogroup FASTA files with padded sequences.

optional arguments:
  -h, --help            show this help message and exit
  --species_file SPECIES_FILE
                        File with expected species list (lines like >species)
  --input_dir INPUT_DIR
                        Directory with cleaned orthogroup files
  --input_ext INPUT_EXT
                        Input file extension to select (default: .fa)
  --output_dir OUTPUT_DIR
                        Directory to write completed orthogroup files
  --output_ext OUTPUT_EXT
                        Output file extension (default: .fasta)
  --line_length LINE_LENGTH
                        Max characters per sequence line (default: 60)
  --missing_report MISSING_REPORT
                        Optional path to write a table of missing species per
                        orthogroup
```

So here, we would run

```
cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline

## for sco all (not needed but to show nothing made in the missing report)
./complete_missing_species.py \
--species_file analysis_lists/species_list.fa \
--input_dir sco_clean_all \
--input_ext .trimmed \
--output_ext .comp \
--output_dir sco_comp_all \
--missing_report analysis_lists/sco_all_missing_species.txt
```

This gives the output of the following 

```
.....
Processing SCOall.txt.OrthoGroup2.trimmed
Processing SCOall.txt.OrthoGroup3.trimmed
Processing SCOall.txt.OrthoGroup4.trimmed
Processing SCOall.txt.OrthoGroup5.trimmed
Processing SCOall.txt.OrthoGroup6.trimmed
Processing SCOall.txt.OrthoGroup7.trimmed
Processing SCOall.txt.OrthoGroup8.trimmed
✅ Missing species report saved to analysis_lists/sco_all_missing_species.txt

✅ All output files passed validation: correct species count and equal lengths.
```

As we used the SCO in all species, there is nothing to be completed, so yay.  

But, if we had used SCO in >95% species, you would see something like this.  

```
.......
Processing SCOgreat95.txt.OrthoGroup98.trimmed
⚠️  Missing species in SCOgreat95.txt.OrthoGroup98.trimmed: fungi26
⚠️  Missing species in SCOgreat95.txt.OrthoGroup98.trimmed: fungi44
Processing SCOgreat95.txt.OrthoGroup99.trimmed
⚠️  Missing species in SCOgreat95.txt.OrthoGroup99.trimmed: fungi26
⚠️  Missing species in SCOgreat95.txt.OrthoGroup99.trimmed: fungi33
✅ Missing species report saved to analysis_lists/sco_g95_missing.txt

✅ All output files passed validation: correct species count and equal lengths.
```

So here you can see that for 
* Orthogroup 98, padding sequences (i.e. ----) were added for fungi26 and fungi44
* Orthogorup 99, padding sequences were added for fungi26 and fungi33

After that, the script checks everything looks good, and you can see we got the ticks for those. If something is of, you get a X.  

## Step 6. Concatenation of the SCOs

The final step before tree generation is combining all the SCO for each species into one long sequence. I have again written a python script called `./combine_seqs.py` which does this for you, and runs some checks.  

```
./combine_seqs.py --help
usage: combine_seqs.py [-h] --expected_file EXPECTED_FILE --output_file
                       OUTPUT_FILE --ortholog_files ORTHOLOG_FILES
                       [ORTHOLOG_FILES ...] [--line_length LINE_LENGTH]

Concatenate ortholog FASTA files across species and validate sequence lengths.

optional arguments:
  -h, --help            show this help message and exit
  --expected_file EXPECTED_FILE
                        File with expected species list (lines like >species)
  --output_file OUTPUT_FILE
                        Path to write concatenated FASTA
  --ortholog_files ORTHOLOG_FILES [ORTHOLOG_FILES ...]
                        Ortholog FASTA files to concatenate
  --line_length LINE_LENGTH
                        Wrap FASTA sequences at this line length (0 = no
                        wrapping)
```

So for the SCO all we using here, we would do the following. 

```
cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline
./combine_seqs.py \
--expected_file analysis_lists/species_list.fa \
--output_file sco_all.fasta \
--ortholog_files sco_comp_all/* \
--line_length 60
```

And this gives the ouput 

```
✅ All species have the same concatenated length: 5164 bp

Concatenated FASTA file created: sco_all.fasta
```

So we know for all the species, the sequence length is 5,164 base paris, and they are all the same length. Woop.  

Again, just for an example, here is the ouput of SCO >95% species. 

```
✅ All species have the same concatenated length: 217065 bp

Concatenated FASTA file created: sco_ag95.fasta
```

So you can see this one, instead of 5,164 bp in the SCO all, had a alignment length for each species of 217,065 bp. This is the power of using SCO > xx% of species, you can seriously up the sequence length that is used in `RAxML`. 

## Step 7. RAxML Tree Generation

So the final stage is to use your generated aligned and concatenated file in `RAxML`. I am not going to go into this one to much, but the below is the code I used to generate the trees for the SCO all. 

```
#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --qos=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --account=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=60G
#SBATCH --job-name=raxml_100bs
#SBATCH --error=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/raxml_all/raxml_100bs.err
#SBATCH --output=/scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/raxml_all/raxml_100bs.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=2

module purge 
eval "$(conda shell.bash hook)"
conda activate raxml_env

cd /scratch/alpine/beyo2625/laboul_all/laboul_trees/auto_pipeline/raxml_all

raxmlHPC-PTHREADS-AVX2 \
-f a \
-s ../sco_all.fasta \
-n autopipeline_scoall_100bs_root \
-m PROTGAMMAAUTO \
-x 8212 \
-N 100 \
-p 1176 \
-T 50 \
-o fungi29,fungi51 ## samples used as outgroups and rooting in RAxML. Not needed but I like to do. 
```
