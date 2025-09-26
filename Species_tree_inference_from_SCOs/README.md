# Species Tree Inference from Single Copy Orthologs

In this tutorial we’ll build a species tree from single-copy orthologs. These genes are ideal because they avoid problems with duplicates and give us a consistent set of markers across species.

We’ll use two complementary approaches:

**Concatenation** – combine all alignments into one big dataset and infer a tree with maximum likelihood.

**Coalescent methods** – build gene trees, then summarize them into a species tree that accounts for gene tree conflict.

Finally, we’ll look at quartet support values (q1, q2, q3), which tell us how much the gene trees agree with each branch of the species tree versus supporting alternative histories. This gives a clearer picture of which parts of the tree are strongly supported and where there’s conflict.

## Step 0. General Setup 

All my written programs are in the bin directory of this tutorial page.  
Additionally, I have provided the conda environment yaml files for programs I used as I know my written programs work with these.  
Finally, below is the general directroy layout (identified using `tree -d`) for this analysis. You are welcome to use a similar makeup, or rename, whatever floats your ⛵️.  

```
cd /scratch/alpine/beyo2625/species_tree_tut # root project working directory
tree -d
.
├── analysis_lists
├── astral_all
├── bin
├── collapsed_trees
├── iqtree_all
├── iqtree_all_reest
├── iqtree_genetree_all ## There is a directory for every single SCO you using, I have cut out most for viewing ease
│   ├── loop_err_out
│   ├── SCOall.txt.OrthoGroup0
│   ├── SCOall.txt.OrthoGroup1
│   ├── SCOall.txt.OrthoGroup10
    ...............
│   └── SCOall.txt.OrthoGroupNN
├── modified_proteomes
├── proteomes
├── prot_ortho
├── sco_align_all
│   └── loop_err_out
├── sco_all
├── sco_all_clean_ngn
└── sco_clean_all
```


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

```
cd /scratch/alpine/beyo2625/species_tree_tut
mamba activate proteinortho_env

proteinortho_grab_proteins.pl \
-exact \
-cpus=5 \
-tofiles=/scratch/alpine/beyo2625/species_tree_tut/sco_all \
prot_ortho/SCO_all.txt \
modified_proteomes/*.fasta
```

Then we align the SCOs

```
#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --qos=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --job-name=muscle
#SBATCH --error=/scratch/alpine/beyo2625/species_tree_tut/sco_align_all/loop_err_out/muscle.err
#SBATCH --output=/scratch/alpine/beyo2625/species_tree_tut/sco_align_all/loop_err_out/muscle.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX

# ---- USER CONFIG ----
BASE_DIR=/scratch/alpine/beyo2625/species_tree_tut
FASTA_DIR=/scratch/alpine/beyo2625/species_tree_tut/sco_all
ALIGN_DIR=/scratch/alpine/beyo2625/species_tree_tut/sco_align_all
LOG_DIR=/scratch/alpine/beyo2625/species_tree_tut/sco_align_all/loop_err_out
# ---------------------

cd "$FASTA_DIR" || exit 1
PALMATA=$(ls *.fasta | sed 's/\.fasta$//')

echo "Files going through muscle5 alignment:"
echo "$PALMATA"

for PALPAL in $PALMATA; do
    JOB_SCRIPT="$LOG_DIR/${PALPAL}_muscle.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --qos=XXXXX
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem=10G
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --job-name=${PALPAL}_muscle
#SBATCH --error=$LOG_DIR/${PALPAL}_muscle.err
#SBATCH --output=$LOG_DIR/${PALPAL}_muscle.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX

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

Next we trim them 

```
mamba activate trimal_env

cd /scratch/alpine/beyo2625/species_tree_tut/sco_align_all
PALMATA=$(ls *.align | sed 's/\.align$//')

cd /scratch/alpine/beyo2625/species_tree_tut
for PALPAL in $PALMATA
do
trimal \
-in sco_align_all/"$PALPAL".align \
-out sco_clean_all/"$PALPAL".trimmed \
-gappyout
done
```

And finally we remove the `_nn|` (with nn being the gene number) from SCOs so we can concatenate by species and run the rest of the pipeline. Yay 🎉. 

```
cd /scratch/alpine/beyo2625/species_tree_tut
mkdir sco_all_clean_ngn
cd sco_clean_all
for f in *.trimmed; do
    base=$(basename "$f")
    sed 's/_\([0-9]\+\)|//' "$f" > "../sco_all_clean_ngn/$base"
done
```

## Step 2. Generating Concatenated Alignment and Partition File 

We will be using `IQTree3` for alot of the remainder of this tutorial/pipeline. The first step we do is use `iqtree3` to generatea concatenated alignment and partition file

```
mamba activate iqtree_env

cd /scratch/alpine/beyo2625/species_tree_tut/iqtree_all

iqtree3 \
-p ../sco_all_clean_ngn/ \
--out-aln alignments_concat.nex
```
```
Alignment was printed to alignments_concat.nex
Partition information was printed to alignments_concat.nex.nex
Partition information in Raxml format was printed to alignments_concat.nex.partitions
```

These will not be used till later on in the pipeline, but I like to generate them early.  

There will also be two files present in our project directory as we ran this in the terminal and not as a job, these can just be `mv` to the `iqtree_all` directory to maintain organisation.  

## Step 3. Identyfying the Best Model for each Gene Tree and Generating Gene Trees

Now we identify the best model for each SCO we are using (in our case for each of the 69 SCO between all species) and then generate the subsequent gene tree.  

There are **alot** of different ways to do this, but I have settled on using `IQTree3` as it does both of these steps in one command which is super nice.  

```
#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --account=xxxxx
#SBATCH --partition=xxxxx
#SBATCH --qos=xxxxx
#SBATCH --nodes=1
#SBATCH --mem=2G
#SBATCH --job-name=iqtree_model_tree
#SBATCH --error=/scratch/alpine/beyo2625/species_tree_tut/iqtree_genetree_all/loop_err_out/mod_tree.err
#SBATCH --output=/scratch/alpine/beyo2625/species_tree_tut/iqtree_genetree_all/loop_err_out/mod_tree.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxxxx

# ---- USER CONFIG ----
BASE_DIR=/scratch/alpine/beyo2625/species_tree_tut
FASTA_DIR=/scratch/alpine/beyo2625/species_tree_tut/sco_all_clean_ngn
OUTPUT_DIR=/scratch/alpine/beyo2625/species_tree_tut/iqtree_genetree_all
LOG_DIR=/scratch/alpine/beyo2625/species_tree_tut/iqtree_genetree_all/iqtree_genetree_all/loop_err_out

# ---------------------

cd "$FASTA_DIR" || exit 1
PALMATA=$(ls *.trimmed | sed 's/\.trimmed$//')

echo "Files going through iqtree3 analysis:"
echo "$PALMATA"

for PALPAL in $PALMATA; do
    JOB_SCRIPT="$LOG_DIR/${PALPAL}_iqtree.sh"

    cat > "$JOB_SCRIPT" <<EOF
#!/bin/bash
#SBATCH --account=xxxxx
#SBATCH --partition=xxxxx
#SBATCH --qos=xxxxx
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=2
#SBATCH --job-name=${PALPAL}_iqtree
#SBATCH --error=$LOG_DIR/${PALPAL}_iqtree.err
#SBATCH --output=$LOG_DIR/${PALPAL}_iqtree.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxxxx

module purge
eval "\$(conda shell.bash hook)"
conda activate iqtree_env

cd $OUTPUT_DIR
mkdir -p ${PALPAL}
cd ${PALPAL}

iqtree3 \\
    -s $FASTA_DIR/${PALPAL}.trimmed \\
    -mset LG,WAG,JTT,Dayhoff \\
    -merit AIC \\
    --prefix ${PALPAL} \\
    -B 1000 \\
    -alrt 1000 \\
    -nt AUTO
EOF

    sbatch "$JOB_SCRIPT"
done
```
There are a few things to note from this script. 

For the `-mset` flag I restricted it to specific inputs. This is because if you use `AUTO` selection in `IQTree3` it chooses the best models based solely on maths, and these may not be approproate for your organism. For example if I kept the `auto` flag here it chooses viral and HIV models which are not appropriate for fungi. Depnding on your organism you will need to carfully use the `IQTree` documentation to pick the `-mset` values you want to use, or follow how others have done it in your organisms in the future.  

* LG -> a modern empirical model built from a large dataset of globular proteins and generally performs best for fungi.
* WAG -> an earlier empirical model derived from a wide range of proteins and is a solid, widely used alternative to LG.
* JTT -> a classic substitution model from the early 1990s based on a smaller dataset and often used as a baseline
* Dayhoff -> the earliest protein model, based on very limited data, and is mostly retained for historical comparison.

For the `-merit` I picked **AIC** which stands for Akaike Information Criterion. From some general reading and googling this is used most commonly as it 
* balances model fit versus complexity
* penalises overfitting but not as harsh as others
* good for predicitive accuracy and not worried about over complexity.

Other options here you can check the `IQTree3` documentation, but options could be **BIC**, **AICc**, or **DTL**.  


## Step 4. Adding Best Models to Partition File

The first step here is to extract the best model from each gene tree generated in the previous step. I have written `select_model.py` which autoimates this for you. Remember, you can run `select_model.py --help` to check what the inputs, outputs and other flags are doing.  

An important note here is that the `best_models.tsv` and partition file SCOs **MUST** match. The `select_model.py` has a suffix flag in it where you can select a suffix to match the partition file. Here, our partition file SCOs (file is the `.nex.nex` in `iqtree_all`) look like `SCOall.txt.OrthoGroupxx.trimmed` whereas if you ran the best_models with no suffix `SCOall.txt.OrthoGroupxx` would be included. Thus, including the `--suffix .trimmed` in the script below we match everything up. 

```
cd /scratch/alpine/beyo2625/species_tree_tut

## g90
bin/select_model.py \
-i iqtree_genetree_all \
-o iqtree_all/best_models.tsv \
-crit AIC
```
```
✅ Wrote 69 results to iqtree_all/best_models.tsv
```

Now we have that and everything matches, we need to merge it with the partition file so that the partition file is showing us 
* Location of each gene for each species
* Model to use for each gene for downstream steps.

I have again written a python script (`update_nexus.py`) which will do this for you. 

```
cd /scratch/alpine/beyo2625/species_tree_tut 

bin/update_nexus.py \
--input-nexus iqtree_all/alignments_concat.nex.nex \
--models iqtree_all/best_models.tsv \
--output-nexus iqtree_all/alignments_concat_withmods.nex
```
```
✅ Loaded 69 models from iqtree_all/best_models.tsv
✅ New NEXUS with models written to iqtree_all/alignments_concat_withmods.nex
```

Okikoke, and now we have that we are ready to move on to the next step. 

## Step 5. Collapsing Poorly Supported Branches

Now we collapse the poorly supported branches in all of the generated gene trees in  **Step 3**. 

This is pretty straightforward and you can just run a simple for loop to do this. 

```
cd /scratch/alpine/beyo2625/species_tree_tut
mamba activate newickutils_env

for tree in iqtree_genetree_all/SCO*/*.treefile; do
    filename=$(basename "$tree")
    echo "Collapsing support ≤10 in $filename"
    nw_ed "$tree" 'i & b<=10' o > collapsed_trees/"${filename%.treefile}.collapsed.tree"
done
```
```
....
Collapsing support ≤10 in SCOall.txt.OrthoGroup6.treefile
Collapsing support ≤10 in SCOall.txt.OrthoGroup7.treefile
Collapsing support ≤10 in SCOall.txt.OrthoGroup8.treefile
Collapsing support ≤10 in SCOall.txt.OrthoGroup9.treefile
```

Okay woinderful, so now we have our individual gene trees with low support collapsed. We can  now build the species tree using these. 

## Step 6. Building the Species Tree

The first step in generating our species tree, using `astral`, is to combine all the individual gene trees into one file. 

```
cd /scratch/alpine/beyo2625/species_tree_tut
cat collapsed_trees/*.collapsed.tree > astral_all/gene_trees_all_SCO.tree
```

And now we can run `astral` to generate the species tree. 

```
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --account=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --job-name=astral_all
#SBATCH --error=/scratch/alpine/beyo2625/species_tree_tut/astral_all/aster.err
#SBATCH --output=/scratch/alpine/beyo2625/species_tree_tut/astral_all/aster.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=2

module purge 
eval "$(conda shell.bash hook)"
conda activate aster_env

cd /scratch/alpine/beyo2625/species_tree_tut/astral_all

astral4 \
-i gene_trees_all_SCO.tree \
-o species_tree_all.nwk \
-u 3 \
-w 1
```

## Step 7. Quartet Support Analysis

And now finally we get to bring everything together and rereun `IQTree3` with 
* `astral` species tree (Step 6). 
* original concatenated genes for each species (Step 2).
* Updated partition file (Step 3 - 4).

```
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --qos=XXXXX
#SBATCH --partition=XXXXX
#SBATCH --account=XXXXX
#SBATCH --nodes=1
#SBATCH --mem=30G
#SBATCH --job-name=reestimate_tree
#SBATCH --error=/scratch/alpine/beyo2625/species_tree_tut/iqtree_all_reest/iqtree_reest.err
#SBATCH --output=/scratch/alpine/beyo2625/species_tree_tut/iqtree_all_reest/iqtree_reest.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=XXXXX
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=2

module purge 
eval "$(conda shell.bash hook)"
conda activate iqtree_env

cd /scratch/alpine/beyo2625/species_tree_tut/iqtree_all_reest
iqtree3 \
-s ../iqtree_all/alignments_concat.nex \
-p ../iqtree_all/alignments_concat_withmods.nex \
-te ../astral_all/species_tree_all.nwk \
--prefix sco_all_reest_root \
--seqtype AA \
-o sample11,sample59,sample66 \
--seed 73 \
-T AUTO \
-nstop 25

iqtree3 \
-s ../iqtree_all/alignments_concat.nex \
-p ../iqtree_all/alignments_concat_withmods.nex \
-te ../astral_all/species_tree_all.nwk \
--prefix sco_all_reest \
--seqtype AA \
--seed 73 \
-T AUTO \
-nstop 25
```

So here we rerun `IQTree3` twice. The first run is with rooting, the second without. The rooting can often fail (especially with SCO>xx% species) so you would just go for the second option with no rooting. Once this runs you are all done and ready to go into your favourite tree visualisation tool. 

Below I show you how to do this in `r` using `tidytree` and `treeio` as that is my preference. 

## Step 8. Visualising Tree with Quartets

So I do this in `r` and I will give the basic high level code to generate your tree with the quartet pie charts included. 

```{r library loading, include = F}
library(tidyverse)
library(ggtree)
library(treeio)
library(tidytree)
library(ggtreeExtra)
library(ggimage)
sessionInfo()
```

It is a little weird for the `tidytree::collapse` function. You need to do the below chunk to make it work. 

```{r making collapse work, include = F}
nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")
```

Okay lovely so first we read in out tree data that we made in **Step 7**. 

```{r reading in the iqtree astral results, include = F}
treeio::read.newick(file = "~/Desktop/iqtree_all_reest/sco_all_reest_root.treefile") -> astral_root
```

Then, you can merge it with your metadata. Remember, we renamed everything toi `sample01 sample02 ..... samplenn` and this is not very useful. Using the `analysis_lists/name_mapping.tsv` (generated in **Step 1**) you can generate a metadata table that matches the root labels (`sample01 sample02 ..... samplenn`) with species names, locations, treatments etc etc (the world is your oyster here). 

```{r importing metadata and annotating iqtree data, include = F}
## read in the metadata
read.csv(file = "~/Desktop/tree_metadata.csv") %>% ## metadata with the name_mapping.tsv
  mutate(across(where(is.character), ~na_if(., ""))) %>% 
  mutate(Genus_species_good = str_replace_all(Genus_species, "_", " ")) %>% 
  mutate(tiplabel = paste0(sample, " | ", Genus_species_good, " | ", Location)) -> tree_md_gp

## making a variable of just the tip sample names
astral_root %>% 
  as_tibble() %>% View()
  filter(str_detect(label, "^[A-Za-z]")) %>%
  dplyr::select(label) -> samples_all

## joining the tree data with the metadata using the sample names as a filter
## do this as you may have more stuff in metadata than in the tree. 
astral_root %>% 
  as_tibble() %>%
  full_join(tree_md_gp %>% 
              dplyr::filter(sample %in% samples_all$label), 
            join_by(label == sample)) %>%
  as.treedata() -> astral_root
```

Now, we need to extract the Q1, Q2, and Q3 data from the treefile. If you do `astral_root %>% as.tibble() %>% View()` you can see that if you scroll down there is **ALOT** of information saved. We only want the Qxx values but all that other info is there if you need it. 

```{r extracting the q values and making the pie data object, include = F}
astral_root %>%
  as_tibble() %>%
  mutate(
    q1 = str_extract(label, "(?<=q1=)[0-9eE\\.\\-]+") %>% as.numeric(),
    q2 = str_extract(label, "(?<=q2=)[0-9eE\\.\\-]+") %>% as.numeric(),
    q3 = str_extract(label, "(?<=q3=)[0-9eE\\.\\-]+") %>% as.numeric()
  ) %>%
  filter(!is.na(q1) & !is.na(q2) & !is.na(q3)) %>%
  select(node, q1, q2, q3) -> pie_data
```

Now we have the quartet data for the pies, you use the `nodepie` function to make the pie charts for the tree, and then can also assign pretty colours. 

```{r using nodepie to make for plotting, include = F}
pies <- nodepie(pie_data, 
                cols = 2:4, 
                alpha = 0.9)
cols <- c("q1" = "steelblue", "q2" = "orange", "q3" = "green")  # your chosen colors
pies <- lapply(pies, function(g) g + scale_fill_manual(values = cols))
```

Then, you are ready to plot yay. Note here I have stuck with `samplenn` as the tip labels as this data is not published yet (once it is I will add the image in). 

```{r viewing g90 no root tree, fig.width=10, fig.height=14, echo = F, results='hide'}
ggtree(astral_root, ladderize = T, size = 0.5) +
  geom_tiplab(aes(label = label, 
                  fontface = "italic"),
              size = 3,
              align = F) +
  geom_inset(pies, 
             width = 0.05, 
             height = 0.05) +
  theme(plot.margin = margin(1, 
                             7, 
                             1,
                             1, 
                             "cm")) +
  geom_rootedge(rootedge = 0.005) +
  coord_cartesian(clip = "off") +
  geom_treescale(
    fontsize = 3,
    x = 0.08,
    linesize = 1,
    offset = 1
  ) +  ## N143 Europe
  geom_cladelab(
    node = 143,
    label = "Europe",
    fontsize = 6,
    align = F,
    offset = 0.05,
    offset.text = 0.001,
    textcolor = 'black',
    barcolor = 'darkgrey',
    family = "Times"
  ) +
  geom_hilight(
    node = 143,
    fill = "steelblue2",
    type = "gradient",
    gradient.direction = 'rt',
    alpha = .3,
    to.bottom = T
  ) + ## N136 Asia
  geom_cladelab(
    node = 136,
    label = "Asia",
    fontsize = 6,
    align = F,
    offset = 0.055,
    offset.text = 0.001,
    textcolor = 'black',
    barcolor = 'darkgrey',
    family = "Times"
  ) +
  geom_hilight(
    node = 136,
    fill = "goldenrod2",
    type = "gradient",
    gradient.direction = 'rt',
    alpha = .3,
    to.bottom = T
  ) + ## N133 Europe
  geom_cladelab(
    node = 133,
    label = "Europe",
    fontsize = 6,
    align = F,
    offset = 0.05,
    offset.text = 0.001,
    textcolor = 'black',
    barcolor = 'darkgrey',
    family = "Times"
  ) +
  geom_hilight(
    node = 133,
    fill = "steelblue2",
    type = "gradient",
    gradient.direction = 'rt',
    alpha = .3,
    to.bottom = T
  ) + ## N136 Europe
  geom_cladelab(
    node = 156,
    label = "Europe",
    fontsize = 6,
    align = F,
    offset = 0.05,
    offset.text = 0.001,
    textcolor = 'black',
    barcolor = 'darkgrey',
    family = "Times"
  ) +
  geom_hilight(
    node = 156,
    fill = "steelblue2",
    type = "gradient",
    gradient.direction = 'rt',
    alpha = .3,
    to.bottom = T
  ) + ## n156 Europe
  geom_cladelab(
    node = 157,
    label = "Europe",
    fontsize = 6,
    align = F,
    offset = 0.05,
    offset.text = 0.001,
    textcolor = 'black',
    barcolor = 'darkgrey',
    family = "Times"
  ) +
  geom_hilight(
    node = 157,
    fill = "steelblue2",
    type = "gradient",
    gradient.direction = 'rt',
    alpha = .3,
    to.bottom = T
  ) + ## N107 North America
  geom_cladelab(
    node = 107,
    label = "North America",
    fontsize = 6,
    align = F,
    offset = 0.06,
    offset.text = 0.001,
    textcolor = 'black',
    barcolor = 'darkgrey',
    family = "Times"
  ) +
  geom_hilight(
    node = 107,
    fill = "forestgreen",
    type = "gradient",
    gradient.direction = 'rt',
    alpha = .3,
    to.bottom = T
  ) + ## N85 North Ameerica
  geom_cladelab(
    node = 85,
    label = "North America",
    fontsize = 6,
    align = F,
    offset = 0.06,
    offset.text = 0.001,
    textcolor = 'black',
    barcolor = 'darkgrey',
    family = "Times"
  ) +
  geom_hilight(
    node = 85,
    fill = "forestgreen",
    type = "gradient",
    gradient.direction = 'rt',
    alpha = .3,
    to.bottom = T
  ) -> iqtree_all
iqtree_all
```

PUT IMAGE IN ONCE PAPER ACCEPTED 

You can see from this that we have 
* pies at each node
* nice tip labels with sample number, species, country
* higher level continent colouring.

While we have done alot of this in `r`, this is now where I would take it into Adobe Illustrator to make the final scale edits and moving of the pies. You can do this in `r` but it is not that easy. 

If exporting to Illustrator, use the following pdf export as some of the `treeio` options dont work with the normal pdf export ☺️. 

```{r saving as high quality pdf g90 astral tree, include = F}
cairo_pdf("~/Desktop/iqtree_all_root.pdf", 
    width = 10, 
    height = 15)
iqtree_all
dev.off()
```


