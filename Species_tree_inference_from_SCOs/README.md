# Species Tree Inference from Single Copy Orthologs

In this tutorial we’ll build a species tree from single-copy orthologs. These genes are ideal because they avoid problems with duplicates and give us a consistent set of markers across species.

We’ll use two complementary approaches:

Concatenation – combine all alignments into one big dataset and infer a tree with maximum likelihood.

Coalescent methods (e.g. ASTRAL) – first build gene trees, then summarize them into a species tree that accounts for gene tree conflict.

Finally, we’ll look at quartet support values (q1, q2, q3), which tell us how much the gene trees agree with each branch of the species tree versus supporting alternative histories. This gives a clearer picture of which parts of the tree are strongly supported and where there’s conflict.


## Step 1. Generating Cleaned and Aligned Single Copy Orthologs  

This step here is broadly the same as Steps xx - xx in the Phylogenomic tutorial. As such, I do not go into detail here but you can refer to other tutorial if you want more information. 



## Step 2. Generating Concatenated Alignment and Partition File 


## Step 3. Identyfying the Best Model for each Gene Tree and Generating Gene Trees


## Step 4. Adding Best Models to Partition File


## Step 5. Collapsing Poorly Supported Branches


## Step 6. Building the Species Tree


## Step 7. Quartet Support Analysis


## Step 8. Visualising Tree with Quartets
