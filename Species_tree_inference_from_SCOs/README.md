# Species Tree Inference from Single Copy Orthologs

In this tutorial we’ll build a species tree from single-copy orthologs. These genes are ideal because they avoid problems with duplicates and give us a consistent set of markers across species.

We’ll use two complementary approaches:

Concatenation – combine all alignments into one big dataset and infer a tree with maximum likelihood.

Coalescent methods (e.g. ASTRAL) – first build gene trees, then summarize them into a species tree that accounts for gene tree conflict.

Finally, we’ll look at quartet support values (q1, q2, q3), which tell us how much the gene trees agree with each branch of the species tree versus supporting alternative histories. This gives a clearer picture of which parts of the tree are strongly supported and where there’s conflict.
