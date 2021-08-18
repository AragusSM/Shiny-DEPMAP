To Whomever this concerns:

---------------------------------------------------------------------------------
This is the source code file that contains source code and data as well as a few
generated plots in order to build a shiny app in the R program. This app explores
relationships between drug, genes, and their effects on cancer through mass data
analysis. The goal is to understand key genes, drugs, or any combination of both
in predicting cancer treatment efficacy. The methodology involves looking at pairs of
drugs and genes across a few hundred shared cancer cell lines and drawing 
correlations between them.

Included in the current version are:
correlation queries, graph displays, as well as network displays for visualizing 
otherwise difficult to comprehend connections between imperative pairs of genes and drugs.
We hope to further develop this application and widen the scope of our visualization
in the hopes of elucidating hitherto undiscovered biological and medical implications.
---------------------------------------------------------------------------------
In order to run the app, simply run this folder in the shiny server after altering
the filepaths in "app.r". Warning: reading the csv files and generating tables/plots
may take a LONG time. In addition, make sure to have enough storage and memory on
your computer.
---------------------------------------------------------------------------------
This folder contains the comprehensive files of the project, 
although there are a couple of not-yet-generated tables, those being the spearman
correlations for the 2nd and third quadrant of the gene codependency table.
The gene codependency data is originally supposed to be a ~17000 by 17000 matrix.
However, due to memory issues, the data is thus spliced into 4 quadrants.
The third and first quadrant are redundant quadrants, so the files titled "genecodep_pearson.csv"
and "genecodep_spearman.csv" are the 1st/3rd quadrant intersection of roughly 8500*8500 genes,
or the first half of the large matrix by the second half. The files titled "genecodep_pearson1.csv"
and "genecodep_pearson2.csv" are the 2nd quadrant and the 4th quadrant, respectively.
Thus the missing data would be the spearman codeps for the 2nd and the 4th quadrant. 

However, rest assured, for the code for generating the tables are positioned with comments 
inside "app.r". Just be prepared to leave your computer on for the entire night.
----------------------------------------------------------------------------------
With that out of the way, the rest of the files are quite self explanatory,
well not quite, that's what comments are for! Please read the comments in each of
the .r files as they will help you obtain a comprehensive view of this project's inner
workings.

Aside from the gene codependency correlations mentioned earlier, the other correlation
data consists of drug codependency titled "drugcodep_pearson.csv" and "drugdep_spearman.csv" and
drug-gene correlation, titled "gene_drug_dependencies_chronos.csv" and "gene_drug_dep_CHRONOS_spearman.csv".
The other data set, "gene_drug_dependencies.csv" uses the CERES method data set for gene dep
and is an alternative to "gene_drug_dependencies_chronos.csv".
-------------------------------------------------------------------------------------
The source data sets are contained in the "PRISM" and "crispr_dep_map" folders, with the drug response
data in PRISM and gene dependency in crispr_dep_map. "The plots and images" folder are a few 
generated images that display high correlation drug and gene pairs as well images that analyze the 
difference between target and sample pairs.
-------------------------------------------------------------------------------------
The main app is "app.r" and the function files implemented inside are "correlation_plots.R",
"correlation_query.R", "DepCorrGenerator.R", and "NameChanger.R" The functions are commented
in detail so this README will not delve into them, simply open in R and look for yourself.

The "Drug_gene_dependencies.R" file is a preliminary file that contains some of the code in app.R
but also contains the code for generating samples and valuable targets above a specific threshold.
The plots in plots and images are generated using this file, so it may prove to be of use.
--------------------------------------------------------------------------------------
Lastly, the "valuable_targets" files are tables that are comprised of gene drug pairs from
the target treatment table inside "app.r". the Target correlations data are not included in the folder
because they can be efficiently generated inside R using the raw data. the files with spearman
are spearman valuable targets above the threshold (5 and 10 p-value), and the files without 
"spearman" are pearson valuable targets with the same p values.
----------------------------------------------------------------------------------------
Note: The Github only contains the code and not actually the raw data. in order to generate data
please go to  https://depmap.org/portal/ and download the data files.

Welp, that is all that is contained inside the current version of the app, explore to your
heart's content. please feel free to make any and all changes to this folder as you like,
including deleting, renaming, changing code, adding new data, etc. 
Thank you for reading!
