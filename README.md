Comprehensive comparison of large-scale tissue expression datasets ![alt text](http://jensenlab.org/images/tissues_icon.png "TISSUES database")
==============

Code to reproduce the fold enrichment analyses and figures from the article
--------------
**Project structure**
- README.md --> *Markdown file*
- makefile  -->	*main script*
- generate\_files.pl --> *script to generate all the necessary files. This script uses the original datasets files (with no filter for unconfident gene-tissue associations)*
- analyses.R --> *script orchestrating the generation of the figures*
- R/ --> *All the scripts necessary to reproduce the figures*
  - summary\_figure\_drawing.R --> *Initial figure with the tissues and number of genes/proteins provided by each dataset*
   - datasets\_expression\_breadth\_analyses.R --> *Expression breadth distribution analysis*
    - fold\_enrichment\_analysis\_per\_dataset.R --> *Generates the fold-enrichment plots for each dataset*
     - fold\_enrichment\_score\_calibration\_analysis.R --> *Score calibration figure*
     - venn\_diagram_analyses.R --> *Creates all the venn diagrams and calculates the p-values of the overalps when comparing common proteins and tissues*
     - external\_function.R --> *External function to switch coordinates in facet grid (Function from* [stackoverflow](http://stackoverflow.com/questions/6625691/is-it-possible-to-switch-the-side-of-y-axis-breaks-and-labels-on-a-faceted-plot))
- data/ 
 - datasests/ --> *Datasets files need to be stored here*
 - dictionary/ --> *Contains the files necessary to perform the tissue backtracking bassed on BRENDA Ontology*
  - labels.tsv --> *The bto terms corresponding to the 21 tissues of interest: tissues\_code  tissue\_name  BTO*
  - bto\_entities.tsv --> *mapping of bto terms to internal identifiers: internal\_code  tissues\_code  BTO*
  - bto\_groups .tsv --> *the parent-children relationships used to do the backtracking: internal\_code  parent\_internal\_code*
- figures/ --> *Folder where all the figures generated are stored*

**Run the analyses**

1. Download the project
2. Make sure you have a default CRAN repository set, for example by
   putting the following into your `~/.Rprofile`:

     ```r
     local({r <- getOption("repos")
            r["CRAN"] <- "http://cran.us.r-project.org"
            options(repos=r)})
     ```

3. Execute the makefile script from the command line:
  `> make`
4. All the files will be generated in the data folder
5. All the figures from the analyses will be created in the figures/ folder

**Generated files**

./data/
- Fold enrichment analyses result files: *DATASET*\_*GOLDSTANDARD*\_fold\_enrichment\_analysis.tsv (goldstandards: UniProt-KB and mRNA reference set)
- Expression breadth analyses result files: *CUTOFF*\_cutoff\_expression\_breadth.tsv (cutoffs: low, medium, high)
- Consistency analyses result files: *CUTOFF*\_consistency\_analysis.tsv (cutoffs: low, medium, high)
- P-value results for venn diagrams where common proteins and tissues are taken into account: pairwise\_pvalue\_results.txt

./figures/
- Fold enrichment plots: *DATASET*\_fold\_enrichment.png
- Venn diagrams:
 - *COMPARISON*\_venn\_diagram.png
 - all\_datasets\_comparison\_venn\_diagram\_*COMMON/nonCOMMON*.png
- Expression breadth plots: *DATASET*\_datasets\_prots\_num\_tissues.png
- Score calibration plot: datasets_score_calibration.png

**Requirements**

- Perl
- R
- curl
