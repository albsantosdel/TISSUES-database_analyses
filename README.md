Comprehensive comparison of large-scale tissue expression datasets ![alt text](http://jensenlab.org/images/tissues_icon.png "TISSUES database")
==============

Code to reproduce the fold enrichment analyses and figures from the article
--------------
**Project structure**
- analyses/
 - README.md `Markdown file`
 - makefile  main `script`
 - generate\_files.pl `script` to generate all the necessary files. This script uses the original datasets files (with no filter for unconfident gene-tissue associations) 
 - analyses.R `script` orchestrating the generation of the figures
 - R/ `All the scripts necessary to reproduce the figures`
   - summary\_figure\_drawing.R `Initial figure with the tissues and number of genes/proteins provided by each dataset`
    - datasets\_expression\_breadth\_analyses.R `Expression breadth distribution analysis`
     - fold\_enrichment\_analysis\_per\_dataset.R `Generates the fold-enrichment plots for each dataset`
      - fold\_enrichment\_score\_calibration\_analysis.R `Score calibration figure`
      - venn\_diagram_analyses.R `Creates all the venn diagrams and calculates the p-values of the overalps when comparing common proteins and tissues`
      - external\_function.R `External function to switch coordinates in facet grid (Function from` [`stackoverflow`](http://stackoverflow.com/questions/6625691/is-it-possible-to-switch-the-side-of-y-axis-breaks-and-labels-on-a-faceted-plot)`)`
 - data/ 
  - datasests/ `Datasets files need to be stored here`
  - dictionary/ `Contains the files necessary to perform the tissue backtracking bassed on BRENDA Ontology`
   - labels.tsv `The bto terms corresponding to the 21 tissues of interest: tissues_code  tissue_name  BTO`
   - bto\_entities.tsv `mapping of bto terms to internal identifiers: internal_code  tissues_code  BTO`
   - bto\_groups .tsv `the parent-children relationships used to do the backtracking: internal_code  parent_internal_code`
 - figures/ `Folder where all the figures generated are stored`

**Run the analyses**

1. Download the project
2. Download the datasets from FigShare and store them into data/datasets folder: `[`datasets`](http://figshare.com/s/cb788d0ef4bd11e4b5ea06ec4b8d1f61)`
3. Modify the analyses.R:
  - Change the setwd(dir=""), to the directory where you downloaded the project, i.e. `setwd(dir="~/albsantosdel/reproduce_analyses")` 
4. Execute the makefile script from the command line:
  `> make`
5. All the files will be generated in the data folder
6. All the figures from the analyses will be created in the 'figures/' folder
