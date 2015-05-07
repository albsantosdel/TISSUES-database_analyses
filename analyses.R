setwd('./')
#Install required libraries
if("RColorBrewer" %in% rownames(installed.packages()) == FALSE) {install.packages("RColorBrewer")}
if("ggplot2" %in% rownames(installed.packages()) == FALSE) {install.packages("ggplot2")}
if("scales" %in% rownames(installed.packages()) == FALSE) {install.packages("scales")}
if("grid" %in% rownames(installed.packages()) == FALSE) {install.packages("grid")}
if("gtable" %in% rownames(installed.packages()) == FALSE) {install.packages("gtable")}
if("gridExtra" %in% rownames(installed.packages()) == FALSE) {install.packages("gridExtra")}
if("shape" %in% rownames(installed.packages()) == FALSE) {install.packages("shape")}
if("VennDiagram" %in% rownames(installed.packages()) == FALSE) {install.packages("VennDiagram")}
if("gridBase" %in% rownames(installed.packages()) == FALSE) {install.packages("gridBase")}
if("lattice" %in% rownames(installed.packages()) == FALSE) {install.packages("lattice")}

#Load libraries
require(RColorBrewer)
require(ggplot2)
require(scales)
require(grid)
require(gtable)
require(gridExtra)
require(shape)
require(VennDiagram)
require(gridBase)
require(lattice)

#External function
source("R/external_function.R")

####################################
#         Summary figure          #
###################################
source("R/summary_figure_drawing.R")

####################################
# Expression breadth distribution #
###################################
source("../R/datasets_expression_breadth_analyses.R")

####################################
#         Fold-enrichment         #
###################################
source("../R/fold_enrichment_analysis_per_dataset.R")
source("../R/fold_enrichment_score_calibration_analysis.R")

####################################
#           Venn-diagrams         #
###################################
source("../R/venn_diagram_analyses.R")
