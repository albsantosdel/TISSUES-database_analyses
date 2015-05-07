.PHONY : all
.DEFAULT: all
PACKAGE      = package
VERSION      = ` date "+%Y.%m%d%" `
RELEASE_DIR  = ..
RELEASE_FILE = $(PACKAGE)-$(VERSION)
SILENT = yes

all: generate_files generate_figures

#Generates all the files necessary for the analyses
generate_files:
	perl generate_files.pl;
#Generates all the figures from the analyses described in the article: Comprehensive comparison of large-scale tissue expression datasets
generate_figures:
	Rscript analyses.R




