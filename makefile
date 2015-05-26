.PHONY : all
.DEFAULT: all
PACKAGE      = package
VERSION      = ` date "+%Y.%m%d%" `
RELEASE_DIR  = ..
RELEASE_FILE = $(PACKAGE)-$(VERSION)
SILENT = yes
DATA_DIR     = data/datasets
DATA_FILES   = 2062920/hpa_rna.tsv 2062917/hpm.tsv 2062918/unigene.tsv 2062919/hpa_ihc.tsv 2062922/text_mining.tsv 2062921/uniprot.tsv 2062923/exon.tsv 2062924/gnf.tsv 2071678/rna_seq.tsv

all: generate_files generate_figures

$(DATA_FILES): 
	curl -L -o $(DATA_DIR)/$(@F) --time-cond $(DATA_DIR)/$(@F) http://files.figshare.com/$@

#Generates all the files necessary for the analyses
generate_files: $(DATA_FILES)
	perl generate_files.pl;
#Generates all the figures from the analyses described in the article: Comprehensive comparison of large-scale tissue expression datasets
generate_figures:
	Rscript analyses.R




