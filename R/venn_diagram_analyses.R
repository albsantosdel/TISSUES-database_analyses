splitvalues<-function(row){
  x<-as.character(row)
  substr(x,1,regexpr("-",x)-1)
}

pair_wise_venn_diagram <- function(l1, l2, val1, val2, val3, val4, val5, c1, c2,title){  
  #initialize plot
  emptyplot(c(-8, 8),main= title)
  
  #plot circles
  plotellipse(mid = c(2, 0),  rx = 4, ry = 4, arrow = FALSE, lwd=2.5, lc=c2)
  plotellipse(mid = c(-2, 0), rx = 4, ry = 4, arrow = FALSE, lwd=2.5, lc=c1)
  plotellipse(mid = c(-6, 0), rx = 4, ry = 4, arrow = FALSE, from=pi/2, to=1.5*pi, lwd=2.5, lty=2, lc=c1)
  plotellipse(mid = c(6, 0),  rx = 4, ry = 4, arrow = FALSE, from=1.5*pi, to=pi/2, lwd=2.5, lty=2, lc=c2)
  
  #add connecting lines
  lines(c(2, 6), c(4, 4), lwd=2.5, lty=2, col=c2)
  lines(c(2, 6), c(-4, -4), lwd=2.5, lty=2, col=c2)
  lines(c(-2, -6), c(4, 4), lwd=2.5, lty=2, col=c1)
  lines(c(-2, -6), c(-4, -4), lwd=2.5, lty=2, col=c1)
  
  #plot numbers
  text(-6.5, 0, val1, adj=c(1.2, 0.5),cex=1)
  text(-2.5, 0, val2, adj=c(1.2, 0.5),cex=1)
  text(0, 0, val3, adj=c(0.5, 0.5),cex=1)
  text(2.5, 0, val4, adj=c(-0.2, 0.5),cex=1)
  text(6.5, 0, val5, adj=c(-0.2, 0.5),cex=1)
  
  #labels
  text(-0.75, -4.5, l1, adj=c(1.2, 1), col=c1,cex=1)
  text(4, -4.5, l2, adj=c(-0.2, 1), col=c2,cex=1)
  
}

get_data<-function(cutoff){
  in_file<-paste(cutoff,"consistency_analysis.tsv",sep="_")
  
  protein_tissues<-read.table(in_file,sep="\t",header=FALSE)
  names(protein_tissues)<-c("dataset","score","tissue","proteins")
  
  return(protein_tissues)
}

filter_dataset_tissue<-function(data,dataset_name,tissues_list,not_in=FALSE,columns){
  if(not_in & !is.null(dataset_name)){
    filtered_data<- data[data$dataset == dataset_name & !(data$tissue %in% tissues_list),columns]
  }
  else if(not_in & is.null(dataset_name)){
    filtered_data<- data[!(data$tissue %in% tissues_list),columns]
  }
  else if(!not_in & !is.null(dataset_name)){
    filtered_data<- data[data$dataset == dataset_name & data$tissue %in% tissues_list,columns] 
  }
  else{
    filtered_data<- data[data$tissue %in% tissues_list,columns]
  }
  
  return(filtered_data)
}


plot_venns<-function(list,cols,filename,title,areas,text_color=NULL){
  size_title<-length(list)-1
  if(is.null(text_color)){
    text_color=cols
  }
  venn.plot <-venn.diagram(x = list,
                           filename= filename,
                           main = title,
                           sub = NULL,
                           main.fontfamily="Helvetica",
                           main.cex = size_title,
                           scaled = F, 
                           lty = "solid",
                           lwd = 2.5,
                           col =  cols,
                           cex = 1,
                           fontface=areas,
                           fontfamily="Helvetica",
                           cat.col=text_color,
                           cat.pos=0,
                           cat.fontfamily="Helvetica")  
  return(venn.plot)
}


print_venn_to_pdf<-function(data_list,colors, title,output_file,height,width,filename=NULL,text_cols=NULL){
  png(file = output_file,
      width = width, height = height)
  areas<-rep("plain",15)
  plot<-plot_venns(data_list, colors, filename, title,areas,text_color=text_cols)
  grid.draw(plot)
  dev.off()
}

calc_pairwise_pvalue <- function(labelA, A, labelB, B, population, cutoff, extra_text = "", append = FALSE){
  common <- length(intersect(x = A, y = B))
  singleA <- length(unique(A)) - common
  singleB <- length(unique(B)) - common
  
  notInAB <- nrow(population) - length(intersect(x = union(x = A, y = B), y = population))
  
  nolabelA <- paste("notIn", labelA, sep ='')
  nolabelB <- paste("notIn", labelB, sep='')
  
  contigency.matrix <- matrix(c(common, singleA, singleB, notInAB),
           nrow = 2, dimnames = list(c(labelB, nolabelB), c(labelA, nolabelA)))
  
  print(contigency.matrix)
  
  p <- fisher.test(contigency.matrix, alternative="greater", simulate.p.value=TRUE)
  
  text <- paste("P-value calculated for ", labelA, " vs ", labelB, " comparison at cutoff ", cutoff," (",extra_text,")","\n", sep='')
  text<-paste(text,"method = ", p$method, "\n", sep='')
  text<-paste(text,"p-value = ", p$p.value, "\n", sep='')
  text<-paste(text,"Confidence interval = ", p$conf.int, "\n", sep='')
  
  write(text, file = "pairwise_pvalue_results.txt", append = append)
  #write.table(contigency.matrix, file = "pvalue_results.txt", sep = "\t", append = TRUE)
  print(p)
}

get_venn_diagrams<-function(protein_tissues,cutoff,uniprot=F,col, append = apppend){
  
  uncommon_tissues<-c("bone")
  
  exon_data<- filter_dataset_tissue(protein_tissues,"exon",c(),not_in=TRUE,c(2,3,4))
  gnf_data<- filter_dataset_tissue(protein_tissues,"gnf",c(),not_in=TRUE,c(2,3,4))
  unigene_data<- filter_dataset_tissue(protein_tissues,"unigene",c(),not_in=TRUE,c(2,3,4))
  rna_data<- filter_dataset_tissue(protein_tissues,"rna",c(),not_in=TRUE,c(2,3,4))
  hparna_data<- filter_dataset_tissue(protein_tissues,"hpa_rna",c(),not_in=TRUE,c(2,3,4))
  
  #Protein data sets
  common_tissues<-c("heart", "liver", "intestine", "nervous system", "kidney","lung","pancreas","adrenal gland","urine","gall bladder")
  #all_data<-filter_dataset_tissue(protein_tissues,NULL,common_tissues, not_in=FALSE,c(3,4))
  #population <- unique(paste(all_data$proteins,all_data$tissue,sep="-"))
  
  
  hpa_data<- filter_dataset_tissue(protein_tissues,"hpa",common_tissues,not_in=FALSE,c(3,4))
  hpm_data<- filter_dataset_tissue(protein_tissues,"hpm",common_tissues,not_in=FALSE,c(3,4))
    
  prot_common_proteins<-intersect(x=hpm_data$proteins,y=hpa_data$proteins)
  
  hpm_data$pairs<-paste(hpm_data$proteins,hpm_data$tissue,sep="-")
  hpm_total_pairs<-length(hpm_data$pairs)
  
  hpa_data$pairs<-paste(hpa_data$proteins,hpa_data$tissue,sep="-")
  hpa_total_pairs<-length(hpa_data$pairs)
    
  hpm_data<-hpm_data[hpm_data$proteins %in% prot_common_proteins,]
  
  hpa_data<-hpa_data[hpa_data$proteins %in% prot_common_proteins,]
  
  
  hpm_unique_pairs<-hpm_total_pairs-length(hpm_data$pairs)
  hpa_unique_pairs<-hpa_total_pairs-length(hpa_data$pairs)
  common_pairs<-length(intersect(hpm_data$pairs,hpa_data$pairs))
  hpm_uncommon<-length(hpm_data$pairs)-common_pairs    
  hpa_uncommon<-length(hpa_data$pairs)-common_pairs
  
  #pairwise p-value calculation for hpa vs hpm
  population <- expand.grid(prot_common_proteins, common_tissues)
    
  pvalue <-calc_pairwise_pvalue("HPA", hpa_data$pairs, "HPM", hpm_data$pairs, population, cutoff, "common proteins", append = append)
  
  prot_title<-paste("../figures/",cutoff,"_","hpm_hpa_venn_diagram.png",sep="")
  png(file = prot_title,
      width = 550, height = 450)
  title<-""
  pair_wise_venn_diagram( "HPA","HPM",hpa_unique_pairs , hpa_uncommon, common_pairs, hpm_uncommon, hpm_unique_pairs,col[6],col[9],title)
  dev.off()
  
  if(cutoff == "medium"){
    hpa_data<- filter_dataset_tissue(protein_tissues,"hpa",c(),not_in=TRUE,c(2,3,4))
    hpm_data<- filter_dataset_tissue(protein_tissues,"hpm",c(),not_in=TRUE,c(2,3,4))
    
    tm_data<- filter_dataset_tissue(protein_tissues,"tm",c(),not_in=TRUE,c(3,4))
    uniprot_data<-filter_dataset_tissue(protein_tissues,"uniprot",c(),not_in=TRUE,c(3,4))
    
    mRNA_data<-rbind(exon_data[exon_data$score>=250,],gnf_data[gnf_data$score>=250,],hparna_data[hparna_data$score>=20,],unigene_data[unigene_data$score>=20,],rna_data[rna_data$score>=5,])
    prot_data<-rbind(hpa_data[hpa_data$score>=10.45,],hpm_data[hpm_data$score>=10,])
    
    prot_data$pairs<-paste(prot_data$proteins,prot_data$tissue,sep="-")
    tm_data$pairs<-paste(tm_data$proteins,tm_data$tissue,sep="-")
    uniprot_data$pairs<-paste(uniprot_data$proteins,uniprot_data$tissue,sep="-")
    mRNA_data$pairs<-paste(mRNA_data$proteins,mRNA_data$tissue,sep="-")
    
    all_list<- list("transcriptomics set"= mRNA_data$pairs,"proteomics set"=prot_data$pairs," "=uniprot_data$pairs," "= tm_data$pairs)
    
    colors<-c("#b6dba4",col[7],col[8],"#99dbf8")
    text_colors<-c("#b6dba4","#99dbf8",col[7],col[8])
    output_file<-"../figures/all_datasets_comparison_venn_diagram_all_proteins_and_tissues.png"
    title<-""
    
    print_venn_to_pdf(all_list,colors, title,output_file,450,550,filename=NULL,text_colors)

    exon_data<-filter_dataset_tissue(exon_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2,3))
    gnf_data<-filter_dataset_tissue(gnf_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2,3))
    unigene_data<-filter_dataset_tissue(unigene_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2,3))
    rna_data<-filter_dataset_tissue(rna_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2,3))
    hparna_data<-filter_dataset_tissue(hparna_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2,3))
    hpa_data<- filter_dataset_tissue(hpa_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2,3))
    hpm_data<- filter_dataset_tissue(hpm_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2,3))
    
    tm_data<- filter_dataset_tissue(tm_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2))
    uniprot_data<-filter_dataset_tissue(uniprot_data,NULL,uncommon_tissues,not_in=TRUE,c(1,2))
    
    mRNA_data<-rbind(exon_data[exon_data$score>=250,],gnf_data[gnf_data$score>=250,],hparna_data[hparna_data$score>=20,],unigene_data[unigene_data$score>=20,],rna_data[rna_data$score>=5,])
    prot_data<-rbind(hpa_data[hpa_data$score>=10.45,],hpm_data[hpm_data$score>=10,])
    
    #common proteins
    common_proteins_med<-intersect(x=intersect(x=intersect(x=tm_data$proteins,y=uniprot_data$proteins),y=prot_data$proteins),y=mRNA_data$proteins)
    #common tissues
    common_tissues <- intersect(x=intersect(x=intersect(x=tm_data$tissue,y=uniprot_data$tissue),y=prot_data$tissue),y=mRNA_data$tissue)
          
    prot_data<-prot_data[prot_data$proteins %in% common_proteins_med,]
    prot_data$pairs<-paste(prot_data$proteins,prot_data$tissue,sep="-")
    tm_data<-tm_data[tm_data$proteins %in% common_proteins_med,]
    tm_data$pairs<-paste(tm_data$proteins,tm_data$tissue,sep="-")
    uniprot_data<-uniprot_data[uniprot_data$proteins %in% common_proteins_med,]
    uniprot_data$pairs<-paste(uniprot_data$proteins,uniprot_data$tissue,sep="-")
    mRNA_data<-mRNA_data[mRNA_data$proteins %in% common_proteins_med,]
    mRNA_data$pairs<-paste(mRNA_data$proteins,mRNA_data$tissue,sep="-")
    
    all_data<-filter_dataset_tissue(protein_tissues,NULL,common_tissues, not_in=FALSE,c(3,4))
    all_data<-all_data[all_data$proteins %in% common_proteins_med,]
    
    population <- expand.grid(common_proteins_med, common_tissues)
    
    pvalue <-calc_pairwise_pvalue("transcriptomics set", mRNA_data$pairs, "proteomics set", prot_data$pairs, population, cutoff, "common_proteins",append = TRUE)
    pvalue <-calc_pairwise_pvalue("transcriptomics set", mRNA_data$pairs, "uniprot", uniprot_data$pairs, population, cutoff, "common_proteins",append = TRUE)
    pvalue <-calc_pairwise_pvalue("transcriptomics set", mRNA_data$pairs, "text mining", tm_data$pairs, population, cutoff, "common_proteins",append = TRUE)
    pvalue <-calc_pairwise_pvalue("proteomics set", prot_data$pairs, "uniprot", uniprot_data$pairs, population, cutoff, "common_proteins",append = TRUE)
    pvalue <-calc_pairwise_pvalue("proteomics set", prot_data$pairs, "text mining", tm_data$pairs, population, cutoff, "common_proteins",append = TRUE)
    pvalue <-calc_pairwise_pvalue("uniprot", uniprot_data$pairs, "text mining", tm_data$pairs, population, cutoff, "common_proteins",append = TRUE)
    
    all_list<- list("transcriptomics set"= mRNA_data$pairs,"proteomics set"=prot_data$pairs," "=uniprot_data$pairs," "= tm_data$pairs)
    
    inAll <- length(intersect(x= mRNA_data$pairs,y= intersect(x = prot_data$pairs, y = intersect(x = uniprot_data$pairs, y = tm_data$pairs))))
    inSingles <- c(length(unique(mRNA_data$pairs)),length(unique(prot_data$pairs)),length(unique(uniprot_data$pairs)),length(unique(tm_data$pairs)))
    pvalue_title <- paste("Pvalues Comparing all datasets (only common proteins) at cutoff ", cutoff, sep='')
    
    ##create supplementary files
    common_prot_mRNA <- intersect(x = mRNA_data$pairs, y=prot_data$pairs)
    common_all <- intersect(x = mRNA_data$pairs, y = intersect( x = prot_data$pairs, y = intersect(x = uniprot_data$pairs, y = tm_data$pairs)))
    unique <- c(setdiff(x = mRNA_data$pairs, y = common_all), setdiff(x = prot_data$pairs, y = common_all),setdiff(x = uniprot_data$pairs, y = common_all),setdiff(x = tm_data$pairs, y = common_all))
    
    
    write.table(common_prot_mRNA, file="Common_gene_tissue_associations_transcriptomic_protomics.tsv", row.names=FALSE, quote = FALSE, col.names = c("Gene-tissue associations"))
    write.table(common_all, file="Common_all_datasets.tsv", row.names=FALSE, quote = FALSE, col.names = c("Gene-tissue associations"))
    write.table(unique, file="gene_tissue_associations_unique.tsv", row.names=FALSE, quote = FALSE, col.names = c("Gene-tissue associations"))
    
    
    
    colors<-c("#b6dba4",col[7],col[8],"#99dbf8")
    output_file<-"../figures/all_datasets_comparison_venn_diagram_common_proteins_and_tissues.png"
    title<-""
    
    print_venn_to_pdf(all_list,colors, title,output_file,450,550,filename=NULL,text_colors)
  }
  common_tissues<-c("heart", "liver", "intestine", "nervous system", "kidney")
  exon_data<-filter_dataset_tissue(exon_data,NULL,common_tissues,not_in=FALSE,)
  gnf_data<-filter_dataset_tissue(gnf_data,NULL,common_tissues,not_in=FALSE,)
  unigene_data<-filter_dataset_tissue(unigene_data,NULL,common_tissues,not_in=FALSE,)
  rna_data<-filter_dataset_tissue(rna_data,NULL,common_tissues,not_in=FALSE,)
  
  hparna_data<-filter_dataset_tissue(hparna_data,NULL,common_tissues,not_in=FALSE,)
  
  #filter only common proteins to all 
  common_proteins<-intersect(intersect(x=intersect(x=intersect(x=exon_data$proteins,y=gnf_data$proteins),y=unigene_data$proteins),y=rna_data$proteins),y=hparna_data$proteins)
  #write.table(file="~/Desktop/high_common_proteins.tsv",common_proteins) 
  
  exon_data<-exon_data[exon_data$proteins %in% common_proteins,]
  exon_data$pairs<-paste(exon_data$proteins,exon_data$tissue,sep="-")
  
  gnf_data<-gnf_data[gnf_data$proteins %in% common_proteins,]
  gnf_data$pairs<-paste(gnf_data$proteins,gnf_data$tissue,sep="-")
  
  unigene_data<-unigene_data[unigene_data$proteins %in% common_proteins,]
  unigene_data$pairs<-paste(unigene_data$proteins,unigene_data$tissue,sep="-")
  rna_data<-rna_data[rna_data$proteins %in% common_proteins,]
  rna_data$pairs<-paste(rna_data$proteins,rna_data$tissue,sep="-")
  hparna_data<-hparna_data[hparna_data$proteins %in% common_proteins,]
  hparna_data$pairs<-paste(hparna_data$proteins,hparna_data$tissue,sep="-")
  
  
  common_pairs<-intersect(intersect(x=intersect(x=intersect(x=exon_data$pairs,y=gnf_data$pairs),y=unigene_data$pairs),y=rna_data$pairs),y=hparna_data$pairs)
  
  population <- expand.grid(common_proteins, common_tissues)
  
  pvalue <-calc_pairwise_pvalue("ExonArray", exon_data$pairs, "GNF", gnf_data$pairs, population, cutoff,"common_proteins", append = TRUE)
  pvalue <-calc_pairwise_pvalue("ExonArray", exon_data$pairs, "RNA-seq", rna_data$pairs, population, cutoff,"common_proteins", append = TRUE)
  pvalue <-calc_pairwise_pvalue("ExonArray", exon_data$pairs, "HPA RNA-seq", hparna_data$pairs, population, cutoff,"common_proteins", append = TRUE)
  pvalue <-calc_pairwise_pvalue("ExonArray", exon_data$pairs, "UniGene", unigene_data$pairs, population, cutoff, "common_proteins",append = TRUE)
  pvalue <-calc_pairwise_pvalue("GNF", gnf_data$pairs, "RNA-seq", rna_data$pairs, population, cutoff, "common_proteins",append = TRUE)
  pvalue <-calc_pairwise_pvalue("GNF", gnf_data$pairs, "HPA RNA-seq", hparna_data$pairs, population, cutoff, "common_proteins",append = TRUE)
  pvalue <-calc_pairwise_pvalue("GNF", gnf_data$pairs, "UniGene", unigene_data$pairs, population, cutoff, "common_proteins",append = TRUE)
  pvalue <-calc_pairwise_pvalue("RNA-seq", rna_data$pairs, "HPA RNA-seq", hparna_data$pairs, population, cutoff, "common_proteins", append = TRUE)
  pvalue <-calc_pairwise_pvalue("RNA-seq", rna_data$pairs, "UniGene", unigene_data$pairs, population, cutoff, "common_proteins", append = TRUE)
  pvalue <-calc_pairwise_pvalue( "HPA RNA-seq", hparna_data$pairs, "UniGene", unigene_data$pairs, population, cutoff, "common_proteins",append = TRUE)
  
  
  
  mRNA_list<- list(" "= exon_data$pairs," "= gnf_data$pairs," "=rna_data$pairs," "=hparna_data$pairs," "=unigene_data$pairs)
  
  inAll <- length(intersect(x= exon_data$pairs,y= intersect(x = gnf_data$pairs, y = intersect(x = rna_data$pairs, y = intersect(x = hparna_data$pairs, y = unigene_data$pairs)))))
  inSingles <- c(length(unique(exon_data$pairs)),length(unique(gnf_data$pairs)),length(unique(rna_data$pairs)),length(unique(hparna_data$pairs)),length(unique(unigene_data$pairs)))
  pvalue_title <- paste("Pvalues Comparing mRNA datasets (only common proteins) at cutoff ", cutoff, sep='')
  
  
  filename<-NULL
  title<-""
  areas<-rep("plain",31)
  
  plot<-plot_venns(mRNA_list, col[1:5], filename, title,areas)
  if(uniprot){
    mRNA_reference_set<- read.table(file="mRNA_reference_set.tsv",header=F,sep="\t")
    names(mRNA_reference_set)<-c("proteins","tissue")
    mRNA_reference_set$pairs<-paste(mRNA_reference_set$proteins,mRNA_reference_set$tissue,sep="-")
    mRNA_total_pairs<-length(mRNA_reference_set$pairs)
    
    uniprot_data<-filter_dataset_tissue(protein_tissues,"uniprot",common_tissues,not_in=FALSE,c(3,4))
    uniprot_data$pairs<-paste(uniprot_data$proteins,uniprot_data$tissue,sep="-")
    uniprot_total_pairs<-length(uniprot_data$pairs)
    
    common_proteins<-intersect(x=mRNA_reference_set$proteins,y=uniprot_data$proteins)
    
    uniprot_data<-uniprot_data[uniprot_data$proteins %in% common_proteins,]
    
    mRNA_reference_set<-mRNA_reference_set[mRNA_reference_set$proteins %in% common_proteins,]
    
    mRNA_unique_pairs<-mRNA_total_pairs-length(mRNA_reference_set$pairs)
    uniprot_unique_pairs<-uniprot_total_pairs-length(uniprot_data$pairs)
    common_pairs<-length(intersect(mRNA_reference_set$pairs,uniprot_data$pairs))
    mRNA_uncommon<-length(mRNA_reference_set$pairs)-common_pairs    
    uniprot_uncommon<-length(uniprot_data$pairs)-common_pairs
    
    population <- expand.grid(common_proteins, common_tissues)
    
    pvalue <-calc_pairwise_pvalue("UniProt-Kb", uniprot_data$pairs, "mRNA reference", mRNA_reference_set$pairs, population, cutoff,"common_proteins", append = TRUE)
    
    
    png(file = "../figures/uniprot_mRNA_refeset_venn_diagram.png",
        width = 550, height = 450)
    title<-""
    pair_wise_venn_diagram("UniProt", "mRNA reference set",uniprot_unique_pairs , uniprot_uncommon, common_pairs, mRNA_uncommon, mRNA_unique_pairs,col[7],"black",title)
    dev.off()
  }
  
  return(plot)
  
}

#Single dataset vs integrated set: quality and coverage analysis
complementarity_analysis<-function(protein_tissues,col){
  #GNF data
  gnf_test_data<- filter_dataset_tissue(protein_tissues,"gnf",c(),not_in=TRUE,c(2,3,4))
  #common_tissues<-unique(gnf_test_data$tissue)
  common_tissues<-unique(gnf_test_data$tissue)
  
  #UniProt Data
  uniprot_data<-filter_dataset_tissue(protein_tissues,"uniprot",common_tissues,not_in=FALSE,c(3,4))
  uniprot_data$pairs<-paste(uniprot_data$proteins,uniprot_data$tissue,sep="-")
  uniprot_total_pairs<-length(uniprot_data$pairs)
  #GNF-UniProt common proteins
  common_proteins_test<-intersect(x=gnf_test_data$proteins,y=uniprot_data$proteins)
  gnf_test_data<-gnf_test_data[gnf_test_data$proteins %in% common_proteins_test,]
  
  
  
  exon_data<- filter_dataset_tissue(protein_tissues,"exon",common_tissues,not_in=FALSE,c(2,3,4))
  gnf_data<- filter_dataset_tissue(protein_tissues,"gnf",common_tissues,not_in=FALSE,c(2,3,4))
  unigene_data<- filter_dataset_tissue(protein_tissues,"unigene",common_tissues,not_in=FALSE,c(2,3,4))
  rna_data<- filter_dataset_tissue(protein_tissues,"rna",common_tissues,not_in=FALSE,c(2,3,4))
  hparna_data<- filter_dataset_tissue(protein_tissues,"hpa_rna",common_tissues,not_in=FALSE,c(2,3,4))
  hpa_data<- filter_dataset_tissue(protein_tissues,"hpa",common_tissues,not_in=FALSE,c(2,3,4))
  hpm_data<- filter_dataset_tissue(protein_tissues,"hpm",common_tissues,not_in=FALSE,c(2,3,4))
  
  exon_data<- exon_data[exon_data$score >= 250 & exon_data$proteins %in% common_proteins_test,c(1,2,3)]
  gnf_data<- gnf_data[gnf_data$score >= 250 & gnf_data$proteins %in% common_proteins_test,c(1,2,3)]
  unigene_data<- unigene_data[unigene_data$score >= 20 & unigene_data$proteins %in% common_proteins_test,c(1,2,3)]
  rna_data<- rna_data[rna_data$score >= 5 & rna_data$proteins %in% common_proteins_test,c(1,2,3)]
  hparna_data<- hparna_data[hparna_data$score >= 20 & hparna_data$proteins %in% common_proteins_test,c(1,2,3)]
  hpa_data<- hpa_data[hpa_data$score >= 10.45 & hpa_data$proteins %in% common_proteins_test,c(1,2,3)]
  hpm_data<- hpm_data[hpm_data$score >= 10 & hpm_data$proteins %in% common_proteins_test,c(1,2,3)]
  
  gnf_test_data$pairs<-paste(gnf_test_data$proteins,gnf_test_data$tissue,sep="-")
  
  exon_data$pairs<-paste(exon_data$proteins,exon_data$tissue,sep="-")
  gnf_data$pairs<-paste(gnf_data$proteins,gnf_data$tissue,sep="-")
  unigene_data$pairs<-paste(unigene_data$proteins,unigene_data$tissue,sep="-")
  rna_data$pairs<-paste(rna_data$proteins,rna_data$tissue,sep="-")
  hparna_data$pairs<-paste(hparna_data$proteins,hparna_data$tissue,sep="-")
  hpa_data$pairs<-paste(hpa_data$proteins,hpa_data$tissue,sep="-")
  hpm_data$pairs<-paste(hpm_data$proteins,hpm_data$tissue,sep="-")
  
  #Union of mRNA high-confidence pairs
  mRNA_pairs<-union(x=union(x=union(union(x=union(x=union(x=exon_data$pairs,y=gnf_data$pairs),y=unigene_data$pairs),y=rna_data$pairs),y=hparna_data$pairs),y=hpa_data$pairs),y=hpm_data$pairs)
  mRNA_proteins<-union(x=union(x=union(union(x=union(x=union(x=exon_data$proteins,y=gnf_data$proteins),y=unigene_data$proteins),y=rna_data$proteins),y=hparna_data$proteins),y=hpa_data$proteins),y=hpm_data$proteins)
  mRNA_total_pairs<-length(mRNA_pairs)
  
  #Order GNF data to get the same number of protein-tissue pairs
  gnf_test_data<-gnf_test_data[with(gnf_test_data, order(-score)), ]
  gnf_test_data<-gnf_test_data[1:mRNA_total_pairs,]
  gnf_test_data<-gnf_test_data[!is.na(gnf_test_data$pairs),]
  gnf_test_total_pairs<-length(gnf_test_data$pairs)
  
  
  
  #mRNA-UniProt common proteins
  common_proteins<-intersect(x=mRNA_proteins,y=uniprot_data$proteins)
  
    
  #mRNA common protein-tissue pairs for common proteins
  mRNA_common_protein_pairs<-union(x=union(x=union(union(x=union(x=union(x=exon_data$pairs[exon_data$proteins %in% common_proteins],y=gnf_data$pairs[gnf_data$proteins %in% common_proteins]),
                                                                 y=unigene_data$pairs[unigene_data$proteins %in% common_proteins]),y=rna_data$pairs[rna_data$proteins %in% common_proteins]),
                                                   y=hparna_data$pairs[hparna_data$proteins %in% common_proteins]),
                                           y=hpa_data$pairs[hpa_data$proteins %in% common_proteins]),
                                   y=hpm_data$pairs[hpm_data$proteins %in% common_proteins])
  
  #UniProt common protein-tissue pairs for common proteins
  uniprot_data<-uniprot_data[uniprot_data$proteins %in% common_proteins,]
  
  #Unique protein-tissue pairs for uncommon proteins
  mRNA_unique_pairs<-mRNA_total_pairs-length(mRNA_common_protein_pairs)
  uniprot_unique_pairs<-uniprot_total_pairs-length(uniprot_data$pairs)
  
  #mRNA-UniProt protein-tissue pairs
  common_pairs<-length(intersect(x=mRNA_pairs,y=uniprot_data$pairs))

  #mRNA-UniProt uncommon protein-tissue pairs
  mRNA_uncommon<-length(mRNA_common_protein_pairs)-common_pairs    
  uniprot_uncommon<-length(uniprot_data$pairs)-common_pairs
  
  #Plot
  png(file = "../figures/uniprot_trans-proteomics_set_venn_diagram.png",
      width = 550, height = 450)
  title<-""
  pair_wise_venn_diagram("UniProt", "transcriptomics\nand proteomics set",uniprot_unique_pairs , uniprot_uncommon, common_pairs, mRNA_uncommon, mRNA_unique_pairs,col[7],"grey",title)
  dev.off()
  
  
  #GNF protein-tissue pairs for common proteins
  gnf_test_pairs<-gnf_test_data$pairs[gnf_test_data$proteins %in% common_proteins_test]
  
  #UniProt protein-tissue pairs for common proteins
  uniprot_data_test<-uniprot_data[uniprot_data$proteins %in% common_proteins_test,]
  #Unique protein-tissue pairs = total number of pairs - protein-tissue pairs for common tissues
  gnf_test_unique_pairs<-gnf_test_total_pairs-length(gnf_test_pairs)
  uniprot_unique_test_pairs<-uniprot_total_pairs-length(uniprot_data_test$pairs)
  #GNF-UniProt common tissue-proteins pairs
  common_pairs_test<-length(intersect(x=gnf_test_pairs,y=uniprot_data$pairs))
  #GNF-UniProt uncommon pairs
  gnf_test_uncommon<-length(gnf_test_pairs)-common_pairs_test    
  uniprot_test_uncommon<-length(uniprot_data$pairs)-common_pairs_test
  
  ##create supplementary files
  write.table(intersect(x=gnf_test_pairs,y=uniprot_data$pairs), file="Common_GNF_UniProtKB.tsv", row.names=FALSE, quote = FALSE, col.names = c("Gene-tissue associations"))
  write.table(intersect(x=mRNA_pairs,y=uniprot_data$pairs), file="Common_Integrated_set_UniProtKB.tsv", row.names=FALSE, quote = FALSE, col.names = c("Gene-tissue associations"))
  write.table(gnf_test_pairs, file="GNF_coverage_associations.tsv", row.names=FALSE, quote = FALSE, col.names = c("Gene-tissue associations"))
  write.table(mRNA_pairs, file="Integrated_set.tsv", row.names=FALSE, quote = FALSE, col.names = c("Gene-tissue associations"))
  
  png(file = "../figures/uniprot_gnf_set_venn_diagram.png",
      width = 550, height = 450)
  title<-""
  pair_wise_venn_diagram("UniProt", "GNF set",uniprot_unique_test_pairs , uniprot_test_uncommon, common_pairs_test, gnf_test_uncommon, gnf_test_unique_pairs,col[7],col[2],title)
  dev.off()
}

###############################
#           Analysis         #
##############################
 #Colors
col<-get_colors()
colScale <- scale_fill_manual(name = "",values = col,labels= c("Exon Array","RNA-seq atlas","HPA RNA-seq","GNF","UniGene","HPA","UniProt","Text-mining", "HPM"))

protein_tissues<-get_data("low")
low_venn<-get_venn_diagrams(protein_tissues,"low",uniprot=F,col, append =FALSE)
complementarity_analysis(protein_tissues,col)

protein_tissues<-get_data("medium")
medium_venn<-get_venn_diagrams(protein_tissues,"medium",uniprot=F,col, append =TRUE)

protein_tissues<-get_data("high")
high_venn<-get_venn_diagrams(protein_tissues,"high",uniprot=T,col, append =TRUE)


png(file = "../figures/All_mRNA_venn_diagrams.png",
    width = 450, height = 850)
#Piece of code borrowed from http://stackoverflow.com/questions/14243609/problems-with-venndiagram
gl <- grid.layout(nrow=3, ncol=1)
# setup viewports
vp.1 <- viewport(layout.pos.col=1, layout.pos.row=1) 
vp.2 <- viewport(layout.pos.col=1, layout.pos.row=2) 
vp.3 <- viewport(layout.pos.col=1, layout.pos.row=3) 
# init layout
pushViewport(viewport(layout=gl))
# access the first position
pushViewport(vp.1)
# start new base graphics in first viewport
grid.draw(low_venn)
# done with the first viewport
popViewport()
# move to the next viewport
pushViewport(vp.2)
grid.draw(medium_venn)
# done with this viewport
popViewport(1)
# move to the next viewport
pushViewport(vp.3)
grid.draw(high_venn)
# done with this viewport
popViewport(1)
dev.off()