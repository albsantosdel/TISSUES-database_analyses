plot_multi_fold_enrichment<-function(data,out_file,title,labels){
  png(out_file,height=850,width=1200)
  g <- gridExtra::borderGrob(type=9, colour="black", lwd=2)
  plot<- ggplot(data, aes(mean_score, fold_enrichment,color=set)) + 
    colScale+
    facet_grid(set~.,space="free",scales="free_x") +
    #ggtitle(title) +
    xlab(label="Mean score")+
    ylab(label="Fold enrichment compared with UniProt")+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    annotation_custom(g)+
    scale_x_log10(breaks=c(1,100,1000,2000,5000))+
    annotation_logticks(sides="b",, short = unit(0.02, "cm"), 
                        mid = unit(0.05, "cm"), long = unit(0.1, "cm"),color="black")+
    theme(strip.background=element_rect(fill='white',color='white'))+
    theme(strip.text.y = element_blank())+
    theme(axis.text.x = element_text(family= "Helvetica",colour = "black"))+
    theme(axis.text.y = element_text(family= "Helvetica",colour = "black"))+
    theme(axis.text = element_text(family= "Helvetica",size=rel(1)))+
    theme(axis.title = element_text(family= "Helvetica",size=rel(1)))+
    theme(legend.title = element_blank())+
    #theme(plot.title = element_text(family= "Helvetica",size = rel(2.5), colour = "#003366",vjust=0.7))+
    theme(legend.position="bottom")+
    theme(legend.text = element_text(family= "Helvetica",size = rel(1)))+
    theme(legend.key.width=unit(1,"line"))+
    theme(legend.key = element_rect(fill = "white",colour = "white")) + 
    guides(colour = guide_legend(override.aes = list(size=3)))+
    theme(legend.background = element_blank())
  
  for(j in 1:length(labels)){
    dataset<-labels[j]
    susbset_data<-subset(data,set = dataset)
    #method="gam",formula=y~s(x,bs='tp')
    plot<- plot + geom_point(data=susbset_data,aes(x=mean_score,y=fold_enrichment),size=1)+
       geom_line(aes(x=low_cutoff))+
       geom_line(aes(x=medium_cutoff))+
       geom_line(aes(x=high_cutoff))
  }
  
  print(plot)
  dev.off()
}
  
plot_single_fold_enrichment<-function(data,out_file,title,show_legend=F,log=F,labels){
  png(out_file,height=650,width=700)
  #Necessary to plot x-axis in all facets
  g <- gridExtra::borderGrob(type=2, colour="black", lwd=2)
  plot<- ggplot(data=data,aes(x=mean_score,y=fold_enrichment)) +
    geom_point(aes(x=mean_score,y=fold_enrichment,shape=factor(shape),color=color,size = size))+
    theme(panel.margin=unit(1,"lines"))+
    scale_size_area(breaks = unique(data$size), 
                    labels = unique(data$shape), 
                    guide = "legend",max_size=max(data$size))+
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_color_identity(breaks = unique(data$color), 
                         labels = unique(data$shape), 
                         guide = "legend")+
    scale_shape_manual(breaks = factor(data$shape), 
                       values = c(19, 19,18))+ 
    scale_x_continuous(breaks = seq(0, max(data$mean_score), by = 1)) +
    scale_y_continuous(breaks = seq(0, max(data$fold_enrichment), by = 0.5))+
    guides(color = guide_legend("title"),shape = guide_legend("title"),size = guide_legend("title"))+
    facet_grid(gold~set,scales="free",space="free")+
    annotation_custom(g)+#Necessary to plot x-axis in all facets
    annotation_logticks(sides="b",, short = unit(0, "cm"), 
                        mid = unit(0, "cm"), long = unit(0.1, "cm"),color="black")+
    xlab(label="Mean score")+
    ylab(label=" Fold enrichment")+
    theme(strip.background=element_rect(fill='white',color='white'))+
    theme(strip.text.y = element_text(family= "Helvetica",size = rel(1), colour = "black"))+
    theme(axis.text.y = element_text(family= "Helvetica",colour = "black"))+
    theme(axis.text = element_text(family= "Helvetica",size=rel(1)))+
    theme(axis.title = element_text(family= "Helvetica",size=rel(1)))+
    theme(legend.title = element_blank())+
    #theme(plot.title = element_text(family= "Helvetica",size = rel(2.5), colour = "#003366",vjust=0.7))+
    theme(legend.position="bottom")+
    theme(legend.text = element_text(family= "Helvetica",size = rel(1)))+
    theme(legend.key.width=unit(1,"line"))+
    theme(legend.key = element_rect(fill = "white",colour = "white")) + 
    theme(legend.background = element_blank())
  if(!show_legend){
    plot<-plot+
      theme(legend.position="none")
  }
  if(!log){
    plot <- plot + geom_line(data=data,aes(x=low_cutoff,color=color))+
      geom_line(aes(x=medium_cutoff,color=color))+
      geom_line(aes(x=high_cutoff,color=color))  
  }
  else{
    plot <- plot + geom_line(data=data,aes(x=log10(low_cutoff)),color=color)+
      geom_line(aes(x=log10(medium_cutoff)),color=color)+
      geom_line(aes(x=log10(high_cutoff)),color=color)
  }
  
  switch_facet_strip(plot,switch = "y")
  
  dev.off()
}
  
get_mRNA_data<-function(levels,labels, low_cutoff,medium_cutoff,high_cutoff){
  dataset<- levels[1]
  file<-paste(dataset,"uniprot_fold_enrichment_analysis.tsv",sep="_")
  
  all_data<-read.csv(file=file,header=FALSE,sep="\t")
  names(all_data)<-c("mean_stars","fold_enrichment","mean_score")
  all_data$set<-dataset
  
  for (i in 2:length(levels) ) {
    dataset<- levels[i]
    file<-paste(dataset,"uniprot_fold_enrichment_analysis.tsv",sep="_")
    data<-read.csv(file=file,header=FALSE,sep="\t")
    names(data)<-c("mean_stars","fold_enrichment","mean_score")
    data$set<-dataset
    all_data<-rbind(all_data,data)
  }
  all_data$set<-factor(all_data$set,levels=levels,labels= labels)
  
  cutoffs<-cbind.data.frame(labels,low_cutoff,medium_cutoff,high_cutoff)
  
  all_data<-merge(all_data,cutoffs,by.x="set",by.y="labels")
  return(all_data)
}

###############################
#           Analysis         #
##############################
#Colors
col<-get_colors(is_label = T)
colScale <- scale_color_manual(name = "",values = col,labels=c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene","HPA","UniProt","Text-mining","HPM"))

#Get mRNA data
levels<-c("exon","gnf","rna","hpa_rna","unigene")
labels<-c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene")
low_cutoff<-c(50,50,0.5,1,1)
medium_cutoff<-c(100,100,1,10,10)
high_cutoff<-c(250,250,5,20,20)
data<-get_mRNA_data(levels,labels,low_cutoff,medium_cutoff,high_cutoff)  

#Plot fold-enrichment
title <- "Datasets fold enrichment per score"
out_file <- "../figures/datasets_fold_enrichment.png"
plot_multi_fold_enrichment(data,out_file,title,labels)

  
  
#Fold enrichment HPA and Text-mining
title <- "HPA fold enrichment vs. UniProt/mRNA reference set"
out_file <- "../figures/hpa_fold_enrichment.png"
color<-col[names(col)=="HPA"][1]

file<-"hpa_uniprot_fold_enrichment_analysis.tsv"
hpa_uniprot_data<-read.csv(file=file,header=FALSE,sep="\t")
names(hpa_uniprot_data)<-c("mean_stars","fold_enrichment","mean_score")
hpa_uniprot_data$set<-"HPA"
hpa_uniprot_data$gold<-"UniProt"

hpa_single<-hpa_uniprot_data[hpa_uniprot_data$mean_score %in% c(1,3,6),]
hpa_multi<-hpa_uniprot_data[!hpa_uniprot_data$mean_score %in% c(1,3,6),]

hpa_single$fold_enrichment[hpa_single$mean_score ==1]<-summary(hpa_single[hpa_single$mean_score ==1,2])[[4]]
hpa_single$fold_enrichment[hpa_single$mean_score ==3]<-summary(hpa_single[hpa_single$mean_score ==3,2])[[4]]
hpa_single$fold_enrichment[hpa_single$mean_score ==6]<-summary(hpa_single[hpa_single$mean_score ==6,2])[[4]]

hpa_single$shape<-"single antibody"
hpa_multi$shape<- "multiple antibodies"


hpa_single$size<- 1.5
hpa_multi$size<- 1


hpa_single$color<-"#f08526"
hpa_multi$color<- color


file<-"hpa_mRNA_reference_combined_fold_enrichment_analysis.tsv"
hpa_mRNA_data<-read.csv(file=file,header=FALSE,sep="\t")
names(hpa_mRNA_data)<-c("mean_stars","fold_enrichment","mean_score")
hpa_mRNA_data$set<-"HPA"
hpa_mRNA_data$gold<-"mRNA reference set"

hpa_single_mRNA<-hpa_mRNA_data[hpa_mRNA_data$mean_score %in% c(1,3,6),]
hpa_multi_mRNA<-hpa_mRNA_data[!hpa_mRNA_data$mean_score %in% c(1,3,6),]

hpa_single_mRNA$shape<-"single antibody"
hpa_multi_mRNA$shape<-"multiple antibodies"

hpa_single_mRNA$size<- 1.5
hpa_multi_mRNA$size<- 1


hpa_single_mRNA$color<-"#f08526"
hpa_multi_mRNA$color<- color

hpa_single_mRNA$fold_enrichment[hpa_single_mRNA$mean_score ==1]<-summary(hpa_single_mRNA[hpa_single_mRNA$mean_score ==1,2])[[4]]
hpa_single_mRNA$fold_enrichment[hpa_single_mRNA$mean_score ==3]<-summary(hpa_single_mRNA[hpa_single_mRNA$mean_score ==3,2])[[4]]
hpa_single_mRNA$fold_enrichment[hpa_single_mRNA$mean_score ==6]<-summary(hpa_single_mRNA[hpa_single_mRNA$mean_score ==6,2])[[4]]

data<-rbind(hpa_multi,hpa_multi_mRNA,hpa_single,hpa_single_mRNA)
data$gold<-factor(data$gold, levels=c("UniProt","mRNA reference set"))

data$low_cutoff<-1
data$medium_cutoff<-6
data$high_cutoff<-10

hpa_data<-data
#plot_single_fold_enrichment(data,out_file,title,show_legend=T)

title <- "Text-mining fold enrichment vs. UniProt/mRNA reference set"
out_file <- "../figures/tm_fold_enrichment.png"
color<-col[names(col)=="Text-mining"][1]

file<-"tm_uniprot_fold_enrichment_analysis.tsv"
tm_uniprot_data<-read.csv(file=file,header=FALSE,sep="\t")
names(tm_uniprot_data)<-c("mean_stars","fold_enrichment","mean_score")
tm_uniprot_data$set<-"Text-mining"
tm_uniprot_data$gold<-"UniProt"
tm_uniprot_data$shape<-"multiple antibodies"

file<-"tm_mRNA_reference_combined_fold_enrichment_analysis.tsv"
tm_mRNA_data<-read.csv(file=file,header=FALSE,sep="\t")
names(tm_mRNA_data)<-c("mean_stars","fold_enrichment","mean_score")
tm_mRNA_data$set<-"Text-mining"
tm_mRNA_data$gold<-"mRNA reference set"
tm_mRNA_data$shape<-"multiple antibodies"

data<-rbind(tm_uniprot_data,tm_mRNA_data)
data$gold<-factor(data$gold, levels=c("UniProt","mRNA reference set"))
data$color<-color
data$size<- 1.1

low_cutoff<-1
medium_cutoff<-2.5
high_cutoff<-3.5
labels<-c("Text-mining")
cutoffs<-cbind.data.frame(labels=labels,low_cutoff,medium_cutoff,high_cutoff)

data<-merge(data,cutoffs,by.x="set",by.y="labels")

plot_single_fold_enrichment(data,out_file,title,show_legend = F,log=F,labels)

#HPM data
title <- "HPM fold enrichment vs. UniProt/mRNA reference set"
out_file <- "../figures/hpm_fold_enrichment.png"
color<-col[names(col)=="HPM"][1]

file<-"hpm_uniprot_fold_enrichment_analysis.tsv"
hpm_uniprot_data<-read.csv(file=file,header=FALSE,sep="\t")
names(hpm_uniprot_data)<-c("mean_stars","fold_enrichment","mean_score")
hpm_uniprot_data$set<-"HPM"
hpm_uniprot_data$gold<-"UniProt"
hpm_uniprot_data$shape<-"not needed"

file<-"hpm_mRNA_reference_combined_fold_enrichment_analysis.tsv"
hpm_mRNA_data<-read.csv(file=file,header=FALSE,sep="\t")
names(hpm_mRNA_data)<-c("mean_stars","fold_enrichment","mean_score")
hpm_mRNA_data$set<-"HPM"
hpm_mRNA_data$gold<-"mRNA reference set"
hpm_mRNA_data$shape<-"not needed"

data<-rbind(hpm_uniprot_data,hpm_mRNA_data)
data$gold<-factor(data$gold, levels=c("UniProt","mRNA reference set"))
data$color<-color
data$size<- 1.2

low_cutoff<-1
medium_cutoff<-5
high_cutoff<-10

cutoffs<-cbind.data.frame(labels=c("HPM"),low_cutoff,medium_cutoff,high_cutoff)

data<-merge(data,cutoffs,by.x="set",by.y="labels")
hpm_data<-data

#plot_single_fold_enrichment(data,out_file,title,show_legend = F)

#Plot HPA and HPM together
hpm_data$mean_score<-log10(hpm_data$mean_score)
hpa_data$mean_score<-log10(hpa_data$mean_score)

all_data <- hpa_data

all_data<-rbind(all_data,hpm_data)
out_file<- "../figures/hpm_hpa_fold_enrichment.png"
labels<-c("HPA","HPM")
plot_single_fold_enrichment(all_data,out_file,title,show_legend = F,log=T,labels)
