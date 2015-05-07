plot_prots_num_tissues<-function(data,title,file_name,sc="free",multi=F,color=NULL,stacked=F){
  
  image_name<- file_name
  
  #Necessary to plot x-axis in all facets
  g <- gridExtra::borderGrob(type=2, colour="black", lwd=2)
  plot<-ggplot(data,aes(x=factor(num_tissues)))+
    facet_grid(cutoff~dataset,scales=sc,space = "free")+
    #ggtitle(title)+
    annotation_custom(g)+#Necessary to plot x-axis in all facets
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())+
    scale_x_discrete(breaks=1:max(data$num),labels=1:max(data$num))+
    scale_y_continuous(labels = percent_format(),expand = c(0,0))+
    theme(strip.background=element_blank())+
    theme(strip.text.x = element_blank())+
    theme(strip.text.y = element_text(family= "Helvetica",size = rel(1), colour = "black",hjust=0.5,vjust=1))+
    theme(axis.text.x = element_text(family= "Helvetica",colour = "black"))+
    theme(axis.text.y = element_text(family= "Helvetica",colour = "black"))+
    theme(axis.text = element_text(family= "Helvetica",size=rel(1)))+
    theme(axis.title = element_text(family= "Helvetica",size=rel(1)))+
    xlab("Number of tissues \n ")+ylab("Proteins")+
    theme(plot.title = element_text(family= "Helvetica",size = rel(1), colour = "#003366",vjust=0.7))+
    theme(legend.position="bottom")+
    theme(legend.text = element_text(family= "Helvetica",size = rel(1)))+
    theme(legend.key.width=unit(1,"line"))
  
  if(multi){
    png(image_name,height=450,width=950)
    plot<-complete_multi_plot(plot)+theme(panel.margin=unit(1,"lines"))
  }
  else if(!stacked){
    png(image_name,height=450,width=450)
    plot<-complete_single_plot(plot,color)+theme(panel.margin=unit(4,"lines"))
  }
  else{
    png(image_name,height=450,width=450)
    plot<-complete_stacked_plot(plot,color)+theme(panel.margin=unit(4,"lines"))
  }
  
  switch_facet_strip(plot,switch = "y")
  dev.off()
}  
#Plot multiple distributions in a single figure
complete_multi_plot<-function(plot){
  plot<-plot+
    colScale+
    geom_bar(aes(num=num,y=..count../num,fill=dataset),width=0.5)
  return(plot)
}
#Plot distribution in a single plot
complete_single_plot<-function(plot,color){
  plot<-plot+
    geom_bar(aes(num=num,y=..count../num),fill=color,width=0.5)+
    colScale
  return(plot)
}
#HPA IHC stacked plot single/multiple antibodies
complete_stacked_plot<-function(plot,color){
  color_single<-"#f08526"
  color_multiple<-color[[1]]
  plot<-plot+
    geom_bar(aes(num=num,y=..count../num,fill=factor(single),order=-single),width=0.5)+
    scale_fill_manual(labels=c("multiple antibodies","single antibody"),values=c(color_multiple,color_single))
  return(plot)
}

###############################
#           Analysis         #
##############################

#Colors
col<-get_colors(is_label = T)
colScale <- scale_fill_manual(name = "",values = col,labels= c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene","HPA","UniProt","Text-mining", "HPM"))

#Load data
high<-read.table(file="high_cutoff_expression_breadth.tsv",header=F,sep="\t")
medium<-read.table(file="medium_cutoff_expression_breadth.tsv",header=F,sep="\t")
low<-read.table(file="low_cutoff_expression_breadth.tsv",header=F,sep="\t")
names(high)<-c("dataset","num_tissues","prot")
names(medium)<-c("dataset","num_tissues","prot")
names(low)<-c("dataset","num_tissues","prot")

high$cutoff<-"High"
medium$cutoff<-"Medium"
low$cutoff<-"Low"
data<-rbind(high,medium,low)
data$cutoff<-factor(data$cutoff, levels=c("Low","Medium","High"))

data$num<-0

data$num[data$dataset == "exon"]<-length(unique(data$prot[data$dataset == "exon"]))
data$num[data$dataset == "gnf"]<-length(unique(data$prot[data$dataset == "gnf"]))
data$num[data$dataset == "rna"]<-length(unique(data$prot[data$dataset == "rna"]))
data$num[data$dataset == "hpa_rna"]<-length(unique(data$prot[data$dataset == "hpa_rna"]))
data$num[data$dataset == "unigene"]<-length(unique(data$prot[data$dataset == "unigene"]))
data$num[data$dataset == "hpa"]<-length(unique(data$prot[data$dataset == "hpa"]))
data$num[data$dataset == "uniprot"]<-length(unique(data$prot[data$dataset == "uniprot"]))
data$num[data$dataset == "tm"]<-length(unique(data$prot[data$dataset == "tm"]))
data$num[data$dataset == "hpm"]<-length(unique(data$prot[data$dataset == "hpm"]))

hpa_data<-data[data$dataset == "hpa",]
hpa_data$dataset<-"HPA"

#Get information about single and multiple antibodies HPA IHC
single_antibody_proteins<-read.table(file="hpa_single_antibody_proteins.tsv",header=F)
hpa_data$single<- 0
hpa_data$single[hpa_data$prot %in% single_antibody_proteins$V1]<-1


uniprot_data<-data[data$dataset == "uniprot" & data$cutoff == "High",]
uniprot_data$dataset<-"UniProt"

tm_data<-data[data$dataset == "tm",]
tm_data$dataset<-"Text-mining"

hpm_data<-data[data$dataset == "hpm",]
hpm_data$dataset<-"HPM"


data<-data[data$dataset %in% c("exon","gnf","rna","hpa_rna","unigene"),]
data$dataset<-factor(data$dataset,levels=c("exon","gnf","rna","hpa_rna","unigene"),labels= c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene"))

#Plot the expression breadth distribution for each dataset
title<- paste("Number of proteins in number of tissues per data set", sep="")
plot_prots_num_tissues(data,title,"../figures/mRNA_datasets_prots_num_tissues.png",sc="free",multi=TRUE,stacked=F)

title<- paste("HPA number of proteins in number of tissues per data set", sep="")
color<-col[names(col)=="HPA"][1]
plot_prots_num_tissues(hpa_data,title,"../figures/hpa_datasets_prots_num_tissues.png",sc="fixed",FALSE,color,stacked=T)

title<- paste("UniProt number of proteins in number of tissues per data set", sep="")
color<-col[names(col)=="UniProt"][1]
plot_prots_num_tissues(uniprot_data,title,"../figures/uniprot_datasets_prots_num_tissues.png",sc="free",multi=FALSE,color,stacked=F)

title<- paste("Text-mining number of proteins in number of tissues per data set", sep="")
color<-col[names(col)=="Text-mining"][1]
plot_prots_num_tissues(tm_data,title,"../figures/tm_datasets_prots_num_tissues.png",sc="fixed",FALSE,color,stacked=F)

title<- paste("HPM number of proteins in number of tissues per data set", sep="")
color<-col[names(col)=="HPM"][1]
plot_prots_num_tissues(hpm_data,title,"../figures/hpm_datasets_prots_num_tissues.png",sc="fixed",FALSE,color,stacked=F)