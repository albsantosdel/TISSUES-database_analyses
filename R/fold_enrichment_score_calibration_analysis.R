datasets<-c("exon","gnf","rna","hpa_rna","unigene","hpa","tm","hpm")
dataset_labels<- c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene","HPA","Text-mining","HPM")

title <- "Datasets fold enrichment per confidence score"
out_file <- "../figures/datasets_score_calibration.png"
png(out_file,height=750,width=750)#,useDingbats=FALSE)


  
#assign colors to each dataset
col<-get_colors()
colScale <- scale_fill_manual(name = "",values = col,labels= c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene","HPA","Text-mining","HPM"))


plot<- ggplot()+colScale


y <- 1.5 
#y <- 1.5 
j<- 0.2
for (i in 1:length(datasets) ) {
  dataset<- datasets[i]
  label<- dataset_labels[i]
  file<-paste(dataset,"uniprot_fold_enrichment_analysis.tsv",sep="_")
  color<-col[names(col)==dataset][1]
    
  data<-read.csv(file=file,header=FALSE,sep="\t")
  names(data)<-c("mean_stars","fold_enrichment","score")
  x <- 3.6
  #x <- 4.5
  y<- y - j
  data2.labels <- data.frame(x = c(x),y = c(y),label = c(label))
  
  
  plot<-plot + geom_point(data=data,aes(x=mean_stars,y=fold_enrichment),size=1,color= color, shape = 20) #+
    #geom_text(data= data2.labels, aes(x=x,y=y,label= label), color=color,size = 7,fontface="bold")
}
plot <- plot  +
  #ggtitle(title) +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank())+
  theme(strip.background=element_rect(fill='white',color='white'))+
  theme(strip.text.y = element_blank())+
  theme(plot.title = element_text(size = rel(1), colour = "#003366"))+
  theme(axis.text=element_text(size=rel(1)))+
  xlab(label="mean confidence score")+
  ylab(label="Fold enrichment compared with UniProt")+
  theme(axis.title=element_text(size=rel(1)))+
  theme(legend.text = element_text(size = 5))+
  theme(axis.text.x = element_text(colour = "black"))+
  theme(axis.text.y = element_text(colour = "black"))+
  theme(legend.title = element_text(size = 5,face="bold"))

print(plot)
dev.off()