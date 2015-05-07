setwd(dir="data")

levels<-c("uniprot","tm","hpa","hpa_rna","gnf","unigene","hpm","rna","exon")
labels<-c("UniProt","Text-mining","HPA","HPA RNA-seq","GNF","UniGene","HPM","RNA-seq atlas","Exon Array")

tis<-c("heart","intestine","kidney","liver","nervous system","muscle","lung","spleen","lymph node","pancreas","skin",
       "stomach","thyroid gland","bone marrow","saliva","adrenal gland","urine","blood","eye","gall bladder","bone")
tis_labels<-c("Heart","Intestine","Kidney","Liver","Nervous\nsystem","Muscle","Lung","Spleen","Lymph\nnode","Pancreas","Skin",
       "Stomach","Thyroid\ngland","Bone\nmarrow","Saliva","Adrenal\ngland","Urine","Blood","Eye","Gall\nbladder","Bone")

data<-read.table(file="datasets_major_tissues.tsv",sep="\t",header=F)
names(data)<-c("dataset","tissues","num_prots")

#Update the number of proteins
data<-data[,c(1,2)]

data$dataset<-factor(data$dataset,levels=levels,labels= labels)
data$tissues<-factor(data$tissues,levels=tis,labels=tis_labels)

gnf_prots<-11015
exon_prots<-12945
unigene_prots<-14045
uniprot_prots<-14722
tm_prots<-14760
rna_prots<-15332
hpa_prots<-15517
hpa_rna_prots<-16082
hpm_prots<-17038

data$num_prots[data$dataset == "GNF"]<- gnf_prots
data$num_prots[data$dataset == "Exon Array"]<- exon_prots
data$num_prots[data$dataset == "RNA-seq atlas"]<- rna_prots
data$num_prots[data$dataset == "HPA RNA-seq"]<- hpa_rna_prots
data$num_prots[data$dataset == "HPA"]<- hpa_prots
data$num_prots[data$dataset == "HPM"]<- hpm_prots
data$num_prots[data$dataset == "UniGene"]<- unigene_prots
data$num_prots[data$dataset == "Text-mining"]<- tm_prots
data$num_prots[data$dataset == "UniProt"]<- uniprot_prots

datasets<-unique(data$dataset)

#Colors
col<-get_colors(is_label = T)
fillcolScale <- scale_fill_manual(name = "",values = col,labels= c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene","HPA","UniProt","Text-mining","HPM"))
colScale <- scale_color_manual(name = "",values = col,labels= c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene","HPA","UniProt","Text-mining","HPM"))

png("../figures/Summary_figure_tissues_num_proteins.png",height=750,width=1100)

plot1<-ggplot(data=data,aes(dataset,tissues,colour=dataset))+geom_point(shape=15,size=7.5)+
  coord_flip()+
  colScale+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  #xlab(label="Tissues and number of proteins per dataset")+
  theme(axis.text.x = element_text(angle = 90,hjust=0))+
  theme(axis.text.y = element_text(angle = 90,hjust = 0.5,vjust=0.5))+
  theme(legend.position="none")+
  theme(line = element_blank(),
        line = element_blank(),
        axis.title.y=element_text(color="white"))+
  theme(axis.text.x = element_text(family= "Helvetica",colour = "black"))+
  theme(axis.text.y = element_text(family= "Helvetica",colour = "black"))+
  theme(axis.text = element_text(family= "Helvetica",size=rel(1)))+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())
  #theme(axis.title.y = element_text(family= "Helvetica",size = rel(3), colour = "#003366",vjust=0.5))

plot2<-ggplot(data=data,aes(dataset,num_prots,fill=dataset))+geom_bar(stat="identity",position = "identity",width=0.2)+
  coord_flip()+
  fillcolScale+
  theme_bw()+
  theme(axis.line = element_line(colour = "white"),
        panel.border = element_blank())+
  ylab(label="Number of proteins")+
  scale_y_reverse()+
  #theme(axis.text.y = element_text(hjust = 0.5,vjust=0.5))+  
  theme(axis.text.x = element_text(angle = 90, family= "Helvetica",colour = "black",size=rel(1)))+  
  theme(axis.title.x = element_text(angle = 180,vjust=1,family= "Helvetica",size = rel(1)))+  
  theme(legend.position="none")+
  theme(axis.ticks = element_line(color = "white"))+
  theme(panel.grid.major.x = element_line(colour = "grey"))+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank())

gA <- ggplotGrob(plot1)
gB <- ggplotGrob(plot2)
maxHeight = grid::unit.pmax(gA$heights[2:5], gB$heights[2:5])
gA$heights[2:5] <- as.list(maxHeight)
gB$heights[2:5] <- as.list(maxHeight)
grid.arrange(gA, gB, ncol=2)
dev.off()