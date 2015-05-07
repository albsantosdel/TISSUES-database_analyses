#Function from http://stackoverflow.com/questions/6625691/is-it-possible-to-switch-the-side-of-y-axis-breaks-and-labels-on-a-faceted-plot
switch_facet_strip <- function(p, switch = c("x", "y")) {
  rbind_gtable <- gtable:::rbind_gtable
  cbind_gtable <- gtable:::cbind_gtable
  
  if ("y" %in% switch)
    p <- p + theme(strip.text.y = element_text(vjust = 0.5, angle = 90))
  
  g <- ggplotGrob(p)
  
  
  gdim <- as.numeric(g$layout[g$layout$name == "background", c("b", "r")])
  tpos <- g$layout[g$layout$name == "strip-top", "b"][1]
  rpos <- g$layout[g$layout$name == "strip-right", "r"][1]
  new_tpos <- g$layout[g$layout$name == "axis-b", "b"][1] + 1
  new_rpos <- g$layout[g$layout$name == "axis-l", "r"][1] - 1
  
  if ("x" %in% switch) {
    g <- rbind_gtable(
      rbind_gtable(
        gtable_add_rows(
          rbind_gtable(g[1:tpos-1, ] , g[(tpos+1):(new_tpos-1), ], "first"),
          unit(5, units = "mm")),
        g[tpos, ], "first"),
      g[new_tpos:gdim[1], ], "first")
  }
  
  if ("y" %in% switch) {
    g <- cbind_gtable(
      cbind_gtable(
        gtable_add_cols(
          cbind_gtable(g[, 1:new_rpos], g[, rpos], "first"),
          unit(1, units = "mm")),
        g[, (new_rpos+2):rpos-1], "first"),
      g[, (rpos+1):gdim[2]], "first")
  }
  
  grid.newpage()
  grid.draw(g)
}

#Define the colors to be used for each dataset
get_colors<-function(is_label = F){
  #Names
  ds<-c("exon","rna","hpa_rna","gnf","unigene","hpa","uniprot","tm", "hpm")
  if(is_label){
    ds<-c("Exon Array","GNF","RNA-seq atlas","HPA RNA-seq","UniGene","HPA","UniProt","Text-mining", "HPM")
  }
  #assign colors to each dataset
  col<-c("#ef4023","#f7931d","#86bf86","#2b948b","#a7a9ab","#0077be","#d4a7cd","#9957a2","#52b2e4")
  names(col) <- ds
  
  return(col)
}


## calculate p-value for N circle venn diagram 
## code inspired from https://bioinfomagician.wordpress.com/2014/03/04/p-value-for-intersection-of-multiple-circle-venn-diagram/
## adapted by Anders Boeck Jensen
calc_venn_pvalue<-function(N,inAll,inSingles, nIterations, title,  append = FALSE){
  #N=7000 # total number of genes for selection
  #inAll=3 # number of overlaps in all three sets
  #inSingles = c(200,250,300) = number of up-regulated genes in each sample
  
  #create a dataframe to store frequency table
  d<-data.frame(0:min(inSingles),rep(0,min(inSingles)+1))
  #To loop over (just to make it once)
  loopValues = c(2:length(inSingles))
  for (i in 1:nIterations){
    A=sample(1:N,size=inSingles[1],replace=FALSE)
    sampledInAll = rep(TRUE, length(A))
    for(i in loopValues){
      B=sample(1:N,size=inSingles[i],replace=FALSE)
      sampledInAll = sampledInAll&(A %in% B)
      if(sum(sampledInAll)==0){
        break
      }
    }
    e = sum(sampledInAll)
    d[e+1,2]<-d[e+1,2]+1
  }
  
  colnames(d)<-c("Intersect","p-value")
  # convert counts to ratios
  d[,2]<-d[,2]/10000
  # calculate p-value
  p<-sum(d[(inAll+1):(min(inSingles)+1),2])
  
  text <- paste("P-value ", title, "\n", sep='')
  text<-paste(text,"p-value = ", p, "\n", sep='')
  
  write(text, file = "multiple_ways_pvalue_results.txt", append = append)
  #print(d)
  return(p)
}
