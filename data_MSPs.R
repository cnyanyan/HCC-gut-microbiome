## MSPs
metadata_DEN <- as.data.frame(read.csv("~/Desktop/Project/data/metadata_DEN.csv",header=T,sep=","))
metadata_DEN$labels = paste(metadata_DEN$Experiment, metadata_DEN$Feedin, sep='_')
index_DEN <- match(metadata_DEN$SampleID,colnames(mgs.med.vect))

metadata_nonDEN <- as.data.frame(read.csv("~/Desktop/Project/data/metadata_nonDEN.csv",header=T,sep=","))
metadata_nonDEN$labels = paste(metadata_nonDEN$Experiment, metadata_nonDEN$Feedin, sep='_')
index_nonDEN <- match(metadata_nonDEN$SampleID,colnames(mgs.med.vect))



### DEN_Fasted - DEN_Fed
n.DEN <- 0
p.DEN <- c()
for (i in 1:168){
  metadata_DEN$mgs <- as.matrix(unlist(mgs.med.vect[i,index_DEN]))
  res.DEN <- wilcox.test(mgs ~ labels, data = metadata_DEN, paired = FALSE)
  p.DEN <- c(p.DEN,res.DEN$p.value)
  if(sum(metadata_DEN$mgs>0)){
  if(res.DEN$p.value<0.01){
    print(paste(rownames(mgs.med.vect)[i],mgs_list[i,2],mgs_list[i,3]));
    n.DEN = n.DEN + 1;
    filename.DEN=paste("DEN",rownames(mgs.med.vect)[i],mgs_list[i,2],".jpeg",sep="-")
    jpeg(file=filename.DEN)
    print(ggboxplot(metadata_DEN, x = "labels", y = "mgs", 
                    color = "labels",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
            stat_compare_means(method = "wilcox.test", label.x = 1.3, label= "p") + 
            labs(x = "DEN_Fasted - DEN_Fed", y = paste(mgs_list[i,2])))
    dev.off()
  } else{}
    }else{}
}
print(n.DEN)
p.DEN <- as.matrix(p.DEN)
rownames(p.DEN) <- rownames(mgs.med.vect)


### nonDEN_Fasted - nonDEN_Fed
n.nonDEN <- 0
p.nonDEN <- c()
for (i in 1:168){
  metadata_nonDEN$mgs <- as.matrix(unlist(mgs.med.vect[i,index_nonDEN]))
  res.nonDEN <- wilcox.test(mgs ~ labels, data = metadata_nonDEN, paired = FALSE)
  p.nonDEN <- c(p.nonDEN,res.nonDEN$p.value)
  if(sum(metadata_nonDEN$mgs>0)){
    if(res.nonDEN$p.value<0.01){
      print(paste(rownames(mgs.med.vect)[i],mgs_list[i,2],mgs_list[i,3]));
      n.nonDEN = n.nonDEN + 1;
      filename.nonDEN=paste("nonDEN",rownames(mgs.med.vect)[i],mgs_list[i,2],".jpeg",sep="-")
      jpeg(file=filename.nonDEN)
      print(ggboxplot(metadata_nonDEN, x = "labels", y = "mgs", 
                      color = "labels",palette =c("#00AFBB", "#E7B800"), add = "jitter") +
              stat_compare_means(method = "wilcox.test", label.x = 1.3, label= "p") + 
              labs(x = "nonDEN_Fasted - nonDEN_Fed", y = paste(mgs_list[i,2])))
      dev.off()
    } else{}
  }else{}
}
print(n.nonDEN)
p.nonDEN <- as.matrix(p.nonDEN)
rownames(p.nonDEN) <- rownames(mgs.med.vect)


