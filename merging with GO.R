#merging the MK test results with Pfil annotations and GO terms.
library(plyr)
library(seqinr)
library(stringr)



mk<- read.csv("../MKT_V5_PfvMc.csv") #13826
blast<- read.csv("manacus_pfil_v1_blast_gene_results.csv") #14497
blast[which(blast$query_gene=="CELF4"),]


pfil<- read.csv("../../Manakins/Data/analysis_code/Gene Ontology/Pipra Annotations/Maggies_annotations_modifiedR.csv")  #547887
pfil<- rename(pfil, replace=c("GeneID"="pfil_gene"))
pfil<- pfil[,c("pfil_gene", "best_anno")]

blast<- rename(blast, replace=c("query_gene"="gene_name"))
blast<- subset(blast, select=c("gene_name","pfil_gene"))
mk$gene_name<- sub("_circ", "", mk$gene_name)
mk2<- merge(blast, mk, by="gene_name")

mk3<- merge(mk2, pfil, by="pfil_gene", all.x=TRUE)

write.csv(mk3, file="mkt_v5_PfvMc_annot.csv", row.names=FALSE)

go<- read.csv("../../Manakins/Data/analysis_code/Gene Ontology/Pipra Annotations/pfil_GO_key_BP.csv")
go<- rename(go, replace=c("gene"="pfil_gene"))

mk4<- merge(mk3, go, by="pfil_gene", all.x=TRUE)
mk4<- mk4[,-24]
write.csv(mk4, file="mkt_v5_PfvMc_annot_go.csv", row.names = FALSE)





