#### --- Manacus MKtest experiment -- ####

# here we are annotating the Manacus V3 genome with Pipra V1 names,

# First I BLASTed Manacus_cds.fasta using Pipra_BLAST_annotation.sh
library(seqinr)
library(stringr)
# Then we need to read in the BLAST resutls
blastres<- read.table("Manacus_pfil_blast_v1")
colnames(blastres)<- c("qseqid","sseqid", "pident","length","mismatch","gapopen","qstart","qend","sstart","send", "evalue","bitscore")



#generate a list of Pfil IDs
pfil<-read.fasta("GCF_003945595.1_ASM394559v1_cds_from_genomic.fna")

pfil_cds<- names(pfil)
pfil_gene<- word(getAnnot(pfil),2)
pfil_gene<- sub("(\\[gene=)((.)+)(\\])", "\\2", pfil_gene)
pfil_protein<- getAnnot(pfil)
pfil_protein<- as.vector(do.call("rbind",pfil_protein))
pfil_protein<- sub("(.+ \\[protein=)(.+\\])","\\2", pfil_protein)
pfil_protein<- sub("\\].+","", pfil_protein)
pfil_gene<- as.data.frame(cbind(pfil_cds, pfil_gene, pfil_protein))
pfil_gene$pfil_cds<- sub("lcl\\|", "", pfil_gene$pfil_cds)
pfil_gene$XP<- sub("(.+_XP_)(.+)(\\..+)","\\2", pfil_gene$pfil_cds)

#manakin query

#make a list of the manakin query.
manakin<- read.fasta("../genome/GCF_001715985.3_ASM171598v3_cds_from_genomic.fna")
manakin_query<- names(manakin)
manakin_gene<- word(getAnnot(manakin),2)
manakin_gene<- sub("(\\[gene=)((.)+)(\\])", "\\2", manakin_gene)
manakin_protein<- getAnnot(manakin)
manakin_protein<- as.vector(do.call("rbind",manakin_protein))
manakin_protein<- sub("(.+ \\[protein=)(.+\\])","\\2", manakin_protein)
manakin_protein<- sub("\\].+","", manakin_protein)
manakin_gene<- as.data.frame(cbind(manakin_query, manakin_gene, manakin_protein))
manakin_gene$manakin_cds<- sub("lcl\\|", "", manakin_gene$manakin_query)
manakin_gene$XP<- sub("(.+_XP_)(.+)(\\..+)","\\2", manakin_gene$manakin_cds)





#let's compare the query against what is actually in the fasta
length(unique(manakin_query)) #24624
length(unique(blastres$qseqid)) #24524

#what are the queries that did not make the e-value cut off?
gross_queries<- manakin_query[!manakin_query %in% unique(blastres$qseqid)]
remainder<- manakin[gross_queries]
#write.fasta(sequences=remainder, names=names(remainder), nbchar = 100, file.out="Manacus_query_remainders.fasta")



res_per_query<- as.data.frame(table(blastres$qseqid))
summary(res_per_query[,2])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    4.00   12.00   38.46   31.25 4608.00 

#let's extract the hit with the smallest e-value for all queries.
ques<- unique(blastres$qseqid)
mine<- list()
for(i in ques){
  sub<- blastres[blastres$qseqid==i,]
  gene<- manakin_gene$manakin_gene[manakin_gene$manakin_query==i]
  sub$query_gene<- gene
  qprotein<- manakin_gene$manakin_protein[manakin_gene$manakin_query==i]
  sub$query_protein<- qprotein
  min<- min(sub$evalue)
  sub<- sub[sub$evalue==min,]
  mine[[i]]<- sub
}
mine<- do.call("rbind", mine)
#took the number of hits from 1886967 to 314356

## annotate the blast with the Pipra annotations & Manacus annotations
out<- list()
for(i in mine$sseqid){
  sub<- mine[mine$sseqid==i,]
  gene<- pfil_gene$pfil_gene[pfil_gene$pfil_cds==i]
  sub$pfil_gene<- gene
  sprotein<- pfil_gene$pfil_protein[pfil_gene$pfil_cds==i]
  sub$pfil_protein<- sprotein
  out[[i]]<- sub
}
mine2<- do.call("rbind", out)



res_per_query<- as.data.frame(table(mine2$qseqid))
summary(res_per_query[,2])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    1.00    4.00    7.89    8.00  232.00  

#how many genes does this represent in pfil ?
genes_per_query<- aggregate(.~ qseqid, mine2[,c("qseqid","pfil_gene")], FUN=function(x) length(unique(x)) )
summary(genes_per_query[,2])

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   2.304   1.000 134.000 
#most of the results only have 1 gene as the hit. Which is good. 
nrow(genes_per_query[genes_per_query$pfil_gene>1,])


## now let's take the best hit using a match of the gene name, and failing that the bit score and percent identity. 

multiples<- res_per_query$Var1[res_per_query$Freq>1] #vector of query sequences with multiple hits.
multiple_hits<- mine2[mine2$qseqid %in% multiples,] #subset the blast results database for multiple hits only.
single_hits<- mine2[!mine2$qseqid %in% multiples,]

#then we are going to take the hit that matches the gene name of the query, and if there are none, we will take the max bit score.

#i<- "AM-NW.020893960.1.cds.XP.027588113.1.177.p1"
out<- list()
for(i in multiples){
  sub<- multiple_hits[multiple_hits$qseqid==i,]
  sub$xp_match<- sub$query_gene==sub$pfil_gene
  if(is.element(TRUE, sub$xp_match)){
    sub<- sub[sub$xp_match==TRUE,]
    sub$xp_match<- NULL
    out[[i]]<- sub
  }else{
    sub$xp_match<- NULL
    max<- max(sub$bitscore)
    sub<- sub[sub$bitscore==max,]
    if(nrow(sub)>1){
      max<- max(sub$pident)
      sub<- sub[sub$pident==max,]
      out[[i]]<- sub
    }else{
      out[[i]]<- sub
    }
  }
}
best_hits<- do.call("rbind", out)

length(unique(best_hits$qseqid)) #check to make sure the number of query sequences is the same as the input
length(unique(multiple_hits$qseqid)) 



res_per_query<- as.data.frame(table(best_hits$qseqid))
summary(res_per_query[,2])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   2.000   4.000   7.495   8.000 158.000 
genes_per_query<- aggregate(.~ qseqid, best_hits[,c("qseqid","pfil_gene")], FUN=function(x) length(unique(x)) )
summary(genes_per_query[,2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   1.001   1.000   7.000 

#extract the single hits from the above process and add them to the pile. 

multiples<- res_per_query$Var1[res_per_query$Freq>1] #vector of query sequences with multiple hits.
multiple_hits<- best_hits[best_hits$qseqid %in% multiples,] #subset the blast results database for multiple hits only.
single_hits2<- best_hits[!best_hits$qseqid %in% multiples,]
#nrow(single_hits)+nrow(single_hits2) + length(multiples) #sanity check.
single_hits<- rbind(single_hits, single_hits2)
remove(list= c("single_hits2"))

#i<- "AM-NW.020896602.1.cds.XP.027602090.1.30221.p1"


out<- list()
for(i in multiples){
  sub<- multiple_hits[multiple_hits$qseqid==i,]
  if(is.element(TRUE, duplicated(sub$pfil_gene))){
    sub<- sub[!duplicated(sub$pfil_gene),]
    out[[i]]<- sub
  }else{
    out[[i]]<- sub
  }
}
best_hits2<- do.call("rbind", out)
length(unique(best_hits2$qseqid)) 

res_per_query<- as.data.frame(table(best_hits2$qseqid))
summary(res_per_query[,2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   1.005   1.000   7.000 
genes_per_query<- aggregate(.~ qseqid, best_hits2[,c("qseqid","pfil_gene")], FUN=function(x) length(unique(x)) )
summary(genes_per_query[,2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   1.001   1.000   7.000 

#extract the single hits from the above process and add them to the pile. 

multiples<- res_per_query$Var1[res_per_query$Freq>1] #vector of query sequences with multiple hits.
multiple_hits<- best_hits2[best_hits2$qseqid %in% multiples,] #subset the blast results database for multiple hits only.
single_hits2<- best_hits2[!best_hits2$qseqid %in% multiples,]
#nrow(single_hits)+nrow(single_hits2) + length(multiples) #sanity check.
single_hits<- rbind(single_hits, single_hits2)

#i<- "lcl|NW_024103432.1_cds_XP_027571377.2_4"
out<- list()
for(i in multiples){
  sub<- multiple_hits[multiple_hits$qseqid==i,]
  sub$desc_match<- sub$query_protein==sub$pfil_protein
  if(is.element(TRUE, sub$desc_match)){
    sub<- sub[sub$desc_match==TRUE,]
    sub<- sub[!duplicated(sub$pfil_protein),]
    sub$desc_match<- NULL
    out[[i]]<- sub
  }else{
    sub$desc_match<- NULL
    sub<- sub[!duplicated(sub$pfil_protein),]
    out[[i]]<- sub
  }
}
best_hits3<- do.call("rbind", out)
length(unique(best_hits3$qseqid)) 



res_per_query<- as.data.frame(table(best_hits3$qseqid))
summary(res_per_query[,2])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.000   1.000   1.000   1.308   1.500   5.000 
genes_per_query<- aggregate(.~ qseqid, best_hits3[,c("qseqid","pfil_gene")], FUN=function(x) length(unique(x)) )
summary(genes_per_query[,2])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   1.308   1.500   5.000 

#extract the single hits from the above process and add them to the pile. 

multiples<- res_per_query$Var1[res_per_query$Freq>1] #vector of query sequences with multiple hits.
multiple_hits<- best_hits3[best_hits3$qseqid %in% multiples,] #subset the blast results database for multiple hits only.
single_hits2<- best_hits3[!best_hits3$qseqid %in% multiples,]
#nrow(single_hits)+nrow(single_hits2) + length(multiples) #sanity check.
single_hits<- rbind(single_hits, single_hits2)

#there are 23 queries that have multiples. However, these seem roughly equivalent in many ways. Therefore I am going to just take the 2nd hit...

out<- list()
for(i in multiples){
  sub<- multiple_hits[multiple_hits$qseqid==i,]
  sub<- sub[2,]
  out[[i]]<- sub
}
best_hits4<- do.call("rbind", out)
length(unique(best_hits4$qseqid)) 

single_hits<- rbind(single_hits, best_hits4)
length(unique(blastres$qseqid)) #sanity check the nrow of object above matches the query length

#NOW let's add on those 148 genes that were not in the BLAST results file. 

remainder_queries<- as.data.frame(cbind(gross_queries, matrix(data=NA, nrow=length(gross_queries), ncol=length(colnames(single_hits))-1)))
colnames(remainder_queries)<- colnames(single_hits)

out<- list()
for(i in  gross_queries){
  sub<- remainder_queries[remainder_queries$qseqid==i,]
  gene<- manakin_gene$manakin_gene[manakin_gene$manakin_query==i]
  sub$query_gene<- gene
  sprotein<- manakin_gene$manakin_protein[manakin_gene$manakin_query==i]
  sub$query_protein<- sprotein
  out[[i]]<- sub
}
remainder_queries2<- do.call("rbind", out)



final_manakin_annotation<- rbind(remainder_queries2, single_hits)
write.csv(final_manakin_annotation, file = "manacus_pfil_v1_blast_1q1r.csv", row.names=FALSE)

#now let's take the gene names

res_per_query<- as.data.frame(table(final_manakin_annotation$query_gene))
summary(res_per_query[,2])

genes<- unique(final_manakin_annotation$query_gene)
out<- list()
for(i in genes){
  sub<- final_manakin_annotation[final_manakin_annotation$query_gene==i,]
  if(nrow(sub)==1){
  out[[i]]<- sub
  } else if (length(unique(sub$pfil_gene))==1){
  sub<- sub[!duplicated(sub$query_gene),]
  out[[i]]<- sub
  }else{
  sub<- sub[!is.na(sub$pfil_gene),]
  sub<- sub[!duplicated(sub$query_gene),]
  out[[i]]<- sub
  }
}
single_genes<- do.call("rbind", out)

write.csv(single_genes, file="manacus_pfil_v1_blast_gene_results.csv", row.names=FALSE)