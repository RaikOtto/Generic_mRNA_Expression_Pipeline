library(fdrtool)
options(max.print=200)

filename = commandArgs(TRUE)[1];
t=read.table(file=filename,header=T,sep="\t")

p_values	=	phyper(q=t$q-1,m=t$m,n=t$n,k=t$k,lower.tail=F)
q_values	=	round(p.adjust(p_values,"BH"),4)
p_values	=	round(p_values, 4)

d = with(t,data.frame(p_values,q_values,q,m,n,k,gene_source,pathway,genes))
d = d[order(p_values,decreasing=F),]

write.table(file=filename, d[d$p_values<=0.05,], sep="\t", quote=FALSE,row.names=FALSE,col.names=T)
