cmd_args = commandArgs(trailingOnly=TRUE);
options(max.print=200)

input = unlist(strsplit(cmd_args," ")[1])
output = unlist(strsplit(cmd_args," ")[2])

t = read.table(sep="\t",header=T,file=input)
t = t[order(t$k,decreasing = TRUE),]

#if (k>0){
#val = phyper( k-1, N, M-N, n, lower.tail=FALSE)
#}else{
#val = sum(dhyper( seq(0,100), N, M-N, n))
#}
#cat(as.numeric(val)); invisible(val)

write.table( t, file=output, sep="\t",row.names=F, col.names=F, quote=F)
