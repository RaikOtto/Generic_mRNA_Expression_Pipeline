library(VariantAnnotation)
library("stringr")
require(graphics)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library("devtools")
library("xlsx")
#library("GMD")
source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")

vcf_filfe = "~/Dropbox/PhD/Klinghammer/upload/VCF/all.vcf"

split_fun = function(x){ 
  
  split = as.character( unlist( str_split( x, "_"  ) ) )
  
  if ( length(split) >= 2 ){
    
    res = tail( unlist( split ), 1 )
    
  } else {
    
    res = ""
  }
  
  return( res )
}

vcf = readVcf(file =  vcf_filfe,"hg19")

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
loc  = locateVariants(query = rowRanges(vcf), subject = txdb,region = CodingVariants())
splt = split(mcols(loc)$GENEID, mcols(loc)$QUERYID)
table(sapply(splt, function(x) length(unique(x)) > 1))
splt <- split(mcols(loc)$QUERYID, mcols(loc)$GENEID)
res =  sapply(splt, function(x) length(unique(x)))
gene_ids = names( res )

library( "biomaRt" )

ensembl     = useMart("ensembl",dataset="hsapiens_gene_ensembl")
entrez_ids  = getBM( attributes = c( "entrezgene", "hgnc_symbol" ), values = unique( gene_ids), filters = "entrezgene" , mart = ensembl)
mapping = match( names( res ), entrez_ids$entrezgene  )
names(res)[mapping] = entrez_ids$hgnc_symbol

#g = summarizeVariants(subject = vcf, query = txdb,  mode = CodingVariants())

vcf_rewrite <- setClass(
 
  "vcf_rewrite",
  
  slots = c(
    data = "matrix",
    type = "list",
    gene_symbols = "list",
    aa = "list"
  ),
  
  prototype = list(
    
    data         = matrix( as.character()  ),
    type         = list(),
    aa           = list(),
    gene_symbols = list()
  )
)

create_overview_vcf = function( vcf_obj, exon_only = T ){
  
  GT           = geno( vcf_obj )$GT
  GT           = apply( GT, FUN=function(x){ x[x=="."]=0; return(x)  }, MARGIN = 1 )
  GT           = t(GT)
  nr_variants  = dim( GT )[ 1 ]
  nr_samples   = dim( GT )[ 2 ]
  FC_list      = as.character( unlist( info( vcf_obj )$FC ) )
  type_list    = unlist( lapply( FC_list, FUN = function( x ){ return( head( unlist( str_split( x, "_"  ) ), 1 ) ) } ) )
  aa_list      = unlist( lapply( FC_list, FUN = split_fun ) )
  re           = vcf_rewrite( data = matrix( as.character(), ncol = nr_samples ) )
  colnames( re@data ) = colnames( vcf_obj )
  rownms       = c()
  
  vcf_genes = unique( as.data.frame( info( vcf )$GI ) )[,3]
  
  for ( i in 1:nr_variants  ){
    
    #exon_member = exon_only == unlist( info( vcf )$EXON )[ i ]
    exon_member = T
    #interesting_alteration = !( FC_list[i] %in% c("Silent","Synonymous","Noncoding") )
    interesting_alteration = T
    
    if (exon_member && interesting_alteration){
      
      gt_row = str_replace_all( as.character( GT[ i, ] ), "\n", "" )
      gt_row[ gt_row != "0" ] = 1
      
      re@data = rbind( re@data, matrix( gt_row, ncol = nr_samples  )  )
      re@aa   = c( re@aa, aa_list[ i ] )
      re@type = c( re@type, type_list[ i ] )
      re@gene_symbols = c( re@gene_symbols, as.character( vcf_genes[ i ] ) )
      rownms = c(rownms,rownames(geno( vcf )$GT)[i])
    }
  }
  
  rownames(re@data) = rownms
  return( re )
}

re = create_overview_vcf( vcf )
data = matrix(as.numeric( re@data  ), ncol = dim(re@data)[2] )
colnames(data) = colnames(re@data)

gene_list = as.character( unique( unlist( re@gene_symbols) ) )
mapping = match( unlist( re@gene_symbols ), gene_list)
gene_overview = matrix( as.numeric(), ncol = dim(data)[2] )

for ( i in 1:length(gene_list) ){
  
  sub_map = which( mapping == i  )
  sub_data = matrix( data[ sub_map, ], ncol = dim(data)[2] )
  res_gene = apply( sub_data, FUN = sum, MARGIN = 2 )
  
  gene_overview = rbind( gene_overview, res_gene  )
}

colnames( gene_overview ) = colnames(data)
rownames( gene_overview ) = gene_list
gene_overview_vis = apply( gene_overview, FUN = function(x){ x[x!=0] = 1; return(x)  }, MARGIN=2)

mycols = colorRampPalette(colors = c("white","yellow","orange","red"))

setwd("~/Dropbox/PhD/Klinghammer/upload/Cel_files/")
cohorts_file = "cohorts.tab" 
main_drug_type = "Docetaxel"
drug_data = read.table( file =  cohorts_file, sep = "\t", header = T)

mapping = match( colnames(data), as.character( drug_data$VCF_id_2), nomatch = 0)

colnames( gene_overview ) = str_replace( string = unlist(as.character(drug_data$ID))[ mapping ] , ".CEL", "")

drugs = c("X5.FU","Carboplatin","Docetaxel","Cetuximab","Everolimus")

if (exists("mycols_drug_mat"))
  remove(mycols_drug_mat)

for (drug_type in drugs ){
  
  pre_select = unlist(as.character(drug_data[mapping,which(colnames(drug_data)==drug_type)]))
  drug_response = as.numeric(str_replace(string = pre_select, ",","."))
  
  mycols_drug = colorRampPalette(colors = c("green","blue"))( length(drug_response) )
  mycols_drug[ is.na(drug_response)  ] = "white"
  
  if (drug_type == main_drug_type){
    
    gene_overview = gene_overview[ ,order( drug_response ) ]
    
  } else{
    
    mycols_drug = mycols_drug[ order(order( drug_response ) )  ]
  }
  
  if (!exists("mycols_drug_mat")){
    
    mycols_drug_mat = matrix( mycols_drug, nrow = length(mapping) )
    colnames(mycols_drug_mat) = c(drug_type)
    
  } else {
    
    mycols_drug_mat       = cbind( mycols_drug_mat, mycols_drug )
    colnames(mycols_drug_mat)[dim(mycols_drug_mat)[2]] = drug_type 
  }
}
dev.off()
heatmap.3( 
  main= "Altertation frequency HDNSC",
  #gene_overview,
  gene_overview_vis,
  col = mycols,
  #srtCol = 40,
  #adjCol = c(1,0),
  trace = "none",
  ColSideColors = mycols_drug_mat,
  ColSideColorsSize = dim(mycols_drug_mat)[1],
  Colv = F,
  dendrogram = "row"
)

### cor tests 

find_cor_drug = function( geno_vec, drug_resp ){
  
  #return( round( cor( geno_vec, drug_resp, use = "pairwise.complete.obs", method = "spearman"  ), 2 )  )
  return( round( cor( geno_vec, drug_resp, use = "pairwise.complete.obs", method = "pearson"  ), 2 )  )
}

if ( exists("cors") )
  remove(cors)

for (drug_type in drugs ){
  
  pre_select = unlist(as.character(drug_data[mapping,which(colnames(drug_data)==drug_type)]))
  drug_resp = as.numeric( str_replace( pre_select, ",", "." ) )
  #drug_resp = scale(drug_resp)

  res = apply( gene_overview, MARGIN = 1, FUN = find_cor_drug, drug_resp)

  if (!exists("cors")){
    
    cors = matrix( res, nrow = dim(gene_overview)[1]  )
    colnames(cors) = c(drug_type)
    
  } else {
    
    cors = cbind( cors, res )
    colnames(cors)[ dim( cors )[ 2 ] ] = drug_type
  }
}

rownames(cors)[ is.na(rownames(cors)) ] = ""

#write.xlsx( x= cors, file =  "/Users/raikotto//Dropbox/PhD/Klinghammer/results/drug_response_correlation_pearson.xlsx")
write.xlsx( x= cors, file =  "/Users/raikotto//Dropbox/PhD/Klinghammer/results/drug_response_correlation_spearman.xlsx")
