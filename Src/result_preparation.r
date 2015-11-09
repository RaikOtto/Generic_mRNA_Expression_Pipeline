library("WriteXLS")

if ( chip_type == "hgu133plus2" ){
  
  library("hgu133plus2.db")
  
  probe_ids = row.names( topall )
  
  hgnc_genes  = as.character(mget( rownames( topall ) ,hgu133plus2SYMBOL))
  hgnc_names  = as.character(mget( rownames( topall ) ,hgu133plus2GENENAME))
  entrez_genes= as.character(mget( rownames( topall ) ,hgu133plus2ENTREZID))
  pathway     = as.character(mget( rownames( topall  ) ,hgu133plus2PATH))
  
  topall_res = data.frame(
    
    "logFC"               = topall$logFC,
    "P_Value"             = topall$P.Val,
    "HGNC_symb"           = hgnc_genes,
    "HGNC_name"           = str_replace_all(hgnc_names,",",";"),
    "entrez"              = entrez_genes,
    "pathway"             = pathway
  )
  
  row.names( topall_res ) = probe_ids
  
  topall_res = topall_res[ order(topall_res$logFC, decreasing = T)  ,]
  
} else if ( chip_type == "hgu133a" ) {

  library("hgu133a.db")
  
  hgnc_genes  = as.character(mget( rownames( topall ) ,hgu133aSYMBOL))
  hgnc_names  = as.character(mget( rownames( topall ) ,hgu133aGENENAME))
  entrez_genes= as.character(mget( rownames( topall ) ,hgu133aENTREZID))
  pathway     = as.character(mget( rownames( topall  ) ,hgu133aPATH))
  
  topall_res = data.frame(
    
    "logFC"               = topall$logFC,
    "P_Value"             = topall$P.Val,
    "HGNC_symb"           = hgnc_genes,
    "HGNC_name"           = str_replace_all(hgnc_names,",",";"),
    "entrez"              = entrez_genes,
    "pathway"             = pathway
  )
  
  topall_res = topall_res[ order(topall_res$logFC, decreasing = T)  ,]
  
} else if ( chip_type %in% c( "pd.hugene.2.0.st", "pd.huex.1.0.st.v2" ) ){
  
  if ( ! exists("index_case"))
    source("Src/annotation.r")
  
  probe_ids = rownames( topall )
  
  index_probes = match( rownames( topall ), rownames(eset), nomatch = 0 )
  exprs_case = rowMeans( exprs( eset )[ index_probes, index_case ] )
  exprs_ctrl = rowMeans( exprs( eset )[ index_probes, index_ctrl ] )
  
  split_fun = function( entry, pos ){ res = unlist( str_split( entry, " // " ) ); if (length(res) > 1){ return( res[pos] ) } else{ return( "" ) } }
  #map=which( rownames(eset) %in% topall$probesetid )
  hgnc_symbols = str_trim( unlist( lapply( topall$geneassignment, FUN=split_fun, 2 ) ) )
  hgnc_names = str_trim( unlist( lapply( topall$geneassignment, FUN=split_fun, 3 ) ) )
  
  #library("biomaRt")

  #ensembl     = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  #entrez  = getBM( attributes = c( "entrezgene", "hgnc_symbol" ), values = unique( hgnc_symbols), filters = "hgnc_symbol" , mart = ensembl)[1]

  #hgnc_map    = match( hgnc_symbols, entrez_ids$hgnc_symbol  , nomatch = 0 )
  #entrez      = hgnc_symbols
  #entrez[ entrez != ""  ] = ""
  #entrez[ which( hgnc_symbols %in%  entrez_ids$hgnc_symbol ) ] = entrez_ids$entrezgene[ hgnc_map ]
  #entrez[ is.na(entrez) ] = ""
    
  topall_res = data.frame(
    "logFC"               = round( topall$logFC,2 ),
    "expr_ctrl"           = round( exprs_ctrl, 2  ),
    "expr_case"           = round( exprs_case, 2  ),
    "P_Value"             = topall$P.Value,
    "HGNC_symb"           = hgnc_symbols,
    "HGNC_names"          = hgnc_names,
    "Probe_ids"           = probe_ids,
    #"entrez"              = entrez,
    #"gene_information"    = topall$geneassignment
  )
  
  row.names( topall_res ) = probe_ids
  topall_res = topall_res[ order( topall_res$logFC, decreasing = T )  ,]
  
  if ( filter_topall_res ){
    topall_res = topall_res[which( topall_res$HGNC_symb != "" ), ]
    topall_res = topall_res[-which( grepl( "microRNA", topall_res$HGNC_names ) ),]
  }

} else if ( chip_type == "HumanHT-12.v4" ){
  
  probe_ids = rownames( topall )
  
  index_probes = match( rownames( topall ), rownames(eset), nomatch = 0 )
  exprs_case = rowMeans( exprs( eset )[ index_probes, index_case ] )
  exprs_ctrl = rowMeans( exprs( eset )[ index_probes, index_ctrl ] )
  
  hgnc_symbols = topall$Symbol
  
  topall_res = data.frame(
    
    "ID"                  = probe_ids
    "logFC"               = round( topall$logFC,2 ),
    "expr_ctrl"           = round( exprs_ctrl, 2  ),
    "expr_case"           = round( exprs_case, 2  ),
    "P_Value"             = topall$P.Value,
    "adj.P.Value"         = topall$adj.P.Value,
    "HGNC_symb"           = hgnc_symbols,
    "entrez"              = topall$Entrez_Gene_ID,
    "AveExpr"             = topall$AveExpr,
    "t"                   = topall$t,
    "B"                   = topall$B
  )
  
  row.names( topall_res ) = probe_ids
  topall_res = topall_res[ order( topall_res$logFC, decreasing = T )  ,]

}
  
dir.create( results_file_path, showWarnings = F)
write.xlsx( topall_res, str_replace(str_replace(name_res_file,"~",user_folder),".csv",".xls"), row.names=F )

# 
message( c( "Amount genes higher in Case cohort:", sum( topall_res$logFC >= lfc_exp ) ) )
message( c( "Amount genes lower in Case cohort:" , sum( topall_res$logFC <= lfc_exp ) ) )

