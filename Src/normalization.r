if ( chip_type %in% c( "hgu133a", "hgu133plus2" ) ){
  
  eset = threestep( raw_data, background.method = "GCRMA", normalize.method="quantile", summary.method="median.polish")
  #eset = threestep( raw_data, background.method = "RMA.2", normalize.method="quantile", summary.method="median.polish")
  
} else if ( chip_type %in% c( "pd.huex.1.0.st.v2" ) ){
  
  eset = oligo::rma( raw_data, target = 'extended', normalize = T)
  message("Normalization finished, loading annotation information")
  featureData( eset ) = getNetAffx( eset, 'transcript') 
  
} else if ( chip_type %in% c( "pd.hugene.2.0.st" ) ){
  
  eset = oligo::rma( raw_data, target = "core", normalize = T)
  message("Normalization finished, loading annotation information")
  featureData( eset ) = getNetAffx( eset, 'transcript') 
  
} else {
  
  message(c("Unknown Chip Type: ",chip_type))
  stop()
}

if (! quality_control_only  ){
  
  if (! zipped ){
    eset = eset[  , match( base::gsub( c(".gz|.CEL|.cel|.GZ"), "", names(cohorts_vec) ), base::gsub( c(".gz|.CEL|.cel|.GZ"), "", colnames(eset) ), nomatch = 0 ) ] # development
  } else {
    eset = eset[  , match( base::gsub( c(".gz|.CEL|.cel|.GZ"), "", names(cohorts_vec) ), base::gsub( c(".gz|.CEL|.cel|.GZ"), "", colnames(eset) ), nomatch = 0 ) ] # development
  }
} else {
  eset = eset
}

# unify the names
rownames(pData(eset)) = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset) ) ) # the reason is to assure, that no .gz ending is present when we add it
mapping_cohort_p = match( base::gsub( c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset) ) ), base::gsub( c(".gz|.CEL|.cel|.GZ"), "", names(cohorts_vec) ), nomatch = 0 )
mapping_cohort_c = match( base::gsub( c(".gz|.CEL|.cel|.GZ"), "", names(cohorts_vec) ) , base::gsub( c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset) ) ) , nomatch = 0 )

eset = eset[,rownames(pData(eset)) %in% base::gsub( c(".gz|.CEL|.cel|.GZ"), "", names(cohorts_vec) )]
pData(eset)$Cohort = cohorts_vec[ mapping_cohort_p ]


if( export_eset ){
  
  ensembl_genes = as.character(mget( rownames(eset), hgu133plus2ENSEMBL ))
  uniprot = as.character(mget( rownames(eset), hgu133plus2UNIPROT ))
  omim = as.character(mget( rownames(eset), hgu133plus2OMIM ))
  enzyme = as.character(mget( rownames(eset), hgu133plus2ENZYME ))
  hgnc_genes  = as.character(mget( rownames( eset ) ,hgu133plus2SYMBOL))
  entrez_genes= as.character(mget( rownames( eset ) ,hgu133plus2ENTREZID))
  fdata = data.frame(fData(eset))
  
  exprs_values = data.frame(exprs(eset))
  
  expressionSet_out = data.frame(
    
    "Affy-ID" = fdata,
    "HGNC" = hgnc_genes,
    "Entrez" = entrez_genes,
    "Ensembl" = ensembl_genes,
    "Uniprot" = uniprot,
    "Omim" = omim,
    "Enzyme" = enzyme
    
    )
  
  expressionSet_out = cbind(expressionSet_out, exprs_values)
  
  write.table(expressionSet_out, file = paste( output_path, "Output/ExpressionSet.txt", sep = "/"), sep = "\t")
}
