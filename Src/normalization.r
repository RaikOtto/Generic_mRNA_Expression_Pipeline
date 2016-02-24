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
 
} else if ( chip_type %in% c( "HumanHT-12.v4" ) ){
  
  eset = raw_data
  message("Normalization finished, loading annotation information")
  #featureData( eset ) = getNetAffx( eset, 'transcript') 
  
} else {
  
  stop(c("Unknown Chip Type: ",chip_type))

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
  
  print("Exporting eset data")
  source("Src/annotation.r")
  
  expressionSet_out = data.frame(
    
    "Exprs" = as.character(exprs(eset)),
    "HGNC" = hgnc_genes
    
  )
  
  expressionSet_out = expressionSet_out[ expressionSet_out$HGNC!=""  ,]
  
  write.table(
    expressionSet_out, 
    file = paste( 
      output_path, 
      paste0(
        paste0
          ("Output/ExpressionSet_",
          project_name),".txt"
      ), 
    sep = "/"),
  sep = "\t", row.names = F, quote = F)
  
}
