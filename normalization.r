print("Running step 4: Normalizing"  )

if ( chip_type %in% c( "hgu133a", "hgu133plus2" ) ){
  
  eset = threestep( raw_data, background.method = "GCRMA", normalize.method="quantile", summary.method="median.polish")
  #eset = threestep( raw_data, background.method = "RMA.2", normalize.method="quantile", summary.method="median.polish")
  
} else if ( chip_type %in% c( "pd.huex.1.0.st.v2" ) ){

  eset = oligo::rma( raw_data, target = 'extended', normalize = T)
  print("Normalization finished, loading annotation information")
  featureData( eset ) = getNetAffx( eset, 'transcript') 
  
} else if ( chip_type %in% c( "pd.hugene.2.0.st" ) ){
  
  eset = oligo::rma( raw_data, target='extended', normalize = T)
  print("Normalization finished, loading annotation information")
  featureData( eset ) = getNetAffx( eset, 'transcript') 
  
} else {
  
  print(c("Unknown Chip Type: ",chip_type))
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

#pData(eset) = phenodata
if (! zipped ){
  mapping_cohort_p = match( rownames(pData(eset)), names(cohorts_vec) , nomatch = 0 )
  mapping_cohort_c = match( names(cohorts_vec), rownames(pData(eset)) , nomatch = 0 )
  
} else{
  
  rownames(pData(eset)) = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset) ) ) # the reason is to assure, that no .gz ending is present when we add it
  mapping_cohort_p = match( base::gsub( c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset) ) ), base::gsub( c(".gz|.CEL|.cel|.GZ"), "", names(cohorts_vec) ), nomatch = 0 )
  mapping_cohort_c = match( base::gsub( c(".gz|.CEL|.cel|.GZ"), "", names(cohorts_vec) ) , base::gsub( c(".gz|.CEL|.cel|.GZ"), "", rownames(pData(eset) ) ) , nomatch = 0 )
  
}

pData(eset)$Cohort = cohorts_vec[ mapping_cohort_p ]
phenodata$ID = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", phenodata$ID )
