# redo this here

print( "Running step 2: Cohort construction"  )

source("Src/set_generic_initial_parameters.r")

if (exists("eset")){
  
  if (zipped){
    p_data = paste0( phenodata$ID, ".gz")
  } else {
    p_data = phenodata$ID
  }
  
}

rownames( pData( raw_data ) ) = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", rownames( pData(raw_data) ) ) # the reason is to assure, that no .gz or .CEL ending is present when we add it
colnames(raw_data) = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", colnames(raw_data) ) 
phenodata$ID = base::gsub( c(".gz|.CEL|.cel|.GZ"), "", phenodata$ID ) # the reason is to assure, that no .gz ending is present when we add it
#print( paste( "Did not find following sampels in raw data, but in cohorts file: ", paste( !( phenodata$ID[phenodata$ID %in% rownames( pData( raw_data ) ) ] ) , sep =" " ) ) ) 
phenodata = phenodata[ phenodata$ID %in% rownames( pData( raw_data ) ), ]
cohorts_vec = as.vector( phenodata[ , which( colnames( phenodata ) == cohorts_type  ) ] )

if ( !( quality_control_only )  ){
  
  cohorts_vec[ cohorts_vec %in% set_ctrl ] = "CTRL"
  cohorts_vec[ cohorts_vec %in% set_case ] = "CASE"
  
  index_cohorts_vec = which( cohorts_vec %in% c("CTRL","CASE") == T  )
  cohorts_vec = cohorts_vec[ cohorts_vec %in% c("CTRL","CASE") ]
  names( cohorts_vec ) = as.character( phenodata$ID[ index_cohorts_vec ] )
  
  nmbr_samples = sum(cohorts_vec %in% c("CTRL","CASE"))
  
  if (stat_design == "intercept") {
    
    design = model.matrix( ~  cohorts_vec )  # intercept
    
  }else{
    
    design = model.matrix( ~ 0 + cohorts_vec ) # contrast
  }
  
  colnames( design )[ colnames( design ) == paste0( "cohorts_vec","CTRL") ] = "CTRL"
  colnames( design )[ colnames( design ) == paste0( "cohorts_vec","CASE") ] = "CASE"
  
}

index_ctrl = as.integer( which( design[ ,colnames(design) == "CTRL" ] == 1 ) )
index_case = as.integer( which( design[ ,colnames(design) == "CASE" ] == 1 ) )

#if (exists("raw_data")){
if( ! ("Group" %in% colnames(pData( raw_data ) ) ) ) {
  
  raw_data_group_vec = rep("",dim( pData(raw_data) )[1] )
  raw_data_group_vec[  index_ctrl ] = "CTRL"
  raw_data_group_vec[  index_case ] = "CASE"
  pData(raw_data) = cbind( pData(raw_data), raw_data_group_vec )
  colnames(pData(raw_data))[-1] = "Group"
  #raw_data = raw_data[,pData(raw_data)$Group!=""]
}
#}