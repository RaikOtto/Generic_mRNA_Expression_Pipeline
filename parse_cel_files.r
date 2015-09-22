print("Running step 3: Parsing raw CEL files"  )

#celFiles = celFiles[ -( excluded_files )  ] # exclude bad samples
#rawdata = ReadAffy( filenames = celFiles )

if ( chip_type %in% c( "hgu133plus2", "hgu133a" ) ){
  
  if ( zipped ){
    celFiles = list.celfiles( cel_files_path, full =T)
  } else {
    celFiles = list.celfiles( cel_files_path, full =T)
  }

  raw_data = read.affybatch( filenames = celFiles )
  #excluded_files = c()
  
  if( project_name == "Plus2_Runx"  ){
    excluded_files = c( 1,2,3,21,94,115,116,117,118,119,120,121,122,123,124,125,129,138,153,154,155,156,157,158,159,160,161,164,165,174,175,176,179)
    raw_data = raw_data[ -( excluded_files )  ] # exclude bad samples
    phenodata = phenodata[ match( phenodata$ID, colnames(raw_data), nomatch = 0 ),]
  }

} else if ( chip_type %in% c( "pd.huex.1.0.st.v2", "pd.hugene.2.0.st" ) ){
  
  library("oligo")
  celFiles = list.celfiles( cel_files_path, full = T)
  raw_data = read.celfiles( filenames = celFiles )

} else {
  
  print(c("Unknown Chip Type: ",chip_type))
  stop()
}
