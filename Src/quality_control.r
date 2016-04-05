suppressMessages(library("arrayQualityMetrics"))
library("affyQCReport")

if ( !( quality_control_only )  ){
  
  mapping = match( colnames(eset),  phenodata$ID )
  
  for ( elem in colnames( phenodata ) ){
    
    if ( (elem != "ID"  ) && ( ! ( elem %in% colnames( pData( eset ) ) )  ) ){ 
      
      index_elem = which( colnames( phenodata ) == elem )
      pData( eset )$elem = phenodata[ mapping, index_elem ]
      colnames( pData( eset ) )[ length( colnames( pData( eset ) ) ) ] = elem
    }
  }
  pData( eset )$Cohort = cohorts_vec[ match( rownames(pData( eset ) ), names(cohorts_vec) ) ]
  
} else {
  
  mapping_eset_pheno = match( as.character( colnames( eset ) ), as.character( phenodata$ID ) )
  pData( eset  )$Group = phenodata$Group[ mapping_eset_pheno  ]
}

pData( raw_data ) = pData( eset )

message( "Running Quality Metrics"  )
QCReport(raw_data, file = "~/Documents/AffyQCReport.pdf")
arrayQualityMetrics( raw_data, intgroup = c("Cohort","Group"), outdir = quality_report_path, force = T, showWarnings = F, do.logtransform = T)
#arrayQualityMetrics( raw_data, intgroup = "Cohort", outdir = quality_report_path, force = T, showWarnings = F)
#arrayQualityMetrics( raw_data, force = T)