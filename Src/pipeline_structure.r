
run_analysis = function(){
  
  if ( run_generic ){
    
    source("Src/set_generic_initial_parameters.r")
  }
  
  ### Part Parse Cel Files
  
  if ( parse_files ){
    
    source("Src/parse_cel_files.r")
  }
  
  ### Part cohort creation
  
  if ( create_cohorts ){
    
    source("Src/cohort_creation.r")
  }
  
  ### Part Normalization
  
  if ( normalize ){
    if ( use_frma_normalization ){
      source("Src/frma_normalization.r")
    }else{
      source("Src/normalization.r")
    }
  }
  
  # Quality Control
  
  if ( qc_control ){
    
    source("Src/quality_control.r")
  }
  
  ### Part annotation
  
  if ( annotate )
    
    source("Src/annotation.r")
  
  ### Part Excluding non-present Genes
  
  if (absent_analysis)
    
    source("Src/absent_genes_analysis.r")
  
  ### Part Differential Expression
  
  if ( dif_exp_ana )
    
    source("Src/differential_expression.r")
  
  ### Part Analysis
  
  #source("")
  
  if ( heatmap_vis  )
    source( "Src/heatmap_probes_of_interest.r" )
  
  ### Part Visualization
  
  if ( export_results )
    
    source( "Src/result_preparation.r" )
  
  if ( create_pathways )
    
    if (! time_series  )
      
      source("Src/create_pathway_maps.r")
  
  if (extract_interest)
    
    source("Src/Extract_interesting_entities.r")
  
  if (annotate_tissue_abbundance)
    
    source("Src/annotate_tissue_abbundance.r")
}