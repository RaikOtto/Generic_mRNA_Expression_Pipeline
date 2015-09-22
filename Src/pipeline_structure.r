
run_analysis = function(){
  
  if ( run_generic ){
    
    source("src/set_generic_initial_parameters.r")
  }
  
  ### Part Parse Cel Files
  
  if ( parse_files ){
    
    source("src/parse_cel_files.r")
  }
  
  ### Part cohort creation
  
  if ( create_cohorts ){
    
    source("src/cohort_creation.r")
  }
  
  ### Part Normalization
  
  if ( normalize ){
    if ( use_frma_normalization ){
      source("src/frma_normalization.r")
    }else{
      source("src/normalization.r")
    }
  }
  
  # Quality Control
  
  if ( qc_control ){
    
    source("src/quality_control.r")
  }
  
  ### Part annotation
  
  if ( annotate )
    
    source("src/annotation.r")
  
  ### Part Excluding non-present Genes
  
  if (absent_analysis)
    
    source("src/absent_genes_analysis.r")
  
  ### Part Differential Expression
  
  if ( dif_exp_ana )
    
    source("src/differential_expression.r")
  
  ### Part Analysis
  
  #source("")
  
  if ( heatmap_vis  )
    source( "src/heatmap_probes_of_interest.r" )
  
  ### Part Visualization
  
  if ( export_results )
    
    source( "src/result_preparation.r" )
  
  if ( create_pathways )
    
    if (! time_series  )
      
      source("src/create_pathway_maps.r")
  
  if (extract_interest)
    
    source("src/Extract_interesting_entities.r")
  
  if (annotate_tissue_abbundance)
    
    source("src/annotate_tissue_abbundance.r")
}