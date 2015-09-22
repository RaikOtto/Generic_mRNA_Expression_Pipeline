
run_analysis = function(){
  
  if ( run_generic ){
    
    source("set_generic_initial_parameters.r")
  }
  
  ### Part cohort creation
  
  if ( create_cohorts ){
    
    source("cohort_creation.r")
  }
  
  ### Part Parse Cel Files
  
  if ( parse_files ){
    
    source("parse_cel_files.r")
  }
  
  ### Part Normalization
  
  if ( normalize ){
    if ( use_frma_normalization ){
      source("frma_normalization.r")
    }else{
      source("normalization.r")
    }
  }
  
  # Quality Control
  
  if ( qc_control ){
    
    source("quality_control.r")
  }
  
  ### Part annotation
  
  if ( annotate )
    
    source("annotation.r")
  
  ### Part Excluding non-present Genes
  
  if (absent_analysis)
    
    source("absent_genes_analysis.r")
  
  ### Part Differential Expression
  
  if ( dif_exp_ana )
    
    source("differential_expression.r")
  
  ### Part Analysis
  
  #source("")
  
  if ( heatmap_vis  )
    source( "heatmap_probes_of_interest.r" )
  
  ### Part Visualization
  
  if ( export_results )
    
    source( "result_preparation.r" )
  
  if ( create_pathways )
    
    if (! time_series  )
      
      source("create_pathway_maps.r")
  
  if (extract_interest)
    
    source("Extract_interesting_entities.r")
  
  if (annotate_tissue_abbundance)
    
    source("annotate_tissue_abbundance.r")
}