kocent = ( system( 'uname -n', intern = T ) == "kocent" )

if ( kocent ){
  pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
} else if ( system( 'uname -n', intern = T ) == "MacBook-Jan-Niklas.local" ){
  pipeline_loc = paste( system( "echo $HOME", intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )
} else{
  pipeline_loc = "/Users/raikotto/Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/"
}

setwd( pipeline_loc ) # Set the path to where the pipeline is located

default_parameters = T
which_project = "test_hgu133plus2"
source("project_files.r")


### generic
stat_design = "contrast"

if ( which_project == "test_pd.huex.1.0.st.v2" ){
  
  message( "Running test-case for platform pd.huex.1.0.st.v2." )
  message( "Checking if all required packages can be loaded:" )
  packages = c( "xlsx", "gdata", "RColorBrewer", "limma", "KEGG.db", "pathview", "stringr", "oligo", "oligoData", "pd.huex.1.0.st.v2", "arrayQualityMetrics", "WriteXLS", "genefilter", "plotly" )
  suppressMessages( sapply( packages, require, character.only=TRUE, quietly = TRUE ) )
   
} else if ( which_project == "test_hgu133plus2" ){
  
  message( "Running test-case for platform hgu133plus2." )
  message( "Checking if all required packages can be loaded:" )
  packages = c( "xlsx", "gdata", "RColorBrewer", "limma", "KEGG.db", "pathview", "stringr", "hgu133plus2.db", "affy", "simpleaffy", "affyPLM", "affycoretools", "affyQCReport", "annaffy", "hgu133a.db", "oligoData", "arrayQualityMetrics", "WriteXLS", "genefilter", "plotly" )
  suppressMessages( sapply( packages, require, character.only=TRUE, quietly = TRUE ) )
  
}

invisible(source( "Src/set_generic_initial_parameters.r" ) )
message( "Set generic initial parameters - succesfull" )

capture.output( suppressMessages( source( "Src/parse_cel_files.r" ) ), file = 'NUL' )
message( "Parsing cel files - successfull")

capture.output( source( "Src/cohort_creation.r" ), file = 'NUL' )
message( "Cohort creation - succesfull")

capture.output( suppressMessages( source( "Src/normalization.r" ) ), file = 'NUL' )
message ( "Normalization - successfull")

#capture.output( suppressMessages( source( "Src/quality_control.r" ) ), file = 'NUL' )
#message( "Quality control - succesfull")

capture.output( source( "Src/annotation.r" ), file = 'NUL' )
message( "Annotation - succesfull")

capture.output( source( "Src/absent_genes_analysis.r" ), file = 'NUL' )
message( "Absent genes analysis - succesfull" )

capture.output( suppressMessages( source( "Src/differential_expression.r" ) ), file = 'NUL' )
message( "Differential expression analysis - succesfull" )

capture.output( suppressMessages( source( "Src/result_preparation.r" ) ), file = 'NUL' )
message( "Result preparation - succesfull" )
reference_results = read.table( reference_path, sep = "\t", header = T)
compare_up = topall_res$Probe_ids[1:100] %in% reference_results$ID[1:100]
compare_down = topall_res$Probe_ids[( length(topall_res$Probe_ids) - 99 ):( length(topall_res$Probe_ids) )] %in% reference_results$ID[( length(reference_results$ID) - 99 ):( length(reference_results$ID) )]
if ( FALSE %in% c(compare_up, compare_down) ){
  #pos_up = which( compare_up %in% FALSE )
  #pos_down = which( compare_down %in% FALSE )
  #pos_down = lapply( pos_down, function(x) x + ( length(topall_res$Probe_ids) - 99 ) ) 
  message( "Calculated results differ from reference results." )
} else{
  message( "Calculated results are matching reference results." )
}

capture.output( suppressMessages( source( "Src/create_pathway_maps.r" )), file = 'NUL' )
message( "Create pathway maps - succesfull" )

capture.output( suppressMessages( source( "Src/Extract_interesting_entities.r" )), file = 'NUL' )
message( "Extract interesting entities - succesfull" )

#capture.output( source( "Src/annotate_tissue_abbundance.r" ), file = 'NUL' )
#message( "Annotate tissue abbundance - succesfull" )

message( "Pipeline is working." )