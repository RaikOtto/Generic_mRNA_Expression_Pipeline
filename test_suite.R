message("Running Test-Suite.")

kocent = ( system( 'uname -n', intern = T ) == "kocent" )

if ( kocent ){
  pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
} else if ( system( 'uname -n', intern = T ) == "MacBook-Jan-Niklas.local" ){
  pipeline_loc = paste( system( "echo $HOME", intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )
} else{
  pipeline_loc = "/Users/raikotto/Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/"
}

default_parameters = T
which_project = "test_case"
source("project_files.r")

setwd( pipeline_loc ) # Set the path to where the pipeline is located

### generic
stat_design = "contrast"

capture.output( source( "Src/set_generic_initial_parameters.r" ), file = 'NUL' )
message( "set generic initial parameters - succesfull" )

capture.output( source( "Src/parse_cel_files.r" ), file = 'NUL' )
message( "Parsing cel files - successfull")

capture.output( source( "Src/cohort_creation.r" ), file = 'NUL' )
message( "Cohort creation - succesfull")

capture.output( source( "Src/normalization.r" ), file = 'NUL' )
message ( "Normalization - successfull")

#capture.output( source( "Src/quality_control.r" ), file = 'NUL' )
#message( "Quality control - succesfull")

capture.output( source( "Src/annotation.r" ), file = 'NUL' )
message( "Annotation - succesfull")

capture.output( source( "Src/absent_genes_analysis.r" ), file = 'NUL' )
message( "Absent genes analysis - succesfull" )

capture.output( source( "Src/differential_expression.r" ), file = 'NUL' )
message( "Differential expression analysis - succesfull" )

capture.output( source( "Src/result_preparation.r" ), file = 'NUL' )
message( "Result preparation - succesfull" )

capture.output( source( "Src/create_pathway_maps.r" ), file = 'NUL' )
message( "Create pathway maps - succesfull" )

capture.outpur( source( "Src/Extract_interesting_entities.r" ), file = 'NUL' )
message( "Extract interesting entities - succesfull" )

capture.output( source( "Src/annotate_tissue_abbundance.r" ), file = 'NUL' )
message( "Annotate tissue abbundance - succesfull" )