kocent = ( system( 'uname -n', intern = T ) == "kocent" )

if ( kocent ){
  
  pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
  
} else {
  
  pipeline_loc = paste( system( "echo $HOME", intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )
  
}

setwd( pipeline_loc )

output_path = paste( pipeline_loc, "Project_files/ag_na_GSEA/Output", sep = "/" )
dir.create(output_path)

GSEA.program.location = paste(pipeline_loc, "GSEA.1.0.R", sep ="/")
source(GSEA.program.location, verbose=T, max.deparse.length=9999)

GSEA.Analyze.Sets(
  directory           = output_path,        # Directory where to store output and results (default: "")
  topgs = 20,                                                           # number of top scoring gene sets used for analysis
  height = 16,
  width = 16
)