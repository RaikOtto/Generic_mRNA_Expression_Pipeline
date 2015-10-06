kocent = ( system('uname -n',intern=T) == "kocent" )

if ( kocent ){
  debug_mode = F
} else {
  debug_mode = T
}

if ( !debug_mode ){
  pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
} else if ( system('uname -n',intern=T) == "MacBook-Jan-Niklas.local" ){
  pipeline_loc = paste( system("echo $HOME",intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )
  pipeline_loc = "/Users/raikotto/Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/"
} else{
  pipeline_loc = "/Users/raikotto/Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/"
}

setwd( pipeline_loc ) # Set the path to where the pipeline is located

source( "Src/pipeline_structure.r" )

default_parameters = T
which_project = "ovarian"
source("project_files.r")

###

create_cohorts    = T # 2
parse_files       = T # 3
normalize         = T # 4
qc_control        = F # 5
annotate          = F # 6
absent_analysis   = F # 7
dif_exp_ana       = F # 8
export_results    = F # 9
create_pathways   = F # 10
extract_interest  = F # 11
annotate_tissue_abbundance = F # 12

### generic
stat_design = "contrast"

run_analysis();print( "Finished" )
