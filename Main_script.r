debug_mode = F
#user = "raik"
user = "janniklas"

if ( !debug_mode ){
  pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
} else if ( user == "raik" ){
  pipeline_loc = "/Users/raikotto/Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/"
} else{
  pipeline_loc = paste( system("echo $HOME",intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )
}

setwd( pipeline_loc ) # Set the path to where the pipeline is located

source( "Src/pipeline_structure.r" )

<<<<<<< HEAD
default_parameters = T
which_project = "immuno_all"
=======
which_project = "hnsc"
>>>>>>> 52530cf1e295546f46ffcc79cfb4e9aec7b8226f
source("project_files.r")

###

create_cohorts    = T # 2
parse_files       = T # 3
normalize         = T # 4
qc_control        = F # 5
annotate          = T # 6
absent_analysis   = F # 7
dif_exp_ana       = T # 8
export_results    = T # 9
create_pathways   = F # 10
extract_interest  = T # 11
annotate_tissue_abbundance = F # 12

### generic
stat_design = "contrast"

run_analysis();print( "Finished" )
