kocent = ( system( 'uname -n', intern = T ) == "kocent" )

debug_mode = F

if ( kocent ){
  
  pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
  
} else if ( system( 'uname -n', intern = T ) == "MacBook-Jan-Niklas.local" ){
  
  pipeline_loc = paste( system( "echo $HOME", intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )
  
} else{
  
  pipeline_loc = "/Users/raikotto/Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/"
}

print(c("Kocent:",kocent ))

setwd( pipeline_loc ) # Set the path to where the pipeline is located

default_parameters = T
which_project = "immuMZ"
source("project_files.r")

source( "Src/pipeline_structure.r" )

###

create_cohorts    = T # 2
parse_files       = T # 3
normalize         = T # 4
<<<<<<< HEAD
qc_control        = T # 5
annotate          = T # 6
absent_analysis   = T # 7
dif_exp_ana       = T # 8
=======
qc_control        = F # 5
annotate          = F # 6
absent_analysis   = F # 7
dif_exp_ana       = F # 8
>>>>>>> 9014c1e20cbd2965055c78d99b98d0baae6e8764
export_results    = F # 9
create_pathways   = F # 10
extract_interest  = F # 11
annotate_tissue_abbundance = F # 12

### generic
stat_design = "contrast"

run_analysis();print( "Finished" )
