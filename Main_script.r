kocent = ( system( 'uname -n', intern = T ) == "kocent" )

debug_mode = F

if ( kocent ){

  pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path

} else if ( system( 'uname -n', intern = T ) == "MacBook-Jan-Niklas.local" ){

  pipeline_loc = paste( system( "echo $HOME", intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )

} else{

  pipeline_loc = "/Users/raikotto/Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/"

}

#pipeline_loc = paste( system( "echo $HOME", intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" )
print(c("Kocent:",kocent ))

setwd( pipeline_loc ) # Set the path to where the pipeline is located

#options(error=traceback)
default_parameters = T

which_project = "bm_vs_blood_all"

source("project_files.r")
options(error=traceback)
source( "Src/pipeline_structure.r" )

### set run parameters

var_filter = T
#p_val = 0.05

create_cohorts    = T # 2
parse_files       = T # 3
normalize         = T # 4
qc_control        = F # 5
annotate          = T # 6
absent_analysis   = T # 7
dif_exp_ana       = T # 8
export_results    = T # 9
create_pathways   = F # 10
extract_interest  = F # 11
create_heatmaps_genes_of_interest = F # 12
annotate_tissue_abbundance = F # 13

### generic
stat_design = "contrast"

run_analysis();print( "Finished" )
