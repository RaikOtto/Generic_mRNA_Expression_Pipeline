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

which_project = "GSE43795"

source("project_files.r")
options(error=traceback)
source( "Src/pipeline_structure.r" )


### MZ_SM_all, MZlike_all, SM_all, SM_MZ_no_sici, MZ_no_sici, SM_no_sici
var_filter = F
<<<<<<< HEAD
p_val = 0.05
create_heatmaps_genes_of_interest = T

=======
>>>>>>> 62f0fce3bbae12e4804aa066fc90cfbf1036c426
create_cohorts    = T # 2
parse_files       = T # 3
normalize         = F # 4
qc_control        = F # 5
annotate          = F # 6
absent_analysis   = F # 7
<<<<<<< HEAD
dif_exp_ana       = F # 8
export_results    = F # 9
=======
dif_exp_ana       = T # 8
export_results    = T # 9
>>>>>>> 62f0fce3bbae12e4804aa066fc90cfbf1036c426
create_pathways   = F # 10
extract_interest  = F # 11
annotate_tissue_abbundance = F # 12

### generic
stat_design = "contrast"

run_analysis();print( "Finished" )
