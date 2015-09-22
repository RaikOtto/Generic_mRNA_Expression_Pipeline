pipeline_loc = "/usr/Generic_mRNA_Expression_Pipeline" #server path
#pipeline_loc = paste( system("echo $HOME",intern = T), "Generic_mRNA_Expression_Pipeline", sep ="/" ) #local path
setwd( paste(pipeline_loc, "Src", sep="/") ) # Set the path to where the pipeline is located

source( "pipeline_structure.r" )

default_parameters = T
which_project = "ovarian"
source("project_files.r")

###

create_cohorts    = T # 2
parse_files       = T # 3
normalize         = T # 4
qc_control        = F # 5
annotate          = T # 6
absent_analysis   = T # 7
dif_exp_ana       = T # 8
export_results    = T # 9
create_pathways   = T # 10
extract_interest  = T # 11
annotate_tissue_abbundance = F # 12

## Misc

# Extra information

expression_data = "~/Dropbox/PhD/NAR_sub_june_2015/expression_all.txt"

### frma test

frma_path = "/media/rayott/Backup/Runx_AML1_Leukemea/"

# generic
stat_design = "contrast"

setwd( paste(pipeline_loc, "Src", sep="/") )

source( "pipeline_structure.r" );run_analysis();print( "Finished" )
#save.image("~/ovarian.Rdata")