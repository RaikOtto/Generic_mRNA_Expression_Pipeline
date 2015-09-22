pipeline_loc = paste( system("echo $HOME",intern = T), "Dropbox/PhD/Generic_Biomarker_mRNA_Pipeline/", sep ="/" )
setwd( pipeline_loc ) # Set the path to where the pipeline is located
source("pipeline_structure.r")

# defaults

time_series           = F # to be changed later; 
quality_control_only  = F; create_heatmaps_genes_of_interest = F
dif_exp_experiment    = T; integrate_vcf_files = F
multi_probe           = F
var_filter            = F # variance based filtering of the expression data jsut before differential expression detection
integrate_drug_data   = F
use_frma_normalization = F
heatmap_vis           = F; zipped = F
use_logFc_only        = F
run_generic           = T # 1
annotage_tissue_abbundance = F
cohorts_type          = "Group"
stat_design           = "contrast"
p_val = 0.1
lfc_exp = 1
absent_genes_file = "absent_genes.tab"
res_file = "dif_exp_res_file.xlsx"
time_series_res_file = "time_series_expression.csv"
genes_of_interest_file = "genes_of_interest.csv"

###

create_cohorts    = T # 2
parse_files       = T # 3
normalize         = T # 4
qc_control        = F # 5
annotate          = T # 6
absent_analysis   = F # 7
dif_exp_ana       = T # 8
export_results    = T # 9
create_pathways   = T # 10
extract_interest  = T # 11
annotate_tissue_abbundance = F # 12

## Misc

# Extra information

expression_data = "~/Dropbox/PhD/NAR_sub_june_2015/expression_all.txt"

### PROJECT FILE

## Characterization of Naive B-Cells

project_name = "Time_series_b_cell"
chip_type = "hgu133plus2"
cel_files_path = "/media/rayott/Backup/HGU133Plus2_HSC"
cohorts_file = "cohorts_time_series.tab"

#set_ctrl = c("HSC")
#set_case = c("ALL")

## Runx

project_name = "Plus2_Runx"
chip_type = "hgu133plus2"
#cel_files_path = "~/Dropbox/PhD/Runx/HGU_plus2/"
cel_files_path = "/media/rayott/Backup/HGU133Plus2_HSC"
cohorts_file = "cohorts.csv"
cohorts_type = "Group"
kegg_file = "pathways.csv"
stat_design   = "contrast"
p_val = 10^(-30)

set_ctrl = c("HSC")
set_case = c("ALL")

## ovarian benign tumor vs. malignant tumor

#project_name = "Ovarian_Gse29156_Malign_Case_vs_Benign_Tumor_Ctrl_cut_off_0_5"
#project_name = "Ovarian_Gse29156_Malign_Stroma_Case_vs_Benign_Stroma_Ctrl_cut_off_0_5"
#project_name = "Ovarian_Gse29156_Malign_Stroma_Tumor_vs_Benign_Stroma_Tumor_cut_off_0_5"
project_name = "test_run"

chip_type = "pd.huex.1.0.st.v2"
cel_files_path = "~/Dropbox/PhD/Ovarian_cancer/GSE29156_RAW/"
cohorts_file = "cohorts.tab"

#set_ctrl = c("Benign Tumor")
#set_ctrl = c("Benign Adjacent Stroma")
set_ctrl = c("Benign Tumor","Benign Adjacent Stroma")

#set_case = c("Malignant Tumor")
#set_case = c("Malignant Adjacent Stroma")
set_case = c("Malignant Tumor","Malignant Adjacent Stroma")

time_series = F # to be changed later
kegg_file = "pathways.csv"
zipped = T
p_val = .05
lfc_exp = 1

# GSE14407

project_name = "Ovarian_GSE14407"
chip_type = "hgu133plus2"
cel_files_path = "~/Dropbox/PhD/Ovarian_cancer/GSE14407_RAW/"
cohorts_file = "cohorts.csv"
set_ctrl = c("Normal")
set_case = c("Tumor")
time_series = F # to be changed later
kegg_file = "pathways.csv"
zipped = T
p_val = .001
lfc_exp = 1

## immunologie MZ

project_name = "Immuno_MZ"
chip_type = "pd.hugene.2.0.st"
cel_files_path = "~/Dropbox/PhD/ag_na_upload/CEL/"
cohorts_file = "cohorts.tab"
set_ctrl = c("Blood")
set_case = c("KM")
time_series = F # to be changed later
kegg_file = "kegg_pathways_of_interest_ag_na.csv"

p_val = 1.0
lfc_exp = 1

## immunologie sm

project_name = "Immuno_sm"
chip_type = "pd.hugene.2.0.st"
cel_files_path = "~/Dropbox/PhD/ag_na_upload/CEL"
cohorts_file = "cohorts.tab"
set_ctrl = c("Blood")
set_case = c("KM")
time_series = F # to be changed later
kegg_file = "kegg_pathways_of_interest_ag_na.csv"

p_val = 1.0
lfc_exp = 1

## immunologie sm bone_marrow vs mz bone_marrow

project_name = "Immuno_sm_vs_mz_no_blood"
chip_type = "pd.hugene.2.0.st"
cel_files_path = "~/Dropbox/PhD/ag_na_upload/CEL"
cohorts_file = "cohorts.tab"
set_ctrl = c("Control")
set_case = c("Case")
time_series = F # to be changed later
kegg_file = "kegg_pathways_of_interest_ag_na.csv"
#use_logFc_only = T
p_val = 1.0
lfc_exp = 1

## immunologie sm blood vs mz blood only

project_name = "Immuno_sm_blood_only_vs_mz_blood_only"
chip_type = "pd.hugene.2.0.st"
cel_files_path = "~/Dropbox/PhD/ag_na_upload/CEL"
cohorts_file = "cohorts.tab"
set_ctrl = c("Control")
set_case = c("Case")
time_series = F # to be changed later
kegg_file = "kegg_pathways_of_interest_ag_na.csv"
#use_logFc_only = T
p_val = 1.0
lfc_exp = 1

### kidney biomarker

project_name = "Kidney_biomarker"
chip_type = "hgu133a"
cel_files_path = "~/Dropbox/PhD/Kidney_Cancer_Biomarker/GSE15641_RAW_GPL96"
cohorts_file = "cohorts.tab"
set_ctrl = c("normal")
set_case = c("cc")
kegg_file = "kegg_pathways_of_interest.csv"
zipped = T
lfc_exp = 2
p_val = .05

## GSE53757

project_name = "Kidney_biomarker_GSE53757"
chip_type = "hgu133plus2"
cel_files_path = "~/Dropbox/PhD/Kidney_Cancer_Biomarker/GSE53757_RAW"
cohorts_file = "cohorts.csv"
set_ctrl = c("Normal")
set_case = c("Tumor")
kegg_file = "pathways.csv"
zipped = T
lfc_exp = .5
p_val = .05

### klinghammer hnsc

#project_name = "Cetuximab_Group_5_high_5_low"
#project_name = "Docetaxel_Group_5_high_5_low"
project_name = "cancer_type_classification"
chip_type = "hgu133plus2"
cel_files_path = "~/Dropbox/PhD/Klinghammer/upload/Cel_files"
cohorts_file = "cohorts.tab"
set_ctrl = c("Responder")
set_case = c("Non_Responder")
kegg_file = "pathways.csv"
vcf_folder = "~/Dropbox/PhD/Klinghammer/upload/VCF/"
zipped = F
lfc_exp = 1
integrate_drug_data = F
p_val = .05
#cohorts_type = "Docetaxel_Group_5_high_5_low"
#cohorts_type = "Cetuximab_Group_5_high_5_low"
#drug_type = "Cetuximab"
drug_type = "Docetaxel"
cohorts_type = "Docetaxel_Group"

### frma test

frma_path = "/media/rayott/Backup/Runx_AML1_Leukemea/"

# generic
stat_design = "contrast"

setwd( pipeline_loc )

source( "pipeline_structure.r" );run_analysis();print( "Finished" )
