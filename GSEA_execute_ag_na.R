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

# GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute 
#
# R script to run GSEA Analysis of the Gender vs C1 example (cut and paste into R console)

GSEA(                                                                      # Input/Output Files :-------------------------------------------
                                                                           input.ds =  paste(cel_files_path, "ExpressionSet.gct", sep ="/"),               # Input gene expression Affy dataset file in RES or GCT format
                                                                           input.cls = paste(cel_files_path, "phenotypes_ag_na_GSEA.cls", sep ="/"),               # Input class vector (phenotype) file in CLS format
                                                                           gs.db =     paste(cel_files_path, "c2.all.v5.0.symbols.gmt", sep ="/"),           # Gene set database in GMT format
                                                                           output.directory      = output_path,            # Directory where to store output and results (default: "")
                                                                           #  Program parameters :----------------------------------------------------------------------------------------------------------------------------
                                                                           doc.string            = " SM_MZ_BM_all_vs_SM_MZ_Blood_all_GSEA",     # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
                                                                           non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
                                                                           reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
                                                                           nperm                 = 1000,            # Number of random permutations (default: 1000)
                                                                           weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
                                                                           nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
                                                                           fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
                                                                           fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
                                                                           topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
                                                                           adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
                                                                           gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
                                                                           gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
                                                                           reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
                                                                           preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
                                                                           random.seed           = 111,             # Random number generator seed. (default: 123456)
                                                                           perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
                                                                           fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
                                                                           replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
                                                                           save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
                                                                           OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
                                                                           use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
)
#--------------------------------------------------------------------------------------------------------------------------------------------------

# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
  directory           = output_path,        # Directory where to store output and results (default: "")
  topgs = 20,                                                           # number of top scoring gene sets used for analysis
  height = 16,
  width = 16
)