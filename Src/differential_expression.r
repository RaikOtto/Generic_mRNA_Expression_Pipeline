library("limma")

if ( chip_type %in% c( "hgu133plus2", "hgu133a" ) ){
  eset = eset[ ! startsWith( rownames(eset), "AFFX-" ), ]
}

eset = eset[ , which( colnames(eset) %in% phenodata$ID  ) ]  
index = which( colnames(phenodata) == cohorts_type )
pData( eset )$Group = phenodata[ match( colnames( eset ), phenodata$ID, nomatch = 0 ), index  ]
eDatSet = eset

if (var_filter ){
  
  exprs( eDatSet ) = exprs( varFilter( eDatSet  , var.cutoff = 0.6 ) )
}

source("Src/cohort_creation.r")

if (stat_design == "contrast"){
    
    fit = lmFit( eDatSet[ , c( index_ctrl, index_case )  ], design )
    cont.matrix = makeContrasts(  contrast = CASE - CTRL,  levels = design)
    fit = contrasts.fit( fit, cont.matrix )
    fit = eBayes( fit )
    volc_all = topTable( fit, coef = "contrast", number  = nrow(eDatSet), adjust  ="BH", p.value = 1, lfc = 0)
    
} else {
    
    fit = lmFit( eDatSet, design )  
    fit = eBayes( fit )
    volc_all = topTable( fit, number  = nrow(eDatSet), adjust  ="BH", p.value = 1, lfc = 0)
}

png( paste( output_path, "Output/logFC_vs_1-PValue.png", sep ="/" ), width = 800, height = 800, res = 150  )
plot( volc_all$logFC,   1-( volc_all$P.Value ))
dev.off()

topall = topTable( fit, coef = "contrast", number  = nrow( eDatSet ), adjust  ="none", p.value = p_val, lfc = lfc_exp)

if ( (dim(topall)[1] == 0) & (dim(topall)[2] == 0) ){
  stop("Topall has dimension zero")
  
}

topall$logFC = round(topall$logFC,2)
topall$AveExpr = round(topall$AveExpr,2)
topall$t = round( topall$t, 2 )
topall$B = round( topall$B, 2 )

topall = topall[ abs(topall$logFC) >= lfc_exp  ,]

message( c( "Amount probes higher in Case cohort:", sum( topall$logFC >= lfc_exp ) ) )
message( c( "Amount probes lower in Case cohort:" , sum( topall$logFC <= lfc_exp ) ) )
