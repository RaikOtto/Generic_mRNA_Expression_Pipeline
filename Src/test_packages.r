pack_bioc = c( "limma", "KEGG.db", "pathview", "stringr", "hgu133plus2.db", "affy", "simpleaffy", "affyPLM", "affycoretools", "affyQCReport", "annaffy", "hgu133a.db", "oligoData", "arrayQualityMetrics", "genefilter", "oligo", "pd.huex.1.0.st.v2", "lumi", "GEOquery", "illuminaHumanv4.db", "biomaRt" )
pack_cran = c( "xlsx", "gdata", "RColorBrewer", "stringr", "WriteXLS")

index_bioc = which( ! pack_bioc %in% installed.packages() )
index_cran = which( ! pack_cran %in% installed.packages() )

if ( length( index_bioc ) > 0){
  current = pack_bioc[index_bioc]
  source("https://bioconductor.org/biocLite.R")
  biocLite(current)
}

if ( length( index_cran ) > 0){
  current = pack_cran[index_cran]
  install.packages(pkgs = current, repos = "http://mirrors.softliste.de/cran/")
}

packages = c(pack_bioc, pack_cran)

print( suppressMessages( sapply( packages, require, character.only=TRUE, quietly = TRUE ) ) )