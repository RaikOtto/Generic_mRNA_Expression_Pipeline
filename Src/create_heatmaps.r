print("Creating heatmaps")
if ( create_heatmaps_genes_of_interest ){
  
  heatmap_eset = topall_res[, which(colnames(topall_res) %in% c("logFC","HGNC_symb")  )  ]
  heatmap_eset = heatmap_eset[ heatmap_eset$HGNC_symb != ""  ,]
  heatmap_eset = heatmap_eset[ order( abs( heatmap_eset$logFC ), decreasing = T ), ]
  
  selected_genes = unique(as.character(heatmap_eset$HGNC_symb))[1:heatmap_list_genes_count]
  eset_selection = exprs(eset)[ which( hgnc_symbols %in% selected_genes)  ,]
  rownames(eset_selection) = hgnc_symbols[which( hgnc_symbols %in% selected_genes)  ]
  
  pdf_name = paste(
    output_path, 
    paste0(
      paste0(
        "Output/ExpressionSet_",
        project_name
      ),
      "_map.pdf"
    ),
    sep = "/"
  )
  
  # heatmap.3 
  
  library("devtools")
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  m=colorRampPalette(colors = c("green","black","red"))( 20 )
  
  pdf(pdf_name)
  library("gplots")
  heatmap.3( eset_selection, col = m)
  dev.off()
}