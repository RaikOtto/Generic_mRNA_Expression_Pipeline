print("Creating heatmaps")
if ( create_heatmaps_genes_of_interest ){
  
  cutoff = sort(abs( topall_res$logFC ), decreasing = T)[heatmap_list_genes_count]
  #probe_ids = rownames(topall_res)[ abs( topall_res$logFC ) >= cutoff ]
  probe_ids = topall_res$ID[ abs( topall_res$logFC ) >= cutoff ]
  hgnc_ids  = topall_res$Gene.symbol[ abs( topall_res$logFC ) >= cutoff ]
  #hgnc_ids  = topall_res$HGNC_symb[ abs( topall_res$logFC ) >= cutoff ]
  eset_selection = exprs(eset)[ match( probe_ids, rownames( exprs( eset ) )  ) , ]
  
  dif = as.double(rowMeans( eset_selection[ , index_ctrl]) - rowMeans( eset_selection[ , index_case]))
  
  rownames( eset_selection) = hgnc_ids
  eset_selection = eset_selection[,  order(order(c(index_ctrl,index_case))) ]
  eset_selection_dif = eset_selection - rowMeans( eset_selection )
  
  eset_selection_dif = eset_selection_dif[ order(dif, decreasing = F), ]
  
  pdf_name = paste(
    output_path, 
    paste0(
      paste0(
        "Output/",
        paste0( "Results_",project_name)
      ),
      "/heatmap.pdf"
    ),
    sep = "/"
  )
  
  # heatmap.3 
  
  logFC_side_bar = t(
    c(
      colorRampPalette(colors = c("red"))( length( dif[ dif < 0 ]) ),
      colorRampPalette(colors = c("green"))( length( dif[ dif > 0 ]) )
    )
  )
  rownames(logFC_side_bar) = c("FoldChange")
  library("devtools")
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  m=colorRampPalette(colors = c("green","black","red"))( 75 )
  
  pdf(pdf_name)
    
    heatmap.3( 
      eset_selection,
      main = paste0( "Sample exprs and logFC ", project_name ),
      #main = paste0( "Dif sample to avg expr all ", project_name ),
      RowSideColors = logFC_side_bar,
      col = m,
      trace = NULL,
      Colv = F,
      Rowv = F,
      dendrogram = "none" ,
      margins = c(9,7),
      cexCol = .75,
      cexRow = .75
    )
  dev.off()
  
  ### dif darstellung
  
  pdf_name = paste(
    output_path, 
    paste0(
      paste0(
        "Output/",
        paste0( "Results_",project_name)
      ),
      "/heatmap_difference_overall_expression.pdf"
    ),
    sep = "/"
  )
  
  # heatmap.3 
  
  logFC_side_bar = t(
    c(
      colorRampPalette(colors = c("red"))( length( dif[ dif < 0 ]) ),
      colorRampPalette(colors = c("green"))( length( dif[ dif > 0 ]) )
    )
  )
  rownames(logFC_side_bar) = c("FoldChange")
  library("devtools")
  source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
  m=colorRampPalette(colors = c("green","black","red"))( 75 )
  
  pdf(pdf_name)
  
  heatmap.3( 
    eset_selection_dif,
    main = paste0( "Sample exprs and logFC ", project_name ),
    #main = paste0( "Dif sample to avg expr all ", project_name ),
    RowSideColors = logFC_side_bar,
    col = m,
    trace = NULL,
    Colv = F,
    Rowv = F,
    dendrogram = "none" ,
    margins = c(9,7),
    cexCol = .75,
    cexRow = .75
  )
  dev.off()
}
