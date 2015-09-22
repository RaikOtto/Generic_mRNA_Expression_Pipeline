library("WriteXLS")
print("Step 11: Extraction of interesting information")
library("stringr")

dir.create( entities_of_interest_path, showWarnings = F)
source("Src/cohort_creation.r")

if ( create_heatmaps_genes_of_interest  ){
  links_heatmaps = c()
  source("Src/create_heatmaps.r")
}


if ( ! exists("hgnc_symbols"))
  source("Src/annotation.r")

topall     = exprs(eset)

exprs_case = rowMeans(topall[,index_case])
exprs_ctrl = rowMeans(topall[,index_ctrl])
dif_exp    = exprs_case - exprs_ctrl

topall = cbind(
  round( dif_exp, 2 ),
  round( exprs_case, 2),
  round( exprs_ctrl, 2),
  exprs(eset)
)

colnames(topall)[1] = "logFC"
colnames(topall)[2] = "Expr_case"
colnames(topall)[3] = "Expr_ctrl"
topall = as.data.frame(topall)

split_fun = function( entry ){ res = unlist( str_split( entry, " // " ) ); if (length(res) > 1){ return( res[2] ) } else{ return( "" ) } }
split_fun_ensembl = function( entry ){ res = unlist( str_split( entry, " // " ) ); if (length(res) > 1){ return( res[1] ) } else{ return( "" ) } }

if ( chip_type == "hgu133plus2" ){
  
  hgnc_symbols = mget( rownames(eset), hgu133plus2SYMBOL ); hgnc_symbols[ is.na(hgnc_symbols)  ] = ""
  
} else if ( chip_type == "hgu133a" ){
  
  hgnc_symbols = mget( rownames(eset), hgu133aSYMBOL ); hgnc_symbols[ is.na(hgnc_symbols)  ] = ""
  
} else {
  
  hgnc_symbols = str_trim( unlist( lapply( fData(eset)$geneassignment , FUN = split_fun ) ) )
  ensembl_list = str_trim( unlist( lapply( fData(eset)$geneassignment , FUN = split_fun_ensembl ) ) )
}

kegg_t = read.table( kegg_file_path, header =T , sep ="," )
cpdb_t = read.table( cpdb_file_path, header =T , sep ="\t", fill = T )
cpdb_ident = str_replace(cpdb_t$external_id, "path:", "")

interesting_pathways_mapping = match( str_trim( kegg_t$Kegg_id ), str_trim( cpdb_ident ) )
interesting_pathways_table_kegg_id = cpdb_t[ interesting_pathways_mapping ,]

### genes

if ( ! is.null( kegg_t$Gene_id_hgnc ) ){
  
  res_gene = c()
  res_ctrl = c()
  res_case = c()
  res_val  = c()
  
  for (gene in unique( kegg_t$Gene_id_hgnc )){
    
    if ( gene != ""  ){ 
      
      mapping = which( hgnc_symbols %in% gene )
      
      for (map in mapping){
        
        exprs_case = round( mean( exprs(eset)[ map, index_case]),2 )
        exprs_ctrl = round( mean( exprs(eset)[ map, index_ctrl]),2 )
        dif_exp = round( exprs_case - exprs_ctrl, 2)
        
        res_ctrl = c( res_ctrl, exprs_ctrl )
        res_case = c( res_case, exprs_case )
        res_gene = c(res_gene, gene)
        res_val = c( res_val, dif_exp)
      }
    }
    
  }
  
  genes = kegg_t$Gene_id_hgnc[ kegg_t$Gene_id_hgnc != ""  ]
  mapping = which( hgnc_symbols %in% genes )
  
  sample_expression = as.matrix( exprs(eset)[ mapping ,])
  
  if (dim(sample_expression)[2] == 1){
    
    sample_expression = t(sample_expression)
  }
  
  res_int = cbind(res_gene,res_val,res_ctrl, res_case, round( sample_expression, 2))
  colnames(res_int) =  c( "hgnc_symbol", "logFC" , "expr_ctrl", "expr_case", colnames(eset))
  res_int = res_int[ order( as.double(res_int[,2]), decreasing = T  ),]
  
  if ( integrate_drug_data  ){
    
    drug_vec = rep( "", length( res_int[1,] )  )
    drug_vec[1] = drug_type
    names( drug_vec ) = colnames( res_int )
    
    if ( ! exists("drug_data") ){
      
      print("Did not find drug data, but drug_response parameter was set")
      exit()
    }
    
    mapping = match( names( drug_vec ), names(drug_data), nomatch = 0  )
    
    if ( ! (T %in% ( mapping != 0 ) ) ){
      
      print("Problem mapping samplenames of genes of interest to drug data from cohorts file")
      exit()
    }
    
    drug_vec[ which( mapping != 0  ) ] = as.double( drug_data[ mapping  ] )
    res_ctrl  = as.double( drug_data[ index_ctrl ] )
    res_case = as.double( drug_data[ index_case ] )
    dif = round( mean( res_case ) - mean(as.double( res_ctrl  )), 1 )
    
    # drug response
    
    drug_vec[ match( "expr_case", names(drug_vec) ) ] = mean( res_case )
    drug_vec[ match( "expr_ctrl", names(drug_vec) ) ] = mean( res_ctrl )
    drug_vec[ match( "logFC", names(drug_vec) ) ]     = dif
    
    # correlations
    index = seq( 5, dim(res_int)[2]  )
    cor_pearson  = c()
    cor_spearman = c()
    
    for (n in seq( 1,dim(res_int)[1]  ) ){
      
      cor_pearson  = c( cor_pearson, cor( as.double( res_int[n,  index ] ), as.double( drug_vec[index] ), method = "pearson", use = "pairwise.complete.obs"))
      cor_spearman = c( cor_spearman, cor( as.double( res_int[n,  index ] ), as.double( drug_vec[index] ), method = "spearman", use = "pairwise.complete.obs"))
    }
    
    cor_pearson = round(cor_pearson, 2)
    cor_spearman = round(cor_spearman, 2)
    
    # cohorts
    
    case_ctrl = phenodata [ match( names(cohorts_vec), phenodata$ID  ) , which( colnames(phenodata) == cohorts_type  )  ]
    case_ctrl = c( "Cohort", "" , "" , "", as.character(case_ctrl)  )
    
    res_int = rbind( res_int, drug_vec  )
    res_int = rbind( res_int, case_ctrl  )
    
    res_int = cbind( res_int[ ,1:4], cor_pearson, cor_spearman, res_int[ , 5:dim(res_int)[2]  ] )
  }
  
  if (exists("links_heatmaps"))
    links_heatmaps = c( links_heatmaps, create_heatmap(res_int)  )
  
  write.xlsx( res_int, str_replace(str_replace(genes_of_interest_file_path,"~",user_folder),".csv",".xls"), row.names=F )
}
### pathways

if ( ! is.null( kegg_t$Kegg_id ) ){
  
  mapping = match( kegg_t$Kegg_id  ,cpdb_ident, nomatch = 0 )
  print(c("Not machted Kegg Pathways:", as.character(kegg_t$Kegg_id[mapping==0]) )  )
  
  mapping = match( cpdb_ident, kegg_t$Kegg_id, nomatch = 0 )
  mapping = mapping[mapping!=0]
  genes_of_interest = cpdb_t$hgnc_symbol_ids[mapping]
  
  for ( i  in  mapping ){
    
    pathway_id    = as.character( kegg_t$Kegg_id[ i ] )
    
    if (pathway_id != ""){
      
      pathway_name  = as.character( kegg_t$Kegg_name[ i ] )
      pathway_name = str_replace(pathway_name,"/","and")
      genes         = cpdb_t$hgnc_symbol_ids[ i ]
      gene_list     = unlist( str_split( genes, "," ) )
      mapping_gene  = match( gene_list, as.character( hgnc_symbols), nomatch = 0 )
      
      exprs_genes = exprs(eset)[mapping_gene,]
      exprs_case = round( rowMeans( exprs_genes[,index_case]),2 )
      exprs_ctrl = round( rowMeans( exprs_genes[,index_ctrl]),2 )
      dif_exp = round( exprs_case - exprs_ctrl, 2)
      
      exprs_genes = cbind( as.double( dif_exp ), as.double( exprs_ctrl), as.double( exprs_case), hgnc_symbols[mapping_gene] , round(exprs(eset)[mapping_gene,],2) )
      colnames(exprs_genes) =  c( "logFC" , "expr_ctrl", "expr_case", "hgnc_symbol", colnames(eset))
      exprs_genes = exprs_genes[ order( as.double(exprs_genes[,1]), decreasing = T  ),]
      
      if ( integrate_drug_data  ){
        
        drug_vec = rep( "", length( exprs_genes[1,] )  )
        drug_vec[1] = drug_type
        names( drug_vec ) = colnames( exprs_genes )
        
        if ( ! exists("drug_data") ){
          
          print("Did not find drug data, but drug_response parameter was set")
          exit()
        }
        
        mapping = match( names( drug_vec ), names(drug_data), nomatch = 0  )
        
        if ( ! (T %in% ( mapping != 0 ) ) ){
          
          print("Problem mapping samplenames of genes of interest to drug data from cohorts file")
          exit()
        }
        
        drug_vec[ which( mapping != 0  ) ] = as.double( drug_data[ mapping  ] )
        res_ctrl  = as.double( drug_data[ index_ctrl ] )
        res_case = as.double( drug_data[ index_case ] )
        dif = round( mean( res_case ) - mean(as.double( res_ctrl  )), 1 )
        
        # drugs
        
        drug_vec[ match( "expr_case", names(drug_vec) ) ] = mean( res_case )
        drug_vec[ match( "expr_ctrl", names(drug_vec) ) ] = mean( res_ctrl )
        drug_vec[ match( "logFC", names(drug_vec) ) ]     = dif
        drug_vec[ match( "hgnc_symbol", names(drug_vec) ) ] = drug_type
        
        # corrs
        
        index = seq( 5, dim(exprs_genes)[2]  )
        cor_pearson  = c()
        cor_spearman = c()
        
        for (n in seq( 1,dim(exprs_genes)[1]  ) ){
          
          cor_pearson  = c( cor_pearson, cor( as.double( exprs_genes[n,  index ] ), as.double( drug_vec[index] ), method = "pearson", use = "pairwise.complete.obs"))
          cor_spearman = c( cor_spearman, cor( as.double( exprs_genes[n,  index ] ), as.double( drug_vec[index] ), method = "spearman", use = "pairwise.complete.obs"))
        }
        
        cor_pearson = round(cor_pearson, 2)
        cor_spearman = round(cor_spearman, 2)
        
        # cohorts
        
        case_ctrl = phenodata [, which( colnames(phenodata) == cohorts_type  )  ]
        case_ctrl = c( "Cohort", "" , "" , "", as.character(case_ctrl)  )
        
        exprs_genes = rbind( exprs_genes, drug_vec  )
        exprs_genes = rbind( exprs_genes, case_ctrl  )
        
        exprs_genes = cbind( exprs_genes[ ,1:4], cor_pearson, cor_spearman, exprs_genes[ , 5:dim(exprs_genes)[2]  ] )
      }
      
      file_name = str_replace( genes_of_interest_file_path, "genes_of_interest", paste( pathway_id, pathway_name, sep ="_")  )
      print(c(i,file_name))
      
      write.xlsx( exprs_genes, str_replace(str_replace(file_name,"~",user_folder),".csv",".xls"), row.names=F )
      
    }
  }
}

#if ( ! is.null( kegg_t$Go_id ) ){
if ( F ){
  
  print( "Go ids"  )
  ensembl     = useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  for ( i  in 1:length(kegg_t$Go_id)  ){
    
    go_id = str_trim(kegg_t$Go_id[i])
    
    if( go_id != ""){
      
      go_name  = as.character( kegg_t$Go_name[ i ] )
      go_name = str_replace(go_name,"/","and")
      
      genes = getBM( attributes = c( "hgnc_symbol" ), values = go_id, filters = "go_id" , mart = ensembl)
      gene_list = as.character(unlist(genes))
      
      mapping_gene  = match( gene_list, as.character( hgnc_symbols), nomatch = 0 )
      
      exprs_genes = exprs(eset)[ mapping_gene,]
      exprs_case = round( rowMeans( exprs_genes[,index_case]),2 )
      exprs_ctrl = round( rowMeans( exprs_genes[,index_ctrl]),2 )
      dif_exp = round( exprs_case - exprs_ctrl, 2)
      
      exprs_genes = cbind( as.double( dif_exp ), as.double( exprs_ctrl), as.double( exprs_case), hgnc_symbols[mapping_gene] , round(exprs(eset)[mapping_gene,],2) )
      colnames(exprs_genes) =  c( "logFC" , "expr_ctrl", "expr_case", "hgnc_symbol", colnames(eset))
      exprs_genes = exprs_genes[ order( as.double(exprs_genes[,1]), decreasing = T  ),]
      
      file_name = str_replace( genes_of_interest_file_path, "genes_of_interest", go_name  )
      print(c(i,file_name))
      
      write.xlsx( exprs_genes, str_replace(str_replace(file_name,"~",user_folder),".csv",".xls"), row.names=F )
    }
  }
}