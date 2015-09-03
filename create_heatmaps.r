print("Creating heatmaps")

create_heatmap = function( res_int ){
  
  library("plotly")
  res_heat = res_int[,2:dim(res_int)[2]]
  
  #data[ dim(res_int)[1] , data[dim(res_int)[1],] == "Non_Responder" ] = "0"
  res_heat[ dim(res_int)[1] , res_heat[dim(res_int)[1],] == "Non_Responder" ] = "0"
  res_heat[ dim(res_int)[1] , res_heat[dim(res_int)[1],] == "Responder" ] = "1"
  
  info_row = as.character( res_int[,1 ] )
  info_col = colnames(res_int)[2:dim(res_int)[2]]
  
  res_heat = matrix( as.double( res_heat ), ncol = dim(res_int)[2]-1, nrow = dim(res_int)[1])
  
  trace1 = list(
  
    z = res_heat,
    x = info_col,
    y = info_row,
    type = "heatmap",
    showscale = F
  )
  
  data = list( trace1 )
  
  layout = list(
    
    title = "test",
    showlegend = T,
    autosize = F,
    showGrid = T,
    
    xaxis = list(
      title = "x"
    ),
    yaxis = list(
      title = "y"
    ),
    zaxis = list(
      title = "z"
    )
  )
  
  py <- plotly(username='bioinf', key='azet3sehsg')
  response <- py$plotly(data, kwargs=list(layout=layout))
  
  url <- response$url
  print(url)
  system(paste("firefox", url))
  remove(url)
}