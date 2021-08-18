


##########################source functions for generating correlation matrices#######################################



#function to shrink the datasets to only intersecting cell lines given second dataset
shrink <- function(data1, data2){
  names1 <- row.names(data1)
  names2 <- row.names(data2)
  rows_to_remove <- NULL
  for(i in length(data1):1){
    if(!is.element(names1[i], names2)){
      rows_to_remove <- c(rows_to_remove, i)
    }
  }
  data1 <- data1[-rows_to_remove,]
  return(data1)
}



#generates the matrix from 2 tables sharing 1 dimension length
# createMatrix <- function(table1, table2){
#   correlations = matrix(nrow = ncol(table1), ncol = ncol(table2))
#   for (i in 1:ncol(table1)){
#     essentiality <- as.numeric(table1[,i])
#     for(j in 1:ncol(table2)){
#       drug_response <- as.numeric(table2[,j])
#       correlations[i, j] <- cor(essentiality, drug_response, method = "pearson")
#     }
#   }
#   return(correlations)
# }



#data1 will become the rows of the newly created data table, data2 will become the columns
correlations <- function(data1, data2){
  library(data.table)
  shrunk <- shrink(data1, data2)
  shrunk2 <- shrink(data2, data1)
  correlation_matrix <- cor(shrunk, shrunk2, use = "pairwise.complete.obs")
  correlation_data <- data.table(correlation_matrix)
  correlation_data <- cbind(gene = colnames(shrunk), correlation_data)
  return(correlation_data)
}

#same as above function but with spearman correlation
correlations2 <- function(data1, data2){
  library(data.table)
  shrunk <- shrink(data1, data2)
  shrunk2 <- shrink(data2, data1)
  correlation_matrix <- cor(shrunk, shrunk2, use = "pairwise.complete.obs", method = "spearman")
  correlation_data <- data.table(correlation_matrix)
  correlation_data <- cbind(gene = colnames(shrunk), correlation_data)
  return(correlation_data)
}

#function to return a correlation matrix based solely on drugs
codependency <- function(data, met = "pearson"){
  library(data.table)
  correlation_matrix <- cor(data, use = "pairwise.complete.obs", method = met)
  correlation_data <- data.table(correlation_matrix)
  correlation_data <- cbind(names = colnames(correlation_data), correlation_data)
  return(correlation_data)
}




