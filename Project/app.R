#install and load packages. 

install.packages("shiny")
install.packages("shinyalert")
install.packages("DT")
install.packages("visNetwork")
install.packages("data.table")
install.packages("bslib")
install.packages("tidyverse")
install.packages("ggplot")
install.packages("shinycssloaders")
install.packages("ggpubr")
install.packages("grid")
library(ggpubr)
library(grid)
library(shiny)
library(bslib)
library(shinyalert)
library(DT)
library(visNetwork)
library(tidyverse)
library(shinycssloaders)

#change the directories to where ever you put the project folder. ex: path/[project name]
## load dependency data ## change directory as needed.
setwd("C:/Users/ningz/OneDrive/Desktop/My files/DrugResponse2021/crispr_dep_map")

#CHRONOS version.
CRISPR_dep <- read.delim(file = "CRISPR_gene_effect_Chronos.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# cell-line info
sample_info_dep <- read.delim(file = "sample_info.csv", stringsAsFactors = FALSE, header = TRUE, sep=",", row.names = 1)

## load drug response data (primary set) ## change directory as needed.
setwd("C:/Users/ningz/OneDrive/Desktop/My files/DrugResponse2021/PRISM")

# drug response data
drug_resp_primary <- read.delim(file = "primary-screen-replicate-collapsed-logfold-change.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# cell-line info
sample_info_primary <- read.delim(file = "primary-screen-cell-line-info.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)
# drug info
treatment_info_primary <- read.delim(file = "primary-screen-replicate-collapsed-treatment-info.csv", sep = ",", stringsAsFactors = FALSE, header = TRUE, row.names = 1)

#change directory as needed.
setwd("C:/Users/ningz/OneDrive/Desktop/My files/DrugResponse2021")
#source functions
source("NameChanger.R")
source("DepCorrGenerator.R")
source("correlation_query.R")
source("correlation_plots.R")

#change the row and column names of the CRISPR data
update_gene <- cell_name(CRISPR_dep, sample_info_dep, name_col = 2)
update_gene <- gene_name(update_gene)

#change the row and column names of the drug data
update_drug <- cell_name(drug_resp_primary, sample_info_dep, name_col = 2)
update_drug <- drug_name(update_drug, treatment_info_primary, name_col = 2)

#generate the large correlations table as well as target pairs from the treatment table
correlations_table <- correlations(update_gene, update_drug)


#generate the large correlations table for spearman correlation
correlations_table2 <- correlations2(update_gene, update_drug)


#generate codependency matrices
drugcodep_pearson <- codependency(update_drug)
drugcodep_spearman <- codependency(update_drug, met = "spearman")

#split the gene dependency table into 2 halves and do 3 pairwise as we cannot allocate the memory
update_gene1 <- update_gene[, 1:8696]
update_gene2 <- update_gene[, 8697:17393]
genecodep_pearson <- correlations(update_gene1, update_gene2)
genecodep_spearman <- correlations2(update_gene1, update_gene2)
genecodep_pearson1 <- codependency(update_gene1)
genecodep_spearman1 <- codependency(update_gene2, met = "spearman")
genecodep_pearson2 <- codependency(update_gene2)
genecodep_spearman2 <- codependency(update_gene2, met = "spearman")


#alternatively, load the tables from the folder
correlations_table <- data.table(read.csv("gene_drug_dependencies_chronos.csv", row.names = 1))
correlations_table2 <- data.table(read.csv("gene_drug_dep_CHRONOS_spearman.csv", row.names = 1))
drugcodep_pearson <- data.table(read.csv("drugcodep_pearson.csv", row.names = 1))
drugcodep_spearman <- data.table(read.csv("drugcodep_spearman.csv", row.names = 1))
genecodep_pearson <- data.table(read.csv("genecodep_pearson.csv", row.names = 1))
genecodep_spearman <- data.table(read.csv("genecodep_spearman.csv", row.names = 1))
genecodep_pearson1 <- data.table(read.csv("genecodep_pearson1.csv", row.names = 1))
genecodep_spearman1 <- data.table(read.csv("genecodep_spearman1.csv", row.names = 1))
genecodep_pearson2 <- data.table(read.csv("genecodep_pearson2.csv", row.names = 1))
genecodep_spearman2 <- data.table(read.csv("genecodep_spearman2.csv", row.names = 1))

#remove the diagonals, don't remove gencodep because that is intersection between halfs of the gene table
col <- drugcodep_pearson[,1]
drugcodep_pearson <- drugcodep_pearson[, -1]
diag(drugcodep_pearson) <- NA
drugcodep_pearson <- cbind(col, drugcodep_pearson)
col <- genecodep_pearson1[,1]
genecodep_pearson1 <- genecodep_pearson1[, -1]
diag(genecodep_pearson1) <- NA
genecodep_pearson1 <- cbind(col, genecodep_pearson1)
col <- genecodep_pearson2[,1]
genecodep_pearson2 <- genecodep_pearson2[, -1]
diag(genecodep_pearson2) <- NA
genecodep_pearson2 <- cbind(col, genecodep_pearson2)


#get treatment pairs based on the tables. target_correlations should be pearson and target_correlations2 should be spearman.
target_correlations <- targetgene(correlations_table, treatment_info_primary)
target_correlations <- target_correlations[order(-correlation)]
target_correlations2 <- targetgene(correlations_table2, treatment_info_primary)
target_correlations2 <- target_correlations2[order(-correlation)]

#load precalculated targets that are above a certain correlation threshold
valuable_targets <- data.table(read.csv(row.names = 1, "valuable_targets_5_percent.csv"))
valuable_targets <- valuable_targets[order(-correlation)]
valuable_targets2 <- data.table(read.csv(row.names = 1, "valuable_targets_10_percent.csv"))
valuable_targets2 <- valuable_targets2[order(-correlation)]

valuable_targets_spearman <- data.table(read.csv(row.names = 1, "valuable_targets_5_percent_spearman.csv"))
valuable_targets_spearman <- valuable_targets_spearman[order(-correlation)]
valuable_targets2_spearman <- data.table(read.csv(row.names = 1, "valuable_targets_10_percent_spearman.csv"))
valuable_targets2_spearman <- valuable_targets2_spearman[order(-correlation)]

#intersection between updated drug and gene tables.
drugsplit <- shrink(update_drug, update_gene)
genesplit <- shrink(update_gene, update_drug)

#get all valuable pairs above a 0.3 absolute value correlation
corlist <- target(correlations_table, corr = 0.7)
corlist <- cbind(corlist, val1type = "drug", val2type = "gene")
corlist1 <- target(drugcodep_pearson, corr = 0.7)
corlist1 <- cbind(corlist1, val1type = "drug", val2type = "drug")
corlist2 <- target(genecodep_pearson, corr = 0.7)
corlist2 <- cbind(corlist2, val1type = "gene", val2type = "gene")
corlist3 <- target(genecodep_pearson1, corr = 0.7)
corlist3 <- cbind(corlist3, val1type = "gene", val2type = "gene")
corlist4 <- target(genecodep_pearson2, corr = 0.7)
corlist4 <- cbind(corlist4, val1type = "gene", val2type = "gene")
corlist_0.7 <- rbind(corlist, corlist1, corlist2, corlist3, corlist4)
colnames(corlist_0.7) <- c("column", "val1", "row", "val2", "correlation", "val1type", "val2type")

##################################Shiny APP code##########################################

#generate network between genes and drugs
returnNet <- function(targets){
  drugs <- targets %>%
    distinct(drug) %>%
    rename(label = drug)
  drugs <- cbind(drugs, shape = "square")
  genes <- targets %>%
    distinct(gene) %>%
    rename(label = gene)
  genes <- cbind(genes, shape = "circle")
  nodes <- full_join(drugs, genes, by = c("label", "shape"))
  nodes <- nodes %>% rowid_to_column("id")
  
  per_route <- data.table(drug = targets$drug, gene = targets$gene,
                          weight = (round(targets$correlation*100)))
  
  edges <- per_route %>% 
    left_join(nodes, by = c("drug" = "label")) %>% 
    rename(from = id)
  
  edges <- edges %>% 
    left_join(nodes, by = c("gene" = "label")) %>% 
    rename(to = id)
  
  edges <- select(edges, from, to, weight)
  
  
  edges <- mutate(edges, width = abs(weight)/5 + 1, color = ifelse(weight > 0, "blue", "red"))
  edges$weight <- abs(edges$weight)
  
  lnodes <- data.frame(label = c("Drug", "Gene"),
                       shape = c( "triangle", "square"), color = list(background = "#4cfc7b",border = "#02b832"),
                       title = "Information", id = 1:2)
  
  ledges <- data.frame(color = c("blue", "red"),
                       label = c("r>0", "r<0"),
                       width = 4,
                       arrows.scaleFactor = 0)
  
  network <- visNetwork(nodes, edges) %>%
    visNodes(size = 30, color = list(background = "#4cfc7b",border = "#02b832"))%>%
    visLegend(width = 0.1, addEdges = ledges, addNodes = lnodes, useGroups = FALSE)%>%
    visInteraction(multiselect = TRUE)
  return(network)
}


#generate network for all data including codependency data
returnNet2 <- function(targets){
  values1 <- targets %>%
    distinct(val1, val1type) %>%
    rename(label = val1)
  values1 <- cbind(values1, shape = ifelse(values1$val1type=="drug", "square", "circle"))
  values2 <- targets %>%
    distinct(val2, val2type) %>%
    rename(label = val2)
  values2 <- cbind(values2, shape = ifelse(values2$val2type=="drug", "square", "circle"))
  nodes <- full_join(values1, values2, by = c("label", "shape"))
  nodes <- nodes[, c(1,3)]
  nodes <- nodes %>% rowid_to_column("id")
  
  per_route <- data.table(value1 = targets$val1, value2 = targets$val2,
                          weight = (round(targets$correlation*100)))
  
  edges <- per_route %>% 
    left_join(nodes, by = c("value1" = "label")) %>% 
    rename(from = id)
  
  edges <- edges %>% 
    left_join(nodes, by = c("value2" = "label")) %>% 
    rename(to = id)
  
  edges <- select(edges, from, to, weight)
  
  
  edges <- mutate(edges, width = abs(weight)/10 + 1, color = ifelse(weight > 0, "blue", "red"))
  edges$weight <- abs(edges$weight)
  
  lnodes <- data.frame(label = c("Drug", "Gene"),
                       shape = c( "triangle", "square"), color = list(background = "#4cfc7b",border = "#02b832"),
                       title = "Information", id = 1:2)
  
  ledges <- data.frame(color = c("blue", "red"),
                       label = c("r>0", "r<0"),
                       width = 4,
                       arrows.scaleFactor = 0)
  
  network <- visNetwork(nodes, edges) %>%
    visNodes(size = 30, color = list(background = "#4cfc7b",border = "#02b832"))%>%
    visLegend(width = 0.1, addEdges = ledges, addNodes = lnodes, useGroups = FALSE)%>%
    visInteraction(multiselect = TRUE)
  return(network)
}

#used to generate selected node ids for the select all option in the Gene and Drug networks tab (panel 1)
returnNodeIds <- function(targets, type = "drug"){
  drugs <- targets %>%
    distinct(drug) %>%
    rename(label = drug)
  drugs <- cbind(drugs, shape = "square")
  genes <- targets %>%
    distinct(gene) %>%
    rename(label = gene)
  genes <- cbind(genes, shape = "circle")
  nodes <- full_join(drugs, genes, by = c("label", "shape"))
  nodes <- nodes %>% rowid_to_column("id")
  
  if(type == "drug"){
    selected <- as.numeric(which(nodes$shape == "square"))
  }else if(type == "gene"){
    selected <- as.numeric(which(nodes$shape == "circle"))
  }else{
    selected <- numeric()
  }
  return(selected)
}

#used to generate selected node ids for the select all option in the Gene and Drug networks tab (panel 2)
returnNodeIds2 <- function(targets, type = "drug"){
  values1 <- targets %>%
    distinct(val1, val1type) %>%
    rename(label = val1)
  values1 <- cbind(values1, shape = ifelse(values1$val1type=="drug", "square", "circle"))
  values2 <- targets %>%
    distinct(val2, val2type) %>%
    rename(label = val2)
  values2 <- cbind(values2, shape = ifelse(values2$val2type=="drug", "square", "circle"))
  nodes <- full_join(values1, values2, by = c("label", "shape"))
  nodes <- nodes[, c(1,3)]
  nodes <- nodes %>% rowid_to_column("id")
  
  if(type == "drug"){
    selected <- as.numeric(which(nodes$shape == "square"))
  }else if(type == "gene"){
    selected <- as.numeric(which(nodes$shape == "circle"))
  }else{
    selected <- numeric()
  }
  return(selected)
}

#generates the table for the cell line panel
returnPlot <- function(data, drugtable, genetable, bins){
  library(ggplot2)
  if(data[2] == "drug response"){
    row <- which(rownames(drugtable) == data[1])
    if(length(row) <= 0){
      df <- data.frame()
      plot <- ggplot(df) + geom_histogram() + xlim(0, 10) + ylim(0, 100) + theme_gray() + 
        labs(title = "The drug data for this cell line is not yet available",
             x = "efficacy", y = "frequency") + theme(plot.title = element_text(hjust = 0.5))
    }
    else{
      row <- as.numeric(row[1])
      values <- as.numeric(drugtable[row, ])
      plot <- ggplot()+aes(values)+geom_histogram(fill = "#0fbf35", color = "black", bins = bins, alpha = 0.7)+ theme_gray()+
        labs(title = paste("efficacy of all drugs for the", rownames(drugtable)[row[1]], "cell line", sep = " "),
             x = "efficacy", y = "frequency") + theme(plot.title = element_text(hjust = 0.5))
    }
  } else {
    row <- which(rownames(genetable) == data[1])
    if(length(row) <= 0){
      df <- data.frame()
      plot <- ggplot(df) + geom_histogram() + xlim(0, 10) + ylim(0, 100) + theme_gray() + 
        labs(title = "The gene data for this cell line is not yet available",
             x = "efficacy", y = "frequency") + theme(plot.title = element_text(hjust = 0.5))
    }
    else{
      row <- as.numeric(row[1])
      values <- as.numeric(genetable[row, ])
      plot <- ggplot()+aes(values)+geom_histogram(fill = "#ed5807", color = "black", bins = bins, alpha = 0.7)+ theme_gray()+
        labs(title = paste("efficacy of all genes for the", rownames(genetable)[row[1]], "cell line", sep = " "),
             x = "efficacy", y = "frequency") + theme(plot.title = element_text(hjust = 0.5))
    }
  }
  return(plot)
}

#function to generate correlation plots for associations
returnPlot2 <- function(data, drug_table, gene_table, met = "pearson"){
  
  num1 <- which(colnames(drug_table)==data[1])
  num2 <- which(colnames(gene_table)==data[2])
  
  drugvector <- as.numeric(drug_table[, num1])
  genevector <- as.numeric(gene_table[, num2])
  
  corr <- cor(x = drugvector, y = genevector, use = "pairwise.complete.obs", method = met)
  
  grob <- grobTree(textGrob(paste("correlation (", met, "):", round(corr, digits = 6)), 
                            x=0.1,  y=0.9, hjust=0,
                            gp=gpar(col="black", fontsize=13, fontface="italic")))
  
  plot <- ggplot() + aes(x = drugvector, y = genevector)+ annotation_custom(grob) + geom_point()+ geom_smooth(method=lm) +
  labs(title = paste("efficacy between", data[1], "and", data[2], "over shared cell lines",  sep = " "),
       x = paste("efficacy of", data[1], sep = " "), y = paste("efficacy of", data[2], sep = " ")) +
    theme(plot.title = element_text(hjust = 0.5)) 
  return(plot)
  
}


#function to generate a plot for the associations panel
returnPlot3 <- function(table, value, type = "drug", bins){
  fillcol <- ifelse(type=="drug", "#0fbf35", "#ed5807")
  v1 <- as.numeric(table[, as.character(value)])
  plot <- ggplot()+aes(v1)+geom_histogram(fill = fillcol, bins = bins, color = "black", alpha = 0.7)+ theme_gray()+
    labs(title = paste("efficacy of", value, "over all cell lines", sep = " "),
         x = "efficacy", y = "frequency") + theme(plot.title = element_text(hjust = 0.5))
  return(plot)
}

#returns a subset of the target correlations table. A potential would be to extend this query to the entire table
#but it might overload shiny.
returnTable <- function(data, targets){
  table <- data.table()
  if(data[1] == "None" && data[2] == "None" && data[4] == ">"){
    table <- targets[correlation > as.numeric(data[3])]
  }else if (data[1] == "None" && data[2] == "None" && data[4] == "<"){
    table <- targets[correlation < as.numeric(data[3])]
  }else if (data[1] == "None" && data[4] == ">"){
    table <- targets[correlation > as.numeric(data[3]) & gene == data[2]]
  }else if (data[1] == "None" && data[4] == "<"){
    table <- targets[correlation < as.numeric(data[3]) & gene == data[2]]
  }else if (data[2] == "None" && data[4] == ">"){
    table <- targets[correlation > as.numeric(data[3]) & drug == data[1]]
  }else if (data[2] == "None" && data[4] == "<"){
    table <- targets[correlation < as.numeric(data[3]) & drug == data[1]]
  }else if (data[4] == ">"){
    table <- targets[correlation > as.numeric(data[3]) & drug == data[1] & gene == data[2]]
  }else if (data[4] == "<"){
    table <- targets[correlation < as.numeric(data[3]) & drug == data[1] & gene == data[2]]
  }
  return(table)
}

#cell line panel
panel1 <- tabPanel(title = "Cell lines",
               fluidRow(column(width = 4, offset = 4, align = "center", actionLink(inputId = "learnmore1", label = "click to learn more about this function"))),
               fluidRow(column(width = 6, align = "center",
                               selectInput(inputId = "cellselect", label = "select a cell line",
                                           choices = unique(c(rownames(update_drug), rownames(update_gene))))), 
                        column(width = 6, align = "center", radioButtons(inputId = "genesordrugs", label = "treatment", choices = c("drug response", "gene knockout")))
               ),
               fluidRow(column(width = 4, offset = 4, align = "center", actionButton(inputId = "update1", label = "Show results"))),
               fluidRow(column(width = 12, tags$p("                         "))),
               fluidRow(column(width = 12, tags$p("                         "))),
               fluidRow(column(width = 12, align = "center", plotOutput("cellhist"))),
               fluidRow(column(width = 12, tags$p("                         "))),
               fluidRow(column(width = 12, align = "center", sliderInput(inputId = "bins", label = "Bins", min = 1, max = 100, value = 30, width = "100%")))
          )

#panels 4 5 and 6 are used in the Associations panel
panel4 <- tabPanel(title = "Genes & Drugs", 
                   fluidRow(column(width = 4, offset = 4, align = "center", actionLink(inputId = "learnmore2", label = "click to learn more about this function"))),
                    fluidRow(column(width = 6, align = "center",
                                    selectInput(inputId = "drugselect1", label = "select a  drug",
                                                choices = unique(as.character(c(colnames(correlations_table)[-1]))))), 
                            column(width = 6, align = "center", selectInput(inputId = "geneselect1", 
                              label = "select a gene", choices = unique(as.character(c(correlations_table$gene)))))
                   ),
                   fluidRow(column(width = 6, align = "center", 
                                   actionButton(inputId = "druggraph", label = "Drug results")),
                            column(width = 6, align = "center",
                                   actionButton(inputId = "genegraph", label = "Gene results"))),
                   br(),
                   br(),
                   fluidRow(column(width = 6, align = "center", plotOutput("drugplot")),
                            column(width = 6, align = "center", plotOutput("geneplot"))),
                   br(),
                   fluidRow(column(width = 6, align = "center", 
                            sliderInput(inputId = "bins1", label = "Bins", min = 1, max = 100, value = 30, width = "100%")),
                            column(width = 6, align = "center",
                            sliderInput(inputId = "bins2", label = "Bins", min = 1, max = 100, value = 30, width = "100%"))),
                   br(),
                   fluidRow(column(width = 4, offset = 4, align ="center", actionButton(inputId = "showcorr", label = "Show Correlation"))),
                   br(),
                   fluidRow(column(width = 6, align = "center", plotOutput("corrplot")),
                            column(width = 6, align = "center", plotOutput("corrplot2"))),
                  )

panel5 <- tabPanel(title = "Drugs Only", 
                   fluidRow(column(width = 4, offset = 4, align = "center", actionLink(inputId = "learnmore3", label = "click to learn more about this function"))),
                   fluidRow(column(width = 6, align = "center",
                                   selectInput(inputId = "drugselect2", label = "select drug one",
                                               choices = unique(as.character(c(colnames(correlations_table)[-1]))))), 
                            column(width = 6, align = "center", selectInput(inputId = "drugselect3", 
                                              label = "select drug two", choices = unique(as.character(c(colnames(correlations_table)[-1])))))
                   ),
                   fluidRow(column(width = 6, align = "center", 
                                   actionButton(inputId = "druggraph1", label = "Drug1 results")),
                            column(width = 6, align = "center",
                                   actionButton(inputId = "druggraph2", label = "Drug2 results"))),
                   br(),
                   br(),
                   fluidRow(column(width = 6, align = "center", plotOutput("drugplot1")),
                            column(width = 6, align = "center", plotOutput("drugplot2"))),
                   br(),
                   fluidRow(column(width = 6, align = "center", 
                                   sliderInput(inputId = "bins3", label = "Bins", min = 1, max = 100, value = 30, width = "100%")),
                            column(width = 6, align = "center",
                                   sliderInput(inputId = "bins4", label = "Bins", min = 1, max = 100, value = 30, width = "100%"))),
                   br(),
                   fluidRow(column(width = 4, offset = 4, align ="center", actionButton(inputId = "showcorr1", label = "Show Correlation"))),
                   br(),
                   fluidRow(column(width = 6, align = "center", plotOutput("corrplot3")),
                            column(width = 6, align = "center", plotOutput("corrplot4"))),
              )

panel6 <- tabPanel(title = "Genes Only", 
                   fluidRow(column(width = 4, offset = 4, align = "center", actionLink(inputId = "learnmore4", label = "click to learn more about this function"))),
                   fluidRow(column(width = 6, align = "center",
                                   selectInput(inputId = "geneselect2", label = "select gene one",
                                               choices = unique(as.character(c(correlations_table$gene))))), 
                            column(width = 6, align = "center", selectInput(inputId = "geneselect3", 
                                                                            label = "select gene two", choices = unique(as.character(c(correlations_table$gene)))))
                   ),
                   fluidRow(column(width = 6, align = "center", 
                                   actionButton(inputId = "genegraph1", label = "Gene1 results")),
                            column(width = 6, align = "center",
                                   actionButton(inputId = "genegraph2", label = "Gene2 results"))),
                   br(),
                   br(),
                   fluidRow(column(width = 6, align = "center", plotOutput("geneplot1")),
                            column(width = 6, align = "center", plotOutput("geneplot2"))),
                   br(),
                   fluidRow(column(width = 6, align = "center", 
                                   sliderInput(inputId = "bins5", label = "Bins", min = 1, max = 100, value = 30, width = "100%")),
                            column(width = 6, align = "center",
                                   sliderInput(inputId = "bins6", label = "Bins", min = 1, max = 100, value = 30, width = "100%"))),
                   br(),
                   fluidRow(column(width = 4, offset = 4, align ="center", actionButton(inputId = "showcorr2", label = "Show Correlation"))),
                   br(),
                   fluidRow(column(width = 6, align = "center", plotOutput("corrplot5")),
                            column(width = 6, align = "center", plotOutput("corrplot6"))),
          )

#pearson correlation panel on correlation query
panel2 <- tabPanel(title = "Pearson Correlations",
                   fluidRow(column(width = 4, offset = 4, align = "center", actionLink(inputId = "learnmore", label = "click to learn more about this function"))),
                   fluidRow(column(width = 4, align = "center",
                                   selectInput(inputId = "drugselect", label = "select a drug",
                                               choices = c("None", unique(target_correlations$drug)))), 
                            column(width = 4, align = "center",
                                   selectInput(inputId = "geneselect", label = "select a gene",
                                               choices = c("None", unique(target_correlations$gene)))), 
                            column(width = 3, align = "center",
                                   numericInput(inputId = "corrselect", label = "select a correlation",
                                                value = 0, min = -1, max = 1, step = 0.05)),
                            column(width = 1, align = "center",
                                   radioButtons(inputId = "greaterorless", label = "region", choices = c(">", "<")))
                            ),
                   fluidRow(column(width = 4, offset = 4, align = "center", actionButton(inputId = "update", label = "Update"))),
                   br(),
                   fluidRow(column(width = 12, align = "center", DT::dataTableOutput("targettable"))),
                   br(),
                   fluidRow(column(width = 4, offset = 4, align = "center", 
                                   actionButton(inputId = "showspecnet", label = "generate network for these pairs"))),
                   br(),
                   fluidRow(column(width = 12, align = "center", visNetworkOutput("specnet")))
                   )

#spearman correlation panel on correlation query
panel7 <- tabPanel(title = "Spearman Correlations",
                   fluidRow(column(width = 4, offset = 4, align = "center", actionLink(inputId = "learnmore5", label = "click to learn more about this function"))),
                   fluidRow(column(width = 4, align = "center",
                                   selectInput(inputId = "drugselect4", label = "select a drug",
                                               choices = c("None", unique(target_correlations2$drug)))), 
                            column(width = 4, align = "center",
                                   selectInput(inputId = "geneselect4", label = "select a gene",
                                               choices = c("None", unique(target_correlations2$gene)))), 
                            column(width = 3, align = "center",
                                   numericInput(inputId = "corrselect1", label = "select a correlation",
                                                value = 0, min = -1, max = 1, step = 0.05)),
                            column(width = 1, align = "center",
                            radioButtons(inputId = "greaterorless1", label = "region", choices = c(">", "<")))
                   ),
                   fluidRow(column(width = 4, offset = 4, align = "center", actionButton(inputId = "updatesp", label = "Update"))),
                   br(),
                   fluidRow(column(width = 12, align = "center", DT::dataTableOutput("targettable1"))),
                   br(),
                   fluidRow(column(width = 6, offset = 3, align = "center", 
                                   actionButton(inputId = "showspecnet1", label = "generate network for these pairs"))),
                   br(),
                   fluidRow(column(width = 12, align = "center", visNetworkOutput("specnet1")))
)

#target pairs panel on drug and gene networks
panel3 <- tabPanel(title = "Target Pairs", 
                   fluidRow(column(width = 12, align = "center", actionLink(inputId = "networklm", 
                                  label = "Displaying connections between Drugs and Genes. Click to learn more."))),
                   br(),
                   fluidRow(column(width = 6, offset = 3, align = "center", radioButtons(inputId = "selectNodes", inline = TRUE,
                                   label = "Select All: ", choices = (c("none","drug", "gene"))))),
                   br(),
                   fluidRow(
                     column(width = 12, visNetworkOutput("network", height = "900px"))
                   ))
#important pairs panel on drug and gene networks
panel8 <- tabPanel(title = "Important Pairs", 
                   fluidRow(column(width = 12, align = "center", actionLink(inputId = "networklm1", 
                                                                            label = "Displaying connections between important items. Click to learn more."))),
                   br(),
                   fluidRow(column(width = 6, offset = 3, align = "center", radioButtons(inputId = "selectNodes2", inline = TRUE,
                                                                                         label = "Select All: ", choices = (c("none","drug", "gene"))))),
                   br(),
                   fluidRow(
                     column(width = 12, visNetworkOutput("network1", height = "1000px"))
                   ))

#main page layout
ui <- navbarPage(title = "DEPMAP",
                 useShinyalert(),
                 theme = bs_theme(version = 4, bootswatch = "minty"),
                 panel1,
                 tabPanel(title = "Associations",
                   tabsetPanel(
                    panel4,
                    panel5,
                    panel6
                   )
                 ),
                 tabPanel(title = "Correlation Query",
                    tabsetPanel(
                       panel2,
                       panel7
                    )         
                 ),
                 tabPanel(title = "Drug and Gene networks",
                    tabsetPanel(
                      panel3,
                      panel8
                    )         
                 )
                 )

#info to display on popup modals
info <- "Both CRISPR knockouts and drug administrations show various responses in terms of inhibitory effects on cancer cell lines.
The following function correlates the effect of different pairs of drug applications and gene knockouts across a few hundred
shared cell lines. The goal is to understand which drugs behave like which CRISPR knockouts and vice versa for future medical applications."

info1 <- "The following shows the distributions of efficacies versus a specific cancer cell line over all drugs and genes.
negative values indicate an inhibitory effect and thus are considered effective."

info2 <- "This panel displays an interactive network between drugs and genes. The data used are pairs above a 5 percent 
significance level from the target treatments. The type of the nodes and the values of the 
edges(correlation between a gene and drug) are indicated in the key, with positive correlations indicated by blue edges and
negative by red. For a specific interactive network based on querying, please refer to 'correlaton query'."

info3 <- "Like cell lines, this panel displays the efficacy of a specific drug and gene across all cell lines. 
Clicking show correlation will display a scatterplot of all the common cell lines and their corresponding drug and gene
efficacy."

info4 <- "This tab is similar to the first but shows codependency between two drugs."
info5 <- "This tab is similar to the first but shows codependency between two genes."

info6 <- "This panel is similar to the first panel but shows the spearman correlation tables"

info7 <- "Shows the network of high value pairs from the entire data including drug and gene codependencies. 
The correlation cutoff is 0.7, indicating high correlation between the pairs. Some of the pairs from Target pairs
may be included in this network as well."


#server functions and variables.
server <- function(input, output){
  observeEvent(input$learnmore, shinyalert(
    title = "CRISPR_gene_knockout vs Drug response",
    text = info,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  observeEvent(input$learnmore1, shinyalert(
    title = "cell line efficacies",
    text = info1,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  observeEvent(input$networklm, shinyalert(
    title = "Target Pair networks",
    text = info2,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  observeEvent(input$networklm1, shinyalert(
    title = "High value networks",
    text = info7,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  observeEvent(input$learnmore2, shinyalert(
    title = "Gene and Drug efficacies across Cell Lines",
    text = info3,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  observeEvent(input$learnmore3, shinyalert(
    title = "Drug efficacies across Cell Lines",
    text = info4,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  observeEvent(input$learnmore4, shinyalert(
    title = "Gene efficacies across Cell Lines",
    text = info5,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  observeEvent(input$learnmore5, shinyalert(
    title = "Spearman Correlations",
    text = info6,
    size = "m", 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = FALSE,
    type = "success",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "OK",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ))
  #select all nodes of either drug or gene on the first panel in gene and drug networks
  observeEvent(input$selectNodes,
    {
      type <- input$selectNodes
      selected <- returnNodeIds(valuable_targets, type = type)
      visNetworkProxy("network")%>%
        visSelectNodes(id = selected)
    }
  )
  #select all nodes on the second panel
  observeEvent(input$selectNodes2,
    {
      type <- input$selectNodes2
      selected <- returnNodeIds2(corlist_0.7, type = type)
      visNetworkProxy("network1")%>%
      visSelectNodes(id = selected)
    }
  )
  #used mainly by the first and second panels to pick genes and drugs
  data1 <- eventReactive(input$update, {c(input$drugselect, input$geneselect, input$corrselect, input$greaterorless)})
  data6 <- eventReactive(input$updatesp, {c(input$drugselect4, input$geneselect4, input$corrselect1, input$greaterorless1)})
  data7 <- eventReactive(input$showspecnet, {c(input$drugselect, input$geneselect, input$corrselect, input$greaterorless)})
  data8 <- eventReactive(input$showspecnet1, {c(input$drugselect, input$geneselect, input$corrselect, input$greaterorless)})
  data2 <- eventReactive(input$update1, {c(input$cellselect, input$genesordrugs)})
  d1 <- eventReactive(input$druggraph, {input$drugselect1})
  g1 <- eventReactive(input$genegraph, {input$geneselect1})
  d2 <- eventReactive(input$druggraph1, {input$drugselect2})
  d3 <- eventReactive(input$druggraph2, {input$drugselect3})
  g2 <- eventReactive(input$genegraph1, {input$geneselect2})
  g3 <- eventReactive(input$genegraph2, {input$geneselect3})
  data3 <- eventReactive(input$showcorr, {c(input$drugselect1, input$geneselect1)})
  data4 <- eventReactive(input$showcorr1, {c(input$drugselect2, input$drugselect3)})
  data5 <- eventReactive(input$showcorr2, {c(input$geneselect2, input$geneselect3)})
  #network outputs
  output$network <- renderVisNetwork({returnNet(valuable_targets)})
  output$network1 <- renderVisNetwork({returnNet2(corlist_0.7)})
  output$specnet <- renderVisNetwork({returnNet(returnTable(data7(), target_correlations))})
  output$specnet1 <- renderVisNetwork({returnNet(returnTable(data8(), target_correlations2))})
  #correlation query outputs
  output$targettable <- DT::renderDataTable({returnTable(data1(), target_correlations)})
  output$targettable1 <- DT::renderDataTable({returnTable(data6(), target_correlations2)})
  #panel one and two outputs
  output$cellhist <- renderPlot({returnPlot(data2(), update_drug, update_gene, input$bins)})
  output$drugplot <- renderPlot({returnPlot3(update_drug, d1(), type = "drug", bins = input$bins1)})
  output$geneplot <- renderPlot({returnPlot3(update_gene, g1(), type = "gene", bins = input$bins2)})
  output$corrplot <- renderPlot({returnPlot2(data3(), drugsplit, genesplit)})
  output$corrplot2 <- renderPlot({returnPlot2(data3(), drugsplit, genesplit, met = "spearman")})
  output$drugplot1 <- renderPlot({returnPlot3(update_drug, d2(), type = "drug", bins = input$bins3)})
  output$drugplot2 <- renderPlot({returnPlot3(update_drug, d3(), type = "drug", bins = input$bins4)})
  output$geneplot1 <- renderPlot({returnPlot3(update_gene, g2(), type = "gene", bins = input$bins5)})
  output$geneplot2 <- renderPlot({returnPlot3(update_gene, g3(), type = "gene", bins = input$bins6)})
  output$corrplot3 <- renderPlot({returnPlot2(data4(), update_drug, update_drug)})
  output$corrplot4 <- renderPlot({returnPlot2(data4(), update_drug, update_drug, met = "spearman")})
  output$corrplot5 <- renderPlot({returnPlot2(data5(), update_gene, update_gene)})
  output$corrplot6 <- renderPlot({returnPlot2(data5(), update_gene, update_gene, met = "spearman")})
}





#run the shiny app
shinyApp(ui = ui, server = server)


