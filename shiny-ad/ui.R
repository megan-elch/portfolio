library(shiny)
library(bslib)

##########################################################################
# UI: single-page layout with a sidebar for user inputs
##########################################################################
ui <- page_sidebar(
  title = "AD Proteomics Explorer",
  
  ### inputs  #############################################################
  sidebar = sidebar(
    selectInput("cell_type", "Cell Type", choices = NULL),     
    selectInput("protein", "Protein", choices = NULL),          
    radioButtons("volcano_type", "Volcano",
                 choices = c("AD vs Non-AD" = "all",
                             "Braak 6 vs Braak 0" = "b16")),  
    selectInput("cell_type_2", "Compare Cell Type", choices = NULL) 
  ),
  
  ### row 1: abundance trajectory and volcano plot  #######################
  layout_columns(
    card(card_header("Abundance Across Braak Stages"),
         plotOutput("path_plot", height = "400px")),
    card(card_header("Volcano Plot"),
         plotOutput("volcano_plot", height = "400px"))
  ),
  
  ### row 2: cross-cell-type FC comparison and per-patient abundance  #######
  layout_columns(
    card(card_header("log2FC: Cell Type Comparison"),
         plotOutput("fc_comparison_plot", height = "400px")),
    card(card_header("Abundance by Patient"),
         plotOutput("patient_plot", height = "400px"))
  )
)