library(shiny)
library(dplyr)
library(ggplot2)

##########################################################################
# Load data upon startup
##########################################################################
abundances          <- readRDS("data/abundances.rds")
da_all              <- readRDS("data/da_all.rds")       # AD vs Non-AD limma results
da_b16              <- readRDS("data/da_b16.rds")       # Braak VI vs Braak 0 limma results
patient_medians_all <- readRDS("data/patient_medians.rds")  # precomputed per patient x protein x cell type
braak_medians_all   <- readRDS("data/braak_medians.rds")    # precomputed per braak stage x protein x cell type

server <- function(input, output, session) {
  
  ##########################################################################
  # Populate main cell type dropdown 
  ##########################################################################
  observe({
    cell_types <- sort(unique(abundances$celltype_assign))
    updateSelectInput(session, "cell_type", choices = cell_types)
  })
  
  ##########################################################################
  # Populate protein dropdown when cell type is selected
  ##########################################################################
  observe({
    req(input$cell_type)
    proteins <- abundances %>%
      filter(celltype_assign == input$cell_type) %>%
      pull(SYMBOL) %>% unique() %>% sort()
    updateSelectInput(session, "protein", choices = proteins)
  })
  
  ##########################################################################
  # Filter to selected cell type and protein
  # Used by path_plot and patient_plot 
  ##########################################################################
  cell_data <- reactive({
    req(input$cell_type, input$protein)
    abundances %>%
      filter(celltype_assign == input$cell_type,
             SYMBOL == input$protein)
  })
  
  ##########################################################################
  # Plot 1: Violin/path plot of abundance across braak stages
  # braak stage summaries are precomputed to reduce lag
  ##########################################################################
  output$path_plot <- renderPlot({
    req(input$cell_type, input$protein)
    
    df <- cell_data()
    
    braak_colors <- c(
      "0" = "#2166AC", "I" = "#6BAED6", "II" = "#9ECAE1",
      "III" = "#CCCCCC", "IV" = "#FC8D59", "V" = "#D7301F", "VI" = "#990000"
    )
    
    stage_means <- braak_medians_all %>%
      filter(celltype_assign == input$cell_type, SYMBOL == input$protein)
    
    ggplot(df, aes(x = factor(brk), y = obs_log2, fill = factor(brk))) +
      geom_violin(alpha = 0.5, color = NA) +
      geom_jitter(aes(color = factor(brk)), width = 0.15, alpha = 0.15, size = 0.5) +
      scale_fill_manual(values = braak_colors, name = "Braak Stage") +
      scale_color_manual(values = braak_colors, guide = "none") +
      geom_line(data = stage_means, aes(x = factor(brk), y = mean_abund, group = 1),
                color = "black", linewidth = 1, inherit.aes = FALSE) +
      geom_point(data = stage_means, aes(x = factor(brk), y = mean_abund),
                 color = "black", size = 3, inherit.aes = FALSE) +
      scale_x_discrete(limits = c("0", "II", "V", "VI")) +
      labs(x = "Braak Stage", y = "log2 Abundance",
           title = paste(input$protein, "—", input$cell_type)) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom",
            axis.text = element_text(size = 14))
  })
  
  ##########################################################################
  # Plot 2: Volcano plot of differential abundance results
  # Can select AD/Non-AD (braak stages 0-3, 4-6) or Braak6/Braak 0
  # Selected protein highlighted in volcano plot
  ##########################################################################
  output$volcano_plot <- renderPlot({
    da <- if (input$volcano_type == "all") da_all else da_b16
    da <- da %>% filter(celltype_assign == input$cell_type)
    selected <- da %>% filter(SYMBOL == input$protein)
    
    x_label <- if (input$volcano_type == "all") {
      "log2 Fold Change (AD / Non-AD)"
    } else {
      "log2 Fold Change (Braak 6 / Braak 0)"
    }
    
    # orange = upregulated, blue = downregulated
    selected_color <- ifelse(nrow(selected) > 0 && selected$log2fc >= 0, "#E86B1A", "#4A90D9")
    
    ggplot(da, aes(x = log2fc, y = -log10(fdr), color = sig)) +
      geom_point(alpha = 0.5, size = 1) +
      scale_color_manual(values = c("NS" = "gray70", "sig" = "#4A90D9"),
                         guide = "none") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray40") +
      geom_point(data = selected, color = selected_color, size = 4) +
      geom_text(data = selected, aes(label = SYMBOL),
                nudge_y = 0.5, color = selected_color, fontface = "bold", size = 6) +
      labs(x = x_label, y = "-log10(FDR)") +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_text(size = 14))
  })
  
  ##########################################################################
  # Populate second cell type dropdown
  ##########################################################################
  observeEvent(input$cell_type, {
    req(input$cell_type)
    cell_types <- sort(unique(abundances$celltype_assign))
    updateSelectInput(session, "cell_type_2", choices = cell_types,
                      selected = cell_types[which(cell_types != input$cell_type)[1]])
  })
  
  ##########################################################################
  # Plot 3: log2FC comparison across two cell types
  # Each point is a protein; selected protein highlighted
  # Useful for assessing whether fold change patterns are cell-type specific
  ##########################################################################
  output$fc_comparison_plot <- renderPlot({
    req(input$cell_type, input$cell_type_2)
    da <- if (input$volcano_type == "all") da_all else da_b16
    
    ct1 <- da %>% filter(celltype_assign == input$cell_type) %>%
      select(SYMBOL, log2fc) %>% rename(fc_1 = log2fc)
    ct2 <- da %>% filter(celltype_assign == input$cell_type_2) %>%
      select(SYMBOL, log2fc) %>% rename(fc_2 = log2fc)
    
    joined   <- inner_join(ct1, ct2, by = "SYMBOL")
    selected <- joined %>% filter(SYMBOL == input$protein)
    
    selected_color <- ifelse(nrow(selected) > 0 && selected$fc_1 >= 0, "#E86B1A", "#4A90D9")
    
    ggplot(joined, aes(x = fc_1, y = fc_2)) +
      geom_point(alpha = 0.3, size = 1, color = "gray60") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
      geom_point(data = selected, color = selected_color, size = 4) +
      geom_text(data = selected, aes(label = SYMBOL),
                nudge_y = 0.03, color = selected_color, fontface = "bold", size = 6) +
      labs(x = paste("log2FC —", input$cell_type),
           y = paste("log2FC —", input$cell_type_2)) +
      theme_minimal(base_size = 14) +
      theme(axis.text = element_text(size = 14))
  })
  
  ##########################################################################
  # Plot 4: Violin/path plot of abundance across patients
  # Patient violins ordered and colored by Braak stage
  # medians are also pre-computed
  ##########################################################################
  output$patient_plot <- renderPlot({
    req(input$cell_type, input$protein)
    
    df <- cell_data()
    
    braak_colors <- c(
      "0" = "#2166AC", "I" = "#6BAED6", "II" = "#9ECAE1",
      "III" = "#CCCCCC", "IV" = "#FC8D59", "V" = "#D7301F", "VI" = "#990000"
    )
    braak_order <- c("0", "I", "II", "III", "IV", "V", "VI")
    
    patient_medians <- patient_medians_all %>%
      filter(celltype_assign == input$cell_type, SYMBOL == input$protein)
    
    # order patients by Braak stage so x-axis reflects disease progression
    patient_order <- patient_medians %>%
      distinct(patient, brk) %>%
      mutate(brk = factor(brk, levels = braak_order)) %>%
      arrange(brk) %>%
      pull(patient)
    
    df              <- df %>% mutate(patient = factor(patient, levels = patient_order))
    patient_medians <- patient_medians %>%
      mutate(patient = factor(patient, levels = patient_order))
    
    ggplot(df, aes(x = patient, y = obs_log2, fill = factor(brk))) +
      geom_violin(alpha = 0.5, color = NA) +
      geom_point(data = patient_medians, aes(x = patient, y = med, color = factor(brk)),
                 size = 2.5, inherit.aes = FALSE) +
      geom_line(data = patient_medians, aes(x = patient, y = med, group = 1),
                color = "black", linewidth = 0.5, inherit.aes = FALSE) +
      scale_fill_manual(values = braak_colors, name = "Braak Stage") +
      scale_color_manual(values = braak_colors, guide = "none") +
      labs(x = "Patient", y = "log2 Abundance",
           title = paste(input$protein, "—", input$cell_type)) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
            legend.position = "bottom")
  })
}