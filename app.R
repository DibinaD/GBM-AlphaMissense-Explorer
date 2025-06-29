library(shiny)
library(DT)
library(dplyr)
library(readr)
library(plotly)
library(shinythemes)

# Load the data
mut_data <- read_csv("C:/Users/dibin/Downloads/My_project/annotated_with_alphamissense.csv")

ui <- navbarPage(
  "🧬 GBM Mutation Explorer",
  theme = shinytheme("flatly"),
  
  tabPanel("Explore",
           sidebarLayout(
             sidebarPanel(
               helpText("Filter GBM mutations by gene, AI pathogenicity, and databases."),
               selectInput("gene", "Select Gene(s):", 
                           choices = sort(unique(mut_data$Hugo_Symbol)),
                           selected = NULL, multiple = TRUE),
               sliderInput("score", "AlphaMissense Pathogenicity Score:",
                           min = 0, max = 1, value = c(0, 1)),
               checkboxInput("only_pathogenic", "Only AlphaMissense Pathogenic", value = FALSE),
               checkboxInput("clinvar_only", "Only ClinVar Confirmed", value = FALSE),
               checkboxInput("cosmic_only", "Only COSMIC Listed", value = FALSE),
               downloadButton("downloadData", "Download Filtered Data")
             ),
             mainPanel(
               verbatimTextOutput("summaryStats"),
               tabsetPanel(
                 tabPanel("📊 Table", 
                          DTOutput("mutation_table"),
                          hr(),
                          h4("🧬 Selected Mutation Details:"),
                          tableOutput("mutationDetails")
                 ),
                 tabPanel("📈 Interactive Plot", plotlyOutput("pathogenicity_plot")),
                 tabPanel("🎯 Gene Scores", plotlyOutput("gene_plot"))
               )
             )
           )
  ),
  
  tabPanel("About",
           fluidRow(
             column(12,
                    h3("🔬 Project Overview: GBM Mutations Annotated with AlphaMissense"),
                    p("This interactive Shiny dashboard enables in-depth exploration of somatic mutations in Glioblastoma Multiforme (GBM), integrating data from The Cancer Genome Atlas (TCGA, PanCancer Atlas) with AI-based pathogenicity predictions from DeepMind's AlphaMissense."),
                    br(),
                    h4("🧠 Why This Project?"),
                    p("Glioblastoma is one of the most aggressive brain cancers, characterized by a highly mutated genome and poor prognosis. Identifying and prioritizing pathogenic mutations is key to advancing therapeutic research."),
                    p("Traditional annotation tools like PolyPhen, SIFT, and CADD offer useful insights but often produce conflicting results. AlphaMissense, a deep learning model by DeepMind, introduces a paradigm shift by offering high-confidence predictions on the pathogenicity of all possible missense mutations in the human genome."),
                    br(),
                    h4("💡 Key Features"),
                    tags$ul(
                      tags$li("✅ AI-based annotation using AlphaMissense predictions"),
                      tags$li("✅ Filters for ClinVar-confirmed pathogenic mutations and COSMIC-listed cancer variants"),
                      tags$li("✅ Real-time interactive visualizations with Plotly"),
                      tags$li("✅ Gene-wise mutation scoring with selection of multiple genes"),
                      tags$li("✅ Clean mutation summaries with download support"),
                      tags$li("✅ Built entirely in R with Shiny, DT, dplyr, and plotly")
                    ),
                    br(),
                    h4("🧬 Who Is This For?"),
                    p("This dashboard is intended for cancer researchers, bioinformaticians, and clinicians who wish to:"),
                    tags$ul(
                      tags$li("Explore GBM mutation landscapes"),
                      tags$li("Prioritize mutations based on AI-predicted pathogenicity"),
                      tags$li("Validate findings with ClinVar and COSMIC databases"),
                      tags$li("Extract publication-ready tables or conduct preliminary gene filtering")
                    ),
                    br(),
                    h4("🛠️ Technical Stack"),
                    tags$ul(
                      tags$li("Frontend: R Shiny, shinythemes"),
                      tags$li("Visualization: plotly, DT"),
                      tags$li("Data Processing: dplyr, readr"),
                      tags$li("Data Sources: TCGA (via cBioPortal), AlphaMissense (DeepMind), ClinVar, COSMIC")
                    ),
                    br(),
                    h4("🧑‍💻 Developed by Dibina Dinakaran"),
                    p("For questions, collaborations, or feedback, please visit the GitHub repository or reach out directly. This project was built as a demonstration of integrating AI with cancer genomics to create useful tools for the research community.")
             )
           )
  )
) 

server <- function(input, output, session) {
  
  filtered_data <- reactive({
    df <- mut_data
    
    if (length(input$gene) > 0) {
      df <- df %>% filter(Hugo_Symbol %in% input$gene)
    }
    
    df <- df %>% filter(am_pathogenicity >= input$score[1],
                        am_pathogenicity <= input$score[2])
    
    if (input$only_pathogenic) {
      df <- df %>% filter(am_class == "pathogenic")
    }
    
    if (input$clinvar_only && "CLIN_SIG" %in% names(df)) {
      df <- df %>% filter(grepl("Pathogenic", CLIN_SIG, ignore.case = TRUE))
    }
    
    if (input$cosmic_only && "COSMIC" %in% names(df)) {
      df <- df %>% filter(!is.na(COSMIC) & COSMIC != "")
    }
    
    df
  })
  
  output$summaryStats <- renderPrint({
    df <- filtered_data()
    cat("Total mutations:", nrow(df), "\n")
    cat("Unique genes:", length(unique(df$Hugo_Symbol)), "\n")
    cat("Pathogenic:", sum(df$am_class == "pathogenic", na.rm = TRUE), "\n")
    cat("Benign:", sum(df$am_class == "benign", na.rm = TRUE), "\n")
  })
  
  output$mutation_table <- renderDT({
    df <- filtered_data()
    df %>% select(Hugo_Symbol, Chromosome, Start_Position, End_Position,
                  Variant_Classification, Variant_Type, Reference_Allele, Tumor_Seq_Allele2,
                  HGVSp_Short, AA_Change, Protein_position,
                  am_pathogenicity, am_class,
                  CLIN_SIG, COSMIC,
                  PolyPhen, SIFT,
                  Tumor_Sample_Barcode, t_ref_count, t_alt_count) %>%
      datatable(options = list(pageLength = 10, scrollX = TRUE))
  })
  
  output$mutationDetails <- renderTable({
    selected <- input$mutation_table_rows_selected
    if (length(selected)) {
      df <- filtered_data()[selected, ]
      df %>% select(Hugo_Symbol, Chromosome, Start_Position, Reference_Allele, 
                    Tumor_Seq_Allele2, HGVSp_Short, am_pathogenicity, am_class, CLIN_SIG, COSMIC)
    } else {
      data.frame(Message = "Select a row in the table to view mutation details.")
    }
  })
  
  output$pathogenicity_plot <- renderPlotly({
    df <- filtered_data()
    plot_ly(df, x = ~am_class, type = "histogram",
            marker = list(color = "rgba(55, 128, 191, 0.7)")) %>%
      layout(title = "AlphaMissense Class Distribution",
             xaxis = list(title = "Class"),
             yaxis = list(title = "Mutation Count"))
  })
  
  output$gene_plot <- renderPlotly({
    df <- filtered_data()
    plot_ly(df, x = ~HGVSp_Short, y = ~am_pathogenicity,
            color = ~am_class, type = "scatter", mode = "markers",
            text = ~paste("Transcript:", Transcript_ID)) %>%
      layout(title = "Mutation Scores by Protein Change",
             xaxis = list(title = "Protein Change"),
             yaxis = list(title = "AlphaMissense Score"))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("Filtered_Mutations_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write_csv(filtered_data(), file)
    }
  )
}

shinyApp(ui, server)
