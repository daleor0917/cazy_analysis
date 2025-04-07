library(shiny)
library(reticulate)
library(DT)
library(ggplot2)
library(reshape2)

use_python("/usr/bin/python3")

# Import Python modules
cazyfetcher <- import_from_path("cazy_fetcher", path="/workspace/cazy_analysis/modules")
blastrunner <- import_from_path("blastrunner", path="/workspace/cazy_analysis/modules")

families <- c("GH8", "GH51", "GH180")

ui <- fluidPage(
    tags$head(
        tags$style(HTML("
            body {
                overflow-y: hidden;
            }
            .dataTables_wrapper .dataTables_paginate {
                float: right;
            }
            .shiny-output-error {
                color: red;
            }
        "))
    ),
    sidebarLayout(
        sidebarPanel(
            selectInput("familia", "Choose a CAZy Family", choices = families),
            actionButton("fetch_sequences", "Fetch Sequences"),
            textOutput("fetch_message"),
            textAreaInput("external_sequences", "Enter External Sequences (FASTA format)", rows = 5, placeholder = ">seq1\nATGC...\n>seq2\nATGCG..."),
            actionButton("add_external_sequences", "Add External Sequences"),
            fileInput("external_file", "Upload External Sequences File (FASTA format)", accept = c(".fasta", ".fa")),
            actionButton("run_blast", "Run BLAST")
        ),
        mainPanel(
            fluidRow(
                column(12, DT::dataTableOutput("blast_results_table")),
                column(12, plotOutput("heatmap", height = "300px"))
            )
        )
    )
)

server <- function(input, output, session) {
    downloader <- NULL
    cazy_sequences_file <- "/workspace/cazy_analysis/output/cazy_sequences/sequences.fasta"
    combined_sequences_file <- "/workspace/cazy_analysis/output/cazy_sequences/updated_sequences.fasta"
    blast_results <- reactiveVal(NULL)

    observeEvent(input$fetch_sequences, {
        output$fetch_message <- renderText({
            tryCatch({
                downloader <<- cazyfetcher$CazyDownloader(email = "your_email@example.com", familia = input$familia, verbose = TRUE)
                downloader$obtener_genbank_ids()
                downloader$obtener_fasta_ncbi()
                sequences <- readLines(cazy_sequences_file)
                if (length(sequences) == 0) stop("Error: No sequences were retrieved from the file.")
                writeLines(sequences, cazy_sequences_file)
                paste("CAZy sequences successfully fetched for the family:", input$familia)
            }, error = function(e) {
                paste("Error:", e$message)
            })
        })
    })

    observeEvent(input$add_external_sequences, {
        req(downloader)
        external_sequences <- input$external_sequences
        if (nzchar(external_sequences)) {
            if (grepl("^>.*\\n([ACDEFGHIKLMNPQRSTVWY\n]+)$", external_sequences)) {
                temp_external_file <- tempfile(fileext = ".fasta")
                writeLines(external_sequences, temp_external_file)
                downloader$agregar_secuencia_externa(external_fasta_file = temp_external_file)
                output$fetch_message <- renderText({"External sequences added successfully."})
            } else {
                output$fetch_message <- renderText({"Error: The input is not in valid FASTA format."})
            }
        } else {
            output$fetch_message <- renderText({"Error: No sequences provided."})
        }
    })

    observeEvent(input$external_file, {
        req(downloader)
        external_file_path <- input$external_file$datapath
        downloader$agregar_secuencia_externa(external_fasta_file = external_file_path)
        output$fetch_message <- renderText({"External sequences from the file added successfully."})
    })

    observeEvent(input$run_blast, {
        req(downloader)
        input_file <- if (is.null(input$external_file)) {
            cazy_sequences_file
        } else {
            combined_sequences_file
        }

        blast_instance <- blastrunner$BLASTrunner(verbose = TRUE)
        blast_instance$make_blast_db(input_file)
        output_file <- "/workspace/cazy_analysis/output/blast_output/blast_results.txt"
        database_prefix <- "/workspace/cazy_analysis/output/data_base/mi_base_de_datos"
        blast_results_data <- blast_instance$run_blastp(input_file, database_prefix, output_file)

        if (!is.null(blast_results_data)) {
            output$fetch_message <- renderText({"BLAST completed successfully. Results saved."})
            blast_results(blast_results_data)
            output$blast_results_table <- DT::renderDataTable({
                DT::datatable(blast_results(), options = list(pageLength = 5, scrollY = FALSE))
            })
        } else {
            output$fetch_message <- renderText({"Error running BLAST."})
        }
    })

    output$heatmap <- renderPlot({
        req(blast_results())
        df <- as.data.frame(blast_results())
        colnames(df) <- gsub(" ", "_", tolower(colnames(df)))
        if (!all(c("query_id", "subject_id", "identity") %in% colnames(df))) return(NULL)
        df_identity <- df[, c("query_id", "subject_id", "identity")]
        identity_matrix <- reshape2::dcast(df_identity, query_id ~ subject_id, value.var = "identity", fill = 0)
        identity_long <- reshape2::melt(identity_matrix, id.vars = "query_id")

        ggplot(identity_long, aes(x = variable, y = query_id, fill = value)) +
            geom_tile(color = "white") +
            geom_text(aes(label = sprintf("%.1f", value)), size = 3) +
            scale_fill_gradient(low = "#e5f5e0", high = "#238b45", name = "Porcentaje de Identidad") +
            labs(title = "BLAST Identity Heatmap", x = "Subject", y = "Query") +
            theme_minimal(base_size = 12) +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, face = "italic"),
                axis.text.y = element_text(face = "italic"),
                axis.ticks = element_blank(),
                panel.grid = element_blank(),
                plot.title = element_text(hjust = 0.5)
            ) +
            scale_x_discrete(position = "top")
    })
}


shinyApp(ui = ui, server = server)
