library(shiny)
library(reticulate)
library(DT)  # For displaying DataFrames as tables

# Configure Python (adjust the path if using a virtual environment)
use_python("/usr/bin/python3")

# Import Python modules
cazyfetcher <- import_from_path("cazy_fetcher", path="/workspace/cazy_analysis/modules")


# Define available families
families <- c("GH8", "GH51", "GH180")

# UI: User interface
ui <- fluidPage(
    titlePanel("Cazy Family BLAST"),
    sidebarLayout(
        sidebarPanel(
            selectInput("familia", "Choose a CAZy Family", choices = families),
            actionButton("fetch_sequences", "Fetch Sequences"),
            textAreaInput("external_sequences", "Enter External Sequences (FASTA format)", rows = 5, placeholder = ">seq1\nATGC...\n>seq2\nATGCG..."),
            actionButton("add_external_sequences", "Add External Sequences"),
            fileInput("external_file", "Upload External Sequences File (FASTA format)", accept = c(".fasta", ".fa")),
            textOutput("fetch_message")  # Success message
        ),
        mainPanel(
            verbatimTextOutput("fetch_output"),
        )
    )
)

# Server: Application logic
server <- function(input, output, session) {
    downloader <- NULL
    cazy_sequences_file <- "/workspace/cazy_analysis/output/cazy_sequences/sequences.fasta"  # Updated path
    combined_sequences_file <- "/workspace/cazy_analysis/output/cazy_sequences/updated_sequences.fasta"  # Updated path

    observeEvent(input$fetch_sequences, {
        output$fetch_message <- renderText({  # Success message
            tryCatch({
                downloader <<- cazyfetcher$CazyDownloader(email = "your_email@example.com", familia = input$familia, verbose = TRUE)
                downloader$obtener_genbank_ids()
                downloader$obtener_fasta_ncbi()  # This writes to the file but may return NULL

                # Read the sequences from the file
                sequences <- readLines(cazy_sequences_file)

                # Check if sequences were read successfully
                if (length(sequences) == 0) {
                    stop("Error: No sequences were retrieved from the file.")
                }

                # Write the sequences to the output file
                writeLines(sequences, cazy_sequences_file)  # Save to the specified path
                paste("CAZy sequences successfully fetched for the family:", input$familia)
            }, error = function(e) {
                paste("Error:", e$message)
            })
        })
    })

    
    observeEvent(input$add_external_sequences, {
        req(downloader)
        external_sequences <- input$external_sequences
        
        # Debugging: Print the external sequences to check the input
        print(external_sequences)
        
        if (nzchar(external_sequences)) {
            # Check if the input is in the correct FASTA format for amino acids
            if (grepl("^>.*\\n([ACDEFGHIKLMNPQRSTVWY\n]+)$", external_sequences)) {
                # Write the external sequences to a temporary file
                temp_external_file <- tempfile(fileext = ".fasta")
                writeLines(external_sequences, temp_external_file)
                
                # Call the method to add external sequences
                downloader$agregar_secuencia_externa(external_fasta_file = temp_external_file)
                output$fetch_message <- renderText({"External sequences added successfully."})  # Success message
            } else {
                output$fetch_message <- renderText({"Error: The input is not in valid FASTA format."})  # Error message
            }
        } else {
            output$fetch_message <- renderText({"Error: No sequences provided."})  # Error message
        }
    })
    observeEvent(input$external_file, {
        req(downloader)
        external_file_path <- input$external_file$datapath
        downloader$agregar_secuencia_externa(external_fasta_file = external_file_path)
        output$fetch_message <- renderText({"External sequences from the file added successfully."})  # Success message
    })
    

}

# Run the app
shinyApp(ui = ui, server = server)
