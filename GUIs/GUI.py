#!/usr/bin/env python

############
# Libraries


import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(project_root, "modules"))

import dash
from dash import Dash, html, dash_table, Input, Output, State, dcc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import base64
from io import BytesIO
from blastrunner import BLASTrunner
from cazy_fetcher import CazyDownloader

# Constants

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")
TQ_CACHE = os.getenv("TQ_CACHE", "F").lower().startswith("t")

# Defining inputs

# Initially, the input_file is left empty, it will be updated later
input_file = ""
database_prefix = "/workspace/cazy_analysis/output/data_base/mi_base_de_datos"
output_file = "/workspace/cazy_analysis/output/blast_output/blast_results.txt"

# Functions


# Global variable to store blast results
# df_blast = pd.DataFrame()
# Define the function to run BLAST
def run_blast(input_file):
    BLAST = BLASTrunner(verbose=TQ_VERBOSE)
    BLAST.make_blast_db(input_file)
    result = BLAST.run_blastp(input_file, database_prefix, output_file)

    # print("BLAST Result:", result)  # Verifica la salida
    # print("Type:", type(result))  # ¿Es un DataFrame?

    return result


# Definición de la aplicación Dash
app = Dash(__name__)

# Layout de la aplicación
app.layout = html.Div(
    style={"display": "flex", "height": "100vh"},
    children=[
        # Sidebar
        html.Div(
            [
                html.H3("CAZy Family"),
                dcc.Dropdown(
                    id="cazy_family_dropdown",
                    options=[
                        {"label": "GH180", "value": "GH180"},
                        {"label": "GH51", "value": "GH51"},
                        # Add more families as needed
                    ],
                    placeholder="Select a CAZy family",
                ),
                html.Br(),
                html.Br(),
                dcc.Textarea(
                    id="input_sequence_text",
                    placeholder="Enter your sequence here...",
                    style={"width": "100%", "height": 100},
                ),
                dcc.Upload(
                    id="upload_external_fasta",
                    children=html.Button("Upload Sequence File"),
                    multiple=False,
                ),
                html.Div(id="upload_status"),  # New div for upload status
                html.Br(),
                html.Br(),
                html.Button("Load Sequences", id="load_sequences_button"),
                html.Div(id="output_fasta_status"),  # To show the download status
                html.Br(),
                html.Br(),
                html.Button("Blast", id="blast_button"),
                html.Div(id="blast_status"),  # New div for BLAST status
                html.H3("Filter Options"),
                dcc.Dropdown(
                    id="seq_dropdown",
                    options=[],  # Initially empty, will be filled later
                    placeholder="Choose a protein",
                    multi=True,
                ),
            ],
            style={
                "width": "25%",
                "padding": "20px",
                "backgroundColor": "#f8f9fa",
                "height": "100vh",
                "boxShadow": "2px 0px 5px rgba(0,0,0,0.1)",
            },
        ),
        # Table + Heatmap Container
        html.Div(
            [
                html.H3("BLAST Results"),
                dash_table.DataTable(
                    id="table",
                    columns=[],
                    data=[],  # Initially empty, will be filled with the callback
                    style_table={'overflowX': 'auto'},
                    page_size=10,
                ),
                html.H3("Heatmap of Identity Percentages"),
                html.Img(id="heatmap", style={"width": "100%"}),
            ],
            style={
                "width": "75%",
                "padding": "20px",
                "overflowY": "auto",
            },
        ),
    ],
)


# Callback to show upload status
@app.callback(
    Output("upload_status", "children"),
    Input("upload_external_fasta", "contents"),
    prevent_initial_call=True,
)
def show_upload_status(contents):
    if contents:
        return "uploaded external sequence"
    return ""


# Callback to load sequences based on user input
@app.callback(
    [
        Output("output_fasta_status", "children"),
        Output("seq_dropdown", "options"),  # Update dropdown options
    ],
    Input("load_sequences_button", "n_clicks"),  # Button to load sequences
    State("cazy_family_dropdown", "value"),  # Selected CAZy family
    State("upload_external_fasta", "contents"),  # Uploaded external sequence
    State("input_sequence_text", "value"),  # Text input for sequence
    prevent_initial_call=True,
)
def load_sequences(
    n_clicks, selected_family, external_fasta_contents, input_sequence_text
):

    if selected_family is None:
        return "Please select a CAZy family.", []

    email = "your_email@example.com"  # Change this to your real email
    downloader = CazyDownloader(email, selected_family, verbose=TQ_VERBOSE)

    # Step 1: Obtain FASTA sequences for the selected family
    genbank_ids = downloader.obtener_genbank_ids()
    if not genbank_ids:
        return "No GenBank IDs found for the selected family.", []

    # Save the FASTA sequences to the temp directory
    fasta_file_temp = "/workspace/cazy_analysis/input/temp/family_sequences.fasta"
    downloader.obtener_fasta_ncbi(output_file=fasta_file_temp)

    # Step 2: Handle the upload of the external file if provided
    external_fasta_file = None  # Initialize the variable

    if external_fasta_contents:
        content_type, content_string = external_fasta_contents.split(",")
        decoded = base64.b64decode(content_string)
        external_fasta_file = (
            "/workspace/cazy_analysis/input/temp/external_sequences.fasta"
        )

        with open(external_fasta_file, "wb") as f:
            f.write(decoded)

    # Step 3: Handle the text input for sequences
    if input_sequence_text:
        input_sequence_file = "/workspace/cazy_analysis/input/temp/input_sequence.fasta"
        with open(input_sequence_file, "w") as f:
            f.write(input_sequence_text)

        # Set the external fasta file to the input sequence file
        external_fasta_file = input_sequence_file

    # Step 4: Combine the FASTA files
    combined_fasta_file = "/workspace/cazy_analysis/input/temp/updated_sequences.fasta"

    # Ensure the temp directory exists
    temp_dir = os.path.dirname(combined_fasta_file)
    os.makedirs(
        temp_dir, exist_ok=True
    )  # Create the temp directory if it doesn't exist

    with open(combined_fasta_file, "w") as outfile:
        # Write the family sequences
        with open(fasta_file_temp) as infile:
            outfile.write(infile.read())
            outfile.write("\n")

        # If an external fasta file exists, write its content
        if external_fasta_file and os.path.exists(external_fasta_file):
            with open(external_fasta_file) as infile:
                outfile.write(infile.read())
                outfile.write("\n")

    # Update the downloader with the combined sequences
    downloader.agregar_secuencia_externa(
        external_fasta_file, output_file=combined_fasta_file
    )

    # Step 5: Update dropdown options with new sequences
    loaded_sequences = downloader.get_loaded_sequences()

    # Debug: Print loaded sequences to check their format
    # print("Loaded sequences:", loaded_sequences)

    # Extract only the sequence IDs from the loaded sequences
    options = []
    for seq in loaded_sequences:
        # Split the sequence on spaces and take the first part (the ID)
        seq_id = seq.split()[0]  # This should give us the ID including '>'
        if seq_id.startswith(">"):
            seq_id = seq_id[1:]  # Remove the '>' character
        options.append({"label": seq_id, "value": seq_id})

    # Show different messages depending on whether an external sequence was uploaded
    status_message = "File/sequence uploaded and combined successfully."

    return status_message, options

# Callback para ejecutar BLAST cuando se hace clic en el botón
@app.callback(
    [
        Output("blast_status", "children"),
        Output("table", "data"),  # Update the DataTable data
        Output("table", "columns"),  # Update the DataTable columns
    ],
    [Input("blast_button", "n_clicks"),
     Input("seq_dropdown", "value")],  # Add dropdown as an input
    State("upload_external_fasta", "contents"),
    State("input_sequence_text", "value"),
    prevent_initial_call=True,
)
def run_blast_callback(n_clicks, selected_sequences, external_fasta_contents, input_sequence_text):
    input_file = "/workspace/cazy_analysis/input/temp/updated_sequences.fasta"
    df_blast = run_blast(input_file)

    if df_blast is not None and not df_blast.empty:
        # Filter the DataFrame based on selected sequences and identity percentage
        if selected_sequences:
            # Assuming 'query_id' is the column that contains the sequence IDs
            filtered_df = df_blast[df_blast['query_id'].isin(selected_sequences) & (df_blast['identity'] > 60)]
        else:
            filtered_df = df_blast  # No filtering if no sequences are selected

        columns = [{"name": col, "id": col} for col in filtered_df.columns]
        return "BLAST executed successfully.", filtered_df.to_dict("records"), columns
    else:
        return "Error executing BLAST.", [], []

# Run the app
if __name__ == "__main__":
    app.run(debug=True)
