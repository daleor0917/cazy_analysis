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
database_prefix = (
    "/workspace/Test_Python_Jupyter/mi_entorno/docs/data_base/mi_base_de_datos"
)
output_file = (
    "/workspace/Test_Python_Jupyter/mi_entorno/docs/blast_output/blast_results.txt"
)

# Global variable to store blast results
df_blast = pd.DataFrame()

# Define the function to run BLAST
def run_blast(input_file):
    BLAST = BLASTrunner(verbose=TQ_VERBOSE, tq_cache=TQ_CACHE)
    BLAST.make_blast_db(input_file)
    return BLAST.run_blastp(input_file, database_prefix, output_file)

# Defining app

app = Dash()

# App layout
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
                    id='input_sequence_text',
                    placeholder='Enter your sequence here...',
                    style={'width': '100%', 'height': 100}
                ),
                dcc.Upload(
                    id='upload_external_fasta',
                    children=html.Button('Upload Sequence File'),
                    multiple=False
                ),
                html.Div(id='upload_status'),  # New div for upload status
                html.Br(),
                html.Br(),
                html.Button("Load Sequences", id="load_sequences_button"),
                html.Div(id='output_fasta_status'),  # To show the download status
                html.Br(),
                html.Br(),
                html.Button("Blast", id="blast_button"),
                html.Div(id='blast_status'),  # New div for BLAST status
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
                    id="blast_table",
                    columns=[],  # Initially empty, will be filled later
                    page_size=10,
                    style_table={"overflowX": "auto"},
                    data=[],
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
    Output('upload_status', 'children'),
    Input('upload_external_fasta', 'contents'),
    prevent_initial_call=True,
)
def show_upload_status(contents):
    if contents:
        return "uploaded external sequence"
    return ""

# Callback to load sequences based on user input
@app.callback(
    [
        Output('output_fasta_status', 'children'),
        Output('seq_dropdown', 'options'),  # Update dropdown options
        Output('blast_table', 'columns'),  # Update table columns
        #Output('blast_table', 'data'),  # Update table data
    ],
    Input('load_sequences_button', 'n_clicks'),  # Button to load sequences
    State('cazy_family_dropdown', 'value'),  # Selected CAZy family
    State('upload_external_fasta', 'contents'),  # Uploaded external sequence
    State('input_sequence_text', 'value'),  # Text input for sequence
    prevent_initial_call=True,
)
def load_sequences(n_clicks, selected_family, external_fasta_contents, input_sequence_text):
    # Check if a CAZy family is selected
    if selected_family is None:
        return "Please select a CAZy family.", [], []#, []

    # Initialize the downloader with the selected family and verbosity
    email = "your_email@example.com"  # Change this to your real email
    downloader = CazyDownloader(email, selected_family, TQ_VERBOSE)

    # Step 1: Obtain FASTA sequences for the selected family
    genbank_ids = downloader.obtener_genbank_ids()
    if not genbank_ids:
        return "No GenBank IDs found for the selected family.", [], []#, []

    downloader.obtener_fasta_ncbi()

    # Step 2: Handle the upload of the external file if provided
    external_fasta_file = None  # Initialize the variable
    external_sequence_uploaded = False  # Flag to check if an external sequence was uploaded

    if external_fasta_contents:
        content_type, content_string = external_fasta_contents.split(',')
        decoded = base64.b64decode(content_string)
        external_fasta_file = "/workspace/Test_Python_Jupyter/mi_entorno/docs/external_sequences.fasta"
        
        with open(external_fasta_file, "wb") as f:
            f.write(decoded)

        external_sequence_uploaded = True  # Set flag to True

    # Step 3: Handle the text input for sequences
    if input_sequence_text:
        input_sequence_file = "/workspace/Test_Python_Jupyter/mi_entorno/docs/input_sequence.fasta"
        with open(input_sequence_file, "w") as f:
            f.write(input_sequence_text)

        external_fasta_file = input_sequence_file  # Use the input sequence file as the external file
        external_sequence_uploaded = True  # Set flag to True

    # Step 4: Combine the FASTA files
    combined_fasta_file = "/workspace/Test_Python_Jupyter/mi_entorno/docs/updated_sequences.fasta"
    with open(combined_fasta_file, "w") as outfile:
        with open("/workspace/Test_Python_Jupyter/mi_entorno/docs/sequences.fasta") as infile:
            outfile.write(infile.read())
            outfile.write("\n")

        if external_fasta_file and os.path.exists(external_fasta_file):
            with open(external_fasta_file) as infile:
                outfile.write(infile.read())
                outfile.write("\n")

    downloader.agregar_secuencia_externa(external_fasta_file, output_file=combined_fasta_file)

    # Step 5: Update dropdown options with new sequences
    options = [{"label": seq, "value": seq} for seq in downloader.get_loaded_sequences()]

    # Show different messages depending on whether an external sequence was uploaded
    status_message = "File/sequence uploaded successfully."

    return status_message, options, []#, []

# Callback to run BLAST when the button is clicked
@app.callback(
    Output('blast_status', 'children'),  # New output for BLAST status
    Input('blast_button', 'n_clicks'),  # Button to run BLAST
    State('upload_external_fasta', 'contents'),  # Uploaded external sequence
    State('input_sequence_text', 'value'),  # Text input for sequence
    prevent_initial_call=True,
)
def run_blast_callback(n_clicks, external_fasta_contents, input_sequence_text):
    # Determine the input file to use for BLAST
    input_file = "/workspace/Test_Python_Jupyter/mi_entorno/docs/updated_sequences.fasta"

    # Run BLAST
    df_blast = run_blast(input_file)
    print(df_blast.head())

    if df_blast is not None:
        return "BLAST executed successfully."
    else:
        return "Error executing BLAST."

# Callback para filtrar y actualizar la tabla de resultados de BLAST
@app.callback(
    Output('blast_table', 'data'),  # Actualiza los datos de la tabla
    Input('seq_dropdown', 'value'),  # Secuencias seleccionadas en el dropdown
    prevent_initial_call=True
)
def filter_blast_results(selected_sequences):
    global df_blast  # Usa el DataFrame global con los resultados de BLAST

    if df_blast.empty:
        return []

    # Verificar si hay selección
    if not selected_sequences:
        return df_blast.to_dict("records")  # Devuelve todo si no hay selección

    # Filtrar por secuencia y porcentaje de identidad > 60%
    filtered_df = df_blast[
        (df_blast["query_id"].isin(selected_sequences)) & (df_blast["identity"] > 60)
    ]

    return filtered_df.to_dict("records")

# Run the app
if __name__ == "__main__":
    app.run(debug=True)