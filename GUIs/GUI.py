#!/usr/bin/env python

############
# Libraries


import sys
import os

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(project_root, "modules"))

import dash
from dash import Dash, html, dash_table, Input, Output, State, dcc
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import base64
from io import BytesIO

from network_builder import NetworkBuilder
from blastrunner import BLASTrunner
from cazy_fetcher import CazyDownloader

from layout_components.sidebar import sidebar_col
from layout_components.main_components import main_components_col


# Constants

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")
TQ_CACHE = os.getenv("TQ_CACHE", "F").lower().startswith("t")

# Defining inputs

# Initially, the input_file is left empty, it will be updated later
input_file = ""
database_prefix = "/workspace/cazy_analysis/output/data_base/mi_base_de_datos"
output_file = "/workspace/cazy_analysis/output/blast_output/blast_results.txt"

# Functions

# Define the function to run BLAST
def run_blast(input_file):
    BLAST = BLASTrunner(verbose=TQ_VERBOSE)
    BLAST.make_blast_db(input_file)
    result = BLAST.run_blastp(input_file, database_prefix, output_file)

    return result

# Definición de la aplicación Dash
app = Dash(
    __name__,
    external_stylesheets=[dbc.themes.SANDSTONE]
)

# Layout de la aplicación
app.layout = dbc.Container(
    fluid=True,
    style={"height": "100vh", "fontSize": "12px"},
    children=[
        dcc.Store(id='blast_results_store'),
        dbc.Row(
            [
                sidebar_col,
                main_components_col,
            ]
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
        Output("seq_dropdown", "options"),
    ],
    Input("load_sequences_button", "n_clicks"),
    State("cazy_family_dropdown", "value"),
    State("upload_external_fasta", "contents"),
    State("input_sequence_text", "value"),
    prevent_initial_call=True,
)
def load_sequences(n_clicks, selected_family, external_fasta_contents, input_sequence_text):
    if not selected_family:
        return "Please select a CAZy family.", []

    email = "your_email@example.com"  # Cambiar por un correo válido
    downloader = CazyDownloader(email, selected_family, verbose=TQ_VERBOSE)

    # Rutas de archivos temporales
    temp_dir = "/workspace/cazy_analysis/input/temp/"
    os.makedirs(temp_dir, exist_ok=True)
    fasta_family = os.path.join(temp_dir, "family_sequences.fasta")
    fasta_external = os.path.join(temp_dir, "external_sequences.fasta")
    fasta_input = os.path.join(temp_dir, "input_sequence.fasta")
    fasta_combined = os.path.join(temp_dir, "updated_sequences.fasta")

    # Paso 1: Descargar secuencias de la familia
    genbank_ids = downloader.obtener_genbank_ids()
    if not genbank_ids:
        return "No GenBank IDs found for the selected family.", []

    downloader.obtener_fasta_ncbi(output_file=fasta_family)

    # Paso 2: Manejar archivo externo (si existe)
    external_files = []
    if external_fasta_contents:
        content_type, content_string = external_fasta_contents.split(",")
        with open(fasta_external, "wb") as f:
            f.write(base64.b64decode(content_string))
        external_files.append(fasta_external)

    # Paso 3: Manejar secuencia ingresada como texto
    if input_sequence_text:
        with open(fasta_input, "w") as f:
            f.write(input_sequence_text)
        external_files.append(fasta_input)

    # Paso 4: Combinar todos los archivos en uno solo
    with open(fasta_combined, "w") as outfile:
        with open(fasta_family, "r") as infile:
            outfile.write(infile.read())
            outfile.write("\n")
        for ext_file in external_files:
            if os.path.exists(ext_file):
                with open(ext_file, "r") as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")

    # Actualizar secuencias cargadas desde el archivo combinado
    #downloader.agregar_secuencia_externa(fasta_combined, output_file=fasta_combined)
    loaded_sequences = downloader.get_loaded_sequences()

    # Crear opciones para el dropdown
    options = []
    for seq in loaded_sequences:
        seq_id = seq.split()[0]
        if seq_id.startswith(">"):
            seq_id = seq_id[1:]
        options.append({"label": seq_id, "value": seq_id})

    status_message = "Sequences loaded successfully."
    if external_files:
        status_message += " External sequences included."

    return status_message, options

@app.callback(
    [
        Output("blast_status", "children"),
        Output("table", "data"),  # Update the DataTable data
        Output("table", "columns"),  # Update the DataTable columns
        Output("blast_results_store", "data"),  # Store the DataFrame
    ],
    [Input("blast_button", "n_clicks"),
     Input("seq_dropdown", "value")],  # Add dropdown as an input
    State("upload_external_fasta", "contents"),
    State("input_sequence_text", "value"),
    prevent_initial_call=True,
)
def run_blast_callback(n_clicks, selected_sequences, external_fasta_contents, input_sequence_text):
    input_file = "/workspace/cazy_analysis/input/temp/updated_sequences.fasta"
    output_file = "/workspace/cazy_analysis/output/blast_output/blast_results.txt"

    # Intenta ejecutar el BLAST
    df_blast = run_blast(input_file)

    # Si el BLAST se ejecutó y hay resultados
    if df_blast is not None and not df_blast.empty:
        # Filtrar el DataFrame basado en las secuencias seleccionadas y el porcentaje de identidad
        if selected_sequences:
            filtered_df = df_blast[df_blast['query_id'].isin(selected_sequences) & (df_blast['identity'] > 50)]
        else:
            filtered_df = df_blast  # No filtrar si no se seleccionan secuencias

        columns = [{"name": col, "id": col} for col in filtered_df.columns]
        return "BLAST executed successfully.", filtered_df.to_dict("records"), columns, df_blast.to_dict('records')

    # Si el BLAST no se ejecutó, intenta cargar los resultados desde el archivo
    try:
        # Cargar el DataFrame desde el archivo de salida de BLAST
        df_blast = pd.read_csv(output_file, sep="\t")

        # Filtrar el DataFrame basado en las secuencias seleccionadas y el porcentaje de identidad
        if selected_sequences:
            filtered_df = df_blast[df_blast['query_id'].isin(selected_sequences) & (df_blast['identity'] > 60)]
        else:
            filtered_df = df_blast  # No filtrar si no se seleccionan secuencias

        columns = [{"name": col, "id": col} for col in filtered_df.columns]
        return "No new sequences to BLAST, displaying previous results.", filtered_df.to_dict("records"), columns, df_blast.to_dict('records')

    except FileNotFoundError:
        return "Error: BLAST results file not found.", [], [], []
    except Exception as e:
        return f"Error loading previous results: {str(e)}", [], [], []



# Callback to update the heatmap
@app.callback(
    Output("heatmap", "src"),
    [Input("table", "data")],  # Use the DataTable data as input
    State("blast_results_store", "data"),  # Get the stored DataFrame
    prevent_initial_call=True,
)
def update_heatmap(filtered_data, blast_results):
    if not filtered_data or blast_results is None or len(blast_results) == 0:
        return ""  # Hide the graph if there is no data

    # Convert to DataFrame
    filtered_df = pd.DataFrame(filtered_data)

    if filtered_df.empty:
        return ""

    # Get unique queries and subjects
    selected_queries = filtered_df["query_id"].unique()
    selected_subjects = filtered_df["subject_id"].unique()

    # Filter the original DataFrame
    heatmap_df = pd.DataFrame(blast_results)[
        pd.DataFrame(blast_results)["query_id"].isin(selected_queries) & 
        pd.DataFrame(blast_results)["subject_id"].isin(selected_subjects)
    ]

    # Create pivot table
    pivot_df = heatmap_df.pivot_table(
        index="query_id", columns="subject_id", values="identity", aggfunc="mean"
    ).fillna(0)

    # Set figure size dynamically based on the number of unique queries and subjects
    plt.figure(figsize=(max(10, len(selected_subjects) * 0.5), max(6, len(selected_queries) * 0.5)))

    # Create the heatmap
    heatmap = sns.heatmap(pivot_df, annot=True, fmt=".1f", cmap="Greens", cbar_kws={'label': 'Porcentaje de Identidad'})

    # Rotate x-axis labels
    heatmap.set_xticklabels(heatmap.get_xticklabels(), rotation=45)

    # Move protein names to the top
    plt.tick_params(axis='x', which='both', bottom=False, top=True, labelbottom=False, labeltop=True)

    # Italicize protein names
    for tick in heatmap.get_yticklabels():
        tick.set_fontstyle("italic")
    for tick in heatmap.get_xticklabels():
        tick.set_fontstyle("italic")

    # Adjust layout
    plt.tight_layout()

    # Save the heatmap to a buffer
    buf = BytesIO()
    plt.savefig(buf, format="png", dpi=300, bbox_inches='tight')
    plt.close()
    buf.seek(0)

    # Encode the image in base64
    encoded_image = base64.b64encode(buf.read()).decode('utf-8')

    return f"data:image/png;base64,{encoded_image}"


# Callback to generate the network
@app.callback(
    Output("network-graph", "elements"),
    Input("generate_network_button", "n_clicks"),
    State("blast_results_store", "data"),
    State("seq_dropdown", "value"),
    prevent_initial_call=True,
)
def generate_network(n_clicks, blast_results, selected_sequences):
    if blast_results is None or len(blast_results) == 0:
        return []  # No results to display

    # Convert to DataFrame
    df_blast = pd.DataFrame(blast_results)

    if selected_sequences:
        # Filter the DataFrame based on selected sequences
        df_blast = df_blast[df_blast['query_id'].isin(selected_sequences)]

    # Create an instance of NetworkBuilder
    nb = NetworkBuilder()

    # Analyze the results
    analyzed_df = nb.analyze_blast_results(df_blast)

    # Build the network
    network = nb.build_network(analyzed_df, threshold=50)

    # Prepare elements for Cytoscape
    elements = []
    selected_color = "#006400"  # Color for selected sequences
    unselected_color = "#98fb98"  # Color for unselected sequences

    for node in network.nodes():
        # Check if the node is in the selected sequences
        color = selected_color if node in selected_sequences else unselected_color
        elements.append({'data': {'id': node, 'label': node, 'color': color}})

    for edge in network.edges(data=True):
        query, subject, data = edge
        weight = data.get("weight", 1)
        elements.append({'data': {'source': query, 'target': subject, 'weight': weight}})

    return elements


# Run the app
if __name__ == "__main__":
    app.run(debug=True)
