#!/usr/bin/env python
# coding: utf-8

import dash_bootstrap_components as dbc
from dash import Dash, html, dash_table, Input, Output, State, dcc
import dash_cytoscape as cyto

# Main Components Column
# Main Components Column
main_components_col = dbc.Col(
    [
        # Row for Heatmap and Significant Alignments Table
        dbc.Row(
            [
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H5("Heatmap of Identity Percentages", className="text-center"),
                                html.Img(id="heatmap", style={"width": "100%", "border": "1px solid #ccc", "padding": "10px"}),
                            ]
                        ),
                    ),
                    width=6,  # Half width for heatmap
                    style={"padding": "2px"},  # Remove padding for alignment
                ),
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H5("Significant Alignments", className="text-left"),
                                dash_table.DataTable(
                                    id="table",
                                    columns=[],
                                    data=[],
                                    style_table={
                                        'overflowX': 'auto',
                                        'maxHeight': '200px',  # Set a max height for scrolling
                                        'overflowY': 'hidden',  # Hide vertical overflow
                                    },
                                    page_size=7,
                                    style_cell={'textAlign': 'left'},  # Align text to the left
                                    style_header={'backgroundColor': 'lightgrey', 'fontWeight': 'bold'},  # Style header
                                ),
                            ]
                        ),
                    ),
                    width=6,  # Half width for table
                    style={"padding": "2px"},  # Remove padding for alignment
                ),
            ],
            style={"padding": "11px", "margin-bottom": "0px"},  # Remove padding and margin for the row
        ),

        # Row for Network Visualization
        dbc.Row(
            dbc.Col(
                dbc.Card(
                    dbc.CardBody(
                        [
                            html.H5("Network Visualization", className="text-center"),
                            cyto.Cytoscape(
                                id='network-graph',
                                layout={'name': 'cose'},
                                style={'width': '100%', 'height': '400px'},  # Set height for network visualization
                                elements=[],
                                stylesheet=[
                                    {
                                        'selector': 'node',
                                        'style': {
                                            'background-color': 'data(color)',
                                            'label': 'data(label)',
                                        }
                                    },
                                ],
                            ),
                        ]
                    ),
                    style={"margin-bottom": "0px"},  # No margin at the bottom
                ),
                width=12,  # Full width for network visualization
            )
        ),
    ],
    width=9,
    style={"display": "flex", "flexDirection": "column", "height": "100vh"},  # Use flexbox for vertical layout
)
