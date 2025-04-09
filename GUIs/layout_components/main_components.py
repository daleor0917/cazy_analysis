#!/usr/bin/env python
# coding: utf-8

import dash_bootstrap_components as dbc
from dash import Dash, html, dash_table, Input, Output, State, dcc
import dash_cytoscape as cyto

#!/usr/bin/env python
# coding: utf-8

import dash_bootstrap_components as dbc
from dash import Dash, html, dash_table, Input, Output, State, dcc
import dash_cytoscape as cyto

main_components_col = dbc.Col(
    [
        # Row for Heatmap and Network Visualization
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
                        style={"margin-bottom": "0px"},  # No margin at the bottom
                    ),
                    width=6,
                ),
                dbc.Col(
                    dbc.Card(
                        dbc.CardBody(
                            [
                                html.H5("Network Visualization", className="text-center"),
                                cyto.Cytoscape(
                                    id='network-graph',
                                    layout={'name': 'cose'},
                                    style={'width': '100%', 'height': '600px'},
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
                    width=6,
                ),
            ],
            style={"padding": "20px", "margin-bottom": "0px"},  
        ),

        # Row for Significant Alignments Table
        dbc.Row(
            dbc.Col(
                dbc.Card(
                    dbc.CardBody(
                        [
                            html.H5("Significant Alignments", className="text-left"),
                            dash_table.DataTable(
                                id="table",
                                columns=[],
                                data=[],
                                style_table={'overflowX': 'auto'},
                                page_size=5,
                                style_cell={'textAlign': 'left'},  # Align text to the left
                                style_header={'backgroundColor': 'lightgrey', 'fontWeight': 'bold'},  # Style header
                            ),
                        ]
                    ),
                ),
                #width=12,  # Full width of the parent column
                style={"padding": "0px"},  # Remove extra padding
            )
        ),
    ],
    width=9,
    style={"overflowY": "auto"},
)
