#!/usr/bin/env python
# coding: utf-8

import dash_bootstrap_components as dbc
from dash import html, dcc


# Sidebar Component
sidebar_col = dbc.Col(
    [
        html.Br(),
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        dbc.CardGroup(
                            [
                                dcc.Dropdown(
                                    id="cazy_family_dropdown",
                                    options=[
                                        {"label": "GH180", "value": "GH180"},
                                        {"label": "GH51", "value": "GH51"},
                                    ],
                                    placeholder="Select a CAZy family",
                                    style={"width": "100%"},
                                ),
                            ]
                        ),
                        html.Br(),
                        dbc.CardGroup(
                            [
                                dbc.Textarea(
                                    id="input_sequence_text",
                                    placeholder="Enter your sequence here...",
                                    style={"height": "60px"},  # Further reduced height
                                    className="mb-2",
                                ),
                            ]
                        ),
                        dbc.CardGroup(
                            [
                                html.Div(
                                    [
                                        dbc.Label(
                                            "or upload file",
                                            html_for="upload_external_fasta",
                                            className="mb-0 me-2",
                                        ),
                                        html.Div(
                                            dcc.Upload(
                                                id="upload_external_fasta",
                                                children=dbc.Button(
                                                    "Choose File", color="secondary", size="sm"  # Smaller button
                                                ),
                                                multiple=False,
                                            ),
                                            style={"display": "inline-block", "marginLeft": "10px"},
                                        ),
                                    ],
                                    className="d-flex align-items-center",
                                ),
                                html.Div(id="upload_status", className="mt-2"),
                            ]
                        ),
                        html.Br(),
                        dbc.Button(
                            "Load Sequences",
                            id="load_sequences_button",
                            color="primary",
                            className="mt-2",  # Reduced margin
                            size="sm"  # Smaller button
                        ),
                        html.Div(id="output_fasta_status", className="mt-2"),
                        html.Hr(),
                        dbc.Button(
                            "Blast", id="blast_button", color="danger", className="mt-2", size="sm"  # Smaller button
                        ),
                        html.Div(id="blast_status", className="mt-2"),
                    ]
                ),
            ],
            className="mb-4",
            style={"boxShadow": "0 0 8px rgba(0, 0, 0, 0.1)", "padding": "10px"},
        ),
        html.Br(),
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H3("Filter Options", className="text-primary mb-3"),
                        dbc.Label("Choose a protein"),
                        dcc.Dropdown(
                            id="seq_dropdown",
                            options=[],
                            placeholder="Choose a protein",
                            multi=True,
                        ),
                        dbc.Button(
                            "Generate Network",
                            id="generate_network_button",
                            color="success",
                            className="mt-2",  # Reduced margin
                            size="sm"  # Smaller button
                        ),
                    ]
                )
            ],
            style={
                "boxShadow": "0 0 8px rgba(0, 0, 0, 0.1)",
                "padding": "10px",
                "backgroundColor": "#ffffff",
            },
        ),
    ],
    width=3,
    style={
        "padding": "10px",  # Reduced padding
        "backgroundColor": "#f8f9fa",
        "resize": "horizontal",
        "min-width": "300px",
        "height": "100vh",
        "boxShadow": "2px 0px 5px rgba(0,0,0,0.1)",
        "overflowY": "auto",
    },
)

