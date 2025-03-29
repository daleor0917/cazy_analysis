#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import networkx as nx
from pyvis.network import Network

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")
TQ_CACHE = os.getenv("TQ_CACHE", "F").lower().startswith("t")


class NetworkBuilder:

    def __init__(
        self,
        base_dir="/workspace/cazy_analysis/output",
        verbose=TQ_VERBOSE,
        tq_cache=TQ_CACHE,
    ):
        self.verbose = verbose
        self.tq_cache = tq_cache
        self.base_dir = base_dir
        self.blast_output_file = os.path.join(
            self.base_dir, "blast_output", "blast_results.txt"
        )

        # Crear el directorio 'network' si no existe
        self.network_dir = os.path.join(self.base_dir, "network")
        os.makedirs(self.network_dir, exist_ok=True)

        self.network_output_file = os.path.join(self.network_dir, "blast_network.html")

    def analyze_blast_results(self, df):
        """
        Extrae las columnas relevantes del DataFrame de resultados de BLAST.

        :param df: DataFrame que contiene los resultados de BLAST.
        :return: DataFrame con las columnas 'query_id', 'subject_id' e 'identity'.
        """
        result_df = df[["query_id", "subject_id", "identity"]]
        return result_df

    def build_network(self, result_df, threshold=75):
        """
        Construye un grafo a partir de los resultados de BLAST.

        :param result_df: DataFrame con los resultados de BLAST.
        :param threshold: Porcentaje mínimo de identidad para establecer un enlace.
        :return: Grafo de NetworkX.
        """
        G = nx.Graph()

        for _, row in result_df.iterrows():
            query = row["query_id"]
            subject = row["subject_id"]
            identity = row["identity"]

            if identity >= threshold:
                G.add_edge(query, subject, weight=identity)

        return G

    def visualize_network(self, G):
        """
        Visualiza el grafo utilizando pyvis.

        :param G: Grafo de NetworkX.
        """
        net = Network(
            notebook=True,
            height="600px",
            width="100%",
            bgcolor="#222222",
            font_color="white",
        )

        for node in G.nodes():
            net.add_node(node, label=node, color="#00ff1e")

        for edge in G.edges(data=True):
            query, subject, data = edge
            weight = data.get("weight", 1)
            net.add_edge(query, subject, value=weight, title=f"Identity: {weight}%")

        net.show(self.network_output_file)
        print(f"Red guardada en {self.network_output_file}")


if __name__ == "__main__":

    nb = NetworkBuilder()

    # Leer el DataFrame desde el archivo de resultados de BLAST
    df_results = pd.read_csv(nb.blast_output_file, sep="\t")

    # Analizar los resultados
    analyzed_df = nb.analyze_blast_results(df_results)

    # Construir la red
    network = nb.build_network(analyzed_df, threshold=75)

    print(f"Número de nodos: {network.number_of_nodes()}")
    print(f"Número de enlaces: {network.number_of_edges()}")

    nb.visualize_network(network)
