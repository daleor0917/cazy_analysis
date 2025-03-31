#!/usr/bin/env python
# coding: utf-8

import os
import requests
from bs4 import BeautifulSoup
from Bio import Entrez

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")


class CazyDownloader:
    def __init__(
        self,
        email,
        familia,
        output_dir="/workspace/cazy_analysis/output/cazy_sequences",
        verbose=TQ_VERBOSE,
    ):
        """
        Inicializa el downloader con el email para NCBI, la familia de CAZy y el modo verbose.
        """
        Entrez.email = email  # Configurar email para NCBI Entrez
        self.familia = familia
        self.genbank_ids = []
        self.verbose = verbose  # Guardar el estado de verbose
        self.loaded_sequences = []  # Lista para almacenar las secuencias cargadas
        self.output_dir = output_dir  # Directorio de salida

        # Crear el directorio de salida si no existe
        os.makedirs(self.output_dir, exist_ok=True)

    def _print(self, message):
        """Imprime el mensaje si verbose está activado."""
        if self.verbose:
            print(message)

    def obtener_genbank_ids(self):
        """
        Obtiene los identificadores GenBank de una página de familia en CAZy.
        """
        url = f"https://www.cazy.org/{self.familia}_all.html"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, "html.parser")

        # Flag to skip the header row
        header_skipped = False

        for row in soup.find_all("tr"):
            cols = row.find_all("td")
            if len(cols) > 3:  # GenBank ID está en la 4ta columna
                genbank_id = cols[3].text.strip()
                if genbank_id:
                    # Skip the header row
                    if not header_skipped and genbank_id == "GenBank":
                        header_skipped = True
                        continue
                    self.genbank_ids.append(genbank_id)

        self._print(f"GenBank IDs obtenidos: {self.genbank_ids}")
        return self.genbank_ids

    def obtener_fasta_ncbi(self, output_file="sequences.fasta"):
        """
        Descarga las secuencias FASTA desde NCBI usando los IDs de GenBank.
        """
        output_path = os.path.join(self.output_dir, output_file)
        with open(output_path, "w") as fasta_file:
            for genbank_id in self.genbank_ids:
                try:
                    handle = Entrez.efetch(
                        db="protein", id=genbank_id, rettype="fasta", retmode="text"
                    )
                    fasta_sequence = handle.read()
                    handle.close()

                    if fasta_sequence.startswith(">"):
                        fasta_file.write(fasta_sequence + "\n")
                        self.loaded_sequences.append(
                            fasta_sequence.splitlines()[0][1:]
                        )  # Store the sequence ID
                    else:
                        self._print(f"No se encontró la secuencia para {genbank_id}")

                except Exception as e:
                    self._print(f"Error con {genbank_id}: {e}")

    def agregar_secuencia_externa(
        self, external_fasta_file=None, output_file="updated_sequences.fasta"
    ):
        """
        Agrega una secuencia externa en formato FASTA a las secuencias descargadas.
        Si no se pasa un archivo externo, solo guarda las secuencias de NCBI.
        """
        try:
            combined_file_path = os.path.join(self.output_dir, output_file)
            with open(combined_file_path, "w") as combined_file:
                # Escribir las secuencias de NCBI
                with open(
                    os.path.join(self.output_dir, "sequences.fasta"), "r"
                ) as fasta_file:
                    combined_file.write(fasta_file.read())

                # Si no se proporciona un archivo externo, no intentamos agregarlo
                if external_fasta_file:
                    if os.path.exists(external_fasta_file):
                        with open(external_fasta_file, "r") as external_file:
                            combined_file.write(external_file.read())
                        self._print(
                            f"Secuencia externa añadida desde {external_fasta_file}."
                        )
                    else:
                        self._print(
                            f"No se encontró el archivo externo: {external_fasta_file}"
                        )
                else:
                    self._print("No se proporcionó un archivo externo.")

            # Mensaje de confirmación siempre se imprime
            self._print(f"Secuencias combinadas guardadas en {combined_file_path}.")

        except Exception as e:
            self._print(f"Error al combinar secuencias: {e}")

    def get_loaded_sequences(self, combined_file_path="/workspace/cazy_analysis/input/temp/updated_sequences.fasta"):
        """
        Devuelve las secuencias cargadas para ser utilizadas en el dropdown.
        Lee las secuencias del archivo combinado y devuelve solo los IDs.
        """
        sequences = []

        # Check if the combined file exists and read from it
        if os.path.exists(combined_file_path):
            with open(combined_file_path, "r") as combined_file:
                for line in combined_file:
                    if line.startswith(">"):  # Only consider header lines
                        sequences.append(line.split()[0][1:])  # Get the sequence ID without '>'
        
        # If no sequences found in combined file, check the individual sequences
        if not sequences:
            for genbank_id in self.genbank_ids:
                sequences.append(genbank_id)  # Add the GenBank IDs

        return sequences


# Ejemplo de uso
if __name__ == "__main__":
    familia = "GH180"  # Cambiar según la familia deseada
    email = "daleor0917@gmail.com"  # Cambiar por tu email real
    verbose = (
        os.getenv("TQ_VERBOSE", "F").lower().startswith("t")
    )  # Obtener el estado de verbose
    downloader = CazyDownloader(email, familia, verbose=verbose)

    downloader.obtener_genbank_ids()
    downloader.obtener_fasta_ncbi()

    # Agregar una secuencia externa
    external_fasta_file = "/workspace/cazy_analysis/input/GH51_short.txt"  # Cambiar por el nombre del archivo externo
    downloader.agregar_secuencia_externa(external_fasta_file)
