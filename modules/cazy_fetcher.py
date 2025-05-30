#!/usr/bin/env python
# coding: utf-8

import os
import re
import requests
import argparse
from bs4 import BeautifulSoup
from Bio import Entrez

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")

class CazyFetcher:
    def __init__(self, email, familias, output_dir="output/cazy_sequences", verbose=TQ_VERBOSE):
        Entrez.email = email
        self.familias = familias  # Ahora es una lista de familias
        self.genbank_ids = []
        self.verbose = verbose
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def _print(self, message):
        if self.verbose:
            print(message)

    def obtener_genbank_ids(self):
        for familia in self.familias:
            url = f"https://www.cazy.org/{familia}_characterized.html"
            response = requests.get(url)
            soup = BeautifulSoup(response.text, "html.parser")

            header_skipped = False

            for row in soup.find_all("tr"):
                cols = row.find_all("td")
                if len(cols) > 4:
                    genbank_text = cols[4].text.strip()
                    if genbank_text:
                        if not header_skipped and genbank_text == "GenBank":
                            header_skipped = True
                            continue

                        valid_ids = re.findall(r'[A-Z_]+\d+\.\d+', genbank_text)
                        if valid_ids:
                            first_valid_id = valid_ids[0]
                            if first_valid_id not in self.genbank_ids:
                                self.genbank_ids.append(first_valid_id)

        self._print(f"GenBank IDs obtenidos ({len(self.genbank_ids)}): {self.genbank_ids}")
        return self.genbank_ids

    def obtener_fasta_ncbi(self, output_file="sequences.fasta"):
        with open(output_file, "w") as fasta_file:
            for genbank_id in self.genbank_ids:
                try:
                    handle = Entrez.efetch(db="protein", id=genbank_id, rettype="fasta", retmode="text")
                    fasta_sequence = handle.read()
                    handle.close()

                    if fasta_sequence.startswith(">"):
                        fasta_file.write(fasta_sequence + "\n")
                    else:
                        self._print(f"No se encontró la secuencia para {genbank_id}")

                except Exception as e:
                    self._print(f"Error con {genbank_id}: {e}")

def main():
    parser = argparse.ArgumentParser(description='Fetch CAZy family sequences')
    parser.add_argument('--family', required=True, help='CAZy family name (e.g., GH62)')
    parser.add_argument('--output', required=True, help='Output file name')
    parser.add_argument('--email', default='daleor0917@gmail.com', help='Email for NCBI')
    
    args = parser.parse_args()
    
    # Crear el fetcher con una sola familia
    familias = [args.family]
    email = args.email
    verbose = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")
    
    fetcher = CazyFetcher(email, familias, verbose=verbose)
    fetcher.obtener_genbank_ids()
    fetcher.obtener_fasta_ncbi(args.output)

if __name__ == "__main__":
    main()