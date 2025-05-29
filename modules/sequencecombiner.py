#!/usr/bin/env python
# coding: utf-8

import os

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")

class SequenceCombiner:
    def __init__(self, cazy_file, external_fasta_file=None, output_dir="output/cazy_sequences", output_file="updated_sequences.fasta", verbose=TQ_VERBOSE):
        self.cazy_file = cazy_file
        self.external_fasta_file = external_fasta_file
        self.output_dir = output_dir
        self.output_file = os.path.join(output_dir, output_file)  # Asegúrate de que el archivo de salida esté en el directorio correcto
        self.verbose = verbose

        # Crear el directorio de salida si no existe
        os.makedirs(self.output_dir, exist_ok=True)

    def _print(self, message):
        if self.verbose:
            print(message)

    def combinar_secuencias(self):
        try:
            with open(self.output_file, "w") as combined_file:
                with open(self.cazy_file, "r") as fasta_file:
                    combined_file.write(fasta_file.read())

                if self.external_fasta_file:
                    if os.path.exists(self.external_fasta_file):
                        with open(self.external_fasta_file, "r") as external_file:
                            combined_file.write(external_file.read())
                        self._print(f"Secuencia externa añadida desde {self.external_fasta_file}.")
                    else:
                        self._print(f"No se encontró el archivo externo: {self.external_fasta_file}")
                else:
                    self._print("No se proporcionó un archivo externo.")

            self._print(f"Secuencias combinadas guardadas en {self.output_file}.")

        except Exception as e:
            self._print(f"Error al combinar secuencias: {e}")

# Ejemplo de uso
if __name__ == "__main__":
    cazy_file = "output/cazy_sequences/sequences.fasta"
    external_fasta_file = "/workspace/cazy_analysis/input/mis_GH51.txt"
    combiner = SequenceCombiner(cazy_file, external_fasta_file)

    combiner.combinar_secuencias()