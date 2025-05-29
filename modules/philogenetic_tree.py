#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from ete3 import Tree, TreeStyle

class PhylogeneticTreeBuilder:
    def __init__(self, input_fasta, output_dir="output/phylogenetic_trees"):
        self.input_fasta = input_fasta
        self.output_dir = output_dir
        self.filtered_fasta = os.path.join(self.output_dir, "filtered_sequences.fasta")
        self.aligned_file = os.path.join(self.output_dir, "aligned_sequences.phy")
        self.cleaned_aligned_file = os.path.join(self.output_dir, "aligned_sequences_clean.phy")
        self.tree_file = os.path.join(self.output_dir, "phylogenetic_tree.nwk")
        os.makedirs(self.output_dir, exist_ok=True)

    def filter_valid_sequences(self):
        """Filtra secuencias vacías o inválidas antes del alineamiento."""
        print("Filtrando secuencias inválidas...")
        count = 0
        with open(self.filtered_fasta, "w") as out_f:
            for record in SeqIO.parse(self.input_fasta, "fasta"):
                seq = str(record.seq).upper()
                if any(base in seq for base in "ACGT"):
                    SeqIO.write(record, out_f, "fasta")
                    count += 1
                else:
                    print(f"Secuencia descartada por ser inválida: {record.id}")
        if count < 3:
            raise ValueError("Se requieren al menos 3 secuencias válidas para construir el árbol.")
        print(f"{count} secuencias válidas conservadas.")

    def align_sequences(self):
        """Alinea las secuencias usando Clustal Omega."""
        print("Ejecutando alineamiento con Clustal Omega...")
        command = [
            "clustalo",
            "-i", self.filtered_fasta,
            "-o", self.aligned_file,
            "--force",
            "--outfmt", "phylip"  # Cambiado a 'phylip'
        ]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Error en Clustal Omega: {result.stderr}")
        print("Alineamiento completado.")
        return self.aligned_file

    def clean_aligned_sequences(self):
        """Elimina secuencias completamente vacías del alineamiento."""
        print("Filtrando secuencias vacías post-alineamiento...")
        alignment = AlignIO.read(self.aligned_file, "phylip")
        filtered = [rec for rec in alignment if set(rec.seq.upper()) - set("-")]
        if len(filtered) < 3:
            raise ValueError("Muy pocas secuencias útiles tras el alineamiento.")
        AlignIO.write(MultipleSeqAlignment(filtered), self.cleaned_aligned_file, "phylip")
        print(f"{len(filtered)} secuencias útiles conservadas tras el alineamiento.")
        return self.cleaned_aligned_file

    def build_tree(self, aligned_file):
        """Construye el árbol filogenético usando FastTree."""
        print("Ejecutando FastTree...")
        command = ["fasttree", "-nt", aligned_file]
        with open(self.tree_file, "w") as output:
            result = subprocess.run(command, stdout=output, stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Error en FastTree: {result.stderr}")
        if not os.path.isfile(self.tree_file):
            raise FileNotFoundError("El archivo del árbol no fue generado.")
        print("Árbol filogenético generado exitosamente.")
        return self.tree_file

    def visualize_tree(self, tree_file):
        """Exporta el árbol filogenético a un archivo de imagen usando ETE3."""
        if not tree_file or not os.path.isfile(tree_file):
            print("El archivo del árbol no existe o es None.")
            return
        print("Exportando árbol filogenético a imagen...")
        
        # Forzar modo sin interfaz gráfica
        os.environ["QT_QPA_PLATFORM"] = "offscreen"

        tree = Tree(tree_file)
        ts = TreeStyle()
        ts.show_leaf_name = True
        output_img = os.path.join(self.output_dir, "phylogenetic_tree.png")
        try:
            tree.render(output_img, tree_style=ts)
            print(f"Árbol filogenético exportado a: {output_img}")
        except Exception as e:
            print(f"Error al exportar el árbol: {e}")


if __name__ == "__main__":
    input_fasta = "output/cazy_sequences/updated_sequences.fasta"
    builder = PhylogeneticTreeBuilder(input_fasta)
    try:
        builder.filter_valid_sequences()
        builder.align_sequences()
        builder.clean_aligned_sequences()
        tree_path = builder.build_tree(builder.cleaned_aligned_file)
        builder.visualize_tree(tree_path)
    except Exception as e:
        print(f"Error durante la ejecución: {e}")