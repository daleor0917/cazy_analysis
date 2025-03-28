#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import pandas as pd

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")

class BLASTrunner: 
                                
    def __init__(self, base_dir="/workspace/cazy_analysis/output", verbose=TQ_VERBOSE):
        self.verbose = verbose
        self.base_dir = base_dir
        self.output_dir_db = os.path.join(self.base_dir, "data_base")
        self.output_dir = os.path.join(self.base_dir, "blast_output")

        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.output_dir_db, exist_ok=True)
    
    def read_input(self, input_file):

        if self.verbose:
            print(f"Reading input file: {input_file}")

        try:
            with open(input_file, "r") as file:  
                lines = file.readlines()  

            sequences = {}
            sequence_name = None 

            for line in lines:
                line = line.strip()  
                if line.startswith(">"):  
                    sequence_name = line[1:] 
                    sequences[sequence_name] = ""  
                else:  
                    sequences[sequence_name] += line  

            return sequences  

        except FileNotFoundError:
            print(f"Error: The file {input_file} was not found.")  
            return None
        except Exception as e:
            print(f"Error: {e}")  
            return None


    def make_blast_db(self, input_file, dbtype='prot', output_db='mi_base_de_datos'):

        input_sequences = self.read_input(input_file)
        
        if input_sequences is None:
            print("Error: No se pudo leer el archivo de entrada.")
            return
        
        # Use the input_file directly instead of trying to get a key from the dictionary
        input_sequences_path = input_file  
        
        command = [
            "makeblastdb",
            "-in", input_sequences_path,
            "-dbtype", dbtype,
            "-out", os.path.join(self.output_dir_db, output_db)
        ]
        
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            if self.verbose:
                print("BLAST database updated successfully.")
        else:
            print("Error updating BLAST database:", result.stderr)


    def run_blastp(self, input_file, database_prefix, output_file):
        if self.verbose:
            print(f"Running BLASTP with input file: {input_file} and database: {database_prefix}")

        command = [
            "blastp",
            "-query", input_file,
            "-db", database_prefix,
            "-out", output_file,
            "-outfmt", "6"  # Formato tabular
        ]
        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            if self.verbose:
                print(f"BLASTP completed successfully. Results saved to {output_file}.")
        else:
            print("Error running BLASTP.", result.stderr)
            return None

        column_names = [
            "query_id", "subject_id", "identity", "alignment_length", "mismatches", "gap_opens",
            "query_start", "query_end", "subject_start", "subject_end", "e_value", "bit_score"
        ]
        # Leer la salida desde el archivo generado
        try:
            blast_df = pd.read_csv(output_file, sep="\t", header=None, names=column_names)
        except Exception as e:
            print(f"Error reading BLASTP output file: {e}")
            return None

        return blast_df

if __name__ == "__main__":
    BLAST = BLASTrunner(verbose=TQ_VERBOSE)
    
    input_file = "/workspace/cazy_analysis/input/GH51_short.txt"
    database_prefix = "/workspace/cazy_analysis/output/data_base/mi_base_de_datos"
    output_file = "/workspace/cazy_analysis/output/blast_output/blast_results.txt"
    
    # Crear la base de datos BLAST
    BLAST.make_blast_db(input_file)

    # Ejecutar BLASTP
    BLAST.run_blastp(input_file, database_prefix, output_file)
