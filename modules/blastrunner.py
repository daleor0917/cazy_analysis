#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import pandas as pd
import hashlib

TQ_VERBOSE = os.getenv("TQ_VERBOSE", "F").lower().startswith("t")
TQ_CACHE = os.getenv("TQ_CACHE", "F").lower().startswith("t")

class BLASTrunner:
    def __init__(
        self,
        input_fasta,
        output_dir="/workspace/cazy_analysis/modules/output/blast",
        verbose=TQ_VERBOSE,
        cache=TQ_CACHE,
    ):
        self.input_fasta = input_fasta
        self.verbose = verbose
        self.cache = cache
        self.output_dir = output_dir
        self.output_dir_db = os.path.join(self.output_dir, "data_base")
        self.output_dir_results = os.path.join(self.output_dir, "blast_results")
        self.cache_dir = os.path.join(self.output_dir, "cache")
        self.md5_cache = {}  # Dictionary to store MD5 hashes

        # Ensure the directories are created correctly
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.output_dir_db, exist_ok=True)
        os.makedirs(self.output_dir_results, exist_ok=True)
        os.makedirs(self.cache_dir, exist_ok=True)

    def read_input(self, input_file):
        if self.verbose:
            print(f"Reading input file: {input_file}")

        try:
            with open(input_file, "r") as file:
                lines = file.readlines()

            sequences = {}
            sequence_name = None
            concatenated_sequences = ""

            for line in lines:
                line = line.strip()
                if line.startswith(">"):
                    sequence_name = line[1:]
                    sequences[sequence_name] = ""
                else:
                    sequences[sequence_name] += line
                    concatenated_sequences += line  # Concatenate for MD5 hash

            # Calculate MD5 hash of the concatenated sequences
            md5_hash = hashlib.md5(concatenated_sequences.encode()).hexdigest()

            return sequences, md5_hash

        except FileNotFoundError:
            print(f"Error: The file {input_file} was not found.")
            return None, None
        except Exception as e:
            print(f"Error: {e}")
            return None, None

    def make_blast_db(self, dbtype="prot", output_db="mi_base_de_datos"):
        input_sequences = self.read_input(self.input_fasta)

        if input_sequences is None:
            print("Error: No se pudo leer el archivo de entrada.")
            return

        command = [
            "makeblastdb",
            "-in",
            self.input_fasta,
            "-dbtype",
            dbtype,
            "-out",
            os.path.join(self.output_dir_db, output_db),
        ]

        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            if self.verbose:
                print("BLAST database updated successfully.")
        else:
            print("Error updating BLAST database:", result.stderr)

    def run_blastp(self, database_prefix, output_file):
        # Read input and get sequences and MD5 hash
        input_sequences, md5_hash = self.read_input(self.input_fasta)
        print(f"Current MD5: {md5_hash}")

        if input_sequences is None:
            print("Error: No se pudo leer el archivo de entrada.")
            return

        # Define the cache file path
        cache_file = os.path.join(self.cache_dir, os.path.basename(self.input_fasta) + ".md5")

        # Check if cache is enabled and if the MD5 hash is the same
        if self.cache:
            # Load existing MD5 hashes from the cache file if it exists
            if os.path.exists(cache_file):
                with open(cache_file, "r") as f:
                    cached_md5 = f.read().strip()
                    self.md5_cache[self.input_fasta] = cached_md5  # Store in the cache dictionary
                if cached_md5 == md5_hash:
                    print("No need to execute BLASTP, file is the same.")
                    return

        if self.verbose:
            print(f"Running BLASTP with input file: {self.input_fasta} and database: {database_prefix}")

        command = [
            "blastp",
            "-query",
            self.input_fasta,
            "-db",
            database_prefix,
            "-out",
            output_file,
            "-outfmt",
            "6",
            "-evalue",
            "1",
        ]

        result = subprocess.run(command, capture_output=True, text=True)

        if result.returncode == 0:
            if self.verbose:
                print(f"BLASTP completed successfully. Results saved to {output_file}.")
            # Save the current MD5 hash to the cache file
            with open(cache_file, "w") as f:
                f.write(md5_hash)
            self.md5_cache[self.input_fasta] = md5_hash  # Update the cache dictionary
        else:
            print("Error running BLASTP.", result.stderr)
            return None

        column_names = [
            "query_id",
            "subject_id",
            "identity",
            "alignment_length",
            "mismatches",
            "gap_opens",
            "query_start",
            "query_end",
            "subject_start",
            "subject_end",
            "e_value",
            "bit_score",
        ]

        # Read the output from the generated file
        try:
            blast_df = pd.read_csv(
                output_file, sep="\t", header=None, names=column_names
            )

            # Filter to keep only the most significant result (lowest e_value) for each query_id and subject_id combination
            best_hits = blast_df.loc[
                blast_df.groupby(["query_id", "subject_id"])["e_value"].idxmin()
            ]

            # Save the filtered DataFrame to the output file
            best_hits.to_csv(
                output_file, sep="\t", index=False
            )  # Save the DataFrame in the same file
        except Exception as e:
            print(f"Error reading BLASTP output file: {e}")
            return None

        return best_hits


if __name__ == "__main__":
    input_fasta = "/workspace/cazy_analysis/input/mis_GH51.txt"
    output_dir = "/workspace/cazy_analysis/modules/output/blast"
    BLAST = BLASTrunner(input_fasta, output_dir, verbose=TQ_VERBOSE, cache=TQ_CACHE)

    database_prefix = os.path.join(output_dir, "data_base/mi_base_de_datos")
    output_file = os.path.join(output_dir, "blast_results/blast_results.txt")

    # Crear la base de datos BLAST
    BLAST.make_blast_db()

    # Ejecutar BLASTP
    BLAST.run_blastp(database_prefix, output_file)