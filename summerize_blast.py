import os
import sys
import pandas as pd
import numpy as np
import argparse
"""
Scipt to summerize the results from the Blast
"""
parser = argparse.ArgumentParser(
    description=
    'Script to summerize the results from the Blast. Check the README file \
    for more information regarding the used arguments')
parser.add_argument("-bf",
                    "--blast-file",
                    help="directory path to the blast file",
                    type=str,
                    required=True)
parser.add_argument("-ff",
                    "--fasta-file",
                    help="directory path to the fasta file",
                    type=str,
                    required=True)
parser.add_argument("-i",
                    "--identity",
                    help="hit identity in percentage to be accepted as a hit, \
                        recommended to use the value of 97",
                    type=float,
                    default=97)
parser.add_argument("-al",
                    "--alignment-length",
                    help="hit alignment length in percentage to be accepted as\
                    a hit, recommeded to use default of 95",
                    default=95,
                    type=int)

no_hits_list = []
rows_dicts = {}
cols_dicts = {}
vtx_name_dict = {}
total = 0

# for fasta sequences
sequences = {}
sample_blast = {}
sample_fasta = {}
fasta_sequence = {}
sequences_total = 0
sequences_match = 0
total_hits = 0
sample_vtx_df = pd.DataFrame()
last_sample_name = None

next_in_sequence = False
sequence_data = ""


def process_blast_line(line):
    global total_hits, sample_vtx_df, total
    line = line.strip()
    col = line.split("\t")
    if len(col) < 17:
        print("Error: invalid blast line")
        return
    # Get the minimum length of the query and subject
    min_len = min(int(col[13]), int(col[14]))
    # Get the alignment length
    align_len = int(col[6])

    # Increase the total counter
    total += 1

    # Print the hit count on every 10000 hits
    if total_hits % 10000 == 0:
        print(f"Hits: {total_hits}/{total}")

    if float(col[4]) >= args.identity and align_len >= min_len * (
            args.alignment_length / 100.0):
        # Get the sample name
        sample = col[0].split("-")[0]
        # Parse vtx name
        vtx_name = col[2].split(" ")[-1]

        # Update the columns for each sample
        if sample not in cols_dicts:
            cols_dicts[sample] = 1
        else:
            cols_dicts[sample] += 1

        # Update the row for each vtx name
        if vtx_name not in rows_dicts:
            rows_dicts[vtx_name] = 1
        else:
            rows_dicts[vtx_name] += 1

        index = "_".join([vtx_name, sample])
        if index in vtx_name_dict:
            vtx_name_dict[index] += 1
        else:
            vtx_name_dict[index] = 1
        total_hits += 1
        sequences[col[0]] = True

    # Print final total hits
    print(f"Total Hits: {total_hits}/{total}")


def parse_blast_file(filename):
    output_filename_vtx = "output/hits_vtx.csv"
    output_filename_vtx_pivot = "output/hits_pivot_table.csv"
    with open(filename, 'r') as file:
        lines = file.readlines()
        # Map the function to the lines
        list(map(process_blast_line, lines))

        # Create vtx name dataframe and save as a csv file
        vtx_name_df = pd.DataFrame.from_dict(vtx_name_dict, orient='index')
        vtx_name_df.to_csv(output_filename_vtx)
        print(f"Dataframe saved as {output_filename_vtx}")

        # Create pivot table
        # Split the raw vtx_name_df into row, column, and value
        vtx_pivot_df = pd.DataFrame()
        vtx_pivot_df[[
            'row', 'column'
        ]] = vtx_name_df.index.to_series().str.extract(r'([^_]+)_([^,]+)')
        vtx_pivot_df['value'] = vtx_name_df[0].astype(int)

        # Convert value to integer
        vtx_pivot_df['value'] = vtx_pivot_df['value'].astype(int)

        # Create a vtx_pivot_df table
        vtx_pivot_table = vtx_pivot_df.pivot_table(index='row',
                                                   columns='column',
                                                   values='value')
        vtx_pivot_table.fillna(0, inplace=True)
        vtx_pivot_table.to_csv(output_filename_vtx_pivot)
        print(f"Pivot table saved as {output_filename_vtx_pivot}")


def process_fasta_line(line):
    global last_sample_name, sequences_total, next_is_sequence, sequence_data
    line = line.strip()
    if line.startswith(">"):
        next_is_sequence = False
        if last_sample_name is not None:
            fasta_sequence[last_sample_name] = sequence_data
        sequence_data = ""
        last_sample_name = None
        sequences_total += 1
        # Print the sequence count on every 10000 sequences
        if sequences_total % 10000 == 0:
            print(f"Sequences: {sequences_match/sequences_total}")

        # Remove > from the line and split the line
        fasta_data = line[1:]

        if fasta_data not in sequences and fasta_data not in sample_fasta:
            sample_fasta[fasta_data] = 1
            last_sample_name = fasta_data
            next_is_sequence = True
		
        elif fasta_data not in sequences:
            sample_fasta[fasta_data] += 1
            last_sample_name = fasta_data
            next_is_sequence = True

    elif next_is_sequence:
        sequence_data += line


def parse_fasta_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        for line in lines:
            process_fasta_line(line)

    # Store the dictionary in the fasta_sequence as a csv file
    fasta_sequence_df = pd.DataFrame.from_dict(fasta_sequence, orient='index')
    fasta_sequence_df.to_csv("output/fasta_sequence.csv")

    # Store the sequence table count and pivot table
    sample_fasta_df = pd.DataFrame.from_dict(sample_fasta, orient='index')
    sample_fasta_df.to_csv("output/sample_fasta.csv")
    print(f"Dataframe saved as sample_fasta.csv")


if __name__ == "__main__":
    print("Parsing Blast File")
    args = parser.parse_args()
    parse_blast_file(args.blast_file)
    parse_fasta_file(args.fasta_file)
