
from collections import Counter
from itertools import combinations
import re
import csv

import decimal
decimal.getcontext().prec = 2

import numpy as np
import pandas as pd
from collections import namedtuple
# example to show how to write named matrices to file with pandas
# ----
# import numpy as np
# import pandas as pd
#
# A = np.random.randint(0, 10, size=36).reshape(6, 6)
# names = [_ for _ in 'abcdef']
# df = pd.DataFrame(A, index=names, columns=names)
# df.to_csv('df.csv', index=True, header=True, sep=' ')

from global_defs import *

 # Utility functions
def is_number(string):
    try:
        float(string)
    except ValueError:
        return False
    else:
        return True

def key_tuple_by_first(tup):
        return tup[0]

def sort_by_first_in_tuple(list_of_tuples):
    return sorted(list_of_tuples, key=key_tuple_by_first)

def print_dict(dict):
    for k, v in dict.iteritems():
        print k, ":", v

def update_average(curr_avg, num_items, added_val):
    return float(curr_avg * num_items + added_val) / (num_items + 1)

# Feature extraction functions
def count_nucleotides(name, sequence):
    header  = ["name"] + [char + "_count" for char in NUCLEOTIDES]
    counter = Counter(sequence)
    counts  = [name] + [counter[char] + counter[char.lower] for char in NUCLEOTIDES]
    return header, counts

def count_dinucleotides(sequence):
    pairs  = [p[0]+p[1] for p in combinations(NUCLEOTIDES, 2)]
    header = [p + "_count" for p in pairs]
    counts = [sequence.count(p) + sequence.count(p.lower) for p in pairs]
    return header, counts

def count_nucleotides_and_dinucleotides(name, sequence):
    header1, counts1 = count_nucleotides(name, sequence)
    header2, counts2 = count_dinucleotides(sequence)
    header = header1 + header2
    counts = counts1 + counts2
    # df = pd.DataFrame(counts, columns=header).set_index("name")
    return {"header": header,"counts": counts}

def update_global_genes_per_RBPs_dict(protein_name, gene_name):
    try:
        all_RBPs[protein_name].append(gene_name)
    except KeyError: # first gene for this RBP, create the list.
        all_RBPs[protein_name] = [gene_name]


# This function aggregates the RBP motifs found in a single gene by RbpMap
def get_RBP_motifs_for_gene(motif_stats_per_gene_dict, gene_name, input_lines):

    motif_stats_per_gene_dict[gene_name] = {} #init a dict of dicts to hold RBPs and their data
    gene_dict = motif_stats_per_gene_dict[gene_name]
    protein = ""

    for line in input_lines:
        if len(line) < 1:
            continue # skip empty lines
        elif not is_number(line[0]):
            if line[0] == "Protein:":
                # Not the first protein in the file
                update_global_genes_per_RBPs_dict(protein, gene_name)

                # For all including the first protein in the file
                protein = line[1]
                gene_dict[protein] = \
                    {
                        motif: "",
                        average_pvalue: 0,
                        occurances: 0,
                        positions: [],
                        p_values: [],
                        k_mers: []
                    }
                protein_dict = gene_dict[protein]

            continue # Not a data line, nothing more to do here


        # A data line, populate the dict with the data
        avg = update_average(protein_dict[average_pvalue], protein_dict[occurances], float(line[PVALUE]))
        avg = round(decimal.Decimal(avg),4)
        protein_dict[average_pvalue] = avg

        protein_dict[motif] = line[MOTIF]
        protein_dict[occurances] += 1
        protein_dict[positions].append(line[POSITION])
        protein_dict[k_mers].append(line[KMER])
        protein_dict[p_values].append(line[PVALUE])
    return



# This function parses the output from RbpMap,
# collects all protein motifs that appear in every gene (one gene at a time)
# and aggregates them in a dict-of-dicts (motif_stats_per_gene_dict)
def get_RBP_motifs_all_genes(rbp_output_file = rbp_data_folder + "All_Predictions.txt"):
    collect_lines = False
    lines = []
    extrated_features = []
    gene_name = ""
    header = ""
    motif_stats_per_gene_dict = {}

    with open(rbp_output_file) as file:
        for line in file:
            line = line.split()
            if len(line) < 1 :
                continue

            # finish a gene section
            if re.match(r'\*+', line[0]) and collect_lines:
                collect_lines = False
                get_RBP_motifs_for_gene(motif_stats_per_gene_dict, gene_name, lines)
                continue

            # add line to be processed by get_RBP_motifs
            if collect_lines:
                lines.append(line)
                continue

            # this line is a gene name
            if re.match(r'ENSMUST.*', line[0]):
                gene_name = re.match(r'ENSMUST.*', line[0]).group(0)
                continue

            # start a new gene section
            if re.match(r'=+', line[0]):
                collect_lines = True
                continue

    genes = []
    frames = []

    # motif_stats_per_gene_dict is a dict of dicts.
    # Reshape it into a list of genes and a list of frames holding the info about the RBPs
    for gene, d in motif_stats_per_gene_dict.iteritems():
        genes.append(gene)
        frames.append(pd.DataFrame.from_dict(d, orient='index'))

    df = pd.concat(frames, keys=genes)
    df_occurances_only = df[occurances].unstack()
    df_occurances_only.fillna(0, inplace=True)
    # print df_occurances_only.head(100)
    df_occurances_only.to_csv(rbp_occurances_per_gene, sep ='\t')

    cols = df_occurances_only.columns
    df_bool = df_occurances_only > 0
    df_motifs_per_gene = df_bool.apply(lambda x: list(cols[x.values]), axis=1)

    # mojo to create an actual dataFrame from a series where every cell is a list of motifs...
    df_motifs_per_gene = df_motifs_per_gene.to_frame()
    df_motifs_per_gene = df_motifs_per_gene[0].apply(pd.Series)
    # print df_motifs_per_gene.head(100)
    df_motifs_per_gene.to_csv(rbp_motifs_per_gene_file, sep = '\t')

def get_kmers_all_genes(all_genes_fasta = utr_output_fasta):
    K = 6
    gene_name = "ERROR! Sequence appears before first gene name!"
    dKmers_per_gene = {}
    with open(all_genes_fasta) as file:
        for line in file:
            line = line.split()
            if len(line) < 1 :
                continue #skip empty lines

            # this line is a gene name
            if re.match(r'>ENSMUST.*', line[0]):
                gene_name = re.match(r'>(ENSMUST.*)', line[0]).group(1)
                continue

            else: # this is a sequence line
                sequence = line[0]
                kmers = [sequence[start:start+K] for start in range(len(sequence) - K)]
                dKmers_per_gene[gene_name] = pd.Series(len(kmers)*[1], index = kmers)
    df_kmers = pd.DataFrame.from_dict(dKmers_per_gene, 'index', "int")
    df_kmers.fillna(0, inplace=True)
    df_kmers.to_csv(kmer_occurances_per_gene, sep = '\t')

    print df_kmers.head(10)

if __name__ == '__main__':
    get_RBP_motifs_all_genes()
    # get_kmers_all_genes()