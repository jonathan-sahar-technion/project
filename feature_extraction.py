# from cogent.parse.consan import sequence
from collections import Counter
from itertools import combinations
import re
import numpy as np
import csv

import decimal
decimal.getcontext().prec = 2

from blastdb_extract_sequences import fasta_output_file
import pandas as pd
# example to show how to write named matrices to file with pandas
# ----
# import numpy as np
# import pandas as pd
#
# A = np.random.randint(0, 10, size=36).reshape(6, 6)
# names = [_ for _ in 'abcdef']
# df = pd.DataFrame(A, index=names, columns=names)
# df.to_csv('df.csv', index=True, header=True, sep=' ')


# Constants
POSITION = 0
MOTIF = 1
KMER  = 2
ZSCORE  = 3
PVALUE  = 4
NUCLEOTIDES = "ACUG"

# Globals
all_RBPs = {}
motif = "motif"
occurances = "occurances"
positions = "positions"
k_mers = "k_mers"
average_pvalue = "average_pvalue"
p_values = "p_values"



#Paths

server = False
workspace = "/home/jonathan/Documents/data/"

if server:
    workspace = "/srv01/technion/jonathans/data/"

rbp_data_folder = workspace + "rbp_motif_analysis/"
rbp_motifs_per_gene_path = rbp_data_folder + "motifs_by_gene_name.tsv"
genes_per_protein_path = rbp_data_folder + "genes_by_protein_name.tsv"
genes_RBPs_matrix = rbp_data_folder + "genes_RBPs_matrix.tsv"
genes_kmers_matrix = rbp_data_folder + "genes_kmers_matrix.tsv"

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
def count_nucleotides(sequence):
    header = [char + "_count" for char in NUCLEOTIDES]
    counter = Counter(sequence)
    counts = [counter[char] + counter[char.lower] for char in NUCLEOTIDES]
    return header, counts

def count_denucleotides(sequence):
    pairs = [p[0]+p[1] for p in combinations(NUCLEOTIDES, 2)]
    header = [p + "_count" for p in pairs]
    counts = [sequence.count(p) + sequence.count(p.lower) for p in pairs]
    return header, counts


def update_global_genes_per_RBPs_dict(protein_name, gene_name):
    try:
        all_RBPs[protein_name].append(gene_name)
    except KeyError: # first gene for this RBP, create the list.
        all_RBPs[protein_name] = [gene_name]

def get_RBP_motifs_for_gene(motif_stats_per_gene_dict, gene_name, input_lines):

    # all_fields = ["motif", "occurances", "positions", "k_mers", "p_values", "average_pvalue"]
    # all_fields = ["motif", "p_values", "average_pvalue"]


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


    # print_dict(gene_dict)



    #
    # # Stringify the lists of values
    # for protein in motif_stats_per_gene_dict.keys():
    #     # motif_stats_per_gene_dict[protein][positions] = ",".join(motif_stats_per_gene_dict[protein][positions])
    #     # motif_stats_per_gene_dict[protein][k_mers] = ",".join(motif_stats_per_gene_dict[protein][k_mers])
    #     motif_stats_per_gene_dict[protein][p_values] = ",".join(motif_stats_per_gene_dict[protein][p_values])
    #
    # # This is a list of lists (will hold all the lines in the file):
    # #   - the outer list is created from motif_stats_per_gene_dict and sorted by the protein name.
    # #   - the inner list is sorted by field name, and sets the order of the fields in the output file.
    #
    #
    # sorted_proteins_and_counts = sort_by_first_in_tuple(motif_stats_per_gene_dict.items())
    #
    # lines = [
    #     [gene_name, prot] +
    #     [value for key, value in sort_by_first_in_tuple(dFields.items())[:-1]] for prot, dFields in sorted_proteins_and_counts
    #     ]
    #
    # # header = ["gene_name", "protein_name", "protein_motif", "num_occurences", "AUG_postions", "k-mers", "p-values", "average_pvalue" ]
    # header = ["gene_name", "protein_name", "protein_motif","num of occurences","average_pvalue" ]
    # return header, lines
    return "a", "b"

    return

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
                header, values = get_RBP_motifs_for_gene(motif_stats_per_gene_dict, gene_name, lines)
                # header, values = get_RBP_motifs_for_gene(motif_stats_per_gene_dict, gene_name, lines)
                extrated_features += values
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

    for gene, d in motif_stats_per_gene_dict.iteritems():
        genes.append(gene)
        frames.append(pd.DataFrame.from_dict(d, orient='index'))

    df = pd.concat(frames, keys=genes)
    df_occurances_only = df[occurances].unstack()
    df_occurances_only.fillna(0, inplace=True)
    print df_occurances_only.head(100)

    exit()



    with open(rbp_motifs_per_gene_path, 'w+') as file:
        writer = csv.writer(file,  delimiter="\t")
        writer.writerow(header)
        writer.writerows(extrated_features)

    header = ['protein_name', 'genes']
    with open(genes_per_protein_path, 'w+') as file:
        writer = csv.writer(file,  delimiter="\t")
        writer.writerow(header)
        all_proteins = [[protein] + genes for protein, genes in sort_by_first_in_tuple(all_RBPs.items())]
        writer.writerows(all_proteins)

        # print  all_proteins

    for k,v in all_RBPs.items():
        all_RBPs[k] = pd.Series([1] * len(v), index = v)
        # print k, all_RBPs[k]

    df_RBPs = pd.DataFrame.from_dict(all_RBPs, 'index', dtype=int)
    df_RBPs.fillna(0, inplace=True)
    df_RBPs.to_csv(genes_RBPs_matrix, sep = '\t')
    # print df_RBPs

def get_kmers_all_genes(all_genes_fasta = fasta_output_file):
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
    df_kmers.to_csv(genes_kmers_matrix, sep = '\t')

    # print df

if __name__ == '__main__':

    get_RBP_motifs_all_genes()
    # get_kmers_all_genes()