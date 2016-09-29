import numpy as np
import csv
POSITION = 0
MOTIF = 1
KMER  = 2
ZSCORE  = 3
PVALUE  = 4

server = False
workspace = "/home/jonathan/Documents/"

if server:
    workspace = "/home/jonathan/"

rbpMap_summary_file_name = workspace + "data/rbp_motifs/All_Predictions.txt"
parsed_summary_file_name = workspace + "data/rbp_motifs/All_Predictions_parsed"

all_RBPs = {}

def is_number(string):
    try:
        float(string)
    except ValueError:
        return False
    else:
        return True

def get_RBP_motifs(gene_name):
    motif = "a - motif"
    occurances = "b - occurances"
    positions = "c - positions"
    k_mers = "d - k_mers"
    p_values = "e - p_values"
    average_pvalue = "f - average_pvalue"
    all_fields = ["motif", "occurances", "positions", "k_mers", "p_values", "average_pvalue"]


    # get_motifs_from_rbpMap(f_name)
    prev_protein = ""
    protein = ""
    bNew_motif = False
    dMotif_counts = {}

    with open(rbpMap_summary_file_name) as rbp_output_file:
        for line in rbp_output_file:
            line = line.split()
            if len(line) < 1: continue
            elif not is_number(line[0]):
                if line[0] == "Protein:":

                    # Not the first protein in the file
                    if protein:
                        # Add current gene to the list of genes that have a binding site for this RBP in their UTR
                        try:
                            all_RBPs[protein].append(gene_name)
                        except KeyError: # first gene for this RBP, create the list.
                             all_RBPs[protein] = [gene_name]

                        # calculate stats for previous protein
                        dMotif_counts[protein][average_pvalue] = np.average(np.array([float(x) for x in dMotif_counts[protein][p_values]]))

                    # For all including the first protein in the file
                    protein = line[1]
                    dMotif_counts[protein] = {  motif: "",
                                                occurances: 0,
                                                positions: [],
                                                k_mers: [],
                                                p_values: [] }

                continue # Not a data line, nothing more to do here

            # A data line, populate the dict with the data
            dMotif_counts[protein][motif] = line[MOTIF]
            dMotif_counts[protein][occurances] += 1
            dMotif_counts[protein][positions].append(line[POSITION])
            dMotif_counts[protein][k_mers].append(line[KMER])
            dMotif_counts[protein][p_values].append(line[PVALUE])

        # This is a list of lists (will hold all the lines in the file):
        #   - the outer list is created from dMotif_counts and sorted by the protein name.
        #   - the inner list is sorted by field name, and sets the order of the fields in the output file.
        lines = [[prot] + [v for k, v in sorted(dFields.items(), key= lambda tup: tup[0])] for prot, dFields in sorted(dMotif_counts.items(), key= lambda tup: tup[0])]
        header = ["protein_name", "protein_motif", "num_occurences", "AUG_postions", "k-mers", "p-values", "average_pvalue" ]
        with open(parsed_summary_file_name, 'w') as parsed_file:
            writer = csv.writer(parsed_file,  delimiter="\t")
            writer.writerow(header)
            for line in lines:
                writer.writerow(line)
            print("written parsed file to {}".format(parsed_summary_file_name))

        # for l in lines:
        #     print l
        # print all_RBPs

get_RBP_motifs("my_gene")