#!/Local/Anaconda-2.0.1-Linux-x86_64/bin/python

import csv
from corebio.seq import protein
from subprocess import check_output
from feature_extraction import count_nucleotides_and_dinucleotides, get_RBP_motifs_all_genes
from global_defs import *
import pandas as pd
from collections import namedtuple
# exon-starts are 0-based, exon-ends are 1-based
# my_input_file = sys.argv[1]
# my_output_file = sys.argv[2]

# line structure
NAME_FIELD = 1
CHROM_FIELD = 2
STRAND_FIELD = 3
TX_START_FIELD = 4
TX_END_FIELD = 5
CDS_START_FIELD = 6
CDS_END_FIELD = 7
EXON_COUNT_FIELD = 8
EXON_STARTS_FIELD = 9
EXON_ENDS_FIELD = 10
NAME2_FIELD = 12
RGB = 111 #some random value
SCORE = 500
START = 0
END = 1000000000

# run blastdbcmd? RBPmap?
bExtract_fasta = True
bCreate_RBPmap_entries = False

# globals to extract features from positional annotations
g_df_features = pd.DataFrame()
df_counts = pd.DataFrame()

def run_blastdbcmd(entries):
    # write entries to temporary file
    entries += "\n"
    with open(blastdb_entries_file, 'w') as file:
        # print "writing to", entries_file, "..."
        file.write("\n".join(entries))

    blastdb_str = "blastdbcmd -db {db} -entry_batch {entries}".format(db = blastdb_path, entries = blastdb_entries_file)
    grep_str = 'grep -v ">"'
    command_str = blastdb_str + "|" + grep_str

    # execute the command
    command_result = check_output(command_str, shell=True)
    return command_result


def run_RBPmap(entries):
    # write entries to temporary file
    with open(rbpmap_entries_file, 'w') as file:
        # print "writing to", entries_file, "..."
        file.write("\n".join(entries))

    rbpmap_cmd = "rbpMap -input {input_file} -genome 'mouse' -db 'mm9'".format(input_file = rbpmap_entries_file)

    # execute the command
    command_result = check_output(rbpmap_cmd, shell=True)
    return command_result


def proccess_line(line, mode): # 'mode' is one of "blastcmd", "rbpmap"

    # init function wide parameters
    blastcmd_entries = []
    exons_to_keep = []
    rbpmap_entries = []

    # extract info from line
    strand = line[STRAND_FIELD]
    name = line[NAME_FIELD]
    name2 = line[NAME2_FIELD]
    joint_name = name2 + "_" + name
    chrom = line[CHROM_FIELD]
    tx_start = line[TX_START_FIELD]
    tx_end = line[TX_END_FIELD]
    cds_end= int(line[CDS_END_FIELD])
    cds_start = int(line[CDS_START_FIELD])

    # A note on 0/1 based representation: the basic data file we're processing here is called bigGenePred, which is almost
    # the internal representation of the UCSC genome browser. It uses start coordinates that are 0 based and end coordinates that are
    # 1 based. We need coordinates for blastdbcmd, which requires consistent 1 based representation.

    # break the list fields into lists
    # +1 because the bed file start indices are zero-based
    exon_ends =  [int(x) for x in line[EXON_ENDS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list
    exon_starts =  [int(x) +1 for x in line[EXON_STARTS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list

    entries = ""

    if strand == '+':

        # 3' ----------------------------------------  5'
        # 5' -----ex_start-----ex_end-----cds_start--> 3'

        # 3' ----------------------------------------  5'
        # 5' -----ex_start----cds_start-----ex_end---> 3'

        # list of tuples of coordinates, where the the start of the exon is before the start of the CDS
        # exons_to_keep = [(int(exon_start),int(exon_ends[i]))
        #                  for i, exon_start in enumerate(exon_starts)
        #                  if int(exon_start) < cds_start]
        exons_to_keep = []
        for i, exon_start in enumerate(exon_starts):
            exon_end = exon_ends[i]
            if int(exon_end) < cds_start:
                exons_to_keep.append((int(exon_start),int(exon_end)))
            elif exon_start < cds_start and int(exon_end) > cds_start:
                exons_to_keep.append((int(exon_start),int(cds_end)))
            else: continue


        blastcmd_entries = ["{chr} {start}-{end} plus".format(chr=chrom, start=start, end=end) for start, end in exons_to_keep]
        rbpmap_entries = ["{chr}:{start}-{end}:plus".format(chr=chrom, start=start, end=end) for start, end in exons_to_keep]

    elif strand == '-':

        # 3' <---cds_start---ex_start-----ex_end----- 5'
        # 5' ---------------------------------------- 3'

        # 3' <-----ex_start----cds_start-----ex_end-- 5'
        # 5' ---------------------------------------- 3'


        # list of tuples of coordinates, where the the start of the exon is before the start of the CDS
        # exons_to_keep = [(int(exon_starts[i]),int(exon_end))
        #                  for i, exon_end in enumerate(exon_ends)
        #                  if int(exon_end) > cds_end]

        exons_to_keep = []
        for i, exon_end in enumerate(exon_ends):
            exon_start = exon_starts[i]

            if int(exon_start) > cds_start:
                exons_to_keep.append((int(exon_start),int(exon_end)))
            elif exon_start < cds_start and int(exon_end) > cds_start:
                exons_to_keep.append((int(cds_start),int(exon_end)))
            else: continue

        blastcmd_entries = ["{chr} {start}-{end} minus".format(chr=chrom, start=start, end=end) for start, end in exons_to_keep]
        blastcmd_entries.reverse() # since this will give the reverse compliment, in order to get a correct concatenation we want to reverse the order of the exons.

        rbpmap_entries = ["{chr}:{start}-{end}:minus".format(chr=chrom, start=start, end=end) for start, end in exons_to_keep]

    if mode == 'blastcmd' and len(blastcmd_entries) > 1:
        print "mode detected: blastcmd"
        first_kept_start = exons_to_keep[0][0]
        first_kept_end = exons_to_keep[0][1]
        first_kept_size = first_kept_end - first_kept_start
        utr_sequence = run_blastdbcmd(blastcmd_entries).replace('\n', '')

        features_from_line = count_nucleotides_and_dinucleotides(joint_name, utr_sequence)["counts"] + [len(utr_sequence), len(exons_to_keep)]
        sequence_params_from_line = [joint_name, chrom, first_kept_start, first_kept_size, strand, utr_sequence]
        return sequence_params_from_line, features_from_line

    elif mode == 'rbpmap':
        print "mode detected: rbpmap"
        # print rbpmap_entries
        # print "\n"
        return [joint_name] + [e for e in rbpmap_entries]
        # return rbpmap_entries

    return []

if __name__ == '__main__':

    print "Proccessing", input_annotations_file, "..."
    output_lines = []
    rbpmap_entries = []
    features = []
    count = 0
    with open(input_annotations_file, 'r') as tsv:
        reader = csv.reader(tsv,  delimiter="\t")
        reader.next()
        for line in reader:
            if LIMIT_INPUT_LINES and count > LIMIT_INPUT_LINES:
                break
            count += 1

            if line[CDS_START_FIELD] == line[CDS_END_FIELD]:
                continue #skeep non-coding RNAs

            if bExtract_fasta:
                extracted_sequence_line, extracted_features = proccess_line(line, 'blastcmd')
                if len(extracted_sequence_line) > 0:
                    print extracted_sequence_line
                    output_lines.append(extracted_sequence_line)
                    features.append([extracted_sequence_line[i] for i in [0, 6, 7]])
            if bCreate_RBPmap_entries:
                extracted_sequence_line = proccess_line(line, 'rbpmap')
                if len(extracted_sequence_line) > 0:
                    print extracted_sequence_line
                    rbpmap_entries.append(extracted_sequence_line)

    output_header_fasta = ["name", "chr", "first_exon_start", "first_exon_size", "strand", "5'_UTR", "UTR_length", "num_of_exons"]
    output_header_rbp = ["name", "chr", "input sequences to RBP"]

    if bExtract_fasta:
        g_df_features = pd.DataFrame(features, columns=["name", "utr_length", "num_exons"]).set_index("name")
        print g_df_features
        exit()

        with open(utr_output_tsv, 'w') as tsv:
            writer = csv.writer(tsv,  delimiter="\t")
            writer.writerow(output_header_fasta)
            writer.writerows(output_lines)
        print("written tab delimeted file to {}".format(utr_output_tsv))

        with open(utr_output_fasta, 'w') as fasta_file:
            for line in output_lines:
                fasta_file.write(">" + line[0]+"\n")
                fasta_file.write(line[5]+"\n\n")
            print("written fasta file to {}".format(utr_output_fasta))

    if bCreate_RBPmap_entries:
        with open(rbpmap_entries_file, 'w') as rbpmap_file:
            writer = csv.writer(rbpmap_file,  delimiter="\t")
            writer.writerow(output_header_rbp)
            writer.writerows([[line] for line in rbpmap_entries])
            print("written rbpmap file to {}".format(rbpmap_entries_file))


            # get_RBP_motifs_all_genes(rbpmap_file)

