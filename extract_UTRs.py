#!/usr/bin/python

import csv
import pandas as pd
import sys


# exon-starts are 0-based, exon-ends are 1-based
# my_input_file = sys.argv[1]
# my_output_file = sys.argv[2]

bServer = False

data_dir = "/home/jonathan/Documents/data"
if bServer:
    data_dir = "/home/jonathan/data"

input_file = data_dir + "/input/mm9_ensGene_eric.gpe"
output_file_mulitple_exons = data_dir + "/output/extract_utrs_output_multi.bed"
output_file_single_exons = data_dir + "/output/extract_utrs_output_single.bed"

START = 0
END = 1000000000


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
all_lines = list()
all_lines_single_exon = list()

print "Proccessing ", input_file, "..."

count = 0
with open(input_file, 'r') as tsv:
    reader = csv.reader(tsv,  delimiter="\t")
    reader.next()
    for line in reader:
        if line[CDS_START_FIELD] == line[CDS_END_FIELD]:
            count += 1
            continue #skeep non-coding RNAs

        # extract info from line
        strand = line[STRAND_FIELD]
        name = line[NAME_FIELD]
        name2 = line[NAME2_FIELD]
        chrom = line[CHROM_FIELD]
        tx_start = line[TX_START_FIELD]
        tx_end = line[TX_END_FIELD]

        # break the list fields into lists
        exon_starts =  [int(x) for x in line[EXON_STARTS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list
        exon_ends =  [int(x) for x in line[EXON_ENDS_FIELD].split(",")[:-1]] # -1 for removing the empty member at the end of the list

        if strand == '+':
            cds_start = int(line[CDS_START_FIELD])

            # 3' ----------------------------------------  5'
            # 5' -----ex_start-----ex_end-----cds_start--> 3'

            # 3' ----------------------------------------  5'
            # 5' -----ex_start----cds_start-----ex_end---> 3'

            exons_to_keep = [index for index, value in enumerate(exon_starts) if int(value) < cds_start]
            starts_to_keep = [str(ex) for i, ex in enumerate(exon_starts) if i in exons_to_keep]

        elif strand == '-':

            # 3' <---cds_start---ex_start-----ex_end----- 5'
            # 5' ---------------------------------------- 3'

            # 3' <-----ex_start----cds_start-----ex_end-- 5'
            # 5' ---------------------------------------- 3'

            cds_end = int(line[CDS_END_FIELD])
            exons_to_keep = [index for index, value in enumerate(exon_ends) if int(value) > cds_end]
            ends_to_keep = [str(ex) for i, ex in enumerate(exon_ends) if i in exons_to_keep]

        all_block_sizes = [end - start for start, end in zip(exon_starts, exon_ends)] #ends are 1-based, starts are 0-based1

        block_sizes =  ",".join([str(b) for i,b in enumerate(all_block_sizes) if i in exons_to_keep])
        block_starts = ",".join([str(ex) for i, ex in enumerate(exon_starts) if i in exons_to_keep])
        num_blocks = len(exons_to_keep)
        # new_lines = [[chrom, start, end, name + "_" + name2, SCORE, strand, start, end, RGB, len(exons_to_keep), block_sizes, block_starts]\
        #              for start,end in zip(exon_starts, exon_ends)]

        if num_blocks == 1:
           new_line = [[chrom, block_starts[0], block_starts[0] + block_sizes[0], name2+ "_" + name, SCORE, strand]] # list with on elist as a member to make += work instead of append...
           all_lines_single_exon += new_line

        if num_blocks > 1:
            new_line = [[chrom, tx_start, tx_end, name2 + "_" + name, SCORE, strand, line[CDS_START_FIELD], line[CDS_END_FIELD], RGB, num_blocks, block_sizes, block_starts]] # list with on elist as a member to make += work instead of append...
            all_lines += new_line

        count += 1

output_header_single_exon = ["chrom", "start", "end", "name", "score", "strand"]
with open(output_file_single_exons, 'w') as tsv:
    writer = csv.writer(tsv,  delimiter="\t")
    print "writing to ", output_file_single_exons, "..."
    writer.writerow(output_header_single_exon)
    writer.writerows(all_lines_single_exon)



output_header_multi = ["chrom", "start", "end", "name", "score", "strand", "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts "]
with open(output_file_mulitple_exons, 'w') as tsv:
    writer = csv.writer(tsv,  delimiter="\t")
    print "writing to ", output_file_mulitple_exons, "..."
    writer.writerow(output_header_multi)
    writer.writerows(all_lines)



