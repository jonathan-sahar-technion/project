
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
# location = "home"
location = "lab"
# location = "ATLAS"

if location == "home":
    workspace = "C:\\Users\\Jonathan\\Documents\\translation_project_data\\"
    rbp_data_folder = workspace + "rbp_motif_analysis\\"
    input_folder = workspace + "input\\"
    output_folder = workspace + "output\\"

elif location == "lab":
    workspace = "/home/jonathan/Documents/data/"
    blastdb_path = workspace + "blastdb/mm9"

elif location == "ATLAS":
    workspace = "/srv01/technion/jonathans/data/"
    blastdb_path = "/storage/md_reut/footprint/mm9/blastdb/mm9"

if location in ["lab", "server"]:
    rbp_data_folder = workspace + "rbp_motif_analysis/"
    input_folder = workspace + "input/"
    output_folder = workspace + "output/"

rbpMap_summary_file_name = workspace + "All_Predictions.txt"
parsed_summary_file_name = workspace + "All_Predictions_parsed"

rbp_motifs_per_gene_file = rbp_data_folder + "motifs_by_gene_name.tsv"
genes_per_protein_file = rbp_data_folder + "genes_by_protein_name.tsv"
rbp_occurances_per_gene = rbp_data_folder + "RBPs_occurances_per_gene.tsv"
kmer_occurances_per_gene = workspace + "kmer_occurances_per_gene.tsv"

input_annotations_file = input_folder + "mm9_ensGene_eric.gpe"

utr_output_tsv = output_folder + "UTR_sequences.tsv"
utr_output_fasta = output_folder + "UTR_sequences.fa"

blastdb_entries_file = input_folder + "blastdb_entries.txt"
rbpmap_entries_file = input_folder + "RBPmap_entries.txt"

LIMIT_INPUT_LINES = 10