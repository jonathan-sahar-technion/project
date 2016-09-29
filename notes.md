* Connecting to ATLAS:
    - IP: 192.114.101.131
    - Hostname: jonathans@tech-ui01.hep.technion.ac.il
    - password: Pierrem3nard

* getting an interactive shell
   - qsub -l nodes=1:ppn=4 -q N -I

* interacting with blastdb
    - makeblastdb/blastdbcmd:

    - Make database:
      makeblastdb -in toy.fasta -input_type fasta -dbtype nucl -parse_seqids -out toy

    - Read range from a chromosome in database:
      blastdbcmd -db toy.fasta -entry chr_Single_Exon -range 1-10
      (works with database created with formatdb)

    - binary location on ATLAS server:
    /Local/ncbi-blast-2.2.30+/bin/blastdbcmd

