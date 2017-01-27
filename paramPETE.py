#!/usr/bin/python3

in_file = "/home/mayere/Analyse_RNAseq/Genomes/Ensembl_Mmul86/Macaca_mulatta.Mmul_8.0.1.86.gtf"
gene_list_file ="/home/mayere/Analyse_RNAseq/Analysis/TF_discovery/TF_discov/Gene_list_EPI_vs_PETE_50avgdiff.txt"
jasp_motifs_file = "/home/mayere/Analyse_RNAseq/Analysis/TF_discovery/Jaspar_sites/pfm_vertebrates.txt"
fasta_file = "/home/mayere/Analyse_RNAseq/Genomes/Ensembl_Mmul86/Macaca_mulatta.Mmul_8.0.1.dna.chromosome.1.fa"
jsonFile = 'resulstJsonPETE.json'
TSVFilse = 'results_tsvPETE.txt'
limite_info= dict( gff_id = [str(i) for i in range(1,21)]+['X','Y'], gff_source = ['ensembl'], gff_type = ['gene'])
DISTANCE = 4000
NUMPROC = 8
