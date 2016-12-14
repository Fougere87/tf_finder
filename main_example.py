#!/usr/bin/python3

from feature_finder import *

in_file = "/home/mayere/Analyse_RNAseq/Genomes/Ensembl_Mmul86/Macaca_mulatta.Mmul_8.0.1.86.gtf"
gene_list_file ="/home/mayere/Analyse_RNAseq/Analysis/TF_discovery/TF_discov/gene_list.txt"
jasp_motifs_file = "/home/mayere/Analyse_RNAseq/Analysis/TF_discovery/Jaspar_sites/pfm_vertebrates.txt"
fasta_file = "/home/mayere/Analyse_RNAseq/Genomes/Ensembl_Mmul86/Macaca_mulatta.Mmul_8.0.1.dna.chromosome.1.fa"
limite_info= dict( gff_id = [str(i) for i in range(1,21)]+['X','Y'], gff_source = ['ensembl'], gff_type = ['gene'])

gene_dict = get_features_from_gff(gff_file=in_file , limite_info=limite_info)
print('Done extracting features from gff')
genes_positions = positions(gene_list_file, gene_dict=gene_dict)
print('Done genes positions')
sequences_records = finding_sequences(fasta=fasta_file,  features=genes_positions)
print('Sequences extracted')
motifs = motifs_list(jasp_motifs_file)
print('Motifs extracted ', len(motifs), " are to be tested per sequence")

results = {}
for s in sequences_records :
    print('Testing gene %s' % s.id )
    bckgrnd = background(s.seq)
    for m in motifs :
        print('testing motif', m.matrix_id)
        alignment = align_motif(mot=m, sequence=s.seq, background= bckgrnd)
        results[s.id] = alignment

print(results)
