#!/usr/bin/python3

from feature_finder import *
import multiprocessing as mp
from itertools import product

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
print('Motifs extracted', len(motifs), "are to be tested per sequence")


from itertools import product

argum = list(product(motifs, sequences_records))
manager = mp.Manager()
return_dict = manager.dict()
jobs = []
numproc = 8
for chunk in range(numproc) :
    #print(argum[int(chunk*len(argum)/numproc):int((chunk+1)*len(argum)/numproc)])
    p = mp.Process(target=testingAllSeq_monoArg2, args=(argum[int(chunk*len(argum)/numproc):int((chunk+1)*len(argum)/numproc)], return_dict,))
    jobs.append(p)
    p.start()

for proc in jobs :
    proc.join()
print(return_dict)
