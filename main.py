#!/usr/bin/python3

from feature_finder import *
import multiprocessing as mp
from itertools import product
import json
from itertools import product
from param import *
import configparser
import sys
import subprocess
import shlex

confFile = sys.argv[1]
conf = configparser.ConfigParser()
conf.read(confFile)
gff_source = [conf.get('General', 'gff_source')]
gff_type =[conf.get('General', 'gff_type')]
chrom_num = int(conf.get('General', 'chrom_num'))

distance = int(conf.get('General', 'DISTANCE'))
gff_file = conf.get('Inputs', 'gff_file')
gene_list_file = conf.get('Inputs', 'gene_list_file')
fasta_file = conf.get('Inputs', 'fasta_file')
xxMotifs = conf.get('External', "XXmotifs_cmd")
out_fasta = conf.get('Outputs', 'out_fasta')
direction  = conf.get('General', 'direction')

##### essayer de faire une structure try pour voir si un fichier de séquence est disponible
##### ===> si oui l'enregistrer comme sequence_records et sauter les 3 premieres étapes qui sont un peu longue.
##### Voir aussi s'il est possible d'indexer le génome et de récupérer sur le  disque plutôt que via RAM. ==> ok

gene_dict = get_features_from_gff(gff_file=gff_file , limite_info=dict( gff_id = [str(i) for i in range(1,chrom_num)]+['X','Y'], gff_source =gff_source, gff_type = gff_type))
print('Done extracting features from gff')
genes_positions = positions(gene_list_file, gene_dict=gene_dict)
print('Done genes positions')
sequences_records = finding_sequences(fasta=fasta_file,  features=genes_positions, distance = distance, direc = direction )
print('Sequences extracted')

if not os.path.isdir(conf['DEFAULT']['out_path']) :
    print(conf['DEFAULT']['out_path'], "not found, creating it.")
    subprocess.run(['mkdir',conf['DEFAULT']['out_path']])

print('Writing Sequences to', out_fasta)
writeFasta(sequences_records, out_fasta)

print('Running XXmotifs on sequences')
xxMotifsCmd = shlex.split(xxMotifs)
print(xxMotifsCmd)
runXXmot = subprocess.run(xxMotifsCmd)
