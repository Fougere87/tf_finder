#!/usr/bin/python3

from feature_finder import *
import multiprocessing as mp
from itertools import product
import json
from itertools import product
from paramPETE import *




gene_dict = get_features_from_gff(gff_file=in_file , limite_info=limite_info)
print('Done extracting features from gff')
genes_positions = positions(gene_list_file, gene_dict=gene_dict)
print('Done genes positions')
sequences_records = finding_sequences(fasta=fasta_file,  features=genes_positions)
print('Sequences extracted')
motifs = motifs_list(jasp_motifs_file)
print('Motifs extracted', len(motifs), "are to be tested per sequence")




argum = list(product(motifs, sequences_records))
manager = mp.Manager()
return_list = manager.list()
jobs = []

for chunk in range(NUMPROC) :
    #print(argum[int(chunk*len(argum)/numproc):int((chunk+1)*len(argum)/numproc)])
    p = mp.Process(target=testingAllSeq, args=(argum[int(chunk*len(argum)/NUMPROC):int((chunk+1)*len(argum)/NUMPROC)], return_list,))
    jobs.append(p)
    p.start()

for proc in jobs :
    proc.join()


for comb in return_list :
    for al in comb :
        print("{}   {}  {}  {}  {}  {}  {}  {}".format(al.id, al.name, al.gene_id, al.gene_name, al.strand, al.start, al.end, al.score  ))


# for k in return_list :
#     if return_dict[k].features != [] :
#         print(return_list[k].id, '\t', str(return_listg[k].annotations), '\t', return_dict[k].name, '\t', end="")
#         [print("{};{};{};{};{}".format(feat.id, feat.name, feat.score, feat.start, feat.strand), end="\t") for feat in return_dict[k].features]
#         print("\n")



with open('resulstJson.json', 'w') as f :
    json.dump(list(return_list), f, default= encoder)
with open(TSVFilse, 'w') as t :
    for comb in return_list :
        for al in comb :
            t.write("{}   {}  {}  {}  {}  {}  {}  {}\n".format(al.id, al.name, al.gene_id, al.gene_name, al.strand, al.start, al.end, al.score  ))
