#!/usr/bin/python3

import  pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import os sys


def get_features_from_gff(gff_file, li) :
    """This function returns a dict object containing all the features
    (seqRecord objects) selected by featType and chromosomes contained
    in the provided GFF file. The chromosome list must be a list of strings"""
    limite_info= dict( gff_id = chromosomes , gff_source = ['ensembl'], gff_type = [featType])
    handle = open(gff_file)
    gene_dict = {}
    for rec in GFF.parse(handle, limit_info=limite_info) :
        for feat in rec.features :
            feat.qualifiers["chrom"] = rec.id
            gene_dict[feat.qualifiers["gene_id"][0]] = feat
    handle.close()
    return gene_dict

class FeatPosition:
    def __init__(self, id, chrom, start, end, strand):
        self.id = id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand


def positions(gene_list, gene_dict) :
    locations_dict = {}
    for line in gene_list_file :
        try :
            feat = gene_dict[line[0:len(line)-1]]
            location = FeatPosition(id = line[0:len(line)-1],
                                    chrom = feat.qualifiers['chrom'],
                                    start = feat.location.nofuzzy_start,
                                    end = feat.location.nofuzzy_end,
                                    strand = feat.location.strand  )
            locations_dict[line[0:len(line)-1]] = location
            #print(line[0:len(line)-1]+gene_dict[line[0:len(line)-1]].qualifiers['chrom']+gene_dict[line[0:len(line)-1]].location)
        except KeyError :
            print('Gene %s non trouvé', line[0:len(line)-1])
    gene_list_file.close()
    return locations_dict

in_file = "/home/mayere/Analyse_RNAseq/Genomes/Ensembl_Mmul86/Macaca_mulatta.Mmul_8.0.1.86.gtf"
gene_list_file =open("/home/mayere/Analyse_RNAseq/Analysis/TF_discovery/TF_discov/gene_list.txt", 'r')
limite_info= dict( gff_id = [str(i) for i in range(1,21)]+['X','Y'], gff_source = ['ensembl'], gff_type = ['gene'])

get_features_from_gff(gff_file=in_file , limite_info=limite_info)

print(gene_dict.keys())
