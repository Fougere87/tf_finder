
import  pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import os sys


def get_features_from_gff(gff_file, featType ) :
    """This function returns
    limite_info= dict( gff_id = [str(i) for i in range(1,20)], gff_source = ['ensembl'], gff_type = [featType])
    handle = open(gff_file)
    gene_dict = {}
    for rec in GFF.parse(handle, limit_info=limite_info) :
        for feat in rec.features :
            gene_dict[feat.qualifiers["gene_id"][0]] = feat
    return gene_dict
print(gene_dict.keys())
