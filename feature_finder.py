#!/usr/bin/python3

from BCBio.GFF import GFFExaminer
from BCBio import GFF
from Bio import SeqIO, motifs
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
import os, sys, gc

class FeatPosition:
    def __init__(self, id, chrom, start, end, strand):
        self.id = id
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

def get_features_from_gff(gff_file, limite_info) :
    """This function returns a dict object containing all the features
    (as seqRecord objects) selected by featType and chromosomes contained
    in the provided GFF file. The chromosome list must be a list of strings"""
    handle = open(gff_file, 'r')
    gene_dict = {}
    for rec in GFF.parse(handle, limit_info=limite_info) :
        for feat in rec.features :
            feat.qualifiers["chrom"] = rec.id
            gene_dict[feat.qualifiers["gene_id"][0]] = feat
    handle.close()
    return gene_dict

def positions(gene_list_file, gene_dict) :
    """Extract the position of the genes of interest from a gene dict"""
    locations_list = []
    gene_list  = open(gene_list_file, 'r')
    for line in gene_list :
        try :
            feat = gene_dict[line[0:len(line)-1]]
            location = FeatPosition(id = line[0:len(line)-1],
                                    chrom = feat.qualifiers['chrom'],
                                    start = feat.location.nofuzzy_start,
                                    end = feat.location.nofuzzy_end,
                                    strand = feat.location.strand )
            locations_list.append(location)
        except KeyError :
            print('Gene %s non trouvé', line[0:len(line)-1])
    gene_list.close()
    return locations_list

def finding_sequences(fasta, features, distance=10000) :
    """Extracts the ustream sequences from the gnenomic fasta file. Possible to set the distance from start"""
    genome_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    sequences_list = []
    for feat in features :
        if feat.strand == '+1' :
            sequence  = genome_dict[feat.chrom].seq[(feat.start-1)-distance:(feat.start)]
            sequence.alphabet= IUPAC.unambiguous_dna
            seqAnnot = SeqRecord(sequence)
            seqAnnot.id = feat.id
        else :
            sequence  = genome_dict[feat.chrom].seq[(feat.end-1):(feat.end+distance-1)].reverse_complement()
            sequence.alphabet= IUPAC.unambiguous_dna
            seqAnnot = SeqRecord(sequence)
            seqAnnot.id = feat.id
        sequences_list.append(seqAnnot)
    del(genome_dict)
    return sequences_list

def background(sequence) :
    """Calculates the background of a sequence : the rate of ACTG  and returns it as a dictionnary"""
    bckgrnd = {}
    bckgrnd['A'] = sequence.count('A')/(len(sequence)-sequence.count('N'))
    bckgrnd['C'] = sequence.count('C')/(len(sequence)-sequence.count('N'))
    bckgrnd['G'] = sequence.count('G')/(len(sequence)-sequence.count('N'))
    bckgrnd['T'] = sequence.count('T')/(len(sequence)-sequence.count('N'))
    return(bckgrnd)

def motifs_list(jasp_motifs_file) :
    jasp_motifs = open(jasp_motifs_file, 'r')
    motifs_list = []
    for m in motifs.parse(jasp_motifs, "jaspar") :
        motifs_list.append(m)
    jasp_motifs.close()
    return motifs_list

def align_motif(mot, sequence, background, precision = 10**4, balance = 1000, pseudocounts = 0.5) :
    mot.background = background
    mot.pseudocounts = pseudocounts
    distribution = mot.pssm.distribution(background=background, precision=precision)
    threshold = distribution.threshold_balanced(balance)
    ret_list = [(position, score) for position, score in enumerate(mot.pssm.calculate(sequence)) if (position > 0)and(score > threshold)]
    return ret_list
