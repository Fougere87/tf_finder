#!/usr/bin/python3

from BCBio.GFF import GFFExaminer
from BCBio import GFF
from Bio import SeqIO, motifs, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
import os, sys, gc

class FeatPosition:
    def __init__(self, id, chrom, start, end, strand):
        self.id = id
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

class Alignement(SeqRecord):
    def _init_(self,id,name, tf_mot_id, tf_mot_name, strand, start, end, chr, chrstart, chrend):



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
        except KeyError :
            print('Gene %s non trouvé', line[0:len(line)-1])
        try :
            location = FeatPosition(id = line[0:len(line)-1],
                                    chrom = feat.qualifiers['chrom'],
                                    name = feat.qualifiers['gene_name'],
                                    start = feat.location.nofuzzy_start,
                                    end = feat.location.nofuzzy_end,
                                    strand = feat.location.strand )
            locations_list.append(location)
        except KeyError :
            print('Pas de gene_name pour %s.', line[0:len(line)-1])
            location = FeatPosition(id = line[0:len(line)-1],
                                    chrom = feat.qualifiers['chrom'],
                                    name = '',
                                    start = feat.location.nofuzzy_start,
                                    end = feat.location.nofuzzy_end,
                                    strand = feat.location.strand )
            locations_list.append(location)
    gene_list.close()
    return locations_list

def finding_sequences(fasta, features, distance=DISTANCE) :
    """Extracts the ustream sequences from the gnenomic fasta file. Possible to set the distance from start"""
    genome_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    sequences_list = []
    for feat in features :
        if feat.strand == '+1' :
            sequence  = genome_dict[feat.chrom].seq[(feat.start-1)-distance:(feat.start-1)]
            sequence.alphabet= IUPAC.unambiguous_dna
            seqAnnot = SeqRecord(sequence)
            seqAnnot.id = feat.id
            seqAnnot.name = feat.name
            seqAnnot.annotations["start"] =  (feat.start-1)-distance
            seqAnnot.annotations["end"] = feat.start-1
            seqAnnot.annotations["chr"] = feat.chrom
            seqAnnot.annotations["strand"] = feat.strand
        else :
            sequence  = genome_dict[feat.chrom].seq[(feat.end+1):(feat.end+1+distance)].reverse_complement()
            sequence.alphabet= IUPAC.unambiguous_dna
            seqAnnot = SeqRecord(sequence)
            seqAnnot.id = feat.id
            seqAnnot.name = feat.name
            seqAnnot.annotations["start"] =  feat.end+1
            seqAnnot.annotations["end"] = feat.end+1+distance
            seqAnnot.annotations["chr"] = feat.chrom
            seqAnnot.annotations["strand"] = feat.strand
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

def align_motif(mot, sequence, background, precision = 10**4, balance = 100000, pseudocounts = 0.1) :
    mot.background = background
    mot.pseudocounts = pseudocounts
    distribution = mot.pssm.distribution(background=background, precision=precision)
    threshold = distribution.threshold_balanced(balance)
    alig_list = [(position+1, score) for position, score in enumerate(mot.pssm.calculate(sequence)) if (position > 0)and(score > threshold)]
    for alig in alig_list :
        sequence.features.append(SeqFeature(id = mot+"."+str(alig[0]), location = SeqFeature.FeatureLocation(start = alig[0], end = "longueur motif à trouver", ))
    print('Found %i motifs' %len(ret_list))
    return ret_list

def testingAllSeq(motifs, sequences) :
    results = {}
    for s in sequences :
        print('Testing gene %s' % s.id )
        bckgrnd = background(s.seq)
        for m in motifs :
            print('testing motif', m.matrix_id)
            alignment = align_motif(mot=m, sequence=s.seq, background= bckgrnd)
            results[s.id] = alignment
    return results


def testingAllSeq_monoArg2(comb, return_dict) :

    for c in comb :
        print('Testing gene %s with motif %s' % (c[1].id, c[0].matrix_id))
        bckgrnd = background(c[1].seq)
        alignment = align_motif(mot=c[0], sequence=c[1].seq, background= bckgrnd)
        return_dict[str((c[1].id,c[0].matrix_id))] = alignment
