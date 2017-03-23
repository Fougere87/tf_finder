#!/usr/bin/python3

from BCBio.GFF import GFFExaminer
from BCBio import GFF
from Bio import SeqIO, motifs
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature , FeatureLocation
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from collections import Counter
import os, sys, gc
import numpy



class FeatPosition:
    def __init__(self, id, name, chrom, start, end, strand):
        self.id = id
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand

class AligFeat:
    def __init__(self, gene_id, gene_name,  id, name, start, end, strand, score):
        self.id = id
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.score = score
        self.gene_id = gene_id
        self.gene_name = gene_name




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
            print('Gene {} non trouvé'.format(line[0:len(line)-1]))
        try :
            location = FeatPosition(id = line[0:len(line)-1],
                                    chrom = feat.qualifiers['chrom'],
                                    name = feat.qualifiers['gene_name'][0],
                                    start = feat.location.nofuzzy_start,
                                    end = feat.location.nofuzzy_end,
                                    strand = feat.strand )
            locations_list.append(location)
        except KeyError :
            print('Pas de gene_name pour {}.'.format(line[0:len(line)-1]))
            location = FeatPosition(id = line[0:len(line)-1],
                                    chrom = feat.qualifiers['chrom'],
                                    name = '',
                                    start = feat.location.nofuzzy_start,
                                    end = feat.location.nofuzzy_end,
                                    strand = feat.strand )
            locations_list.append(location)
    gene_list.close()
    return locations_list

def finding_sequences(fasta, features, distance, mode = "Cache", direc = "both") :
    """Extracts the up and downstream sequences from the gnenomic fasta file.
    Possible to set the distance from start.
    If you are short in RAM you can use mode = "Index"."""
    print('Indexing genome')
    if mode == "Cache" :
        genome_dict = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    elif mode == "Index"    :
        genome_dict = SeqIO.index(fasta, "fasta")
    print('Done')
    sequences_list = []
    for feat in features :
        if feat.strand == '+1' :
            if direc == "both" :
                sequence  = genome_dict[feat.chrom].seq[(feat.start)-distance:(feat.start+distance)]
                sequence.alphabet= IUPAC.unambiguous_dna
                seqAnnot = SeqRecord(sequence)
                seqAnnot.annotations["start"] =  (feat.start)-distance
                seqAnnot.annotations["end"] = feat.start
            elif direc == "up" :
                sequence  = genome_dict[feat.chrom].seq[(feat.start)-distance:(feat.start)]
                sequence.alphabet= IUPAC.unambiguous_dna
                seqAnnot = SeqRecord(sequence)
                seqAnnot.annotations["start"] =  (feat.start)-distance
                seqAnnot.annotations["end"] = feat.start
            elif direc == "down" :
                sequence  = genome_dict[feat.chrom].seq[feat.start:(feat.start+distance)]
                sequence.alphabet= IUPAC.unambiguous_dna
                seqAnnot = SeqRecord(sequence)
                seqAnnot.annotations["start"] =  (feat.start)
                seqAnnot.annotations["end"] = feat.start+distance
            seqAnnot = SeqRecord(sequence)
            seqAnnot.id = feat.id
            seqAnnot.name = feat.name
            seqAnnot.annotations["chr"] = feat.chrom
            seqAnnot.annotations["strand"] = feat.strand
        else :
            if direc == "both" :
                sequence  = genome_dict[feat.chrom].seq[(feat.end)-distance:(feat.end)+distance].reverse_complement()
                sequence.alphabet= IUPAC.unambiguous_dna
                seqAnnot = SeqRecord(sequence)
                seqAnnot.annotations["start"] =  (feat.end)-distance
                seqAnnot.annotations["end"] = feat.end+distance
            elif direc == "down" :
                sequence  = genome_dict[feat.chrom].seq[(feat.end)-distance:(feat.end)].reverse_complement()
                sequence.alphabet= IUPAC.unambiguous_dna
                seqAnnot = SeqRecord(sequence)
                seqAnnot.annotations["start"] =  (feat.end)-distance
                seqAnnot.annotations["end"] = feat.end
            elif direc == "up" :
                sequence  = genome_dict[feat.chrom].seq[feat.end:(feat.end+distance)].reverse_complement()
                sequence.alphabet= IUPAC.unambiguous_dna
                seqAnnot = SeqRecord(sequence)
                seqAnnot.annotations["start"] =  (feat.end)
                seqAnnot.annotations["end"] = feat.end+distance
            seqAnnot.id = feat.id
            seqAnnot.name = feat.name
            seqAnnot.annotations["chr"] = feat.chrom
            seqAnnot.annotations["strand"] = feat.strand
        sequences_list.append(seqAnnot)
    del(genome_dict)
    return sequences_list


def writeFasta(records, output) :
    with open(output, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")

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

def align_motif(mot, sequenceR, precision = 10**4, balance = 100000, pseudocounts = 0.1) :
    mot.background = background(sequenceR.seq)
    mot.pseudocounts = pseudocounts
    distribution = mot.pssm.distribution(mot.background, precision=precision)
    threshold = distribution.threshold_balanced(balance)
    aligList = [(position+1, score) for position, score in enumerate(mot.pssm.calculate(sequenceR.seq)) if (position > 0)and(score > threshold)]
    for alig in aligList :
        sequenceR.features.append(AligFeat(id = mot.matrix_id, name = mot.name,  start = alig[0], end = alig[0] + len(mot), strand = 1, score =  alig[1], gene_id=sequenceR.id, gene_name=sequenceR.name))

    seqRC = sequenceR.seq.reverse_complement()
    mot.background = background(seqRC)
    distribution = mot.pssm.distribution(background=mot.background, precision=precision)
    threshold = distribution.threshold_balanced(balance)
    aligListNeg = [(position+1, score) for position, score in enumerate(mot.pssm.calculate(seqRC)) if (position > 0)and(score > threshold)]
    for alig in aligListNeg :
        sequenceR.features.append(AligFeat(id = mot.matrix_id, name = mot.name,  start = alig[0], end = alig[0] + len(mot), strand = -1, score =  alig[1], gene_id=sequenceR.id, gene_name=sequenceR.name))

    print('Found %i motifs' %len(sequenceR.features))
    return sequenceR.features



def testingAllSeq(prod, return_list) :
    for c in prod :
        print('Testing gene %s with motif %s' % (c[1].id, c[0].matrix_id))
        alig_seq = align_motif(mot=c[0], sequenceR=c[1])
        if alig_seq != [] :
            return_list.append(alig_seq)


def encoder(obj) :
    if hasattr(obj, "__dict__") :
        return obj.__dict__
    elif isinstance(obj, numpy.floating) :
        return float(obj)
    else:
        raise TypeError("Unserializable object {} of type {}".format(obj, type(obj)))
