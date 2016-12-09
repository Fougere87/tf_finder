from BCBio.GFF import GFFExaminer
import sys pprint


"""Examine the GFF to know it's structure"""

def examine_GFF(gff_file) :
    handle = open(gff_file)
    pprint.pprint(examiner.available_limits(handle))
    handle.close()



#### main

examine_GFF(sys.argv[1])
