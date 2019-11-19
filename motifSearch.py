from Bio.Emboss.Applications import FuzznucCommandline
import os
import argparse
import csv

parser = argparse.ArgumentParser(description='find motifs using EMBOSS fuzznuc')
parser.add_argument("sequenceFasta", type=str, help="path to the reference sequence fasta file")
parser.add_argument("motifFile", type=str, help="text file with motif(s) to search for in IUPAC notation, one motif per line")
parser.add_argument("outfile", type=str, help="path + name of outputfile")
parser.add_argument("--complement", type=str, default="yes", choices=["yes", "no"], help="also search in reverse sequence")

args = parser.parse_args()

def parseMotifFile(pFileNameString):
    motifs = []
    with open(pFileNameString) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            motifs.append(row)
    return motifs


#parse the motifs
directory,fname = os.path.split(args.outfile)
motifList = parseMotifFile(args.motifFile)
fileList =  [os.path.join(directory, fname + motif["name"] + ".fuzznuc") for motif in motifList ]
print("using fuzznuc to search for the following patterns:")
for motif in motifList:
    print(motif["motif"])
print("writing output to " + os.path.join(directory, fname))

#run fuzznuc on each motif
for motif, filename in zip(motifList, fileList):
    print("pattern: " + motif["motif"])
    print("mismatches allowed: " + motif["mismatches"])
    fuzznucResults = FuzznucCommandline(sequence=args.sequenceFasta, \
                                    pattern=motif["motif"], \
                                    complement=args.complement, \
                                    outfile=filename, \
                                    pmismatch=motif["mismatches"], \
                                    rformat="excel")
    stdout, stderr = fuzznucResults()
    print("done")
    

#combine into single outfile, remove temporary outputs
with open(args.outfile, 'w') as outfile:
    for filename in fileList:
        with open(filename) as infile:
            for line in infile:
                outfile.write(line)

for filename in fileList:
    os.remove(filename)