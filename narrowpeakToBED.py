import fileParsers
import argparse
import csv

parser = argparse.ArgumentParser(description="convert narrowpeak file to fasta")
parser.add_argument("narrowpeakfile", type=str, help="narrowpeak formatted text file")
parser.add_argument("outfile", help="fasta file to write")

args = parser.parse_args()

#get the first 6 columns out of the narrowpeak file
peaks = fileParsers.parseNarrowpeak(args.narrowpeakfile)
peaksToWrite = []
for row in peaks:
    dDict = { "Chromosome": row["Chromosome"], \
                "Start": str(row["Start"]), \
                "End": str(row["End"]), \
                "Name": row["Name"], \
                "Score": str(row["Score"]) , \
                "Strand": row["Strand"] }
    peaksToWrite.append(dDict)

#write the bed file without the narrowpeak-specific columns
with open(args.outfile, "w", newline="") as outfile:
    fieldnames = list(peaksToWrite[0].keys())
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
    #writer.writeheader()
    for row in peaksToWrite:
        writer.writerow(row)




