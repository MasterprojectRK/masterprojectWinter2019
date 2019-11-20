import fileParsers
import fileWriters
import argparse
import csv

#convert various file formats to BED format
#supports fimo tsv output, fuzznuc excel output and narrowpeak for now


parser = argparse.ArgumentParser(description="convert various file formats to BED")
parser.add_argument("infile", type=str, help="input text file to convert to BED")
parser.add_argument("format", type=str, help= "format of input file", choices=["tsv", "fuzznuc", "narrowpeak"])
parser.add_argument("outfile", type=str, help="BED file to write")

args = parser.parse_args()

peaksToWrite = []
if args.format == "narrowpeak": 
    peaks = fileParsers.parseNarrowpeak(args.infile)
    for row in peaks:
        dDict = {"chrom": row["Chromosome"], \
                 "chromStart": str(row["Start"]), \
                 "chromEnd": str(row["End"]), \
                 "name": row["Name"], \
                 "score": str(row["Score"]) , \
                 "strand": row["Strand"] }
        peaksToWrite.append(dDict)

if args.format == "tsv":
    peaks = fileParsers.parseFimo(args.infile)
    for row in peaks:
        dDict = { "chrom": row["SeqName"], \
                    "chromStart": str(row["Start"]), \
                    "chromEnd": str(row["End"]), \
                    "name": row["Motif"], \
                    "score": str(row["PVal"]) , \
                    "strand": row["Strand"] }
        peaksToWrite.append(dDict)

if args.format == "fuzznuc":
    peaks = fileParsers.parseFuzznuc(args.infile)
    for row in peaks:
        dDict = { "chrom": row["SeqName"], \
                    "chromStart": str(row["Start"]), \
                    "chromEnd": str(row["End"]), \
                    "name": row["Motif"], \
                    "score": str(row["Mismatch"]) , \
                    "strand": row["Strand"] }
        peaksToWrite.append(dDict)


fileWriters.writeBED(pOutfileName = args.outfile, pListOfDictsToWrite = peaksToWrite)



