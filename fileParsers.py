import csv

def __fuzznucMismatchToNumber(pMismatchString):
    #return numeric value for number of mismatches
    #in the result files, zero mismatch is given as . instead of number 0
    if pMismatchString.isnumeric():
        return int(pMismatchString)
    else:
        return 0 

def parseFuzznuc(pFuzznucResultFileString):
    #parse result file from fuzznuc in excel format
    #which is a tab-delimited text file
    result = []
    with open(pFuzznucResultFileString) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            if not row["SeqName"] == "SeqName": #there are additional headings within the data, ignore these rows
                dDict = {"SeqName": row["SeqName"], "Start": int(row["Start"]), "End": int(row["End"]), \
                        "Strand": row["Strand"], "Motif": row["Pattern"].strip("pattern:"), \
                        "Mismatch": __fuzznucMismatchToNumber(row["Mismatch"]) }
                result.append(dDict)
    return result

def parseFimo(pFimoResultFileString):
    #parse result file in tsv-format from fimo
    #which is a tab-delimited text file
    result = []
    with open(pFimoResultFileString) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            if not row["motif_id"].startswith("#"): #ignore rows with comments
                dDict = { "SeqName": row["sequence_name"], "Start": int(row["start"]), "End": int(row["stop"]), \
                        "Strand": row["strand"], "Motif": row["motif_id"], \
                        "PVal": float(row["p-value"]), "QVal": float(row["q-value"]) }
                result.append(dDict)
    return result

def parseNarrowpeak(pNarrowpeakFileString):
    #parse a narrowpeak file
    #this is a BED text file with the first 6 columns as in standard BED plus 4 special columns 
    result = []
    fieldnames = ["Chromosome", "Start", "End", "Name", "Score", "Strand", "SignalValue", "PVal", "QVal", "Peak"]
    with open(pNarrowpeakFileString) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t", fieldnames=fieldnames)
        for row in reader:
            dDict = { "Chromosome": row["Chromosome"], \
                    "Start": int(row["Start"]), \
                    "End": int(row["End"]), \
                    "Name": row["Name"], \
                    "Score": float(row["Score"]) , \
                    "Strand": row["Strand"] ,\
                    "SignalValue": float(row["SignalValue"]) ,\
                    "PVal": float(row["PVal"]) ,\
                    "QVal": float(row["QVal"]) ,\
                    "Peak": int(row["Peak"])
                    }
            result.append(dDict)
    return result

def parseBedFile(pBedFileString):
    #parse a standard BED file
    #only the first 6 (standard) columns for now
    #ignore comments (#), track, and browser data
    result = []
    fieldnames = ["SeqName", "Start", "End", "FeatureName",	"Score", "Strand", \
                  "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blckStarts"]
    with open(pBedFileString) as csvfile:
        reader = csv.DictReader([row for row in csvfile if not row.startswith('#') \
                                 and not row.startswith("track")
                                 and not row.startswith("browser")], \
                                 fieldnames=fieldnames, 
                                 delimiter="\t")
        for row in reader:
            dDict = { "SeqName": row["SeqName"], \
                        "Start": int(row["Start"]), \
                        "End": int(row["End"]), \
                        "FeatureName": row["FeatureName"], \
                        "Score": float(row["Score"]), \
                        "Strand": row["Strand"]}
                        #remaining fields currently not needed
            result.append(dDict)
    return result
