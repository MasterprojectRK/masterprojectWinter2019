import csv

#write list of dictionaries to BED file 
#the keys specify the fieldnames
def writeBED(pOutfileName, pListOfDictsToWrite, pTrackString="", \
             pBrowserString="", pCommentString="", pHeaderBool = True):
    fieldnames = list(pListOfDictsToWrite[0].keys())
    with open(pOutfileName, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t", quotechar="#")
        header = []
        if(pTrackString != ""):
            header.append("track: " + pTrackString.rstrip())
            header.append("\n")
        if(pBrowserString != ""):
            header.append("browser: " + pBrowserString.rstrip())
            header.append("\n")
        if(pCommentString != ""):
            header.append("# " + pCommentString.rstrip())
            header.append("\n")
        if pHeaderBool:
            commentedHeader = "#"
            for name in fieldnames:
                commentedHeader += (name + "\t")
            header.append(commentedHeader.rstrip())
            header.append("\n")
        if len(header) != 0:
            for row in header:
                outfile.write(row)

        for row in pListOfDictsToWrite:
            writer.writerow(row)

def writeNarrowPeak(pOutfileName, pListOfDictsToWrite):
    fieldnames = list(pListOfDictsToWrite[0].keys())
    if len(fieldnames) != 10:
        print("Warning: wrong number of columns in narrowPeak file")

    with open(pOutfileName, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter="\t")
        for row in pListOfDictsToWrite:
            writer.writerow(row)