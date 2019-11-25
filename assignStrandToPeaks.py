import argparse
import csv
import matplotlib.pyplot as plt
import matplotlib.transforms as tf
import matplotlib.colors as colors
import numpy as np
import fileParsers
import fileWriters
import matplotlib.ticker as ticker
from operator import itemgetter



def strandAssigner(pChipSeqList, pBindingSiteList, pDmax):
    #try to assign a strand for every chip seq peak
    #given putative positions of the forward and reverse binding sites
    #this needs to be done on a per-chromosome basis, since the positions are on a per-chromosome basis, too
    
    chromSet = set([row["Chromosome"] for row in pChipSeqList])

    for chromosome in chromSet:
        print("assigning strands for chromosome {0:s}".format(chromosome))
        chipSeqTmpList = [row for row in pChipSeqList if row["Chromosome"] == chromosome]
        posFwdList = [row for row in pBindingSiteList if row["SeqName"] == chromosome and row["Strand"] == "+"]
        posRevList = [row for row in pBindingSiteList if row["SeqName"] == chromosome and row["Strand"] == "-"]

        for element in chipSeqTmpList:
            #get the peak position 
            peakPos = element["Start"]
            if element["Peak"] >= 0:
                peakPos += element["Peak"] #"Peak" is given as offset from "Start", -1 if not specified
            else:
                peakPos += element["End"]
                peakPos = int( peakPos / 2) #just take the middle if no "Peak" available
            
            i = 1
            bestdiffFwd = abs(posFwdList[0]["Start"] - peakPos)
            while i < len(posFwdList):
                d = abs(posFwdList[i]["Start"] - peakPos)
                if d < bestdiffFwd:
                    bestdiffFwd = d
                i += 1
            
            j = 1
            bestdiffRev = abs(posRevList[0]["Start"] - peakPos)
            while j < len(posRevList):
                d = abs( posRevList[j]["Start"] - peakPos )
                if d < bestdiffRev:
                    bestdiffRev = d
                j += 1

            #both forward and reverse could be within tolerance range
            #for the moment, this is ignored and the closer one (if equal, the negative) strand is chosen
            if bestdiffRev <= bestdiffFwd and bestdiffRev <= pDmax:
                element["Strand"] = "-"
            elif bestdiffRev > bestdiffFwd and bestdiffFwd <= pDmax:
                element["Strand"] = "+"
            else:
                element["Strand"] = min(bestdiffFwd, bestdiffRev)
                #element["Strand"] = "."

def peakSorted(pPeakList, pColumnNameList):
    result = pPeakList
    name1 = pColumnNameList[0]
    name2 = pColumnNameList[1]
    if pPeakList and pPeakList[0][name1] and pPeakList[0][name2]:
        chromSet = set([row[name1] for row in pPeakList])
        chromSetNumeric = chromSet - set(["chrX", "chrY"])
        chromSetAlpha = chromSet - chromSetNumeric
        peakListNumeric = [ [ row, int(str.strip(row[name1],"chr")), row[name2] ] for row in pPeakList if row[name1] in chromSetNumeric]
        peakListAlpha = [ [ row, str.strip(row[name1],"chr"), row[name2] ] for row in pPeakList if row[name1] in chromSetAlpha]
        peakListNumeric.sort(key=itemgetter(1,2))
        peakListAlpha.sort(key=itemgetter(1,2))
        result = [row[0] for row in peakListNumeric] + [row[0] for row in peakListAlpha]

    return result

parser = argparse.ArgumentParser(description="assign strand information to ChIP-seq data")
parser.add_argument("bindingSites", type=str, nargs="+", help="List of BED-files specifying binding sites")
parser.add_argument("chipseq", type=str, help="ChIP-seq narrowpeak file")
parser.add_argument("outfile", type=str, help="Path + Name of narrowPeak output file")
parser.add_argument("--chromosomes", type=str, nargs="+", help="list of chromosomes, e.g. chr1, chr2, chrX; default: all")
parser.add_argument("--printGraphs", type=bool, default=False, help="output graphs")

args = parser.parse_args()

chromosomeSet = None
chipSeqGlobalList = fileParsers.parseNarrowpeak(args.chipseq)
availableChroms = set([row["Chromosome"] for row in chipSeqGlobalList])
if args.chromosomes == None:
    chromosomeSet = availableChroms
else:
    chromosomeSet = set(args.chromosomes) & availableChroms
chipSeqFilteredList = [row for row in chipSeqGlobalList if row["Chromosome"] in chromosomeSet]
chipSeqFilteredList = peakSorted(chipSeqFilteredList, ["Chromosome","Start"])

#extract putative binding sites from BED files
bindingSiteGlobalList = []
for bindingSiteFile in args.bindingSites:
    bindingSiteGlobalList += fileParsers.parseBedFile(bindingSiteFile)

bindingSiteFilteredList = [row for row in bindingSiteGlobalList if row["SeqName"] in chromosomeSet]
bindingSiteFilteredList = peakSorted(bindingSiteFilteredList, ["SeqName", "Start"])

#extract and print out some figures from the data  
peakLengthsList = [row["End"] - row["Start"] for row in chipSeqFilteredList]
meanPeakLengthFloat = np.mean(peakLengthsList)
medianPeakLengthFloat = np.median(peakLengthsList)
maxPeakLengthInt = np.max(peakLengthsList)
minPeakLengthInt = np.min(peakLengthsList)
dmaxInt = int(0.5 * medianPeakLengthFloat)
numberChipseqPeaksInt = len(chipSeqFilteredList)
print("ChIPseq data - globally")
print("number of ChIPseq peaks: {0:d}".format(numberChipseqPeaksInt ))
print("mean ChIPseq peak length: {0:.2f}".format(meanPeakLengthFloat))
print("median ChIPseq peak length: {0:.2f}".format(medianPeakLengthFloat))
print("max ChIPseq peak length: {0:d}".format(maxPeakLengthInt))
print("min ChIPseq peak length: {0:d}".format(minPeakLengthInt))
print("ChIPseq data - per chromosome")
for chromosome in sorted(list(chromosomeSet)):
    tmpChipseqList = [row for row in chipSeqFilteredList if row["Chromosome"] == chromosome]
    tmpNumberPeaksInt = len(tmpChipseqList)
    tmpPeakLengthList = [row["End"] - row["Start"] for row in tmpChipseqList]
    tmpMaxInt = np.max(tmpPeakLengthList)
    tmpMinInt = np.min(tmpPeakLengthList)
    tmpMedianFloat = np.median(tmpPeakLengthList)
    print("{0:s}:".format(chromosome))
    print("number of peaks: {0:d}".format(tmpNumberPeaksInt))
    print("median ChIPseq peak length: {0:.2f}".format(tmpMedianFloat))
    print("max ChIPseq peak length: {0:d}".format(tmpMaxInt))
    print("min ChIPseq peak length: {0:d}".format(minPeakLengthInt))

#try to assign a strand to the ChIP-seq peaks  
print("Found {0:d} putative binding sites in {1:d} bed file(s)".format(len(bindingSiteFilteredList), len(args.bindingSites)))
print("assigning peaks closer than {0:d} bp to putative binding sites".format(dmaxInt))
print("This will take some time...")
strandAssigner(chipSeqFilteredList, bindingSiteFilteredList, dmaxInt)
chipSeqFilteredList = peakSorted(chipSeqFilteredList, ["Chromosome","Start"])
fileWriters.writeNarrowPeak(args.outfile, [row for row in chipSeqFilteredList if row["Strand"] in ("+", "-")] )
          
chipSeqFwd = [row for row in chipSeqFilteredList if row["Strand"] == "+"]
chipSeqRev = [row for row in chipSeqFilteredList if row["Strand"] == "-"]
chipSeqUnass = [row for row in chipSeqFilteredList if row["Strand"] not in ("+", "-")]
numberChipseqFwdInt = len(chipSeqFwd)
numberChipseqRevInt = len(chipSeqRev)
numberChipseqUnassInt = len(chipSeqUnass)
chipSeqValFwd = [row["SignalValue"] for row in chipSeqFwd ]
chipSeqValRev = [row["SignalValue"] for row in chipSeqRev ]
chipseqValUnass = [row["SignalValue"] for row in chipSeqUnass ]
print("assigned to fwd. strand: {0:d} ({1:.1f}%)".format(numberChipseqFwdInt, 100*numberChipseqFwdInt/numberChipseqPeaksInt))
print("assigned to rev. strand: {0:d} ({1:.1f}%)".format(numberChipseqRevInt, 100*numberChipseqRevInt/numberChipseqPeaksInt))
print("unassigned peaks: {0:d} ({1:.1f}%)".format(numberChipseqUnassInt, 100*numberChipseqUnassInt/numberChipseqPeaksInt))
print("median sig. val. unassigned: {0:0.2f}".format(np.median(chipseqValUnass)))
print("median sig. val fwd: {0:.2f}".format(np.median(chipSeqValFwd)))
print("median sig. val rev: {0:.2f}".format(np.median(chipSeqValRev)))

print("narrowPeaks file with assigned strands has been written to {0:s}".format(args.outfile))

if args.printGraphs:
    #Draw everything in a graph
    fig, ax = plt.subplots()

    #x- and y-values for the graph
    chipSeqXFwd = [row["Start"] for row in chipSeqFwd]
    chipSeqXRev = [row["Start"] for row in chipSeqRev]
    chipSeqXUnass = [row["Start"] for row in chipSeqUnass]
    chipSeqYFwd = ["forward" for row in chipSeqXFwd]
    chipSeqYRev = ["reverse" for row in chipSeqXRev]
    chipSeqYUnass = ["unassigned" for row in chipSeqXUnass]


    chipSeqX = chipSeqXFwd + chipSeqXRev + chipSeqXUnass
    chipSeqY = chipSeqYFwd + chipSeqYRev + chipSeqYUnass
    chipSeqVal = chipSeqValFwd + chipSeqValRev + chipseqValUnass

    bindingSiteFwdList = [row for row in bindingSiteFilteredList if row["Strand"] == "+"]
    bindingSiteRevList = [row for row in bindingSiteFilteredList if row["Strand"] == "-"]
    fwdSitesXList = [(row["Start"] + row["End"])/2 for row in bindingSiteFwdList]
    revSitesXList = [(row["Start"] + row["End"])/2 for row in bindingSiteRevList]
    fwdSitesYList = ["forward" for row in fwdSitesXList]
    revSitesYList = ["reverse" for row in revSitesXList]
    bindingSitesXList = fwdSitesXList + revSitesXList
    bindingSitesYList = fwdSitesYList + revSitesYList

    """ #max and min values for the graph. Use max over all chromosomes, but obey pval cutoff
    print("fimo search statistics:")
    maxQFloat = max([row["QVal"] for row in fimoDataList if row["PVal"] < args.pval])
    minQFloat = min([row["QVal"] for row in fimoDataList if row["PVal"] < args.pval])
    print("max. q-value: ", maxQFloat)
    print("min. q-value: ", minQFloat)
    maxPFloat = max([row["PVal"] for row in fimoDataList if row["PVal"] < args.pval])
    minPFloat = min([row["PVal"] for row in fimoDataList if row["PVal"] < args.pval])
    print("max. p-value: ", maxPFloat)
    print("min. p-value: ", minPFloat) """

    #labels for the legend, include the names of the patterns
    """ fuzznucNameString = ""
    for pattern in fuzznucPatternSet:
        fuzznucNameString += pattern
        fuzznucNameString += "\n"
    fuzznucLabelString = "fuzznuc with motif(s)\n" + fuzznucNameString + "max. no. of mismatches: " + str(args.mismatch)
    """
    """ fimoNameString = ""
    for pattern in fimoPatternSet:
        fimoNameString += pattern
        fimoNameString += ",\n"
    fimoLabelString = "fimo with profile(s)\n" + fimoNameString + "p < " + str(args.pval)

    fimoColorValue = fimoPVal """

    #apply a transform so that the datasets do not overlap in the plot
    dx, dy = 0., 6.
    bindingSitesTF = tf.offset_copy(ax.transData, fig=fig, x=dx, y=dy, units='dots')
    chipSeqTF = tf.offset_copy(ax.transData, fig=fig, x=dx, y=-dy, units='dots')


    ax.scatter(x = bindingSitesXList, y = bindingSitesYList, transform = bindingSitesTF)
    ax.scatter(x = chipSeqX, y = chipSeqY, transform = chipSeqTF)


    ax.grid(axis="x")
    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1e6))
    ax.xaxis.set_major_formatter(ticks)
    ax.set_xlabel("genomic position / Mbp")
    ax.set_title("CTCF-binding sites for " + args.chromosomes)

    #ax.legend(loc='lower left', bbox_to_anchor= (1.5, 0), labelspacing=1.0, frameon=False, scatterpoints=3)

    """ #text box with some useful figures
    numberSitesFimoStr= "fimo fwd. " + str(len(fimoFwdList)) + "\n" + "fimo rev. " + str(len(fimoRevList))
    numberSitesFuzznucStr = "fuzznuc fwd. " + str(len(fuzznucFwdList)) + "\n" + "fuzznuc rev. " + str(len(fuzznucRevList))
    numberSitesString = "number of sites: \n" + numberSitesFimoStr + "\n" + numberSitesFuzznucStr + "\n"
    if not args.chipseq == "":
        numberSitesString += ("ChIPseq: " + str(len(chipSeqFilteredList)) + "\n")
    numberSitesString += "Reference genome: hg19"
    ax.text(0.05, 0.5, numberSitesString, transform=ax.transAxes, fontsize=12,
            verticalalignment='center') """

    """ cb = plt.colorbar(fimoScatter)
    cb.set_label("fimo p-Values")
    ticks = np.exp(np.linspace(np.log(minPFloat), np.log(maxPFloat), 5))
    cb.set_ticks(ticks)
    cb.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
    cb.ax.minorticks_off() """

    """ if not args.chipseq == "":
        cb2 = plt.colorbar(chipseqScatter)
        cb2.set_label("ChIPseq Signal Value")
        cb2.ax.minorticks_off() """

    plt.subplots_adjust(right=0.75)

    plt.show()

