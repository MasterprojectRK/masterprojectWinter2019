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



def strandAssigner(pChipSeqList, pPosFwdList, pPosRevList, pDmax):
    #try to assign a strand for every chip seq peak
    #given putative positions of the forward and reverse binding sites
    for element in pChipSeqList:
        #get the peak position 
        peakPos = element["Start"]
        if element["Peak"] >= 0:
            peakPos += element["Peak"] #"Peak" is given as offset from "Start", -1 if not specified
        else:
            peakPos += element["End"]
            peakPos = int( peakPos / 2) #just take the middle if no "Peak" available
        
        i = 1
        bestdiffFwd = abs(pPosFwdList[0]["Start"] - peakPos)
        while i < len(pPosFwdList):
            d = abs(pPosFwdList[i]["Start"] - peakPos)
            if d < bestdiffFwd:
                bestdiffFwd = d
            i += 1
        
        j = 1
        bestdiffRev = abs(pPosRevList[0]["Start"] - peakPos)
        while j < len(pPosRevList):
            d = abs( pPosRevList[j]["Start"] - peakPos )
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


parser = argparse.ArgumentParser(description="assign strand information to ChIP-seq data")
parser.add_argument("bindingSites", type=str, nargs="+", help="List of BED-files specifying binding sites")
parser.add_argument("chipseq", type=str, help="ChIP-seq narrowpeak file")
parser.add_argument("chromosome", type=str, help="chromosome, e.g. chr1, chr2, chrX")

args = parser.parse_args()


#extract putative binding sites from BED files
bindingSiteGlobalList = []
for bindingSiteFile in args.bindingSites:
    bindingSiteGlobalList += fileParsers.parseBedFile(bindingSiteFile)

bindingSiteFilteredList = [row for row in bindingSiteGlobalList if row["SeqName"] == args.chromosome]
bindingSiteFwdList = sorted([row for row in bindingSiteFilteredList if row["Strand"] == "+"], key=itemgetter("Start"))
bindingSiteRevList = sorted([row for row in bindingSiteFilteredList if row["Strand"] == "-"], key=itemgetter("Start"))

#extract ChIP-seq peaks from narrowpeaks file
chipSeqList = fileParsers.parseNarrowpeak(args.chipseq)
chipSeqFilteredList = [row for row in chipSeqList if row["Chromosome"] == args.chromosome]
chipSeqFilteredList = sorted(chipSeqFilteredList, key=itemgetter("Start"))
    
peakLengthsList = [row["End"] - row["Start"] for row in chipSeqFilteredList]
meanPeakLengthFloat = np.mean(peakLengthsList)
medianPeakLengthFloat = np.median(peakLengthsList)
maxPeakLengthInt = np.max(peakLengthsList)
minPeakLengthInt = np.min(peakLengthsList)
dmaxInt = int(0.5 * medianPeakLengthFloat)
numberChipseqPeaksInt = len(chipSeqFilteredList)
print("ChIPseq data:")
print("number of ChIPseq peaks in chromosome {0:s}: {1:d}".format(args.chromosome, numberChipseqPeaksInt ))
print("mean ChIPseq peak length: {0:.2f}".format(meanPeakLengthFloat))
print("median ChIPseq peak length: {0:.2f}".format(medianPeakLengthFloat))
print("max ChIPseq peak length: {0:d}".format(maxPeakLengthInt))
print("min ChIPseq peak length: {0:d}".format(minPeakLengthInt))
print("assigning peaks closer than {0:d} bp to putative binding sites".format(dmaxInt))

#try to assign a strand to the ChIP-seq peaks
strandAssigner(chipSeqFilteredList, bindingSiteFwdList, bindingSiteRevList, dmaxInt)
fileWriters.writeNarrowPeak("test.narrowPeak", [row for row in chipSeqFilteredList if row["Strand"] in ("+", "-")] )
          
chipSeqFwd = [row for row in chipSeqFilteredList if row["Strand"] == "+"]
chipSeqRev = [row for row in chipSeqFilteredList if row["Strand"] == "-"]
chipSeqUnass = [row for row in chipSeqFilteredList if row["Strand"] not in ("+", "-")]
chipSeqXFwd = [row["Start"] for row in chipSeqFwd]
chipSeqXRev = [row["Start"] for row in chipSeqRev]
chipSeqXUnass = [row["Start"] for row in chipSeqUnass]
chipSeqYFwd = ["forward" for row in chipSeqXFwd]
chipSeqYRev = ["reverse" for row in chipSeqXRev]
chipSeqYUnass = ["unassigned" for row in chipSeqXUnass]
chipSeqValFwd = [row["SignalValue"] for row in chipSeqFwd ]
chipSeqValRev = [row["SignalValue"] for row in chipSeqRev ]
chipseqValUnass = [row["SignalValue"] for row in chipSeqUnass ]
numberChipseqFwdInt = len(chipSeqXFwd)
numberChipseqRevInt = len(chipSeqXRev)
numberChipseqUnassInt = len(chipSeqXUnass)
print("assigned to fwd. strand: {0:d} ({1:.1f}%)".format(numberChipseqFwdInt, 100*numberChipseqFwdInt/numberChipseqPeaksInt))
print("assigned to rev. strand: {0:d} ({1:.1f}%)".format(numberChipseqRevInt, 100*numberChipseqRevInt/numberChipseqPeaksInt))
print("unassigned peaks: {0:d} ({1:.1f}%)".format(numberChipseqUnassInt, 100*numberChipseqUnassInt/numberChipseqPeaksInt))
print("median sig. val. unassigned: {0:0.2f}".format(np.median(chipseqValUnass)))
print("median sig. val fwd: {0:.2f}".format(np.median(chipSeqValFwd)))
print("median sig. val rev: {0:.2f}".format(np.median(chipSeqValRev)))

chipSeqX = chipSeqXFwd + chipSeqXRev + chipSeqXUnass
chipSeqY = chipSeqYFwd + chipSeqYRev + chipSeqYUnass
chipSeqVal = chipSeqValFwd + chipSeqValRev + chipseqValUnass



#Draw everything in a graph
fig, ax = plt.subplots()

#x- and y-values for the graph
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
ax.set_title("CTCF-binding sites for " + args.chromosome)

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

