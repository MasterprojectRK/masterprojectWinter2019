import argparse
import csv
import matplotlib.pyplot as plt
import matplotlib.transforms as tf
import matplotlib.colors as colors
import numpy as np
import fileParsers
import matplotlib.ticker as ticker
from operator import itemgetter



def strandAssigner(pChipSeqList, pPosFwdList, pPosRevList, pDmax):
    #try to assign a strand for every chip seq peak
    #given putative positions of the forward and reverse binding sites
    for element in pChipSeqList:
        #get the peak position 
        peakPos = element["Start"]
        if element["Peak"] >= 0:
            peakPos += element["Peak"] #"Peak" is offset from start
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



parser = argparse.ArgumentParser(description="display motif positions found by EMBOSS fuzznuc vs Meme FIMO")
parser.add_argument("fimoTsvString", type=str, help="path to FIMO tsv output file")
parser.add_argument("fuzznucResFileString", type=str, help="path to fuzznuc output file in excel format")
parser.add_argument("chromosome", type=str, help="chromosome, e.g. chr1, chr2, chrX")
parser.add_argument("--chipseq", default = "", type=str, help="path to ChIP-seq narrowpeak file")
parser.add_argument("--pval", type=float, default=1.0, help="max. false positive rate (0...1]")
parser.add_argument("--mismatch", type=int, default = 0, help="allowed number of mismatches for fuzznuc")

args = parser.parse_args()


#extract start positions of deteted motif sites from fuzznuc results file.
#result file must be in fuzznuc "excel" output format
#filter for given chromosome and number of mismatches
fuzznucDataList = fileParsers.parseFuzznuc(args.fuzznucResFileString)
fuzznucFilteredList = [row for row in fuzznucDataList if row["SeqName"] == args.chromosome and row["Mismatch"] <= args.mismatch]
fuzznucFilteredList = sorted(fuzznucFilteredList, key=itemgetter("Start"))
fuzznucFwdList = [row for row in fuzznucFilteredList if row["Strand"] == "+"]
fuzznucRevList = [row for row in fuzznucFilteredList if row["Strand"] == "-"]
fuzznucPatternSet = set()
for row in fuzznucFilteredList:
    fuzznucPatternSet.add(row["Motif"])

#extract start positions of deteced motif sites from fimo results file
#result file must be in fimo TSV, tab separated value, format
#filter for given chromosome and q-value
fimoDataList = fileParsers.parseFimo(args.fimoTsvString)
fimoFilteredList = [row for row in fimoDataList if  row["SeqName"] == args.chromosome and row["PVal"] < args.pval ]
fimoFilteredList = sorted(fimoFilteredList, key=itemgetter("Start"))
fimoFwdList = [row for row in fimoFilteredList if row["Strand"] == "+"]
fimoRevList = [row for row in fimoFilteredList if row["Strand"] == "-"]
fimoPatternSet = set()
for row in fimoFilteredList:
    fimoPatternSet.add(row["Motif"])


#extract ChIP-seq peaks from narrowpeaks file
#and try to assign a strand
chipSeqList = chipSeqFilteredList = chipSeqX = chipSeqY = chipSeqVal = []
if not args.chipseq == "":
    chipSeqList = fileParsers.parseNarrowpeak(args.chipseq)
    chipSeqFilteredList = [row for row in chipSeqList if row["Chromosome"] == args.chromosome]
    chipSeqFilteredList = sorted(chipSeqFilteredList, key=itemgetter("Start"))
    
    peakLengthsList = [row["End"] - row["Start"] for row in chipSeqFilteredList]
    meanPeakLengthFloat = np.mean(peakLengthsList)
    medianPeakLengthFloat = np.median(peakLengthsList)
    maxPeakLengthInt = np.max(peakLengthsList)
    minPeakLengthInt = np.min(peakLengthsList)
    dmaxInt = int(2* meanPeakLengthFloat)
    print("ChIPseq data:")
    print("mean ChIPseq peak length: {0:.2f}".format(meanPeakLengthFloat))
    print("median ChIPseq peak length: {0:.2f}".format(medianPeakLengthFloat))
    print("max ChIPseq peak length: {0:d}".format(maxPeakLengthInt))
    print("min ChIPseq peak length: {0:d}".format(minPeakLengthInt))
    print("assigning peaks closer than {0:d} bp to putative binding sites".format(dmaxInt))
    strandAssigner(chipSeqFilteredList, fimoFwdList + fuzznucFwdList, fimoRevList + fuzznucRevList, dmaxInt)

          
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
    print("assigned to fwd. strand: ", len(chipSeqXFwd))
    print("assigned to rev. strand: ", len(chipSeqXRev))
    print("unassigned peaks: ", len(chipSeqXUnass))
    print("median sig. val. unassigned: {0:0.2f}".format(np.median(chipseqValUnass)))
    print("median sig. val fwd: {0:.2f}".format(np.median(chipSeqValFwd)))
    print("median sig. val rev: {0:.2f}".format(np.median(chipSeqValRev)))

    chipSeqX = chipSeqXFwd + chipSeqXRev + chipSeqXUnass
    chipSeqY = chipSeqYFwd + chipSeqYRev + chipSeqYUnass
    chipSeqVal = chipSeqValFwd + chipSeqValRev + chipseqValUnass



#Draw everything in a graph
fig, ax = plt.subplots()

#x- and y-values for the graph
fuzznucXfwd = [row["Start"] for row in fuzznucFwdList]
fuzznucYfwd = ["forward" for row in fuzznucFwdList]

fuzznucXrev = [row["Start"] for row in fuzznucRevList]
fuzznucYrev = ["reverse" for row in fuzznucRevList]

fuzznucX = fuzznucXfwd + fuzznucXrev
fuzznucY = fuzznucYfwd + fuzznucYrev

fimoXfwd = [row["Start"] for row in fimoFwdList]
fimoYfwd = ["forward" for row in fimoFwdList]
fimoPValfwd = [row["PVal"] for row in fimoFwdList]

fimoXrev = [row["Start"] for row in fimoRevList]
fimoYrev = ["reverse" for row in fimoRevList]
fimoPValrev = [row["PVal"] for row in fimoRevList]

fimoX = fimoXfwd + fimoXrev
fimoY = fimoYfwd + fimoYrev
fimoPVal = fimoPValfwd + fimoPValrev

#max and min values for the graph. Use max over all chromosomes, but obey pval cutoff
print("fimo search statistics:")
maxQFloat = max([row["QVal"] for row in fimoDataList if row["PVal"] < args.pval])
minQFloat = min([row["QVal"] for row in fimoDataList if row["PVal"] < args.pval])
print("max. q-value: ", maxQFloat)
print("min. q-value: ", minQFloat)
maxPFloat = max([row["PVal"] for row in fimoDataList if row["PVal"] < args.pval])
minPFloat = min([row["PVal"] for row in fimoDataList if row["PVal"] < args.pval])
print("max. p-value: ", maxPFloat)
print("min. p-value: ", minPFloat)

#apply a transform so that the datasets do not overlap in the plot
dx, dy = 0., 10.
fimoTF = tf.offset_copy(ax.transData, fig=fig, x=dx, y=dy, units='dots')
fuzznucTF = tf.offset_copy(ax.transData, fig=fig, x=dx, y=-dy, units='dots')


#labels for the legend, include the names of the patterns
fuzznucNameString = ""
for pattern in fuzznucPatternSet:
    fuzznucNameString += pattern
    fuzznucNameString += "\n"
fuzznucLabelString = "fuzznuc with motif(s)\n" + fuzznucNameString + "max. no. of mismatches: " + str(args.mismatch)

fimoNameString = ""
for pattern in fimoPatternSet:
    fimoNameString += pattern
    fimoNameString += ",\n"
fimoLabelString = "fimo with profile(s)\n" + fimoNameString + "p < " + str(args.pval)

fimoColorValue = fimoPVal

fuzznucScatter = ax.scatter(x = fuzznucX, y = fuzznucY, marker="+", color="blue", alpha=0.7, label=fuzznucLabelString, transform=fuzznucTF )
fimoScatter = ax.scatter(x = fimoX, \
                       y = fimoY, \
                       marker="x", 
                       c=fimoColorValue, 
                       label=fimoLabelString, 
                       transform=fimoTF, 
                       cmap = plt.cm.get_cmap("viridis"), 
                       vmin=minPFloat, vmax=maxPFloat, 
                       norm=colors.LogNorm() 
                       )
if not args.chipseq == "":
    chipseqScatter = ax.scatter(x = chipSeqX, \
                                y = chipSeqY, \
                                marker = "o", \
                                c = chipSeqVal, \
                                label = "ChIP-seq data", \
                                cmap =  plt.cm.get_cmap("RdYlGn"),  
                                )

ax.grid(axis="x")
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1e6))
ax.xaxis.set_major_formatter(ticks)
ax.set_xlabel("genomic position / Mbp")
ax.set_title("CTCF-binding sites for " + args.chromosome)

ax.legend(loc='lower left', bbox_to_anchor= (1.5, 0), labelspacing=1.0, frameon=False, scatterpoints=3)

#text box with some useful figures
numberSitesFimoStr= "fimo fwd. " + str(len(fimoFwdList)) + "\n" + "fimo rev. " + str(len(fimoRevList))
numberSitesFuzznucStr = "fuzznuc fwd. " + str(len(fuzznucFwdList)) + "\n" + "fuzznuc rev. " + str(len(fuzznucRevList))
numberSitesString = "number of sites: \n" + numberSitesFimoStr + "\n" + numberSitesFuzznucStr + "\n"
if not args.chipseq == "":
    numberSitesString += ("ChIPseq: " + str(len(chipSeqFilteredList)) + "\n")
numberSitesString += "Reference genome: hg19"
ax.text(0.05, 0.5, numberSitesString, transform=ax.transAxes, fontsize=12,
        verticalalignment='center')

cb = plt.colorbar(fimoScatter)
cb.set_label("fimo p-Values")
ticks = np.exp(np.linspace(np.log(minPFloat), np.log(maxPFloat), 5))
cb.set_ticks(ticks)
cb.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2e'))
cb.ax.minorticks_off()

if not args.chipseq == "":
    cb2 = plt.colorbar(chipseqScatter)
    cb2.set_label("ChIPseq Signal Value")
    cb2.ax.minorticks_off()

plt.subplots_adjust(right=0.75)

plt.show()

