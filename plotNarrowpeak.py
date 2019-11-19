import argparse
import fileParsers
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
from operator import itemgetter

parser = argparse.ArgumentParser(description="plot narrowpeak signals")
parser.add_argument("narrowpeakfile", help="narrowpeak formatted text file")
parser.add_argument("chromosome", help="chromosome to plot data for")

args = parser.parse_args()

chipSeqList = fileParsers.parseNarrowpeak(args.narrowpeakfile)
chipSeqFilteredList = [row for row in chipSeqList if row["Chromosome"] == args.chromosome]
chipSeqFilteredList = sorted(chipSeqFilteredList, key=itemgetter("Start"))
chipSeqX = [row["Start"] + row["Peak"] for row in chipSeqFilteredList]
chipSeqY = [0.5 for row in chipSeqFilteredList]
chipSeqVal = [row["SignalValue"] for row in chipSeqFilteredList]
chipSeqLength = [row["End"] - row["Start"] for row in chipSeqFilteredList]
print("number of chipseq peaks:", len(chipSeqX))
print("max. length of chipseq peaks", max(chipSeqLength))
print("min. length of chipseq peaks", min(chipSeqLength))
print("med. length of chipseq peaks {0:d}".format( int(np.median(chipSeqLength)) ) )

fig, (ax1, ax2) = plt.subplots(nrows=2,ncols=1)
p1 = ax1.scatter(x = chipSeqX, \
                y = chipSeqVal,
                c = chipSeqVal,
                cmap = plt.cm.get_cmap("viridis")
                )
p2 = ax2.hist(x=chipSeqLength, bins=100)
ax1.set_title("ChIPseq peaks " + args.chromosome)
ax2.set_title("Peak length histogram " + args.chromosome)
ax1.set_xlabel("Peak position / Mbp")
ax2.set_xlabel("Peak length / bp")
ax1.set_ylabel("Signal Value")
ax2.set_ylabel("number of peaks")
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1e6))
ax1.xaxis.set_major_formatter(ticks)
plt.tight_layout()
plt.show()
