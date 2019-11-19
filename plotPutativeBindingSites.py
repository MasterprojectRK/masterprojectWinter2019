import argparse
import fileParsers
import matplotlib.pyplot as plt
from matplotlib import ticker

#read in and plot BED file w.r.t given chromosome


parser = argparse.ArgumentParser(description="plot putative binding sites from bed files")
parser.add_argument("bedfile", type=str, help="BED formatted text file")
parser.add_argument("chromosome", type=str, help="chromosome to plot data for")

args = parser.parse_args()

bindingSiteList = fileParsers.parseBedFile(args.bedfile)
bindingSiteFilteredList = [row for row in bindingSiteList if row["SeqName"] == args.chromosome]

bindingSiteFwdList = [row for row in bindingSiteFilteredList if row["Strand"] == "+"]
bindingSiteRevList = [row for row in bindingSiteFilteredList if row["Strand"] == "-"] 

bindingSiteFwdX = [int((row["Start"] + row["End"]) / 2) for row in bindingSiteFwdList]
bindingSiteFwdY = [row["Score"] for row in bindingSiteFwdList]

bindingSiteRevX = [int((row["Start"] + row["End"]) / 2) for row in bindingSiteRevList]
bindingSiteRevY = [-row["Score"] for row in bindingSiteRevList]

fig, ax = plt.subplots()
ax.scatter(x = bindingSiteFwdX, y=bindingSiteFwdY)
ax.scatter(x = bindingSiteRevX, y=bindingSiteRevY)

ax.set_title("putative binding sites for " + args.chromosome)
ax.set_xlabel("Peak position / Mbp")
ax.set_ylabel("Score (+forward / -reverse)")
ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x/1e6))
ax.xaxis.set_major_formatter(ticks)

plt.show()