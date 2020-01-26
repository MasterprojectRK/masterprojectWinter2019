import click
import pybedtools
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm
import seaborn
import os


@click.option('--taddomainfile', '-tdf', type=click.Path(exists=True), required=True, help="bed file with TAD boundaries")
@click.option('--narrowpeakfile', '-npf', type=click.Path(exists=True), required=True, help="ChIP-seq narrowPeak file")
@click.option('--outfolder','-o',type=click.Path(writable=True, exists=True), required=True, help="filename for resulting csv file")
@click.option('--chromosome', '-chr', type=str, required=True, help="chromosome to compute values for, omit 'chr' prefix")
@click.option('--chromsizefile', '-csf', type=click.Path(exists=True), required=True, help="chromosome.size file (required for binning)")
@click.option('--proteinname', '-prot', type=str, required=True, help="name of protein or Histone")
@click.option('--cellline', '-cl', type=str, required=True, help="name of cellline")
@click.option('--numberOfSamples', '-nrs', type=click.IntRange(min=1000), default=10000, help="number of sample TADs to draw")
@click.option('--resolution', '-res', type=click.IntRange(min=5000), default=20000, help="resolution used for binning, must correspond to bin size of matrix from which TADs were created")
@click.command()
def compareTadsToPeaks(taddomainfile, narrowpeakfile, outfolder, chromosome, chromsizefile, proteinname, cellline, numberofsamples, resolution):

    #load and bin the narrowpeak file
    peakDf = getProteinDataFromPeakFile(narrowpeakfile)
    peakDf, numberOfPeaks = binProteinDataFromPeaks(peakDf, chromosome, 'mean', resolution)
    
    #add bins which have no protein peaks
    chromsize = getChromSizes([chromosome],chromsizefile)[chromosome]
    maxBinInt = math.ceil(chromsize / resolution)
    proteinDf = pd.DataFrame(columns=['bin_id'])
    proteinDf['bin_id'] = list(range(0,maxBinInt))
    proteinDf.set_index('bin_id', inplace=True)
    proteinDf = proteinDf.join(peakDf, how='outer')
    proteinDf.fillna(0.0,inplace=True)
    normalizeSignalValue(proteinDf)

    #get the TAD boundaries and assign bins to them
    tadDf = getTadDataFromBedFile(taddomainfile)
    tadDf = binTadsFromTadData(tadDf, chromosome, resolution)

    #get all proteins which are at TAD starts and end
    tadDf['startProteins'] = list(proteinDf['signalValue'][tadDf['bin_id_start']])
    tadDf['endProteins'] = list(proteinDf['signalValue'][tadDf['bin_id_end']])
    #get all proteins which are at TAD starts and TAD ends +/- nr_bins
    nr_bins = 3
    for i in range(1,nr_bins):
        maxStartBinId = int(tadDf.bin_id_start.max())
        maxEndBinId = int(tadDf.bin_id_end.max())
        bin_ids_minus = (tadDf['bin_id_start'] - i).clip(lower=0, upper=maxStartBinId)
        bin_ids_plus = (tadDf['bin_id_start'] + i).clip(lower=0, upper=maxStartBinId)
        tadDf['startProteins-' + str(i)] = list(proteinDf['signalValue'][bin_ids_minus])
        tadDf['startProteins+' + str(i)] = list(proteinDf['signalValue'][bin_ids_plus])
        bin_ids_minus = (tadDf['bin_id_end'] - i).clip(lower=0, upper=maxEndBinId)
        bin_ids_plus = (tadDf['bin_id_end'] + i).clip(lower=0, upper=maxEndBinId)
        tadDf['endProteins-' + str(i)] = list(proteinDf['signalValue'][bin_ids_minus])
        tadDf['endProteins+' + str(i)] = list(proteinDf['signalValue'][bin_ids_plus])
    
    #add window proteins
    tadDf['distance'] = (tadDf['bin_id_end'] - tadDf['bin_id_start']).astype('uint32')
    winSize = tadDf.distance.max() + 1
    windowDf = buildWindowDataset(proteinDf, 'signalValue', winSize, 'mean')
    #get the window proteins into an array and slice it to get all values at once 
    #there might be a more efficient way using pandas
    distWindowArr = windowDf.to_numpy()
    slice1 = list(tadDf['bin_id_end'])
    slice2 = list(tadDf['distance'])
    slice3 = (slice1, slice2)
    windowProteins = np.array(distWindowArr[slice3])
    tadDf['windowProteins'] = np.float32(windowProteins)

    #compute mean/min/max signalValue at boundaries 
    meanSignalValueAtStart = tadDf.startProteins.mean()
    medianSignalValueAtStart = tadDf.startProteins.median()
    maxSignalValueAtStart = tadDf.startProteins.max()
    minSignalValueAtStart = tadDf.startProteins.min()
    meanSignalValueAtEnd = tadDf.endProteins.mean()
    medianSignalValueAtEnd = tadDf.endProteins.median()
    maxSignalValueAtEnd = tadDf.endProteins.max()
    minSignalValueAtEnd = tadDf.endProteins.min()
    meanSignalValueAtBoundaries = pd.concat([tadDf['startProteins'], tadDf['endProteins']],ignore_index=True,sort=False).mean()
    
    #compute mean/min/max signal value for whole chromosome
    meanSignalValueOverall = proteinDf.signalValue.mean()
    medianSignalValueOverall = proteinDf.signalValue.median()
    maxSignalValueOverall = proteinDf.signalValue.max()
    minSignalValueOverall = proteinDf.signalValue.min()

    #compute mean/min/max window protein value for TADs
    meanWindowValueAtTads = tadDf.windowProteins.mean()
    medianWindowValueAtTads = tadDf.windowProteins.median()
    maxWindowValueAtTads = tadDf.windowProteins.max()
    minWindowValueAtTads = tadDf.windowProteins.min()

    #get all TAD bins which have zero signal value at start or end
    zeroAtStartMask = tadDf['startProteins'] == 0.0
    zeroAtStartCount = tadDf[zeroAtStartMask].shape[0]
    zeroAtEndMask = tadDf['endProteins'] == 0.0
    zeroAtEndCount = tadDf[zeroAtEndMask].shape[0]
    zeroAtStartEndCount = tadDf[zeroAtStartMask & zeroAtEndMask].shape[0]
    numberOfTads = tadDf.shape[0]

    #prepare data for plotting TAD distances
    tadBinStepsize = 3
    tadsWithProteinsAtBoundariesDf = tadDf[~zeroAtStartMask & ~zeroAtEndMask]
    distBins = np.arange(0, tadDf.distance.max(), tadBinStepsize)
    binnedTadsWithBoundariesList = list(pd.cut(tadsWithProteinsAtBoundariesDf['distance'], bins=distBins).value_counts(sort=False))
    binnedTadsList = list(pd.cut(tadDf['distance'], bins=distBins).value_counts(sort=False))
    
    #get all TAD bins which have zero signal value at start +/- nr_bins bins, end +/- 2 bins
    for i in range(1,nr_bins):
        zeroPlusMask = tadDf['startProteins+' + str(i)] == 0.0
        zeroMinusMask = tadDf['startProteins-' + str(i)] == 0.0
        startPlusMinusMask = zeroAtStartMask & zeroPlusMask & zeroMinusMask
        zeroPlusMask = tadDf['endProteins+' + str(i)] == 0.0
        zeroMinusMask = tadDf['startProteins-' + str(i)] == 0.0
        endPlusMinusMask = zeroAtEndMask & zeroPlusMask & zeroMinusMask
    zeroAtStartPlusMinusCount = tadDf[startPlusMinusMask].shape[0]
    zeroAtEndPlusMinusCount = tadDf[endPlusMinusMask].shape[0]
    zeroAtStartEndPlusMinusCount = tadDf[startPlusMinusMask & endPlusMinusMask].shape[0]

    #print some figures
    writeout = [cellline + " " + proteinname]
    
    msg = "{0:d} TADs for chr{1:s} found in input {2:s}"
    msg = msg.format(tadDf.shape[0],chromosome, os.path.basename(taddomainfile))
    writeout.append(msg)
    print(msg)
    
    msg = "{0:d} protein peaks for chr{1:s} found in input {2:s}"
    msg = msg.format(numberOfPeaks,chromosome, os.path.basename(narrowpeakfile))
    writeout.append(msg)
    print(msg)
    
    print("after binning:")
    msg = "chr{0:s} comprises {1:d} bins, {2:d} of which have protein peaks"
    msg = msg.format(chromosome, maxBinInt, len(proteinDf.signalValue.to_numpy().nonzero()[0]) )
    writeout.append(msg)
    print(msg)
    
    msg = "mean signal value: {0:.3f} (min {1:.3f}, max {2:.3f}, med {3:.3f})"
    msg = msg.format(meanSignalValueOverall, minSignalValueOverall, maxSignalValueOverall, medianSignalValueOverall)
    writeout.append(msg)
    print(msg)
    
    msg = "mean signal value at TAD boundaries: {0:.3f}"
    msg = msg.format(meanSignalValueAtBoundaries)
    writeout.append(msg)
    print(msg)
    
    msg = "mean signal value at TAD starts: {0:.3f} (min {1:.3f}, max {2:.3f}, med {3:.3f})"
    msg = msg.format(meanSignalValueAtStart, minSignalValueAtStart, maxSignalValueAtStart, medianSignalValueAtStart)
    writeout.append(msg)
    print(msg)
    
    msg = "mean signal value at TAD ends: {0:.3f} (min {1:.3f}, max {2:.3f}, med {3:.3f})"
    msg = msg.format(meanSignalValueAtEnd, minSignalValueAtEnd, maxSignalValueAtEnd, medianSignalValueAtEnd)
    writeout.append(msg)
    print(msg)
    
    msg = "mean signal value for TAD window proteins: {0:.3f} (min {1:.3f}, max {2:.3f}, med {3:.3f})"
    msg = msg.format(meanWindowValueAtTads, minWindowValueAtTads, maxWindowValueAtTads, medianWindowValueAtTads)
    writeout.append(msg)
    print(msg)
    
    msg = "{0:d} of {1:d} TADs have no protein at start bin"
    msg = msg.format(zeroAtStartCount, numberOfTads)
    writeout.append(msg)
    print(msg)
    
    msg = "{0:d} of {1:d} TADs have no protein at start bin +/- {2:d} bins"
    msg = msg.format(zeroAtStartPlusMinusCount, numberOfTads, nr_bins-1)
    writeout.append(msg)
    print(msg)
    
    msg = "{0:d} of {1:d} TADs have no protein at end bin"
    msg = msg.format(zeroAtEndCount, numberOfTads)
    writeout.append(msg)
    print(msg)
    
    msg = "{0:d} of {1:d} TADs have no protein at end bin +/- {2:d} bins"
    msg = msg.format(zeroAtEndPlusMinusCount, numberOfTads, nr_bins-1)
    writeout.append(msg)
    print(msg)
    
    msg = "{0:d} of {1:d} TADs have no protein at start and end bin"
    msg = msg.format(zeroAtStartEndCount, numberOfTads)
    writeout.append(msg)
    print(msg)
    
    msg = "{0:d} of {1:d} TADs have no protein at start and end bin +/- {2:d} bins"
    msg = msg.format(zeroAtStartEndPlusMinusCount, numberOfTads, nr_bins-1)
    writeout.append(msg)
    print(msg)

    #also write the text to a file
    resFileName = outfolder + cellline + "_" + proteinname + "_stats.txt"
    with open(resFileName, "w") as outfile:
        outfile.writelines("\n".join(writeout))


    #build pdfs for mean at TAD boundaries, number of TADs with no protein
    tadBoundaryProtArr = []
    tadBoundaryWithoutProtArr = []
    windowProtArr = []
    tadBoundaryMeanProtArr = []
    windowMeanProtArr = []
    maxWinSize = proteinDf.shape[0]
    indexList = [x for x in range(0, maxWinSize)]
    windowDf = buildWindowDataset(proteinDf, 'signalValue', maxWinSize, 'mean')
    distWindowArr = windowDf.to_numpy()
    for i in tqdm(range(0,numberofsamples), miniters=1000, desc="drawing indices for random TADs and computing proteins"):
        #draw random TAD indices with equal probability and without replacement
        #then get protein data for these indices
        tadIndices = np.random.choice(indexList, numberOfTads + 1, replace=False)
        tadSignalValues = list(proteinDf.loc[tadIndices,'signalValue'])
        tadBoundaryProtArr.extend(tadSignalValues)
        tadBoundaryMeanProtArr.append(np.mean(tadSignalValues))
        tadBoundaryWithoutProtArr.append(tadSignalValues.count(0.0))
        tadIndices = sorted(tadIndices)
        tadStartIndices = tadIndices[0:-2]
        tadEndIndices = tadIndices[1:-1]
        distanceList = list(np.subtract(tadEndIndices, tadStartIndices))
        #slice the window array and get the mean out of all the windows        
        windowSlice = (tadEndIndices, distanceList)
        windowProteins = np.array(distWindowArr[windowSlice])
        windowProtArr.extend(windowProteins)
        windowMeanProtArr.append(np.mean(windowProteins))

    #plots
    fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, constrained_layout=True, figsize = (15,7))
    proteinDf.signalValue.plot(kind='hist', density=True, bins=100, ax=ax1)
    seaborn.distplot(list(tadDf.startProteins), norm_hist=True,  bins=100, ax=ax2, hist=True, kde=False, label="signalValue at starts")
    seaborn.distplot(list(tadDf.endProteins), norm_hist=True,  bins=100, ax=ax2, hist=True, kde=False, label="signalValue at ends")
    windowProts = list(tadDf.windowProteins)
    seaborn.distplot(windowProts, norm_hist=True, ax=ax4, hist=True, kde=False, label="window proteins")
    pieSizes = [numberOfTads - zeroAtStartCount - zeroAtEndCount + zeroAtStartEndCount, 
            zeroAtStartCount - zeroAtStartEndCount,
            zeroAtEndCount - zeroAtStartEndCount,
            zeroAtStartEndCount]
    pieLabels = ["with start- and end proteins", 
                "without start proteins",
                "without end proteins",
                "without start- and end proteins"]
    ax3.pie(pieSizes, labels=pieLabels, startangle=-45 , autopct='%1.1f%%')
    ax3.axis('equal')
    ax1.set_title("probability density for signal value over whole chromosome")
    ax2.set_title("probability density for signal value at TAD boundaries")
    ax4.set_title("probability density for TAD window proteins")
    ax3.set_title("TADs with and without proteins, total = " + str(numberOfTads) + " TADs")
    ax1.set_ylabel("Frequency")
    ax2.set_ylabel("Frequency")
    ax4.set_ylabel("Frequency")
    ax1.set_xlabel("signal value")
    ax2.set_xlabel("signal value")
    ax4.set_xlabel("window signal value")
    ax1.set_xlim([-0.05, 1.05])
    ax2.set_xlim([-0.05, 1.05])
    ax4.set_xlim([-0.05, 1.05])
    ax1.set_ylim([0,10])
    ax2.set_ylim([0,10])
    ax4.set_ylim([0,20])
    ax1.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax2.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax4.set_yticks([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20])
    ax1.set_yticks([0, 2, 4, 6, 8, 10])
    ax2.set_yticks([0, 2, 4, 6, 8, 10])
    
    
    
    fig1.suptitle("Metrics for proteins and actual TADs " + cellline + " chr" + chromosome + " " + proteinname)
    #fig1.show()
    fig1.savefig(outfolder + cellline + "_" + proteinname + "_proteinAndTadStats.png")

    fig2 = plt.figure(constrained_layout=True, figsize=(15,7))
    fig2Grid = fig2.add_gridspec(2,2)
    ax5 = fig2.add_subplot(fig2Grid[0,1])
    ax6 = fig2.add_subplot(fig2Grid[:,0])
    ax7 = fig2.add_subplot(fig2Grid[1,1])
    seaborn.distplot(tadBoundaryProtArr, bins=100, hist=True, norm_hist=True,ax=ax5, kde=False)
    ax5.axvline(x=meanSignalValueAtBoundaries, color='red')
    ax5.axvline(x=np.mean(tadBoundaryProtArr), color='blue')
    seaborn.distplot(tadBoundaryWithoutProtArr, bins=max(tadBoundaryWithoutProtArr)-min(tadBoundaryWithoutProtArr), hist=True, norm_hist=True, ax=ax6, kde=True)
    ax6.axvline(x=(zeroAtStartCount + zeroAtEndCount - zeroAtStartEndCount), color='red')
    ax6.axvline(x=np.mean(tadBoundaryWithoutProtArr), color='blue' )
    seaborn.distplot(windowProtArr, hist=True, norm_hist=True, kde=False, ax=ax7)
    ax7.axvline(x=meanWindowValueAtTads,color='red')
    ax7.axvline(x=np.mean(windowProtArr), color='blue')
    ax5.set_title("probability density for signal value at TAD boundaries")
    ax6.set_title("probability density for number of TAD boundaries without protein")
    ax7.set_title("probability density for TAD window signal value")
    ax5.set_ylabel("Frequency")
    ax6.set_ylabel("Frequency")
    ax7.set_ylabel("Frequency")
    ax5.set_xlabel("signal value")
    ax6.set_xlabel("number of TAD boundaries w/o protein")
    ax7.set_xlabel("window signal value")
    ax5.set_xlim([-0.05, 1.05])
    ax7.set_xlim([-0.05, 1.05])
    ax5.set_ylim([0, 10])
    ax7.set_ylim([0,20])
    ax5.set_yticks([0, 2, 4, 6, 8, 10])
    ax7.set_yticks([0, 2.5, 5, 7.5, 10, 12.5, 15, 17.5, 20])
    ax5.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax7.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    fig2.suptitle("Values from " + str(numberofsamples) + " sampled, random TADs vs. actuals, " + cellline + " chr" + chromosome + " " + proteinname)
    #fig2.show()
    fig2.savefig(outfolder + cellline + "_" + proteinname + "_sampleVsReal.png")

    fig3, (f3ax1, f3ax2) = plt.subplots(nrows=1, ncols=2, constrained_layout=True, figsize = (15,7))
    seaborn.distplot(tadBoundaryMeanProtArr, hist=True, kde=True, norm_hist=True, ax=f3ax1)
    seaborn.distplot(windowMeanProtArr, hist=True,kde=True,norm_hist=True,ax=f3ax2)
    f3ax1.axvline(x=np.mean(tadBoundaryMeanProtArr),color='blue')
    f3ax1.axvline(x=meanSignalValueAtBoundaries, color='red')
    f3ax2.axvline(x=np.mean(windowMeanProtArr),color='blue')
    f3ax2.axvline(x=meanWindowValueAtTads, color='red')
    f3ax1.set_xlim([0, 0.2])
    f3ax2.set_xlim([0, 0.2])
    f3ax1.set_ylim([0, 80])
    f3ax2.set_ylim([0, 120])
    f3ax1.set_xticks([0, 0.05, 0.1, 0.15, 0.2])
    f3ax2.set_xticks([0, 0.05, 0.1, 0.15, 0.2])
    f3ax1.set_yticks(np.linspace(0,80,num=7))
    f3ax2.set_yticks(np.linspace(0,120,num=11))
    f3ax1.set_ylabel("Frequency")
    f3ax2.set_ylabel("Frequency")
    f3ax1.set_xlabel("mean signal value")
    f3ax2.set_xlabel("mean window signal value")
    f3ax1.set_title("probability density for mean signal value at TAD boundaries")
    f3ax2.set_title("probability density for mean TAD window signal value")
    fig3.suptitle("Values from " + str(numberofsamples) + " sampled, random TADs vs. actuals, " + cellline + " chr" + chromosome + " " + proteinname)
    fig3.savefig(outfolder + cellline + "_" + proteinname + "_mean_sampleVsReal.png")

    if not tadsWithProteinsAtBoundariesDf.empty:
        fig4, (f4ax1,f4ax2) = plt.subplots(nrows=2, ncols=1, constrained_layout=True)
        xValsList = [resolution * x for x in distBins][0:-1]
        f4ax1.bar(x=xValsList, height=binnedTadsWithBoundariesList, width=tadBinStepsize*resolution)
        f4ax2.bar(x=xValsList, height=binnedTadsList,width=tadBinStepsize*resolution)
        fig4.suptitle("TADs by distance " + cellline + " chr" + chromosome + " " + proteinname)
        f4ax1.set_title("only TADs with proteins at boundaries")
        f4ax2.set_title("All TADs")
        f4ax1.set_xlim(f4ax2.get_xlim())
        f4ax1.set_ylim(f4ax2.get_ylim())
        f4ax1.set_xticklabels([x/1000000 for x in f4ax1.get_xticks()])
        f4ax2.set_xticklabels([x/1000000 for x in f4ax2.get_xticks()])
        f4ax1.set_xlabel("distance / Mbp")
        f4ax2.set_xlabel("distance / Mbp")
        f4ax1.set_ylabel("Frequency")
        f4ax2.set_ylabel("Frequency")
        fig4.savefig(outfolder + cellline + "_" + proteinname + "_TADsByDistance.png")

def getTadDataFromBedFile(pBedFilePath):
    tadData = None
    try:
        tadData = pd.read_table(pBedFilePath, header=None, names=['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb'], index_col=False)
    except Exception as exc:
        msg = "could not parse bed file {0:s} \n"
        msg = "probably no valid TAD-domains file"
        print("exception:", exc)
        raise ValueError(msg.format(pBedFilePath))
    tadData.drop(columns=['thickStart', 'thickEnd', 'itemRgb'], inplace=True)
    return tadData

def binTadsFromTadData(pTadDf, pChrom, pResolution):
    mask = pTadDf['chrom'] == "chr" + str(pChrom)
    binnedDf = pTadDf[mask]
    binnedDf['bin_id_start'] = (binnedDf['chromStart']/pResolution).astype('uint32')  
    binnedDf['bin_id_end'] = (binnedDf['chromEnd']/pResolution).astype('uint32')  
    binnedDf.drop(columns=['chrom', 'chromStart', 'chromEnd', 'score' , 'strand'], inplace=True)
    return binnedDf

def getProteinDataFromPeakFile(pPeakFilePath):
    try:        
        bedToolFile = pybedtools.BedTool(pPeakFilePath)
        malformedFeatures = [features for features in bedToolFile if len(features.fields) not in [9,10]]
    except:
        msg = "could not parse protein peak file {0:s} \n"
        msg += "probably no valid narrowPeak or broadPeak file"
        raise ValueError(msg.format(pPeakFilePath))
    if malformedFeatures:
            msg = "protein file {0:s} seems to be an invalid narrow- or broadPeak file\n"
            msg += "there are rows with more than 10 or less than 9 columns"
            raise ValueError(msg.format(pPeakFilePath))
    ### compute min and max for the normalization
    protData = bedToolFile.to_dataframe()
    columnNames = ['chrom', 'chromStart', 'chromEnd', 'name', 'score',
                            'strand', 'signalValue', 'pValue', 'qValue', 'peak']
    if protData.shape[1] == 10: #narrowPeak format
        protData.columns = columnNames
        mask = protData['peak'] == -1 
        if mask.any(): #no peak summit called
            protData[mask]['peak'] = ( (protData[mask]['chromEnd']-protData[mask]['chromStart'])/2 ).astype('uint32')
    elif protData.shape[1] == 9: #broadPeak format, generally without peak column in the data
        protData.columns = columnNames[0:9]
        protData['peak'] = ((protData['chromEnd'] - protData['chromStart']) / 2).astype('uint32')
    return protData

def binProteinDataFromPeaks(pProteinDf, pChrom, pMergeParam, pResolution):    
    mask = pProteinDf['chrom'] == "chr" + str(pChrom)
    proteinDf = pProteinDf[mask].copy()
    nrOfPeaks = proteinDf.shape[0]
    proteinDf.drop(columns=['chrom','name', 'score', 'strand'], inplace=True)
    proteinDf['bin_id'] = ((proteinDf['chromStart'] + proteinDf['peak'])/pResolution).astype('uint32')
    #print(proteinDf.head(10))
    if pMergeParam == 'max':
        binnedDf = proteinDf.groupby('bin_id')[['signalValue']].max()
    else:
        binnedDf = proteinDf.groupby('bin_id')[['signalValue']].mean()
    return binnedDf, nrOfPeaks

def getChromSizes(pChromNameList, pChromSizeFile):
    chromSizeDict = dict()
    try:
        chromSizeDf = pd.read_csv(pChromSizeFile,names=['chrom', 'size'],header=None,sep='\t')
    except:
        msg = "Error: could not parse chrom.sizes file {0:s}\n".format(pChromSizeFile)
        msg += "Maybe wrong format, not tab-separated etc.?"
        raise Exception(msg)
    for chromName in pChromNameList:
        sizeMask = chromSizeDf['chrom'] == "chr" + str(chromName)
        if not sizeMask.any():
            msg = "no entry for chromosome {0:s} in chrom.sizes file".format(chromName)
            raise ValueError(msg)
        elif chromSizeDf[sizeMask].shape[0] > 1:
            msg = "multiple entries for chromosome {0:s} in chrom.sizes file".format(chromName)
            raise ValueError(msg)
        else:
            try:
                chromSizeDict[chromName] = int(chromSizeDf[sizeMask]['size'])
            except ValueError:
                msg = "entry for chromosome {0:s} in chrom.sizes is not an integer".format(chromName)
    return chromSizeDict


def buildWindowDataset(pProteinsDf, pProteinNr, pWindowSize, pWindowOperation):
    df = pd.DataFrame()
    proteinIndex = str(pProteinNr)
    for winSize in range(pWindowSize):
        if pWindowOperation == "max":
            windowColumn = pProteinsDf[proteinIndex].rolling(window=winSize+1).max()
        elif pWindowOperation == "sum":
            windowColumn = pProteinsDf[proteinIndex].rolling(window=winSize+1).sum()
        else:
            windowColumn = pProteinsDf[proteinIndex].rolling(window=winSize+1).mean()
        df[str(winSize)] = windowColumn.round(3).astype('float32')
    return df

def normalizeSignalValue(pProteinDf):
    #inplace zero-to-one normalization of signal value
    try:
        maxSignalValue = pProteinDf.signalValue.max()
        minSignalValue = pProteinDf.signalValue.min()
    except:
        msg = "can't normalize Dataframe without signalValue"
        raise Warning(msg)
    if minSignalValue == maxSignalValue:
        msg = "no variance in protein data file"
        pProteinDf['signalValue'] = 0.0
        raise Warning(msg)
    else:
        diff = maxSignalValue - minSignalValue
        pProteinDf['signalValue'] = ((pProteinDf['signalValue'] - minSignalValue) / diff).astype('float32')


if __name__ == "__main__":
    compareTadsToPeaks()