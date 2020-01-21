import click
import pybedtools
import pandas as pd
import math
import numpy as np
from matplotlib import pyplot as plt


@click.option('--tadboundaryfile', '-tbf', type=click.Path(exists=True), required=True, help="bed file with TAD boundaries")
@click.option('--narrowpeakfile', '-npf', type=click.Path(exists=True), required=True, help="ChIP-seq narrowPeak file")
@click.option('--outfile','-o',type=click.Path(writable=True), required=True, help="filename for resulting csv file")
@click.option('--chromosome', '-chr', type=str, required=True, help="chromosome to compute values for, omit 'chr' prefix")
@click.option('--chromsizefile', '-csf', type=click.Path(exists=True), required=True, help="chromosome.size file (required for binning)")
@click.command()
def compareTadsToPeaks(tadboundaryfile, narrowpeakfile, outfile, chromosome, chromsizefile):
    resolution = 20000
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

    #get the TAD boundaries and assign bins to them
    tadDf = getTadDataFromBedFile(tadboundaryfile)
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

    #get all TAD bins which have zero signal value at start or end
    zeroAtStartMask = tadDf['startProteins'] == 0.0
    zeroAtStartCount = tadDf[zeroAtStartMask].shape[0]
    zeroAtEndMask = tadDf['endProteins'] == 0.0
    zeroAtEndCount = tadDf[zeroAtEndMask].shape[0]
    zeroAtStartEndCount = tadDf[zeroAtStartMask & zeroAtEndMask].shape[0]
    tadCount = tadDf.shape[0]
    
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
    msg = "{0:d} TADs for chr{1:s} found in input"
    print(msg.format(tadDf.shape[0],chromosome))
    msg = "{0:d} protein peaks for chr{1:s} found in input"
    print(msg.format(numberOfPeaks,chromosome))
    print("after binning:")
    msg = "chr{0:s} comprises {1:d} bins, {2:d} of which have protein peaks"
    print(msg.format(chromosome, maxBinInt, len(proteinDf.signalValue.to_numpy().nonzero()[0]) ))
    msg = "mean signal value: {0:.3f} (min {1:.3f}, max {2:.3f}, med {3:.3f})"
    print(msg.format(meanSignalValueOverall, minSignalValueOverall, maxSignalValueOverall, medianSignalValueOverall))
    msg = "mean signal value at TAD boundaries: {0:.3f}"
    print(msg.format(meanSignalValueAtBoundaries))
    msg = "mean signal value at TAD starts: {0:.3f} (min {1:.3f}, max {2:.3f}, med {3:.3f})"
    print(msg.format(meanSignalValueAtStart, minSignalValueAtStart, maxSignalValueAtStart, medianSignalValueAtStart))
    msg = "mean signal value at TAD ends: {0:.3f} (min {1:.3f}, max {2:.3f}, med {3:.3f})"
    print(msg.format(meanSignalValueAtEnd, minSignalValueAtEnd, maxSignalValueAtEnd, medianSignalValueAtEnd))
    msg = "{0:d} of {1:d} TADs have no protein at start bin"
    print(msg.format(zeroAtStartCount, tadCount))
    msg = "{0:d} of {1:d} TADs have no protein at start bin +/- {2:d} bins"
    print(msg.format(zeroAtStartPlusMinusCount, tadCount, nr_bins-1))
    msg = "{0:d} of {1:d} TADs have no protein at end bin"
    print(msg.format(zeroAtEndCount, tadCount))
    msg = "{0:d} of {1:d} TADs have no protein at end bin +/- {2:d} bins"
    print(msg.format(zeroAtEndPlusMinusCount, tadCount, nr_bins-1))
    msg = "{0:d} of {1:d} TADs have no protein at start and end bin"
    print(msg.format(zeroAtStartEndCount, tadCount))
    msg = "{0:d} of {1:d} TADs have no protein at start and end bin +/- {2:d} bins"
    print(msg.format(zeroAtStartEndPlusMinusCount, tadCount, nr_bins-1))

    #build a pdf for signal values
    valList = list(proteinDf.signalValue)
    probabList = np.ones(len(valList)) * (1/proteinDf.shape[0])
    
    #build cdf for mean at TAD boundaries, nr TADs with no protein
    meanArr = []
    tadStartArr = []
    for i in range(0,50000):
        tadSignalVals = np.random.choice(valList, tadCount, replace=False,p=probabList)
        meanArr.append(np.mean(tadSignalVals))
        tadStartArr.append(np.count_nonzero(tadSignalVals == 0.0))

    #plot
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1)
    proteinDf.signalValue.plot(kind='hist', density=True, bins=100, ax=ax1)
    ax1.set_title("probability density for signalValue at chr" + chromosome)
    p2 = ax2.hist(x=meanArr, bins=100, density=True)
    ax2.set_title("probability density for mean at TAD boundaries at chr" + chromosome)
    p3 = ax3.hist(x=tadStartArr, bins=60, density=True)
    ax3.set_title("probability density for number of TAD starts/ends without protein at chr" + chromosome)
    plt.show()

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


if __name__ == "__main__":
    compareTadsToPeaks()