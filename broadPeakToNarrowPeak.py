import click
import pandas as pd
import os


@click.option('--broadpeakfile', '-bpf', type=click.Path(exists=True), help="broadPeak file to convert")
@click.option('--outfile', '-o', type=click.Path(writable=True), default=None, help="narroPeak file to write", required=False)
@click.command()
def broadPeakToNarrowPeak(broadpeakfile, outfile):
    ### convert broadPeak files to valid narrowPeak files with
    ### the peak column set to -1, i.e., no peak called
    
    peaksDf = None
    try:
        peaksDf = pd.read_table(broadpeakfile, header=None, sep="\t")
    except Exception as e:
        print(e, "\n")
        msg = "could not read broadPeak file {0:s} \n"
        msg += "maybe this is no proper broadPeak file?"
        msg = msg.format(broadpeakfile)
        raise ValueError(msg)
    if peaksDf.shape[1] != 9:
        msg = "broadPeak file {0:s} should have 9 columns.\n"
        msg += "However, it has {1:d} and thus can't be converted"
        msg = msg.format(broadpeakfile, peaksDf.shape[1])
        raise ValueError(msg)
    
    #set the peak column to -1, which means no peak has been called
    peaksDf['peak'] = -1

    #write out as valid narrowPeak file
    if not outfile:
        outfile = os.path.splitext(broadpeakfile)[0] + ".narrowPeak"
    peaksDf.to_csv(outfile, header=None, sep="\t", index=False, float_format='%.6f')


if __name__ == "__main__":
    broadPeakToNarrowPeak()