# cosmic_flag_rfi_kurtosis

This package is designed to flag heavy RFI frequency bins in data that comes from COSMIC.

## Installation:
This package is currently available on [Test PyPI](https://test.pypi.org/project/cosmic-flag-rfi-kurtosis/). It can be installed with the following command:

```pip install -i https://test.pypi.org/simple/ cosmic-flag-rfi-kurtosis```

## Summary of the Process
As of July 14th, 2023, almost all of the functionality of this package is a work in progress. What follows is a summary of the general experience expected in the coming weeks.

### Analysis Workflow
The user begins by specifying a filterbank file. This file is used to generate a [blimpy](https://github.com/UCBerkeleySETI/blimpy) waterfall object. The power is then averaged over time, and the waterfall is then split into a number of frequency bins that can be specified by the user. The default number of bins is 256. The kurtosis of each bin is then calculated using [scipy](https://github.com/scipy/scipy), and bins with a "high" kurtosis are flagged. (*Note: The minimum threshold to flag high kurtosis bins can also be specified by the user. One might use a higher threshold if they are using this package as a quick and dirty way to flag problematic frequency ranges. A lower threshold would be useful if the user knows there are strong yet short-lived RFI signatures in certain regions of the data*). The high RFI bins are output in a table (.csv file) with the following columns:

[//]: # (TODO: Get the actual names of the columns in here when you're ready!)
- Frequency bin bottoms
- Frequency bin tops
- Kurtosis of corresponding bin
- Time-averaged power of corresponding bin

The user can choose to include all frequency bins in the output file so that they can determine the minimum kurtosis threshold on independently.

### Plotting Functions
The user can choose to generate two types of plots to aid them in their scientific endeavors:
1. plot_mask_kurtosis: Plot the kurtosis of each bin against their corresponding bin bottoms. One can choose to include up to 3 of the following data categories to plot:
   1. Unfiltered data
   2. Clean/low RFI channels
   3. Dirty/high RFI channels
2. plot_tavg_powertemp: Plot the time-averaged power spectrum (i.e., time-averaged power *vs*. frequency). The user can specify whether or not to show the bins that have been flagged as having heavy RFI. Flagged bins will appear as transparent red rectangles that span the height of the graph and are "behind" the main spectrum so as not to steal the show.
    
## Examples
WIP!
