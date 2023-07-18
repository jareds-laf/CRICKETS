# cosmic_flag_rfi_kurtosis

This package is designed to flag heavy RFI frequency bins in data that comes from COSMIC (Tremblay, Varghese, et al. (in preparation)). This is accomplished by generated a time-averaged power spectrum from an input .fil file and analyzing the excess kurtosis (A.K.A. exkurt) of the power in a specified number of frequency bins.

In its current state, this package is **NOT** for differentiating any signals of scientific interest from RFI. Its strength is combing through observations of sources with little to no fine frequency emissions that could be mistaken for RFI. Assuming noise follows a Gaussian distribution, any frequency bins with an excess kurtosis outside of a specifiable range around 0 are likely RFI.


With the limitations of this program in mind, the best use of this package is to flag frequency ranges that are heavy in RFI during the time of observation and masking these frequencies after the data from the primary observations have been collected. The package is also effective for studying the overall RFI environment of a series of observations at different times.

**A note on terminology: With the typical (Pearson) definition of kurtosis, a Gaussian distribution has a kurtosis of 3 and an "excess kurtosis" of 0. This package uses the Pearson definition of kurtosis as of July 18th, 2023. There may be some references to "kurtosis" (not excess kurtosis) throughout the code and documentation. If you come across any of these instances (other than the name of the repo and package itself), please contact sofairj@lafayette.edu :)**

## Installation:
As of July 18th, 2023, the package is still a work in progress and it is not entirely functional. Once completed, there will be two primary ways to install the package.

### Dependencies
The versions of the following required packages shown are simply a snapshot of the versions used in the development of this package. It is very likely that other versions of the dependencies, both past and future, can be used for this package.

- blimpy==2.1.4
- matplotlib==3.7.2
- numpy==1.25.1
- pandas==2.0.3
- scipy==1.11.1

### Install by cloning the repository
Find or create the folder you would like to clone this repository to, then use the following command to clone via HTTPS:

```
git clone https://github.com/jareds-laf/cosmic-flag-rfi-kurtosis.git
```

Alternatively, you can use this command to clone via SSH:

```
git clone git@github.com:jareds-laf/cosmic-flag-rfi-kurtosis.git
```

### Install via pip (WIP)
This package is currently available on [Test PyPI](https://test.pypi.org/project/cosmic-flag-rfi-kurtosis/). It can be installed with the following command:

```
pip install -i https://test.pypi.org/simple/ cosmic-flag-rfi-kurtosis
```


## Summary of the Process
Once again, as of July 18th, 2023, most of the functionality of this package is a work in progress. What follows is a summary of the general process expected in the coming weeks.

### Flow of Analysis
The input filterbank file is used to generate a [blimpy](https://github.com/UCBerkeleySETI/blimpy) waterfall object. The power is then averaged over the time domain, and the waterfall object is split into a specifiable number of frequency bins. The default number of bins is 256. So as to avoid any infinite excess kurtosis, the data is rescaled. Then, the excess kurtosis of each bin is then calculated using [scipy](https://github.com/scipy/scipy), and bins with a high* excess kurtosis are flagged. The data is once again checked for any unwanted infinities. The high RFI bins are output in a .csv file with the following columns:

- **rfi_bin_bots**: High RFI frequency bin bottoms
- **rfi_bin_tops**: High RFI frequency bin tops
- **exkurt**: Excess kurtosis of corresponding bin

*The minimum threshold to flag high excess kurtosis bins can be specified by the user. The threshold constitutes the ends of the range of kurtoses that is flagged. That is, the flagged bins fall within the range

threshold &#8804; excess kurtosis &#8804; threshold ⇔ threshold &#8804; |kurtosis|

The user can also choose to include all frequency bins in the output file so that they can determine the minimum excess kurtosis threshold more effectively. One might use a higher threshold if they are using this package as a quick and dirty way to flag problematic frequency ranges. A lower threshold would be useful if the user knows there are strong yet short-lived RFI signatures at certain frequencies in their data. 

### Plotting Functions
The user can choose to generate two types of plots to aid them in their scientific endeavors using the following options:
1. --plot_kurt: Plot the excess kurtosis of each bin against their corresponding bin bottoms. One can choose to include up to 3 of the following data categories to plot the excess kurtosis of:
   1. Unfiltered data
   2. "Clean"/low RFI channels
   3. "Dirty"/high RFI channels
2. --plot_pwr: Plot the time-averaged power spectrum (i.e., time-averaged power *vs*. frequency). The user can specify whether or not to show the bins that have been flagged as having heavy RFI. Flagged bins will have as transparent red rectangles that span the height of the graph and are "behind" the main spectrum so as not to steal the show.
    
## Usage
All functions of this package can be run from the command line. The general syntax is as follows:

```
python3 <path-to-rfi_flag_bulk.py> [Options]
```

Options:
- ```--input_filename``` (Required) Path to input filterbank file (including file name).
- ```--output_filename``` (Required) Path to output csv file (including file name).
- ```--threshold, -T``` (Required) Minimum value of excess kurtosis used to flag channels with significant RFI. Can be any decimal number.
- ```--ndivs, -N``` (Required) Number of frequency bins to split waterfall object into. Can be any integer.

Plotting options are being improved as of July 18th, 2023!

Use either of the following commands to list all options, as well!

```
python3 <path-to-rfi_flag_bulk.py> -h
``` 

```
python3 <path-to-rfi_flag_bulk.py> --help
```

## Examples
### Running analysis without plots
The template for running analysis without generating any plots is as follows:

```
python3 <path/to/rfi_flag_bulk.py> --input_filename <path/to/filterbank.fil> --output_filename <path/to/output.csv> -T <kurtosis_threshold(float)> -N <number_of_bins(int)>
```

Here is an example with a minimum kurtosis threshold of 5 where the waterfall object gets broken into 256 frequency bins.

```
python3 /home/alice/cosmic-flag-rfi-kurtosis/cosmic_flag_rfi_kurtosis/rfi_flag_bulk.py --input_filename /home/alice/filterbank/filterbank1.fil --output_filename /home/alice/example_output.csv -T 5 -N 256
```

Here is a screenshot of the output table taken in Microsoft Excel:

![Alt text](/examples/example_output_excel.png?raw=true "Example output table")

### Running analysis with plots
WIP!