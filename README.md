# CRICKETS

CRICKETS (Categorization of RFI In COSMIC with Kurtosis for Extraterrestrial Searches) is a packaged designed to flag heavy RFI frequency bins in data that comes from [COSMIC](https://science.nrao.edu/facilities/vla/observing/cosmic-seti). This is accomplished by generated a time-averaged power spectrum from an input .fil file and analyzing the excess kurtosis (exkurt) of the power in a specified number of frequency bins.

In its current state, this package is **NOT** for differentiating any signals of scientific interest from RFI. Its strength is combing through observations of sources with little to no fine frequency emissions that could be mistaken for RFI. Assuming noise follows a Gaussian distribution, any frequency bins with an excess kurtosis outside of a specifiable range around 0 are likely RFI. With the limitations of this program in mind, the best use of this package is to flag frequency ranges that are heavy in RFI and masking these frequencies after the data from the primary observations have been collected. The package is also effective for studying the overall RFI environment of a series of observations at different times.

*A note on terminology: With the typical (Pearson) definition of kurtosis, a Gaussian distribution has a kurtosis of 3 and an "excess kurtosis" of 0. This package uses the Pearson definition of kurtosis as of July 18th, 2023. There may be some outdated references to "kurtosis" (not excess kurtosis) throughout the code and documentation from the pre-July-18th era when the Fisher definition was used. If you come across any of these instances (other than the name of the repo and package itself), please contact sofairj@lafayette.edu :)*

# Installation:
As of July 25th, 2023, the package is still a work in progress and it is not entirely functional. Once completed, there will be two primary ways to install the package.

## Dependencies
The versions of the following required packages shown are simply a snapshot of the versions used in the development of this package. If there are any issues with this package, try installing these specific versions of the code as a troubleshooting step. Otherwise, install any version and the code *shoud* work.

- blimpy==2.1.4
- matplotlib==3.7.2
- numpy==1.25.1
- pandas==2.0.3
- scipy==1.11.1

## Install by cloning the repository
Find or create the folder you would like to clone this repository to, then use the following command to clone via HTTPS:

```
git clone https://github.com/jareds-laf/CRICKETS.git
```

Alternatively, you can use this command to clone via SSH:

```
git clone git@github.com:jareds-laf/CRICKETS.git
```

## Install via pip (WIP)
This package is currently available on [Test PyPI](https://test.pypi.org/project/crickets/). It can be installed with the following command:

```
pip install -i https://test.pypi.org/simple/ crickets
```

# Summary of the Process
As of July 25th, 2023, most of the functionality of this package works! Some of it is still a work in progress. Within the coming weeks, most bugs should be ironed out and every function listed below should be implemented :)

## Flow of Analysis
An input filterbank file is used to generate a [blimpy](https://github.com/UCBerkeleySETI/blimpy) waterfall object. The power is then averaged over the time domain, and the waterfall object is split into a specifiable number of frequency bins. The default number of bins is 256. So as to avoid any infinite excess kurtosis, the data is rescaled. Then, the excess kurtosis of each bin is then calculated using [scipy](https://github.com/scipy/scipy), and bins with a high* excess kurtosis are flagged. The data is once again checked for any unwanted infinities. The high RFI bins are output in a .csv file with the following columns:

- **rfi_bin_bots**: High RFI frequency bin bottoms
- **rfi_bin_tops**: High RFI frequency bin tops
- **exkurt**: Excess kurtosis of corresponding bin

*The minimum threshold to flag high excess kurtosis bins can be specified by the user. The threshold constitutes the ends of the range of kurtoses that is flagged. That is, the flagged bins fall within the range

threshold &#8804; |excess kurtosis|

(WIP) The user can also choose to include all frequency bins in the output file so that they can determine the minimum excess kurtosis threshold by any means they see fit. One might use a higher threshold if they are using this package as a quick and dirty way to flag problematic frequency ranges. A lower threshold would be useful if the user knows there are strong yet short-lived RFI signatures at certain frequencies in their data. 

### Plotting Functions
The user can choose to generate two plots to aid them in their scientific endeavors using the ```-p``` option, or by importing the necessary functions into a Python script or notebook (WIP). The two types of plots are:
1. exkurt: Plot the excess kurtosis of each bin against their corresponding bin bottoms. One can choose to include up to 3 of the following data categories to plot the excess kurtosis of:
   1. Unfiltered data
   2. "Clean"/low RFI channels
   3. "Dirty"/high RFI channels
2. tavg_pwr: Plot the time-averaged power spectrum (i.e., time-averaged power *vs*. frequency). The user can specify whether or not to show the bins that have been flagged as having heavy RFI. Flagged bins will have as transparent red rectangles that span the height of the graph and are "behind" the main spectrum so as not to steal the show.
    
# Usage
## From the Command Line
All functions of this package can be performed from the command line. The general syntax is as follows:

```
python3 <path-to-rfi_flag_bulk.py> [Options]
```

Options:
- ```--input_file``` (Required) Path to input filterbank file (including file name).
- ```--output_file``` (Required) Path to output csv file (optionally including file name).
- ```--threshold, -t``` (Required) Minimum value of excess kurtosis used to flag channels with significant RFI. Can be any decimal number.
- ```--ndivs, -n``` (Required) Number of frequency bins to split waterfall object into. Can be any integer.
- ```--plot, -p``` (Optional) Choose whether or not to generate time-averaged power spectrum and excess kurtosis vs. frequency plots. Give output path for plots here (NOT including file name).
- ```plot_file_types``` (Optional, unless -p is given) Output file types (can be pdf, png, and/or jpg). Specify as many of these as you want! 
- ```--verbose, -v``` (Optional) Print more information about the input variables and the processes currently running.
Use either of the following commands to list all options, as well!

```
python3 <path-to-rfi_flag_bulk.py> -h
``` 

```
python3 <path-to-rfi_flag_bulk.py> --help
```
## As a Python Package
WIP!

# Examples
## Command Line Usage
### Running analysis from the command line without plots

The template for running analysis without generating plots is as follows:

```
python3 <path/to/rfi_flag_bulk.py> --input_file <path/to/filterbank.fil> --output_filename <path/to/output.csv> -t <kurtosis_threshold> -n <number_of_bins>
```

Here is an example with a minimum excess kurtosis threshold of 5 where the waterfall object gets broken into 256 frequency bins.

```
python3 /home/alice/cosmic-flag-rfi-kurtosis/cosmic_flag_rfi_kurtosis/rfi_flag_bulk.py --input_file /home/alice/filterbank/filterbank1.fil --output_file /home/alice/example_output.csv -t 5 -n 256
```

Here is a screenshot of the output table taken in Microsoft Excel:

<img src="https://github.com/jareds-laf/CRICKETS/blob/main/examples/example_output_excel.png" alt="Example table output. This table can be found at example/example_output_excel.png" width="227" height="350" />

### Running analysis from the command line with plots
***WIP!***

The template for running analysis and generating plots is as follows:

```
python3 <path/to/rfi_flag_bulk.py> --input_file <path/to/filterbank.fil> --output_filename <path/to/output.csv> -t <kurtosis_threshold> -n <number_of_bins> -p <path/to/plots.png> --plot_file_types <.png, .jpg, .pdf>
```

Here is the same example as above (minimum excess kurtosis threshold of 5 with 256 frequency bins), this time generating the two plots.

```
python3 /home/alice/cosmic-flag-rfi-kurtosis/cosmic_flag_rfi_kurtosis/rfi_flag_bulk.py --input_file /home/alice/filterbank/filterbank1.fil --output_file /home/alice/example_output.csv -t 5 -n 256 -p /home/alice/plots --plot_file_types png pdf jpg
```

Here are some screenshots of the resulting plots.

.png:

.jpg:

.pdf.

## Usage as a Package
WIP!
```
from crickets.analysis import get_exkurt, write_output_table

from crickets.plotting import plot_tavg_power, plot_exkurt
```