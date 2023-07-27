import os
import numpy as np
from scipy.stats import norm, kurtosis
import numpy.ma as ma
import pandas as pd
import matplotlib.pyplot as plt
from blimpy import Waterfall
from blimpy import calcload
import time
import glob

##### Analysis functions #####

def normalize_path(in_path):
    # A quick function to ensure that any input paths are properly referenced
	return os.path.normpath(os.path.realpath(os.path.expanduser(in_path)))

def create_info_table(wf_file_full, saveloc="./"):
    """
    This function creates a table containing the time-averaged power and frequencies as read from an input filterbank file.
    It does so by initializing a blimpy Waterfall object and then averaging the power.
    It outputs the time-averaged power and frequencies into a .csv file in a specified location.
    
    Inputs:
        wf_file_full: Input filterbank file
        saveloc: Location to save output table (NOT including file name)
    """
    # Normalize save path and path to wf_file_full
    t_init = time.time()
    saveloc = normalize_path(saveloc)
    wf_file_full = normalize_path(wf_file_full)

    # Get file name    
    wf_file = os.path.basename(wf_file_full)[:os.path.basename(wf_file_full).find('.')]

    # Initialize blimpy waterfall object
    t0 = time.time()
    print('\nGenerating waterfall object...')
    ml = calcload.calc_max_load(wf_file_full)
    wf = Waterfall(os.path.normpath(wf_file_full), max_load = ml)
    t1 = time.time()
    print(f'Done. Elapsed time: {t1 - t0}')

    # Get power and frequency in increasing order
    if wf.header['foff'] < 0:
        pows = np.flip(wf.data)
        freqs = wf.get_freqs()[::-1]
    else:
        pows = wf.data
        freqs = wf.get_freqs()

    # So:
    # pows_flipped is all of the powers in increasing order,
    # freqs_flipped is all of the frequencies in increasing order

    # Time-average the power
    pows_mean = np.mean(pows, axis=0)[0]

    # Create table with time-averaged power and frequencies
    table = pd.DataFrame(columns=['freq', 'tavg_power'])
    table['freq'] = freqs
    table['tavg_power'] = pows_mean

    # Save table
    save_path = os.path.join(saveloc, f'info_table_{wf_file}.csv')
    table.to_csv(save_path, index=False)
    
    print(f'\nInfo table saved to {save_path}\n')
    t_final = time.time()
    print(f'\nTotal elapsed time: {t_final - t_init}')

def create_info_table_dir(wf_dir, saveloc):
    file_list = glob.glob(wf_dir + '/*.fil')
    print(f"\nwf_dir from create_info_table_dir: {wf_dir}")
    print(f"\nfile list from create_info_table_dir: {file_list}\n")

    for file in file_list:
        create_info_table(wf_file_full=file, saveloc=saveloc)

    # for root, dirs, files in os.walk(wf_dir):
    #     for file in files:
    #         if file.endswith('.fil'):
    #             wf_file_full = os.path.join(root, file)
    #             create_info_table(wf_file_full=wf_file_full, saveloc=saveloc)

def get_exkurt(info_table, n_divs=256, threshold=50):
    """
    This function grabs the excess kurtosis of channels of a specified size for a blimpy waterfall object (section 1)
    and flags bins with high excess kurtosis (i.e., heavy RFI) and returns the necessary information about these bins (section 2).
    Inputs:
        info_table: Location of the info table containing the frequencies and time-averaged powers (including file name)
        n_divs: Number of divisions to break wf_in into
            32 is the correct number of channels to break a waterfall into assuming the frequency range
            of the waterfall is 32 MHz
        threshold: Minimum excess kurtosis for a channel to be flagged as 'RFI-heavy'
    """

##### Section 1 #####

    # TODO: Update this!
    # Get frequencies and powers from info_table
    info_table = normalize_path(info_table)
    info_table_read = pd.read_csv(info_table)
    freqs = np.array(info_table_read['freq'])
    pows = np.array(info_table_read['tavg_power'])

    # Split frequency and time-averaged power into n_divs channels
    freqs = np.array_split(freqs, n_divs)
    pows = np.array_split(pows, n_divs)
    
    # Get excess kurtosis of all channels
    exkurts_list = []
    
    for division in pows:
        exkurts_list.append(kurtosis(division/(10**9))) # Rescaling data so that excess kurtosis != inf ever (hopefully)
    
    exkurts = np.array(exkurts_list, dtype=np.float64)
    
    # Check to see if any of the kurtoses are infinite
    if np.any(np.isfinite(exkurts)) == False:
        print(f'Infinite kurtoses located at indices {np.where(np.isfinite(exkurts) == False)}')
    else:
        print('No infinite kurtoses found')

    # Binning frequencies such that the labeled frequency is the bottom of the bin
    # i.e., if chnl[0] is 2010 MHz and each channel is 1 MHz, then the bin from 2010 MHz to 2010.99 MHz will have
    # a value of "2010"
    bins = []
    for chnl in freqs:
        bins.append(chnl[0])

##### Section 2 #####
    # This part of the function flags bins with high excess kurtosis.
    
    # masked_kurts is an array that has all channels with |excess kurtosis| > threshold masked out
    masked_kurts = ma.masked_where(np.abs(exkurts) > threshold, exkurts)
    bin_mask = ma.getmask(masked_kurts)
    
    # flagged_bins is an array that has the frequencies of the channels with excess kurtosis > threshold NOT masked out
    # flagged_kurts masks the opposite elements as masked_kurts (i.e., it contains all of the kurtoses of the
    # high RFI channels)
    flagged_bins = ma.masked_array(bins, mask=~bin_mask)
    flagged_kurts = ma.masked_array(exkurts, mask=~bin_mask)
    
    # The reason I am no longer using ma.count(masked_array) is because there can be an error where masked elements
    # are converted to NaNs.
    print(f'{len(np.where(bin_mask == True)[0])} out of {n_divs} channels flagged as having substantial RFI')
    

    # TODO: Make sure this isn't completely broken!

    # Bin width (for the files I am testing this with):
    # There are 7 Hz in between each entry of freqs, and 32 MHz total
    # Given a total of 4194304 elements in freqs, if there are 256 bins, then each bin spans 16384 elements
    # If each bin spans 16384 elements, then it spans 125 kHz (125 kHz * 256 = 32 MHz)
    
    # Grab the bin width in terms of MHz (for convenience if needed in the future)
    full_freq_range = np.array(info_table_read['freq'])[-1] - np.array(info_table_read['freq'])[0]
    # print(np.shape())
    # print(freqs[-1])
    # print()
    # print(freqs[0])
    bin_width = full_freq_range / n_divs


    print(f"\n\n\nBin width in terms of f: {bin_width}\n\n\n")
    
    # Grab the bin width in terms of the number of elements per bin
    bin_width_elements = int(np.floor(len(freqs) / n_divs))
        
    masked_freqs = ma.masked_array(freqs)
    
    for rfi_bin in flagged_bins:
        try:
            # Get the frequency indices of the masked frequency bins and put them in a list
            xmin = np.where(freqs == rfi_bin)[0][0]
            xmax = xmin + bin_width_elements
            masking_indices = np.arange(xmin, xmax)
            
            # Create a masked array that masks the indices of the high RFI bins
            masked_freqs[masking_indices] = ma.masked
            freq_mask = ma.getmask(masked_freqs)
        except:
            pass
    
    # Summary of returned variables:
        # bins: All frequency bins for the frequency range of the waterfall
        # kurts: Excess kurtosis of all frequency bins
        # pows: Time-averaged power of each frequency bin
        # flagged_bins: Channels with high excess kurtosis (i.e., high RFI)
        # flagged_kurts: Excess kurtosis of each channel that was flagged as having high RFI
        # masked_kurts: Excess kurtosis of all 'clean' (low RFI) channels with high RFI channels masked out
        # masked_freqs: List of all frequencies with high RFI channels masked out
        # bin_mask: The mask used to generate masked_kurts -- Masks the frequency *bins*
        # freq_mask: A mask to be used to block out the *actual* frequencies, rather than the frequency *bins*
    return bins, exkurts, pows, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask

def write_output_table(info_table, output_filepath='./', n_divs=256, threshold=50, all=False):
    """This function does as it says: It writes the output table. It does so in a .csv format with columns of:
        (bin_top) Frequency bin tops
        (bin_bot) Frequency bin bottoms
        (exkurt) Excess kurtosis of each bin
    Inputs:
        info_table: Location of the info table containing the frequencies and time-averaged powers (including file name)
        output_filepath: Path to output .csv file
        n_divs: See get_exkurt() function definition
        threshold: See get_exkurt() function definition
        all: If True, include all frequency bins in the output table, even if they are not flagged as RFI. If False, only the flagged bins will be included in the output table."""

    # Assign all the base variables and ensure file export path (export_path) is normalized
    print(f"\nInfo table in write_output_table: {info_table}\n")
    export_path = normalize_path(output_filepath)
    bins, kurts, pows, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask = get_exkurt(info_table, n_divs, threshold)

    # Get frequencies and powers from info_table
    info_table = normalize_path(info_table)
    info_table_read = pd.read_csv(info_table)

    # Get bin tops
    freqs = np.array(info_table_read['freq'])
    max_freq = freqs[-1]
    min_freq = freqs[0]
    bin_width_elements = int(np.floor(len(freqs) / n_divs)) # Grab the bin width in terms of the number of elements per bin
    # bin_tops = flagged_bins + bin_width_elements
    # bin_tops = [freqs[int(b + bin_width_elements)] for b in flagged_bins if not np.ma.is_masked(b)]
    
    full_freq_range = freqs[-1] - freqs[0]
    print(f"full_freq_range: {full_freq_range}")
    bin_width = full_freq_range / n_divs
    print(f"bin_width: {bin_width}")
    bin_tops = []
    for bin_bot in flagged_bins:
        bin_tops.append(bin_bot + bin_width)

    print(f"bin_tops: {bin_tops}")
    # Get bin tops


    print(f"\nFlagged bins: {flagged_bins}, {type(flagged_bins)}")
    print(f"\nBin width elements kurts: {bin_width_elements}, {type(bin_width_elements)}")
    print(f"\nBin tops: {bin_tops}, {type(bin_tops)}, {len(bin_tops)}")

    # Format flagged_bins into a regular (non-masked) numpy array
    flagged_bins = ma.filled(flagged_bins, fill_value=np.NaN)

    # Turns the numpy arrays into pandas dataframes so they can be concatenated and exported
    export_bin_bots = pd.DataFrame(data=flagged_bins, columns=['rfi_bin_bots'])
    export_bin_tops = pd.DataFrame(data=bin_tops, columns=['rfi_bin_tops'])
    export_bin_kurt = pd.DataFrame(data=flagged_kurts, columns=['kurt'])

    # Concatenate dataframes
    export_concat = pd.concat([export_bin_bots,
                               export_bin_tops,
                               export_bin_kurt], axis=1)

    # Sort dataframe by frequency, remove pandas indexing, remove blank lines
    export_df = export_concat.sort_values(by=['rfi_bin_bots']).reset_index(drop=True).dropna(how='all')

    # Write dataframe to csv at export_path
    export_df.to_csv(export_path, index=False)

    if all:
        # Quick shout out to GitHub copilot for helping me write this function! :)
        # Parts of this comment as well as the lines in this part of the function were written by GitHub copilot.
        # It just knew what to type!

        # Get bin tops
        # bin_tops = bins + bin_width_elements
        print(f"Bins: {type(bins)}\n")

        all_bin_bots = pd.DataFrame(data=bins, columns=['all_bin_bots'])
        all_bin_tops = pd.DataFrame(data=bins, columns=['all_bin_tops'])
        all_bin_kurt = pd.DataFrame(data=kurts, columns=['all_kurt'])

        all_concat = pd.concat([all_bin_bots, all_bin_tops, all_bin_kurt], axis=1)
        export_all_df = all_concat.sort_values(by=['all_bin_bots']).reset_index(drop=True).dropna(how='all')
        export_path_all = export_path.replace('.csv', '_all.csv')

        print(f"\nAll table output location: {export_path_all}\n")


        export_all_df.to_csv(export_path_all, index=False)
    else:
        return export_df


##### Plotting functions #####


# TODO: Make it so that these functions access the table that
# was created by the write_output_table() function!
'''^^ This is definitely a stretch goal. It's more of an optimization thing.
It would take some more structural changes I don't have time to deal with this summer!'''

def save_fig(filename, types=['png']):
    """Allows user to save figures as multiple file types at once
    Credit: https://stackoverflow.com/questions/17279651/save-a-figure-with-multiple-extensions"""
    fig = plt.gcf()
    for filetype in types:
        fig.savefig(f'{normalize_path(filename)}.{filetype}')

def plot_tavg_power(info_table,
                    f_start=0, f_stop=6000,
                    p_start=0, p_stop=5*10**10, n_divs=256, threshold=50,
                    show_filtered_bins=True,
                    output_dest='', output_type=['png']):
    """Plot the time-averaged power spectrum for a given blimpy waterfall object
    Inputs:
        info_table: Location of the info table containing the frequencies and time-averaged powers (including file name)
        t: The integration number
        f_start: Lower bound for frequency (horizontal) axis
        f_stop: Upper bound for frequency (horizontal) axis
        p_start: Lower bound for time-averaged power (veritcal) axis
        p_start: Lower bound for time-averaged power (veritcal) axis
        show_filtered: If true, mark high RFI channels with a vertical red box
        output_dest: Location (including filename) to save output file
        output_type: Filetype of output"""

    # Get frequencies and powers from info_table    
    freqs = np.array(info_table['freq'])
    pows = np.array(info_table['tavg_pows'])

    # Plot time-averaged power
    fig, ax = plt.subplots()
    
    ax.set_xlim(f_start, f_stop)
    ax.set_ylim(p_start, p_stop)
    
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Time-Averaged Power (Counts)')

    ax.plot(freqs, pows,
            label='Time-averaged power spectrum',
            c='#1f1f1f')

    # Grab info for RFI masking
    bins, kurts, pows_mean, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask = get_exkurt(info_table, n_divs, threshold)

    # Plot frequency bins that were flagged as RFI
    if show_filtered_bins == True:
        full_freq_range = freqs[-1] - freqs[0]
        print(f"full_freq_range: {full_freq_range}")
        bin_width = full_freq_range / n_divs

        for i, rfi_bin in enumerate(flagged_bins.filled(np.nan)):
            xmin = rfi_bin
            xmax = rfi_bin + bin_width
            flagged_line = plt.axvspan(xmin=xmin, xmax=xmax, ymin=0, ymax=1, color='red', alpha=0.5)

        flagged_line.set_label('Dirty channels')
        ax.legend(fancybox=True,shadow=True, loc='lower center', bbox_to_anchor=(1, 1), ncols=1)
    else:
        ax.legend(fancybox=True,shadow=True, loc='lower center', bbox_to_anchor=(1, 1), ncols=1)

    save_fig(os.path.realpath(os.path.expanduser(output_dest)), types=output_type)

def plot_exkurt(info_table, n_divs=256, threshold=50,
                       unfiltered=True, clean_chnls=True, rfi=False,
                      f_start=2000, f_stop=4000,
                      k_start=-5, k_stop=500,
                      output_dest='', output_type='png'):
    """This function plots the excess kurtosis of each frequency channel for a specified waterfall object.
    Inputs:
        info_table: Location of the info table containing the frequencies and time-averaged powers (including file name)
        n_divs: Number of frequency bins to divide the data into
        threshold: Threshold for flagging a frequency bin as dirty
        unfiltered: If true, plot the data before any RFI filtering has occurred
        clean_chnls: If true, plot the data after RFI has been filtered out
        rfi: If true, plot the channels that have been marked as RFI
        output_dest: Location (including filename) to save output file
        output_type: Filetype of output"""
    
    bins, kurts, pows_mean, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask = get_exkurt(info_table, n_divs, threshold)
    
    # Create the plot
    fig, ax = plt.subplots()
    
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Excess Kurtosis')
    
    # Plot all data
    if unfiltered:
        ax.plot(bins, kurts, 'o', c='black', label='Unfiltered data') # Color is a nice black
    # Plot the low RFI channels
    if clean_chnls:
        # I use the .filled() method here because sometimes a warning pops up saying that a masked element was converted to a nan.
        # I figure removing the masked elements to begin with is a good way to avoid this.
        ax.plot(bins, masked_kurts.filled(np.nan), '.', c='#43cc5c', label='Clean channels') # Color is a nice green
    # Plot the high RFI channels
    if rfi:
        ax.plot(flagged_bins, flagged_kurts, '.', c='red', label='Dirty channels') # Color is a nice red
    
    # TODO: Change this condition... :)
    if np.any([f_start, f_stop, k_start, k_stop]) != 0:
        ax.set_xlim(f_start, f_stop)
        ax.set_ylim(k_start, k_stop)
    
        ax.legend(fancybox=True,shadow=True, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncols=3)
    else:
        ax.legend(fancybox=True,shadow=True, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncols=3)
    
    save_fig(os.path.realpath(os.path.expanduser(output_dest)), types=output_type)