import matplotlib.pyplot as plt
import numpy as np
import os
import numpy.ma as ma
from analysis import get_exkurt


# TODO: Make it so that these functions access the table that
# was created by the write_output_table() function!
'''^^ This is definitely a stretch goal. It's more of an optimization thing.
It would take some more structural changes I don't have time to deal with this summer!'''

# Allows user to save figures as multiple file types at once
# Credit: https://stackoverflow.com/questions/17279651/save-a-figure-with-multiple-extensions
def save_fig(filename, types=['png']):
    fig = plt.gcf()
    for filetype in types:
        fig.savefig(f'{os.path.realpath(os.path.expanduser(filename))}.{filetype}')

def plot_tavg_power(wf_in,
                    f_start=0, f_stop=6000,
                    p_start=0, p_stop=5*10**10, n_divs=256, threshold=50,
                    show_filtered_bins=True,
                    output_dest='', output_type=['png']):
    # Plot the time-averaged power spectrum for a given blimpy waterfall object
    # Inputs:
        # wf: The desired input waterfall object
        # t: The integration number
        # f_start: Lower bound for frequency (horizontal) axis
        # f_stop: Upper bound for frequency (horizontal) axis
        # p_start: Lower bound for time-averaged power (veritcal) axis
        # p_start: Lower bound for time-averaged power (veritcal) axis
        # show_filtered: If true, mark high RFI channels with a vertical red box
        # output_dest: Location (including filename) to save output file
        # output_type: Filetype of output

    # Time average the power
    wf_pwr_mean_arr = np.mean(wf_in.data, axis=0)
    wf_pwr_mean = wf_pwr_mean_arr[0]
   
    # Plot time-averaged power
    fig, ax = plt.subplots()
    
    ax.set_xlim(f_start, f_stop)
    ax.set_ylim(p_start, p_stop)
    
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Time-Averaged Power (Counts)')

    ax.plot(wf_in.get_freqs(), wf_pwr_mean,
            label='Time-averaged power spectrum',
            c='#1f1f1f')

    # Grab info for RFI masking
    bins, kurts, pows_mean, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask = get_exkurt(wf_in, n_divs, threshold)

    # Plot frequency bins that were flagged as RFI
    if show_filtered_bins == True:
        full_freq_range = np.amax(wf_in.get_freqs()) - np.amin(wf_in.get_freqs())
        bin_width = full_freq_range / n_divs

        for i, rfi_bin in enumerate(flagged_bins.filled(np.nan)):
            xmin = rfi_bin
            xmax = rfi_bin + bin_width
            if i == 0:
                print('\n\nI think this is where the problem is happening?')
                print(f'Info about flagged_bins: {np.shape(flagged_bins)}, {type(flagged_bins)}, {ma.count(flagged_bins)}/{(ma.count_masked(flagged_bins) + ma.count(flagged_bins))}\n\n')
                print(f'Here is flagged_bins: {flagged_bins}\n\n')
                print('\n\nAnd here is a comparison of the shapes and counts of flagged_bins before and after using .filled:')
                print(f'flagged_bins\n{np.shape(flagged_bins)}, {ma.count(flagged_bins)}, {ma.count_masked(flagged_bins)}')
                print(f'flagged_bins.filled\n{np.shape(flagged_bins.filled(np.nan))}, {np.count_nonzero(~np.isnan(flagged_bins.filled(np.nan)))}, {np.count_nonzero(np.isnan(flagged_bins.filled(np.nan)))}')

            flagged_line = plt.axvspan(xmin=xmin, xmax=xmax, ymin=0, ymax=1, color='red', alpha=0.5)

        flagged_line.set_label('RFI-heavy channels')
        ax.legend(fancybox=True,shadow=True, loc='lower center', bbox_to_anchor=(1, 1), ncols=1)
    else:
        ax.legend(fancybox=True,shadow=True, loc='lower center', bbox_to_anchor=(1, 1), ncols=1)

    save_fig(os.path.realpath(os.path.expanduser(output_dest)), types=output_type)

def plot_mask_exkurt(wf_in, n_divs=256, threshold=50,
                       unfiltered=True, clean_chnls=True, rfi=False,
                      f_start=2000, f_stop=4000,
                      k_start=-5, k_stop=500,
                      output_dest='', output_type='png'):
    # This function plots the excess kurtosis of each frequency channel for a specified waterfall object.
    # Inputs:
        # wf_in: See get_exkurt() function definition
        # n_divs: See get_exkurt() function definition
        # threshold: See get_exkurt() function definition
        # unfiltered: If true, plot the data before any RFI filtering has occurred
        # clean_chnls: If true, plot the data after RFI has been filtered out
        # rfi: If true, plot the channels that have been marked as RFI
        # output_dest: Location (including filename) to save output file
        # output_type: Filetype of output
    
    bins, kurts, pows_mean, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask = get_exkurt(wf_in, n_divs, threshold)
    
    # Create the plot
    fig, ax = plt.subplots()
    
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Excess Kurtosis')
    
    if unfiltered: # Plot all data
        ax.plot(bins, kurts, 'o', c='black', label='Unfiltered data') # Color is a nice black
    if clean_chnls: # Plot the low RFI channels
        print('\n\nConfirmation that we are running properly\n\n')
        ax.plot(bins, masked_kurts.filled(np.nan), '.', c='#43cc5c', label='Clean channels') # Color is a nice green
    if rfi: # Plot the high RFI channels
        ax.plot(flagged_bins, flagged_kurts, '.', c='red', label='Heavy RFI') # Color is a nice red
    
    # TODO: Change this condition... :)
    if np.any([f_start, f_stop, k_start, k_stop]) != 0:
        ax.set_xlim(f_start, f_stop)
        ax.set_ylim(k_start, k_stop)
    
        ax.legend(fancybox=True,shadow=True, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncols=3)
    else:
        ax.legend(fancybox=True,shadow=True, loc='upper center', bbox_to_anchor=(0.5, 1.05), ncols=3)
    
    save_fig(os.path.realpath(os.path.expanduser(output_dest)), types=output_type)