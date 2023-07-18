import matplotlib
import matplotlib.pyplot as plt
import blimpy
from blimpy import Waterfall
from blimpy import calcload
import os
import glob
import numpy as np
import time
from scipy.stats import norm, kurtosis
import scipy
import numpy.ma as ma
import pandas as pd

def get_kurtosis(wf_in, n_divs=256, threshold=50):
    # This function grabs the kurtosis of channels of a specified size for a blimpy waterfall object (section 1)
    # and flags bins with high kurtosis (i.e., heavy RFI) and returns the necessary information about these bins (section 2).
    # Inputs:
        # wf_in: Specified blimpy waterfall object
        # n_divs: Number of divisions to break wf_in into
            # 32 is the correct number of channels to break a waterfall into assuming the frequency range
            # of the waterfall is 32 MHz
        # threshold: Minimum kurtosis for a channel to be flagged as 'RFI-heavy'

#     np.set_printoptions(threshold=4)

##### Section 1 #####

    # Get power and frequency in increasing order
    if wf_in.header['foff'] < 0:
        pows_flipped = np.flip(wf_in.data)
        freqs_flipped = wf_in.get_freqs()[::-1]
    
    # Time-average the power
    pows_mean_flipped = np.mean(pows_flipped, axis=0)[0]   

    # Split frequency and time-averaged power into n_divs channels
    freqs = np.array_split(freqs_flipped, n_divs)
    pows_mean = np.array_split(pows_mean_flipped, n_divs)
    
    # So:
    # pows_flipped is all of the powers in increasing order,
    # freqs_flipped is all of the frequencies in increasing order
    
    # Get kurtosis of all channels
    kurts_list = []
    
    for i, division in enumerate(pows_mean):
        kurts_list.append(kurtosis(division/(10**9))) # Rescaling data so that kurtosis != inf ever (hopefully)
    
    kurts = np.array(kurts_list, dtype=np.float64)
    
    # Check to see if any of the kurtoses are infinite
    if np.any(np.isfinite(kurts)) == False:
        print(f'Infinite kurtoses located at indices {np.where(np.isfinite(kurts) == False)}')
    else:
        print('No infinite kurtoses found')

    # Binning frequencies such that the labeled frequency is the bottom of the bin
    # i.e., if chnl[0] is 2010 MHz and each channel is 1 MHz, then the bin from 2010 MHz to 2010.99 MHz will have
    # a value of "2010"
    bins = []
    for chnl in freqs:
        bins.append(chnl[0])

##### Section 2 #####
    # This part of the function flags bins with high kurtosis.
    
    # masked_kurts is an array that has all channels with |kurtosis| > threshold masked out
    masked_kurts = ma.masked_where(np.abs(kurts) > threshold, kurts)
    bin_mask = ma.getmask(masked_kurts)
    
    # flagged_bins is an array that has the frequencies of the channels with kurtosis > threshold NOT masked out
    # flagged_kurts masks the opposite elements as masked_kurts (i.e., it contains all of the kurtoses of the
    # high RFI channels)
    flagged_bins = ma.masked_array(bins, mask=~bin_mask)
    flagged_kurts = ma.masked_array(kurts, mask=~bin_mask)
    
    # The reason I am no longer using ma.count(masked_array) is because there can be an error where masked elements
    # are converted to NaNs.
    print(f'{len(np.where(bin_mask == True)[0])} out of {n_divs} channels flagged as having substantial RFI')
    
    # Get frequency in increasing order, first
    if wf_in.header['foff'] < 0:
        freqs_flipped = wf_in.get_freqs()[::-1]
        freqs = ma.masked_array(freqs_flipped)
    
    # Bin width (for the files I am testing this with):
    # There are 7 Hz in between each entry of freqs, and 32 MHz total
    # Given a total of 4194304 elements in freqs, if there are 256 bins, then each bin spans 16384 elements
    # If each bin spans 16384 elements, then it spans 125 kHz (125 kHz * 256 = 32 MHz)
    
    # Grab the bin width in terms of MHz (for convenience if needed in the future)
    full_freq_range = freqs[-1] - freqs[0]
    # bin_width = full_freq_range / n_divs
    
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
        # kurts: Kurtosis of all frequency bins
        # pows_mean: Time-averaged power of each frequency bin
        # flagged_bins: Channels with high kurtosis (i.e., high RFI)
        # flagged_kurts: Kurtosis of each channel that was flagged as having high RFI
        # masked_kurts: Kurtosis of all 'clean' (low RFI) channels with high RFI channels masked out
        # masked_freqs: List of all frequencies with high RFI channels masked out
        # bin_mask: The mask used to generate masked_kurts -- Masks the frequency *bins*
        # freq_mask: A mask to be used to block out the *actual* frequencies, rather than the frequency *bins*
    return bins, kurts, pows_mean, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask

def write_output_table(wf_in, output_filepath='./', n_divs=256, threshold=50):
    # This function does as it says: It writes the output table. It does so in a .csv format with columns of:
        # (bin_top) Frequency bin tops
        # (bin_bot) Frequency bin bottoms
        # (kurtosis) Kurtosis of each bins
        # (tavg_pwr) Time-averaged power of each bin
    # Inputs:
        # output_filepath: Path to output .csv file
        # wf_in: See get_tavg_kurtosis() function definition
        # n_divs: See get_tavg_kurtosis() function definition
        # threshold: See get_mask_kurtosis() function definition

    # Assign all the base variables and ensure file export path (export_path) is normalized
    export_path = os.path.normpath(output_filepath)
    bins, kurts, pows_mean, flagged_bins, flagged_kurts, masked_kurts, masked_freqs, bin_mask, freq_mask = get_kurtosis(wf_in, n_divs, threshold)

    # Get bin tops
    bin_width = (np.amax(wf_in.get_freqs()) - np.amin(wf_in.get_freqs())) / n_divs
    bin_tops = flagged_bins + bin_width

    # Format flagged_bins into a regular (not masked) numpy array
    flagged_bins = ma.filled(flagged_bins, fill_value=np.NaN)

    # Turns the numpy arrays into pandas dataframes so they can be concatenated and exported
    export_bin_bots = pd.DataFrame(data=flagged_bins, columns=['rfi_bin_bots'])
    export_bin_tops = pd.DataFrame(data=bin_tops, columns=['rfi_bin_tops'])
    export_bin_sk = pd.DataFrame(data=flagged_kurts, columns=['sk'])

    # Concatenate dataframes
    export_concat = pd.concat([export_bin_bots,
                               export_bin_tops,
                               export_bin_sk], axis=1)

    # Sort dataframes by frequency
    export_df = export_concat.sort_values(by=['rfi_bin_bots']).reset_index(drop=True)

    export_df.to_csv(export_path, index=False, na_rep='nan')

    return export_df