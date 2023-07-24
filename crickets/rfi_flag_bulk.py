import argparse
import blimpy
from blimpy import calcload, Waterfall
import time
import os
import numpy as np
from crickets.analysis import get_exkurt, write_output_table
from crickets.plotting import plot_mask_exkurt, plot_tavg_power

# def file_checker(str):
# 	if str[len(str) - 4:] != '.fil':
# 		raise argparse.ArgumentTypeError(f'Input file must be a filterbank (.fil) file. Specified input file path: {str}')

parser = argparse.ArgumentParser(
                    description='Flag RFI heavy frequency channels based on the excess kurtosis of each channel.')

parser.add_argument('--input_filename',
		    help='(Required) Path to input filterbank file (including file name).', 
			type=str,
			required=True)
parser.add_argument('--output_filename',
		    help='(Required) Path to output csv file (including file name).',
			type=str,
			required=True)
parser.add_argument('-T', '--threshold',
		    help='(Required) Minimum value of excess kurtosis used to flag channels with significant RFI.',
			type=float,
			default=5,
			required=True)
parser.add_argument('-N', '--ndivs',
		    help='(Required) Number of frequency bins to split waterfall object into.',
			type=int,
		    default=256,
			required=True)


# TODO: Better implementation of the plotting arguments!
parser.add_argument('-P', '--plot',
		    help='(Optional) Choose whether or not to generate time-averaged power spectrum and excess kurtosis vs. frequency plots.',
			action='store_true',
			required=False)
parser.add_argument('--plot_file_types',
			help='(Optional, unless -P is given) Output file types (can be pdf, png, and/or jpg). Specify as many of these as you want!',
			choices=['png', 'pdf', 'jpg'],
			nargs='+',
			required=False)
parser.add_argument('--plot_output_path',
		    help='(Optional, unless -P is given) Output path for plots (NOT including file name)',
			type=str,
			required=False)

# parser.add_argument('-P', '--plot_types',
# 		    help='(Optional) List of plot types. tavg_pwr: Time-averaged power spectrum. exkurt: Excess kurtosis vs. frequency plot.',
# 			choices=['exkurt', 'tavg_pwr'],
# 			nargs='+',
# 			# type=str,
# 			# default=[None, None],
# 			required=False)
# TODO: Figure out how to implement custom plot bounds
# parser.add_argument('--plot_bnds',
# 		    help='(Optional) x and y bounds for plots.',
# 			required=False)

args = parser.parse_args()

# Normalize path to input filterbank file
filfil = os.path.normpath(args.input_filename)

# Check for valid filterbank file and output location
if (filfil[len(filfil) - 4:] != '.fil'):
	parser.error(f'Input file must be a valid filterbank (.fil) file. Specified input file path: {filfil}')
if (os.path.isfile(filfil) == False):
	parser.error(f'Input file does cannot be found. Specified input file path: {filfil}')

# TODO: Figure out how to check for valid output filepath (might not need to do this)
# Note: Data type checking is automatically covered, so long as you specify a type= ! :D

# Generate waterfall object
t0 = time.time()
print('\nGenerating waterfall object...')
ml = blimpy.calcload.calc_max_load(filfil)
wf = Waterfall(os.path.normpath(filfil), max_load = ml)
t1 = time.time()
print(f'Done. Elapsed time: {t1 - t0}')

# Run analysis code and generate the output table
write_output_table(wf_in=wf, output_filepath=args.output_filename, n_divs=args.ndivs, threshold=args.threshold)

# TODO: Check to see if plot output and plot file types are given if -P is specified


# Print argument values in case something isn't working :)
print("\nArgument values:")
print(f"Input filename: {args.input_filename}, {type(args.input_filename)}")
print(f"Output filename: {args.output_filename}, {type(args.output_filename)}")
print(f"Threshold: {args.threshold}, {type(args.threshold)}")
print(f"Number of bins: {args.ndivs}, {type(args.ndivs)}")
print(f"Plot (bool): {args.plot}, {type(args.plot)}")
print(f"Plot output path: {args.plot_output_path}, {type(args.plot_output_path)}")
print(f"Plot file type(s): {args.plot_file_types}, {type(args.plot_file_types)}")

# Plotting code
if args.plot:
		f_min = np.floor(np.amin(wf.get_freqs()))
		f_max = np.ceil(np.amax(wf.get_freqs()))
		p_min = np.floor(np.amin(wf.data))
		p_max = np.ceil(np.amax(wf.data))

		if args.input_filename.rfind('/') != -1:
			name_index_start = args.input_filename.rfind('/') + 1
		else:
			name_index_start = 1
		name_index_end = args.input_filename.rfind('.')
		
		# Excess kurtosis plot
		print("\nGenerating exkurt plot...")
		exkurt_plot_name = f'plot_exkurt_{args.input_filename[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}'
		plot_names = []
		for i in args.plot_file_types:
			exec(f"exkurt_plot_name_{i} = f'plot_exkurt_{args.input_filename[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}.{i}'")
			exec(f"plot_names.append(exkurt_plot_name_{i})")
			# exec(f"print(exkurt_plot_name_{i})")

		plot_mask_exkurt(wf_in=wf, n_divs=args.ndivs, threshold=args.threshold,
				unfiltered=True, clean_chnls=True, rfi=True,
				f_start=f_min, f_stop=f_max,
				output_dest=os.path.join(args.plot_output_path, exkurt_plot_name),
				output_type=args.plot_file_types)
			# TODO: Once you figure out how to do plot boundaries, put in k_start and k_stop! :D
		# for filetype in args.plot_file_types:
		
		for i in plot_names:
			exec(f"print(f'exkurt plot generated at {os.path.join(args.plot_output_path, i)}')")

		# Time-averaged power spectrum plot
		print("\nGenerating tavg_pwr plot...")
		tavg_pwr_plot_name = f'plot_tavg_pwr_{args.input_filename[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}'
		plot_names = []
		for i in args.plot_file_types:
			exec(f"tavg_pwr_plot_name_{i} = f'plot_tavg_pwr_{args.input_filename[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}.{i}'")
			exec(f"plot_names.append(tavg_pwr_plot_name_{i})")

		plot_tavg_power(wf_in=wf, n_divs=args.ndivs, threshold=args.threshold,
				f_start=f_min, f_stop=f_max,
				p_start=p_min, p_stop=p_max, show_filtered_bins=True,
				output_dest=os.path.join(args.plot_output_path, tavg_pwr_plot_name),
				output_type=args.plot_file_types)
		
		for i in plot_names:
			exec(f"print(f'tavg_pwr plot generated at {os.path.join(args.plot_output_path, i)}')")