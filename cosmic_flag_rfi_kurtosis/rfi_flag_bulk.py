import argparse
import blimpy
from blimpy import calcload, Waterfall
import time
import os
import numpy as np
from analysis_functions import get_exkurt, write_output_table
from plotting_functions import plot_mask_exkurt, plot_tavg_power

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
parser.add_argument('-P', '--plot_types',
		    help='(Optional) List of plot types. tavg_pwr: Time-averaged power spectrum. exkurt: Excess kurtosis vs. frequency plot.',
			choices=['exkurt', 'tavg_pwr'],
			nargs='+',
			# type=str,
			# default=[None, None],
			required=False)
parser.add_argument('--plot_output_path',
		    help='(Optional, unless -P is given) Output path for plots (NOT including file name)',
			type=str,
			required=False)
parser.add_argument('--plot_file_types',
			help='(Optional, unless -P is given) Output file types (can be pdf, png, and/or jpg). Specify as many of these as you want!',
			choices=['png', 'pdf', 'jpg'],
			nargs='+',
			required=False)
# TODO: Figure out how to implement custom plot bounds
# parser.add_argument('--plot_bnds',
# 		    help='(Optional) x and y bounds for plots.',
# 			required=False)

args = parser.parse_args()

filfil = os.path.normpath(args.input_filename)

# Check for valid filterbank file and output location
if (filfil[len(filfil) - 4:] != '.fil'):
	parser.error(f'Input file must be a valid filterbank (.fil) file. Specified input file path: {filfil}')
if (os.path.isfile(filfil) == False):
	parser.error(f'Input file does cannot be found. Specified input file path: {filfil}')

# TODO: Figure out how to check for valid output filepath (might not need to do this)
# Note: Data type checking is automatically covered, so long as you specify a type= ! :D

t0 = time.time()
print('Generating waterfall object...')
ml = blimpy.calcload.calc_max_load(filfil)
wf = Waterfall(os.path.normpath(filfil), max_load = ml)
t1 = time.time()
print(f'Done. Elapsed time: {t1 - t0}')

# TODO: Check to make sure the plot_types given are valid
# if (np.any(args.plot_types) != 'exkurt') & (np.any(args.plot_types) != 'tavg_pwr'):
# 	parser.error(f'No valid inputs given for --plot_types. Valid inputs are exkurt and tavg_pwr. Inputs given: {args.plot_types}')
# if len(args.plot_style) > 2:
# 	parser.error(f'Expected 2 arguments for --plot_style, given {len(args.plot_style)}. Arguments given: {args.plot_style}')

write_output_table(wf_in=wf, output_filepath=args.output_filename, n_divs=args.ndivs, threshold=args.threshold)

# TODO: Check to see if plot output and plot file types are given if -P is specified



# Print argument values in case something isn't working :)
print("\nArgument values:")
print(args.input_filename, type(args.input_filename))
print(args.output_filename, type(args.output_filename))
print(args.threshold, type(args.threshold))
print(args.ndivs, type(args.ndivs))
print(args.plot_types, type(args.plot_types))
print(args.plot_output_path, type(args.plot_output_path))
print(args.plot_file_types, type(args.plot_file_types))


if args.plot_types != None:
	f_min = np.floor(np.amin(wf.get_freqs()))
	f_max = np.ceil(np.amax(wf.get_freqs()))
	p_min = np.floor(np.amin(wf.data))
	p_max = np.ceil(np.amax(wf.data))

	if args.input_filename.rfind('/') != -1:
		name_index_start = args.input_filename.rfind('/') + 1
	else:
		name_index_start = 1
	name_index_end = args.input_filename.rfind('.')
	

	# plot_mask_exkurt(wf_in=wf, n_divs=args.ndivs, threshold=args.threshold,
	# 			unfiltered=True, clean_chnls=True, rfi=True,
	# 			f_start=f_min, f_stop=f_max)
	# TODO: Once you figure out how to do plot boundaries, put in k_start and k_stop! :D
	if np.any(np.asarray(args.plot_types) == "exkurt"):
		exkurt_plot_name = f'plot_exkurt_{args.input_filename[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}.{args.plot_file_types}'
		plot_mask_exkurt(wf_in=wf, n_divs=args.ndivs, threshold=args.threshold,
		     unfiltered=True, clean_chnls=True, rfi=True,
			 f_start=f_min, f_stop=f_max,
			 output_dest=os.path.join(args.plot_output_path, exkurt_plot_name),
			 output_type=args.plot_file_types)
			#  k_start=, k_stop=)

		print(f'exkurt plot generated at {os.path.join(args.plot_output_path, exkurt_plot_name)}')
	
	if np.any(np.asarray(args.plot_types) == 'tavg_pwr'):
		tavg_pwr_plot_name = f'plot_tavg_pwr_{args.input_filename[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}.{args.plot_file_types}'
		plot_tavg_power(wf_in=wf, n_divs=args.ndivs, threshold=args.threshold,
		   f_start=f_min, f_stop=f_max,
		   p_start=p_min, p_stop=p_max, show_filtered_bins=True,
		   output_dest=os.path.join(args.plot_output_path, tavg_pwr_plot_name),
		   output_type=args.plot_file_types)
		
		
		print(f'tavg_pwr plot generated at {os.path.join(args.plot_output_path, tavg_pwr_plot_name)}')


# plot_tavg_power(wf, f_start=f_min, f_stop=f_max, n_divs=args.ndivs, threshold=args.threshold)
# plot_mask_exkurt(wf)