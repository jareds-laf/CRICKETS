import argparse
from blimpy import calcload, Waterfall
import time
import os
import numpy as np
from analysis import get_exkurt, write_output_table, normalize_path, plot_exkurt, plot_tavg_power

parser = argparse.ArgumentParser(
                    description='Flag RFI heavy frequency channels based on the excess kurtosis of each channel.')

parser.add_argument('--input_file',
		    help='(Required) Path to input filterbank file (including file name).', 
			type=str,
			required=True)
parser.add_argument('--output_file',
		    help='(Required) Path to output csv file (optionally including file name).',
			type=str,
			required=True)
parser.add_argument('--threshold', '-t',
		    help='(Required) Minimum value of excess kurtosis used to flag channels with significant RFI. Can be any decimal number.',
			type=float,
			default=5,
			required=True)
parser.add_argument('--ndivs', '-n',
		    help='(Required) Number of frequency bins to split waterfall object into. Can be any integer.',
			type=int,
		    default=256,
			required=True)

# TODO: Better implementation of the plotting arguments!
parser.add_argument('--plot', '-p',
		    help='(Optional) Choose whether or not to generate time-averaged power spectrum and excess kurtosis vs. frequency plots. Give output path for plots here (NOT including file name).',
			# action='store_true',
			type=str,
			required=False)
parser.add_argument('--plot_file_types',
			help='(Optional, unless -p is given) Output file types (can be pdf, png, and/or jpg). Specify as many of these as you want!',
			choices=['png', 'pdf', 'jpg'],
			nargs='+',
			required=False)
# parser.add_argument('--plot_output_path',
# 		    help='(Optional, unless -p is given) ',
# 			type=str,
# 			required=False)

parser.add_argument('--verbose', '-v',
		    help='(Optional) Print more information about the input variables and the processes currently running.',
			action='store_true',
			required=False)
# parser.add_argument('-p', '--plot_types',
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

# Normalize path to input filterbank file and path to output csv file
filfil = normalize_path(args.input_file)
out_path = normalize_path(args.output_file)

# Check for valid filterbank file and output location
if (filfil[len(filfil) - 4:] != '.fil'):
	parser.error(f'Input file must be a valid filterbank (.fil) file. Specified input file path: {filfil}')
if (os.path.isfile(filfil) == False):
	parser.error(f'Input file cannot be found. Specified input file path: {filfil}')

# TODO: Figure out how to check for valid output filepath (might not need to do this)
# Note: Data type checking is automatically covered, so long as you specify a type= ! :D

if os.path.isdir(out_path):
	if out_path[-1] != '/':
		out_path += '/'
		print("The output file path is a directory.")
		output_file_loc = out_path + f"crickets_{args.input_file[args.input_file.rfind('/')+1:len(args.input_file)-4]}_{args.ndivs}_{args.threshold}.csv"
	else:
		print("The output file path is a directory.")
		output_file_loc = out_path + f"crickets_{args.input_file[args.input_file.rfind('/')+1:len(args.input_file)-4]}_{args.ndivs}_{args.threshold}.csv"
if (os.path.isfile(out_path)) | (out_path[-4:] == '.csv'):
	print("The output file is a file.")
	output_file_loc = out_path

# Generate waterfall object
t0 = time.time()
print('\nGenerating waterfall object...')
ml = calcload.calc_max_load(filfil)
wf = Waterfall(os.path.normpath(filfil), max_load = ml)
t1 = time.time()
print(f'Done. Elapsed time: {t1 - t0}')

# Run analysis code and generate the output table
write_output_table(wf_in=wf, output_filepath=output_file_loc, n_divs=args.ndivs, threshold=args.threshold)
print(f'Output table generated at {output_file_loc}')
# TODO: Check to see if plot output and plot file types are given if -p is specified

# Plotting code
if args.plot:
		f_min = np.floor(np.amin(wf.get_freqs()))
		f_max = np.ceil(np.amax(wf.get_freqs()))
		p_min = np.floor(np.amin(wf.data))
		p_max = np.ceil(np.amax(wf.data))

		if args.input_file.rfind('/') != -1:
			name_index_start = args.input_file.rfind('/') + 1
		else:
			name_index_start = 1
		name_index_end = args.input_file.rfind('.')
		
		# Excess kurtosis plot
		print("\nGenerating exkurt plot...")
		exkurt_plot_name = f'plot_exkurt_{args.input_file[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}'
		plot_names = []
		for i in args.plot_file_types:
			exec(f"exkurt_plot_name_{i} = f'plot_exkurt_{args.input_file[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}.{i}'")
			exec(f"plot_names.append(exkurt_plot_name_{i})")
			# exec(f"print(exkurt_plot_name_{i})")

		plot_exkurt(wf_in=wf, n_divs=args.ndivs, threshold=args.threshold,
				unfiltered=True, clean_chnls=True, rfi=True,
				f_start=f_min, f_stop=f_max,
				output_dest=os.path.join(args.plot, exkurt_plot_name),
				output_type=args.plot_file_types)
			# TODO: Once you figure out how to do plot boundaries, put in k_start and k_stop! :D
		# for filetype in args.plot_file_types:
		
		for i in plot_names:
			exec(f"print(f'exkurt plot generated at {os.path.join(args.plot, i)}')")

		# Time-averaged power spectrum plot
		print("\nGenerating tavg_pwr plot...")
		tavg_pwr_plot_name = f'plot_tavg_pwr_{args.input_file[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}'
		plot_names = []
		for i in args.plot_file_types:
			exec(f"tavg_pwr_plot_name_{i} = f'plot_tavg_pwr_{args.input_file[name_index_start:name_index_end]}_{args.ndivs}_{args.threshold}.{i}'")
			exec(f"plot_names.append(tavg_pwr_plot_name_{i})")

		plot_tavg_power(wf_in=wf, n_divs=args.ndivs, threshold=args.threshold,
				f_start=f_min, f_stop=f_max,
				p_start=p_min, p_stop=p_max, show_filtered_bins=True,
				output_dest=os.path.join(args.plot, tavg_pwr_plot_name),
				output_type=args.plot_file_types)
		
		for i in plot_names:
			exec(f"print(f'tavg_pwr plot generated at {os.path.join(args.plot, i)}')")

# Verbose outputs
if args.verbose:
	# Print argument values in case something isn't working :)
	print("\nArgument values:")
	print(f"Input filename: {args.input_file}, {type(args.input_file)}")
	print(f"Output filename: {args.output_file}, {type(args.output_file)}")
	print(f"Threshold: {args.threshold}, {type(args.threshold)}")
	print(f"Number of bins: {args.ndivs}, {type(args.ndivs)}")
	print(f"Plot (bool): {args.plot}, {type(args.plot)}")
	# print(f"Plot output path: {args.plot_output_path}, {type(args.plot_output_path)}")
	print(f"Plot file type(s): {args.plot_file_types}, {type(args.plot_file_types)}")