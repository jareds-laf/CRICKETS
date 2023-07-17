import argparse
import blimpy
from blimpy import calcload, Waterfall
import time
import os
import numpy as np
from analysis_functions import get_kurtosis, write_output_table

# def file_checker(str):
# 	if str[len(str) - 4:] != '.fil':
# 		raise argparse.ArgumentTypeError(f'Input file must be a filterbank (.fil) file. Specified input file path: {str}')

parser = argparse.ArgumentParser(
                    description='Flag RFI heavy frequency channels based on the kurtosis of each channel.')

parser.add_argument('input_filename',
		    help='(Required) Path to input filterbank file.', 
			type=str)
parser.add_argument('output_filename',
		    help='(Required) Path to output csv file.',
			type=str)
parser.add_argument('-T', '--threshold',
		    help='(Required) Minimum value of kurtosis used to flag channels with significant RFI.',
			type=float,
			default=5,
			required=True)
parser.add_argument('-N', '--ndivs',
		    help='(Required) Number of frequency bins',
			type=int,
		    default=256,
			required=True)
parser.add_argument('-P', '--plot_types',
		    help='(Optional) List of plot types. tavg_pwr: Time-averaged power spectrum. sk: Spectral kurtosis plot.',
			nargs='*',
			type=str,
			default=[None, None],
			required=False)

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


write_output_table(wf_in=wf, filepath='~/kurtosis_output.csv', n_divs=args.ndivs, threshold=args.threshold)

# Print argument values in case something isn't working :)
print("\nArgument values:")
print(args.input_filename, type(args.input_filename))
print(args.output_filename, type(args.output_filename))
print(args.threshold, type(args.threshold))
print(args.ndivs, type(args.ndivs))
print(args.plot_types, type(args.plot_types))