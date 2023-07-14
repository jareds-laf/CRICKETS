import argparse
#import blimpy
#from blimpy import calcload, Waterfall
import time
import os

t0 = time.time()
print(f"about to import numpy")


import numpy


t1 = time.time()
print(f"Done importing numpy: {t1 - t0} s elapsed")
# from analysis_functions import get_mask_kurtosis

def file_checker(str):
	if str[len(str) - 4:] != '.fil':
		raise argparse.ArgumentTypeError(f'Input file must be a filterbank (.fil) file. Specified input file path: {str}')

parser = argparse.ArgumentParser(
                    description='Flag RFI heavy frequency channels based on the kurtosis of each channel.')

parser.add_argument('filename', help='Path to input filterbank file', type=str)
parser.add_argument('-T', help='Minimum value of kurtosis used to flag RFI', type=float, required=True)
parser.add_argument('-N', help='Number of frequency bins', type=int, required=True)
# parser.add_argument('--ohi', help='Print "Hello world!"!', action='store_true', required=False)
# parser.add_argument('filename', help='Path to input filterbank file')

args = parser.parse_args()

filename = os.path.normpath(args.filename)

if filename[len(filename) - 4:] != '.fil':
	parser.error(f'Input file must be a filterbank (.fil) file. Specified input file path: {filename}')


#t0 = time.time()
#print('Generating waterfall object...')
#ml = blimpy.calcload.calc_max_load(filename)
#wf = Waterfall(os.path.normpath(file), max_load = ml)
#t1 = time.time()
#print(f'Done. Elapsed time: {t1 - t0}')



print("Argument values:")
print(args.filename, type(args.filename))
print(args.T, type(args.T))
print(args.N, type(args.N))




# parser.add_argument('filename')           # positional argument
# parser.add_argument('-c', '--count')      # option that takes a value
# parser.add_argument('-v', '--verbose',
#                     action='store_true')  # on/off flag
