import argparse
# from analysis_functions import get_mask_kurtosis

parser = argparse.ArgumentParser(
                    description='Flag RFI heavy frequency channels based on the kurtosis of each channel.')

parser.add_argument('filename', help='Path to input filterbank file')
parser.add_argument('--ohi', help='Print "Hello world!"!', action='store_true', required=False)
# parser.add_argument('filename', help='Path to input filterbank file')

args = parser.parse_args()



print("Argument values:")
print(args.filename, type(args.filename))
print(args.ohi, type(args.ohi))



# parser.add_argument('filename')           # positional argument
# parser.add_argument('-c', '--count')      # option that takes a value
# parser.add_argument('-v', '--verbose',
#                     action='store_true')  # on/off flag
