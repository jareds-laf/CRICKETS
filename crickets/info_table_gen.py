from analysis import normalize_path, create_info_table
import argparse

parser = argparse.ArgumentParser(
                    description='Create table of frequency and time-averaged power for input filterbank file.')

parser.add_argument('--wf_path', '-w',
		    help='(Required) Path to input filterbank file (including file name).', 
			type=str,
			required=True)
parser.add_argument('--save_loc', '-s',
		    help='(Required) Directory to save info table to.',
			type=str,
			required=True)

args = parser.parse_args()

wf_path = normalize_path(args.wf_path)
saveloc = normalize_path(args.save_loc)

create_info_table(wf_file_full=wf_path, saveloc=saveloc)