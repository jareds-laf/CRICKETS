from analysis import normalize_path, create_info_table_dir, create_info_table
import argparse
import configparser
# import logging
# import sys

# log_format = '%(asctime)s - %(levelname)s - %(message)s'

parser = argparse.ArgumentParser(
                    description='Create table of frequency and time-averaged power for input filterbank file.')

parser.add_argument('--wf_path', '-w',
		    help='(Required) Path to directory containing all input filterbank files.', 
			type=str,
			required=True)
parser.add_argument('--info_table_loc', '-i',
		    help='(Required) Directory to save info tables to.',
			type=str,
			required=True)
# parser.add_argument('--verbose', '-v',
# 		    help='(Optional) Print more information about the input variables and the processes currently running.',
# 			action='store_true',
# 			required=False)

args = parser.parse_args()

# Trying to add a nice verbose feature failed :(
# if args.verbose:
# 	print("HELLO THIS IS VERBOSE")
# 	# logger = logging.getLogger()
# 	# logger.setLevel(logging.DEBUG)

# 	# handler = logging.StreamHandler(sys.stdout)
# 	# handler.setLevel(logging.DEBUG)
# 	# formatter = logging.Formatter(log_format)
# 	# handler.setFormatter(formatter)
# 	logging.basicConfig(level=logging.DEBUG, format=log_format)#, stream=sys.stdout)
# 	logger = logging.getLogger()

# 	# logger.addHandler(handler)
# else:
# 	logging.basicConfig(level=logging.INFO, format=log_format, stream=sys.stdout)
# 	logger = logging.getLogger()

# logger.debug("test1")

wf_path = normalize_path(args.wf_path)
itloc = normalize_path(args.info_table_loc)

# An interesting way to save the info table location to pass to crickets_runner.py...
# config = configparser.ConfigParser()
# config.add_section('info_table')

# config['info_table']['location'] = itloc

# with open('crickets/config.ini', 'w') as configfile:
#     config.write(configfile)

# logger.debug("test2")

# create_info_table(wf_file_full=wf_path, saveloc=itloc)

# logger.debug("test3")
# print(f"\nwf_path: {wf_path}\n")
create_info_table_dir(wf_dir=wf_path, saveloc=itloc)