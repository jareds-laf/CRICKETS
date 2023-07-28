from analysis import normalize_path, create_info_table_dir, create_info_table
import argparse
# import configparser
import logging
import sys

# log_format = '%(asctime)s - %(levelname)s - %(message)s'
logger = logging.getLogger('info_table_gen')

# logger.propagate = False
# stdout = logging.StreamHandler(stream=sys.stdout)

# # Create formatter
# formatter = logging.Formatter('\n%(asctime)s :: %(name)s :: %(levelname)s :: %(message)s\n')

# # Add formatter to ch
# stdout.setFormatter(formatter)

# # Add ch to logger
# logger.addHandler(stdout)

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

args = parser.parse_args()

wf_path = normalize_path(args.wf_path)
itloc = normalize_path(args.info_table_loc)

create_info_table(wf_dir=wf_path, saveloc=itloc)

"""
An interesting way to save the info table location
to pass to crickets_runner.py...

This would be in place instead of specifying the location twice
in the two commands you have to put.
"""

# config = configparser.ConfigParser()
# config.add_section('info_table')

# config['info_table']['location'] = itloc

# with open('crickets/config.ini', 'w') as configfile:
#     config.write(configfile)

# logger.debug("test2")

# create_info_table(wf_file_full=wf_path, saveloc=itloc)

# logger.debug("test3")
# print(f"\nwf_path: {wf_path}\n")