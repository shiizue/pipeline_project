import os
import sys
import argparse


def check_arg(args=None):
    parser = argparse.ArgumentParser(
        description="Function to parse command line arguments"
    )
    parser.add_argument("-i", "--input", help="input file", required=True)
    parser.add_argument("-o", "--output", help="output file", required=True)
    return parser.parse_args(args)


arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

#############################
with open(infile, "r") as f:
    for line in f.readlines():
        accession = line
        fasterq_command = "fasterq-dump " + accession
        os.system(fasterq_command)
