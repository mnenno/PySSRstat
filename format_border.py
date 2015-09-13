#!/usr/bin/env python3
"""Format the file border.txt with spaces or tabs

With the optional parameter -idt or --idtrunc the id (position 0) can be
truncated at the first underscore character

Usage: format_border.py border.txt tab|space [-idt|--idtrunc]

Example of line in border.txt

scaffold756_39.3	41	p2	(CT)15	30	54494	54523
IDX 0               1   2     3     4     5       6


"""
import os
import sys
import time
import argparse


def main(borderfile, delimiter, idtrunc):
    if 'space' == delimiter:
        out_filename = 'border-space.txt'
        line_format = "{0:<30} {1:>4} {2:>7} {3:>7}"
    else:
        out_filename = 'border-tab.txt'
        line_format = "{0}\t{1}\t{2}\t{3}"

    with open(borderfile, 'r') as borderfile:
        with open(out_filename, 'w') as outfile:
            for line in borderfile:
                line = line.rstrip()
                if "\t" in line:
                    # remove double tabs
                    line = line.replace("\t\t", "\t")
                    cols = line.split("\t")
                    seq_id = cols[0]
                    # optionally trunc id, e.g. name_1 to name
                    if idtrunc:
                        seq_id = seq_id.split('_', 1)[0]
                    line_formated = line_format.format(
                        seq_id, cols[3], cols[5], cols[6])
                    outfile.write(line_formated+"\n")
                else:
                    outfile.write(line+"\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("borderfile", help="file with sequences and repeats")
    parser.add_argument("delimiter", help="how to sort the sequences", choices=['space', 'tab'])
    parser.add_argument("-idt", "--idtrunc", help="trunc id after underscore", action="store_true")
    args = parser.parse_args()
    if args.borderfile and args.delimiter:
        main(args.borderfile, args.delimiter, args.idtrunc)