#!/usr/bin/env python3
"""Filter the MISA file by minimum and maximum repeat length

Usage:
filterrepeatsmisa.py <MISA-file> <min> <max> <sortmode>

Output:
filtered-repeats-sequence-list.txt

Author: Mario Nenno
Version: 2015-08-31
Copyright: see file LICENCE.txt

"""

import os
import sys
import time
import argparse

# program version
_version_ = '1.0'

# output file
out_filename = 'filtered-repeats-sequence-list.txt'

meres = ['Di', 'Tri', 'Tetra', 'Penta', 'Hexa']


def getsequencelines(filename, minlength, maxlength):
    validssrtypes = ['p2', 'p3', 'p4', 'p5', 'p6']
    if args.imperfect:
        validssrtypes.append('c')

    ssrtype2mere = {'p2': 'Di', 'p3': 'Tri',  'p4': 'Tetra', 'p5': 'Penta', 'p6': 'Hexa'}
    if args.imperfect:
        ssrtype2mere['c'] = 'combined'

    groupedseqlines = {}
    with open(filename, "r") as infile:
        for line in infile:
            line = line.rstrip()
            items = line.split("\t")
            #        ID	    SSR nr.	SSR type SSR	size	start	end
            #     ['PK14952.1', '1', 'p1', '(T)11', '11', '2418', '2428\n']
            #  Idx  0            1    2      3       4     5
            ssrtype = items[2]
            if ssrtype in validssrtypes:
                mere = ssrtype2mere[ssrtype]
                repeatlen = int(items[4])
                # search SSR in list of longest grouped by mere
                if minlength <= repeatlen <= maxlength:
                    # print(items)
                    if mere in groupedseqlines.keys():
                        listlines = groupedseqlines[mere]
                        listlines.append(items)
                        groupedseqlines[mere] = listlines
                    else:
                        listlines = []
                        listlines.append(items)
                        groupedseqlines[mere] = listlines
        return groupedseqlines


def printgroupedseqlines(grouped, misafile, execution_time, minlength, maxlength, sortmode):
    with open(out_filename, 'w') as outfile:
        outfile.write("Program: {} {}\n".format(os.path.basename(sys.argv[0]), _version_))
        outfile.write("Date: {} Execution time: {:.2f} sec\n".format(time.strftime("%Y-%m-%d %H:%M"), execution_time))
        outfile.write("Input file of sequences: {}\n".format(misafile))
        outfile.write("Repeat length min: {}, max: {}, sort sequence by: {} length\n".format(
            minlength, maxlength, sortmode))
        outfile.write("Output: {} Execution time: {:.2f} sec\n".format(out_filename, execution_time))
        outfile.write("="*70 + "\n")
        # --- item structure
        #     ['PK14952.1', '1', 'p1', '(T)11', '11', '2418', '2428\n']
        #  Idx  0            1    2      3       4     5       6

        # copy into a new list without grouping by mere
        listofall = []
        for mere in meres:
            if mere in grouped.keys():
                for items in grouped[mere]:
                    listofall.append(items)

        # find sequences with multiple SSR repeats
        seqwithmultirepeats = {}
        for items in listofall:
            seqid = items[0]
            if seqid in seqwithmultirepeats.keys():
                num = seqwithmultirepeats[seqid]
                num += 1
                seqwithmultirepeats[seqid] = num
            else:
                seqwithmultirepeats[seqid] = 1

        seqlineformat = "{:12} {}: {}\n"
        if "motif" == sortmode:
            for mere in meres:
                outfile.write("----- {}\n".format(mere))
                if mere in grouped.keys():
                    # print(mere, grouped[mere])
                    for items in grouped[mere]:
                        seqline = seqlineformat.format(items[3], items[4], items[0])
                        outfile.write(seqline)
        elif "repeat" == sortmode:
            # print("contains %d" % len(listofall))
            # to sort by length loop over range
            for l in range(minlength, maxlength+1):
                for items in listofall:
                    repeatlength = int(items[4])
                    if l == repeatlength:
                        seqline = "\t".join(items)+"\n"
                        outfile.write(seqline)
        if len(seqwithmultirepeats) > 0:
            found = 0
            for seq in seqwithmultirepeats:
                if seqwithmultirepeats[seq] > 1:
                    found += 1
                    if found == 1:
                        outfile.write("\n----- Sequences with multiple repeats -----\n")
                    outfile.write("{}: {}\n".format(seqwithmultirepeats[seq], seq))


def main(misafile, min_length, max_length, sortmode):
    # if optional arguments imperfect, add the combined mere to the list of meres
    if args.imperfect:
        meres.append('combined')

    # Start timing
    start_time = time.time()

    # for calculations
    minlength = int(min_length)
    maxlength = int(max_length)

    # filter the file for minimum and maximum length
    groupedseqlines = getsequencelines(misafile, minlength, maxlength)

    # End timing
    execution_time = time.time() - start_time

    # write results to output
    printgroupedseqlines(groupedseqlines, misafile, execution_time, minlength, maxlength, sortmode)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("misafile", help="file with sequences and repeats")
    parser.add_argument("min", help="minimum length of repeat")
    parser.add_argument("max", help="maximum length of repeat")
    parser.add_argument("sortmode", help="how to sort the sequences", choices=['motif', 'repeat'])
    parser.add_argument("-i", "--imperfect", help="include imperfect SSRs", action="store_true")
    args = parser.parse_args()
    if args.misafile and args.min and args.max and args.sortmode:
        main(args.misafile, args.min, args.max, args.sortmode)