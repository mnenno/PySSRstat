#!/usr/bin/env python3
"""Find accessions of the longest repeats

Usage:
statgetlongest.py repeats_analysis.txt <MISA-file>

Output:
longest-sequences-list.txt

Author: Mario Nenno 2015
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
out_filename = 'longest-sequences-list.txt'

meres = ['Mono', 'Di', 'Tri', 'Tetra', 'Penta', 'Hexa', 'Septa', 'Octa', 'Nona', 'Deca']


def getlineslongest(analysisfile):
    with open(analysisfile, "r") as infile:
        is_start = False
        lines_with_longest = []
        for line in infile:
            line = line.rstrip()
            line_length = len(line)
            if ' Longest motives ' in line:
                is_start = True
                continue
            elif is_start is True and line_length > 0:
                lines_with_longest.append(line)
            elif is_start is True and line_length == 0:
                break
    return lines_with_longest


def formatrepeat(motif, motiflength):
    return "({}){}".format(motif, motiflength)


def create_repeat_names(list_lines):
    # Create list of repeats grouped by mere
    # Data structure is a map of lists keyed by mere
    # 'Di': ['(AGA)16', '(ATA)16', '(ATG)16', '(TCT)16']
    groupedrepeats = {}
    for line in list_lines:
        if line.startswith("Mono"):
            # Mono  A: 47, T: 43, C: 0, G: 10
            # split in motifs
            listofrepeats = []
            monos = line.split("Mono  ")
            motif_lengths = monos[1].split(",")
            for motif_len in motif_lengths:
                motif, motiflength = motif_len.split(": ")
                motif = motif.lstrip()
                if not motiflength == "0":
                    listofrepeats.append(formatrepeat(motif, motiflength))
            groupedrepeats["Mono"] = listofrepeats
        else:
            # Tri     16: AGA, ATA, ATG, CAA, TCT
            listofrepeats = []
            mere_length, motifs = line.split(":")
            mere, motiflength = mere_length.split()
            listofmotifs = motifs.split(", ")
            for motif in listofmotifs:
                motif = motif.lstrip()
                #print(motif)
                listofrepeats.append(formatrepeat(motif, motiflength))
            groupedrepeats[mere] = listofrepeats
    return groupedrepeats


def getsequencelines(filename, groupedrepeats):
    # ID	     SSR nr.	SSR type	SSR	       size	start	end
    # PK00768.1	 1	        p3	      (GGA)5	    15	        336	    350
    validssrtypes = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10']
    ssrtype2mere = {'p1': 'Mono', 'p2': 'Di', 'p3': 'Tri',  'p4': 'Tetra', 'p5': 'Penta',
                    'p6': 'Hexa', 'p7': 'Septa', 'p8': 'Octa', 'p9': 'Nona', 'p10': 'Deca'}
    groupedseqlines = {}
    with open(filename, "r") as infile:
        for line in infile:
            line = line.rstrip()
            items = line.split("\t")
            #     ['PK14952.1', '1', 'p1', '(T)11', '11', '2418', '2428\n']
            #  Idx  0            1    2      3
            ssrtype = items[2]
            if ssrtype in validssrtypes:
                mere = ssrtype2mere[ssrtype]
                ssr = items[3]
                # search SSR in list of longest grouped by mere
                if ssr in groupedrepeats[mere]:
                    if mere in groupedseqlines.keys():
                        listlines = groupedseqlines[mere]
                        listlines.append(items)
                        groupedseqlines[mere] = listlines
                    else:
                        listlines = []
                        listlines.append(items)
                        groupedseqlines[mere] = listlines
        return groupedseqlines


def printgroupedseqlines(grouped, analysisfile, misafile, execution_time):
    with open(out_filename, 'w') as outfile:
        outfile.write("Program: {} {}\n".format(os.path.basename(sys.argv[0]), _version_))
        outfile.write("Date: {} Execution time: {:.2f} sec\n".format(time.strftime("%Y-%m-%d %H:%M"), execution_time))
        outfile.write("Input file of repeat analysis: {}\n".format(analysisfile))
        outfile.write("Input file of sequences: {}\n".format(misafile))
        outfile.write("Output: {}\n".format(out_filename))
        outfile.write("="*50 + "\n")
        for mere in meres:
            if mere in grouped:
                outfile.write("----- {}\n".format(mere))
                for items in grouped[mere]:
                    #     ['PK14952.1', '1', 'p1', '(T)11', '11', '2418', '2428\n']
                    #  Idx  0            1    2      3       4
                    seqline = "{} {}: {}\n".format(items[3], items[4], items[0])
                    outfile.write(seqline)


def main(analysisfile, misafile):
    # Start timing
    start_time = time.time()

    lines_with_longest = getlineslongest(analysisfile)
    groupedrepeats = create_repeat_names(lines_with_longest)
    groupedseqlines = getsequencelines(misafile, groupedrepeats)

    # End timing
    execution_time = time.time() - start_time

    # write to output
    printgroupedseqlines(groupedseqlines, analysisfile, misafile, execution_time)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("analysisfile", help="file with output from statistics_misa")
    parser.add_argument("misafile", help="file with sequences and repeats")
    args = parser.parse_args()
    if args.analysisfile and args.misafile:
        main(args.analysisfile, args.misafile)