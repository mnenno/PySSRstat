#!/usr/bin/env python3
"""Stat about imperfect repeats in MISA-file

Usage:
imperfect.py <MISA-file>

Output:
imperfect.txt

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
out_filename = "imperfect.txt"


def get_imperfect(filename):
    """ Read MISA-file as input
    """
    # ID	     SSR nr.	SSR type	SSR	   size   start  end
    # PK00768.1	 1	        p3	      (GGA)5   15    336    350

    # SSR types: p = perfect, c = imperfect, c* = compound
    validssrtypes = ['c', 'c*']
    list_imperfect = []
    list_ssr = []
    with open(filename, "r") as infile:
        for line in infile:
            line = line.rstrip()
            items = line.split("\t")
            #     ['PK14952.1', '1', 'p1', '(T)11', '11', '2418', '2428\n']
            #  Idx  0            1    2      3       4
            ssrtype = items[2]
            ssr = items[3]
            if ssrtype in validssrtypes:
                # print(items)
                if ssr not in list_ssr:
                    list_ssr.append(ssr)
                    list_imperfect.append(items)
    return list_imperfect


def analyze(list_items):
    """ Count the different types of imperfect repeats
    """
    analysis = {}
    count = 0
    num_imperfect = 0
    num_compound = 0
    listoflength = []
    #     ['PK14952.1', '1', 'p1', '(T)11', '11', '2418', '2428\n']
    #  Idx  0            1    2      3       4
    for items in list_items:
        count += 1
        l = int(items[4])
        # Count by ssr type
        ssrtype = items[2]
        if ssrtype == 'c':
            num_imperfect += 1
        elif ssrtype == 'c*':
            num_compound += 1

        if l not in listoflength:
            listoflength.append(l)
    # Note: is an unordered list
    analysis['num_imperfect'] = num_imperfect
    analysis['num_compound'] = num_compound
    analysis['total'] = count
    analysis['listoflength'] = listoflength
    return analysis


def extract(list_imperfect, analysis):
    # Order listoflength desc
    extracted = []
    listoflength = analysis['listoflength']
    listoflength.sort(reverse=True)
    for l in listoflength:
        # print('Length: ', l)
        for items in list_imperfect:
            if l == int(items[4]):
                # print(items)
                extracted.append(items)
    # print('Number of length', len(listoflength))
    return extracted


def save(misafile, extracted, analysis, start_time):
    """ Save into output file
    """
    execution_time = time.time() - start_time
    with open(out_filename, 'w') as outfile:
        outfile.write("Program: {} {}\n".format(os.path.basename(sys.argv[0]), _version_))
        outfile.write("Date: {}, duration: {:.2f} sec\n".format(time.strftime("%Y-%m-%d %H:%M"), execution_time))
        outfile.write("Input, MISA file: {}\n".format(misafile))
        outfile.write("Imperfect: {}, compound: {}, Total: {}\n".format(
            analysis['num_imperfect'], analysis['num_compound'], analysis['total']))
        outfile.write("Order: longest first, shortest last\n")
        outfile.write("Output: {}\n".format(out_filename))
        outfile.write("=" * 70 + "\n")
        for items in extracted:
            seqline = "{}\t{}\t{}\t{}\n".format(items[0], items[2], items[4], items[3])
            outfile.write(seqline)


def main(misafile):
    # Start timing
    start_time = time.time()

    imperfect = get_imperfect(misafile)
    analysis = analyze(imperfect)
    extracted = extract(imperfect, analysis)
    save(misafile, extracted, analysis, start_time)


if __name__ == '__main__':
    # command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("misafile", help="file with sequences and repeats")
    args = parser.parse_args()
    if args.misafile:
        main(args.misafile)
