#!/usr/bin/env python3
"""Extract the previously filtered accessions from database file

With the optional parameter for border-mode filter it extracts only sequences that have a 'border'
of n bp up- and downstream of SSR and this border should not contain Ns

Usage:
getsequences.py repeats-sequence-list.txt <db-file> [-b|--border <border-length-in-bp>]

Output:
index.txt, repeats-sequences.fas

or with optional parameter -b nnn
index.txt, repeats-sequences-border.fas, border.txt


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

# name of output files
out_filename_info = "getsequences-info.txt"
out_filename_seq = "repeats-sequences.fas"
out_filename_seq_border = "repeats-sequences-border.fas"
out_filename_border = "border.txt"


def readseqidfrominfile(filename, border):
    repeats = []
    found_seqids = []
    border_len = 0
    if border:
        border_len = int(border)
    # read file with sequence ids in memory as list
    fhseqids = open(filename, "r")
    lines = fhseqids.readlines()
    fhseqids.close()
    # loop over list and extract sequence ids
    is_start = False
    linenum = 0
    for line in lines:
        linenum += 1
        line = line.rstrip()
        ll = len(line)
        if line.startswith("====="):
            is_start = True
            continue
        elif is_start and ll > 0:
            # Id   SSR nr.  type SSR   size start end
            # PK13324.1	1   p2	(TC)15	30	29	 58
            # 0         1   2   3        4   5    6
            # split by tab and use first element
            items = line.split("\t")
            # do not add duplicate sequence ids
            if items[0] in found_seqids:
                continue
            else:
                found_seqids.append(items[0])
            # optionally, test star poisition (downstream border)
            if border:
                start = int(items[5])
                if start >= border_len:
                    repeats.append(items)
            else:
                repeats.append(items)
        elif is_start and ll == 0:
            # is_start = False
            break
        else:
            continue
    return repeats


def create_load_index(seq):
    indexfile = "index.txt"
    if not os.path.isfile(indexfile):
        i = 0
        with open(indexfile, "w") as file:
            for line in seq:
                i += 1
                if line.startswith(">"):
                    line = line.rstrip()
                    # for compatibility change seq id like misa
                    # replace space with underscore
                    seq_id = line[1:].replace(' ', '_')
                    file.write(seq_id + "\t" + str(i) + "\n")
            # add end of last sequence
            file.write("end"+"\t"+str(i+1)+"\n")
    # load index
    fo = open(indexfile, "r")
    index = fo.readlines()
    fo.close()
    return index


def extract_print_seq(accession_repeats, sequences, border, start_time):
    # db file must be in fasta format
    num = 0
    repeats_border_ok = []
    border_len = 0
    if border:
        border_len = int(border)
    # switch filename with sequences
    if border:
        outfilename = out_filename_seq_border
    else:
        outfilename = out_filename_seq
    with open(outfilename, "w") as seq_file:
        # load/create index file for sequences
        index = create_load_index(sequences)
        # search for: PK13324.1
        # in PK13324.1  1234
        i = 0
        for items in accession_repeats:
            i += 1
            idfrominput = items[0]
            idx = 0
            row_number = 0
            for entry in index:
                id_found, row_number = entry.split("\t")
                if id_found == idfrominput:
                    break
                idx += 1
            start = int(row_number)
            # end is the row_number - 1 of the next entry in the index
            next_entry = index[idx+1].rstrip()
            dummy, row_number = next_entry.split("\t")
            end = int(row_number)-1
            # reformat to one sequence a line
            tmp_list = []
            for l in sequences[start:end]:
                tmp_list.append(l.rstrip())
            sequence = "".join(tmp_list)
            if not border:
                # write to file
                seq_file.write(">{}\n".format(idfrominput))
                seq_file.write(sequence + "\n")
                num += 1
            else:
                # Id   SSR nr.  type SSR   size start end
                # PK13324.1	1   p2	(TC)15	30	29	 58
                # 0         1   2   3        4   5    6
                # --- check for N in border upstream
                start_repeat = int(items[5])
                border_up = sequence[(start_repeat - border_len - 1):(start_repeat-1)]
                # print(items[0], border_up, "\n")
                if not 'N' in border_up:
                    # --- check for boder downstream
                    # 1) size: rest >= border
                    end_repeat = int(items[6])
                    rest_after_end = len(sequence) - end_repeat
                    if rest_after_end >= border_len:
                        border_down = sequence[end_repeat:(end_repeat + border_len)]
                        # print(items[0], end_repeat, border_down, "\n")
                        # 2) downstream border does not contain Ns?
                        if not 'N' in border_down:
                            repeats_border_ok.append(items)
                            seq_file.write(">{}\n".format(idfrominput))
                            seq_file.write(sequence + "\n")
                            num += 1
    if border:
        # write the list of repeats in separate border file
        num_repeats_border_ok = len(repeats_border_ok)
        execution_time = time.time() - start_time
        if num_repeats_border_ok > 0:
            with open(out_filename_border, "w") as border_file:
                border_file.write("Program: {} {}\n".format(os.path.basename(sys.argv[0]), _version_))
                border_file.write("Date: {} , duration: {:.2f} sec\n".format(time.strftime("%Y-%m-%d %H:%M"), execution_time))
                border_file.write("Repeats file: {}\n".format(args.listfile))
                border_file.write("Db file: {}\n".format(args.allseqfile))
                border_file.write("{} repeats with border ({} bp)\n".format(num_repeats_border_ok, border))
                border_file.write("=" * 80 + "\n")
                for items in repeats_border_ok:
                    border_file.write("{}\n".format("\t".join(items)))
    return num


def load_dbfile(filename):
    fh = open(filename, "r")
    seq = fh.readlines()
    fh.close()
    return seq


def main(listfile, allseqfile, border):
        # Start timing
        start_time = time.time()

        # Read the list file and extract accession and repeat info
        accession_repeat_info = readseqidfrominfile(listfile, border)
        numofseqids = len(accession_repeat_info)

        # load db file in memory
        db_sequences = load_dbfile(allseqfile)

        # search and extract accessions of repeats in db file of sequences
        numfound = 0
        if numofseqids > 0:
            numfound = extract_print_seq(accession_repeat_info, db_sequences, border, start_time)

        # write summary to file
        execution_time = time.time() - start_time
        with open(out_filename_info, "w") as info_file:
            info_file.write("Program: {}\n".format(os.path.basename(sys.argv[0]), _version_))
            info_file.write("Date: {}, duration: {:.2f} sec\n".format(time.strftime("%Y-%m-%d %H:%M"), execution_time))
            info_file.write("Input, list of sequences and repeats from: {}\n".format(listfile))
            info_file.write("Input, db file with sequences FASTA: {}\n".format(allseqfile))
            info_file.write("Num lines of sequence in db file: {}\n".format(len(db_sequences)))
            if border:
                info_file.write("Mandatory border of bp: {}\n".format(str(border)))
            info_file.write("Found {} repeats\n".format(numfound, execution_time))
            if border:
                info_file.write("Output, accessions with border in FASTA format: {}\n".format(out_filename_seq_border))
                info_file.write("Output, list of accessions with border: {}\n".format(out_filename_border))
            else:
                info_file.write("Output, accessions in FASTA format: {}\n".format(out_filename_seq))


if __name__ == '__main__':
    # Checking command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("listfile", help="list of sequences and repeats")
    parser.add_argument("allseqfile", help="db file in FASTA format with all sequences")
    parser.add_argument("-b", "--border", help="must have n pb border up- and downstream")
    args = parser.parse_args()
    if args.listfile and args.allseqfile:
        main(args.listfile, args.allseqfile, args.border)
