#!/usr/bin/env python3
"""Extract additional statistical data form MISA-statistics-file and MISA-file

containing data form the MISA statistic file and
- relative abundance
- relative density
- longest motifs by repeat length
- total abundance of repeat types in absolute numbers and percentage
- number of perfect, imperfect and compound SSRs (experimental)

Usage:
statistics_misa.py <MISA-statstic-file> <MISA-file> [-rpc|--repeatclassses]

Output:
repeats_analysis.txt

Author: Mario Nenno
Version: 2015-08-31
Copyright: see file LICENCE.txt

"""
    
import os
import time
import sys
import argparse
import re
from collections import OrderedDict

# program version
_version_ = '1.0'

# output file
out_filename = 'repeats_analysis.txt'


def main(statisticsfile, misafile, repeatclasses):
    """Main function contains all logic
    """
    # Names of the repeat length
    mono = 'Mono'
    di = 'Di'
    tri = 'Tri'
    tetra = 'Tetra'
    penta = 'Penta'
    hexa = 'Hexa'
    septa = 'Septa'
    octa = 'Octa'
    nona = 'Nona'
    deca = 'Deca'

    definement_line = ''

    # Mapping repeat length to name
    repeat_len_name = {1: mono, 2: di, 3: tri, 4: tetra, 5: penta, 6: hexa, 7: septa, 8: octa, 9: nona, 10: deca}

    total_length_examined_misa = 0
    total_number_ssr_misa = 0
    total_ssr_in_compound_misa = 0
    max_repeat_unit_length = 0  # extract from MISA statistic file

    abundances = {}             # dictionary of abundance by repeat length
    longest_motifs = {}         # dictionary of longest, k is motif, value is longest repeat
    motif_total_abundance = {}  # dictionary of total abundance, k is motif
    total_abundance = {}        # dictionary of total abundance of repeat type

    is_definement = False
    is_distribution_repeat_types = False
    is_motifs = False
    is_in_motifs = False
    is_results = False
    is_repeat_types = False
    is_in_repeat_types = False

    # Init monos
    longest_motifs["A"] = 0
    longest_motifs["T"] = 0
    longest_motifs["C"] = 0
    longest_motifs["G"] = 0

    # take start time to calc duration later
    start_time = time.time()

    def percent_abundance(abund, total):
        return (abund*100)/total

    def get_num_right(part):
        text, num = part.split(":")
        num.rstrip().lstrip()
        num = int(num)
        return num

    def percent_repeat(d):
        percents = {}
        sum_of_grp = 0
        for k in d:
            sum_of_grp += int(d[k])
        for k in d:
            abund = int(d[k])
            # sum_of_grp = 100%
            # ab  = x%   = (ab*100)/sum_of_grp
            pc = (abund*100)/sum_of_grp
            percents[k] = pc
        return percents

    def list_as_string(a_list):
        # Sort the list before joining
        return ', '.join([str(x) for x in sorted(a_list)])

    def parse_definement(line_definement):
        # get list of Definement of microsatellites (unit size / minimum number of repeats)
        list_definements = []
        # extract by regex
        pattern = re.compile(r'\((\d+)/(\d+)\)')
        for pair in line_definement.split():
            m = pattern.match(pair)
            # group to tuple
            t = (int(m.group(1)), int(m.group(2)))
            # add tuple to list
            list_definements.append(t)
        return list_definements

    # Read from the MISA statistics file
    with open(statisticsfile, "r") as infile:
        idx_definement = 0
        for line in infile:
            if line.startswith("Definement of microsatellites"):
                is_definement = True
            if line.startswith("RESULTS OF MICROSATELLITE SEARCH"):
                is_results = True
                is_definement = False
            if line.startswith("Distribution to different repeat type"):
                is_distribution_repeat_types = True
            if line.startswith("Frequency of identified SSR motifs"):
                is_motifs = True
                is_distribution_repeat_types = False
            if line.startswith("Frequency of classified repeat types"):
                is_repeat_types = True
                is_motifs = False
                is_in_motifs = False

            if is_definement:
                if idx_definement == 0:
                    idx_definement += 1
                    continue
                elif idx_definement == 1:
                    definement_line = line.rstrip()
                    # parse list_definements to get the max repeat unit length
                    definements = parse_definement(line)
                    max_repeat_unit_length, dummy = definements[-1]
                    is_definement = False
                    continue

            if is_distribution_repeat_types:
                if not line[0].isdigit():
                    continue
                else:
                    s, n = line.rstrip().split("\t")
                    size = int(s)
                    number = int(n)
                    abundances[repeat_len_name[size]] = number

            if is_results:
                if len(line) == 0:
                    continue
                else:
                    if line.startswith("Total size of examined sequences"):
                        total_length_examined_misa = get_num_right(line)
                    elif line.startswith("Total number of identified SSRs"):
                        total_number_ssr_misa = get_num_right(line)
                    elif line.startswith("Number of SSRs present in compound"):
                        total_ssr_in_compound_misa = get_num_right(line)
            if is_motifs:
                is_results = False
                if line.startswith("Repeats"):
                    # in header line
                    is_in_motifs = True
                    header = line.split("\t")
                else:
                    if is_in_motifs:
                        line = line.rstrip()
                        if len(line) == 0:
                            continue
                        else:
                            # numLinesMotifs += 1
                            # split the data line
                            list_of_data = line.split("\t")
                            # get the SSR motif
                            motif = list_of_data[0]
                            # get total
                            if motif != "Repeats":
                                motif_total_abundance[motif] = list_of_data[-1]
                            else:
                                print("motif is ", motif)
                            # column index is used later to get number of repeats
                            repeat_idx = len(header)-1
                            # for longest find first column, not empty, in reverse order
                            for ntimes in reversed(list_of_data[:-1]):
                                repeat_idx -= 1
                                if ntimes == "":
                                    continue
                                else:
                                    longest_motifs[motif] = int(header[repeat_idx])
                                    break
            if is_repeat_types:
                if line.startswith("Repeats"):
                    # in header line
                    is_in_repeat_types = True
                    header = line.split("\t")
                if is_in_repeat_types:
                    line = line.rstrip()
                    if len(line) == 0:
                        continue
                    else:
                        # numLinesRepeatTypes += 1
                        # split the data line
                        list_of_data = line.split("\t")
                        # get the repeat type
                        repeat_type = list_of_data[0]
                        # get the total abundance of motif
                        if repeat_type != "Repeats":
                            total_abundance[repeat_type] = list_of_data[-1]

    # ----------------------  Analysis of longtest Motifs (Di, Tri, etc. ) ------------------------------
    # Longest motif length
    # initialize dict of longest all longest to 0
    longest_repeat = {}
    for i in range(2, max_repeat_unit_length + 1):
            longest_repeat[repeat_len_name[i]] = 0

    # find max of multi
    for key in longest_motifs:
        key_len = len(key)
        count = longest_motifs[key]
        # for Di to Deca
        if 1 < key_len < (max_repeat_unit_length + 1):
            if count > longest_repeat[repeat_len_name[key_len]]:
                longest_repeat[repeat_len_name[key_len]] = count

    # Motifs of longest motif length
    # Collect all motifs having length of longest motif
    # initialize the dict of motifs with empty lists
    longest_repeat_motifs = {}
    for i in range(2, (max_repeat_unit_length + 1)):
        longest_repeat_motifs[repeat_len_name[i]] = []

    for key in longest_motifs:
        key_len = len(key)
        count = longest_motifs[key]
        # A dictionary of repeat_len_name with a list as value
        if 1 < key_len < (max_repeat_unit_length + 1):
            if count == longest_repeat[repeat_len_name[key_len]]:
                longest_repeat_motifs[repeat_len_name[key_len]].append(key)

    # -------------------  Calculate relative abundance and density  --------------
    total_length_ssr = 0
    for i in range(1, (max_repeat_unit_length + 1)):
        total_length_ssr += (i * abundances[repeat_len_name[i]])

    total_length_examined_mb_misa = total_length_examined_misa/1000000
    relatve_abundance = total_number_ssr_misa / total_length_examined_mb_misa
    relative_density = total_length_ssr / total_length_examined_mb_misa

    # --------------------------  Repeat types  ----------------------------------
    dict_type_mono = {}
    dict_type_di = {}
    dict_type_tri = {}
    dict_type_tetra = {}
    dict_type_penta = {}
    dict_type_hexa = {}
    dict_type_septa = {}
    dict_type_octa = {}
    dict_type_nona = {}
    dict_type_deca = {}

    # len   repeat_len_name (repeat_len)
    # 3 > mono (1) = (len-1)/2
    # 5 > di (2) = (len-1)/2
    for repeat_type in total_abundance:
        abundance = total_abundance[repeat_type]
        type_length = len(repeat_type)
        # if len(repeat_type) == 3:
        if 3 == type_length:
            dict_type_mono[repeat_type] = abundance
        elif 5 == type_length:
            dict_type_di[repeat_type] = abundance
        elif 7 == type_length:
            dict_type_tri[repeat_type] = abundance
        elif 9 == type_length:
            dict_type_tetra[repeat_type] = abundance
        elif 11 == type_length:
            dict_type_penta[repeat_type] = abundance
        elif 13 == type_length:
            dict_type_hexa[repeat_type] = abundance
        elif 15 == type_length:
            dict_type_septa[repeat_type] = abundance
        elif 17 == type_length:
            dict_type_octa[repeat_type] = abundance
        elif 19 == type_length:
            dict_type_nona[repeat_type] = abundance
        elif 21 == type_length:
            dict_type_deca[repeat_type] = abundance

    # --------------------------  composition types  ----------------------------------
    # composition types: p = perfect, c = imperfect, c* = compound
    num_perfect_ssr = 0
    num_imperfect_ssr = 0
    num_compound_ssr = 0
    if repeatclasses:
        with open(misafile, "r") as infile:
            for line in infile:
                line = line.rstrip()
                items = line.split("\t")
                #     ['PK14952.1', '1', 'p1', '(T)11', '11', '2418', '2428\n']
                #  Idx  0            1    2      3
                ssrtype = items[2]
                if 'c' == ssrtype:
                    num_imperfect_ssr += 1
                elif 'c*' == ssrtype:
                    num_compound_ssr += 1
                elif ssrtype.startswith("p"):
                    num_perfect_ssr += 1

    # End timing
    execution_time = time.time() - start_time

    # write result into output file
    with open(out_filename, 'w') as outfile:
        outfile.write("Program: {} {}\n".format(os.path.basename(sys.argv[0]), _version_))
        outfile.write("Statistic file analysed: {}\n".format(statisticsfile))
        outfile.write("Definement of microsatellites (MISA): {}\n".format(definement_line))
        outfile.write("Misa file analysed: {}\n".format(misafile))
        outfile.write("Date: {} Execution time: {:.2f} sec\n".format(time.strftime("%Y-%m-%d %H:%M"), execution_time))
        outfile.write("Note: Numbers label with '(MISA)' are not calculated but read from input\n")
        if repeatclasses:
            outfile.write("Option repeatclasses: {}\n".format(repeatclasses))
        outfile.write("="*50 + "\n\n")

        outfile.write("Total length of Sequence examined (bp) (MISA): {}\n".format(total_length_examined_misa))
        outfile.write("Total number of SSR (MISA) : {:>8}\n".format(total_number_ssr_misa))
        outfile.write("Total length of SSR (bp)   : {:>8}\n".format(total_length_ssr))
        outfile.write("Relative abundance (SSR/Mb): {0:11.2f}\n".format(relatve_abundance))
        outfile.write("Relative density (bp/Mb)   : {0:11.2f}\n".format(relative_density))

        outfile.write("\n====== Distribution of motif length  =========\n")
        outfile.write("{0:<5} {1:<4}%   {2:6}  {3:>7}\n".format("Motif", "", "Number", "Length[bp]"))
        outfile.write("-" * 40 + "\n")
        out_format = "{0:<5}: {1:4.1f}% ({2:>6})    {3:>7}\n"
        rep_total_abundance = 0
        rep_total_length = 0
        for i in range(1, (max_repeat_unit_length + 1)):
            ab = abundances[repeat_len_name[i]]
            rep_total_abundance += ab
            rep_total_length += i*ab
            outfile.write(out_format.format(
                repeat_len_name[i], percent_abundance(ab, total_number_ssr_misa), ab, i*ab))
        outfile.write("-" * 40 + "\n")
        outfile.write("{0:<5} {1:<5}   {2:6}     {3:>7}\n".format("Total", "", rep_total_abundance, rep_total_length))

        # Verify total abundance
        if rep_total_abundance != total_number_ssr_misa:
            outfile.write("ERROR, total number is NOT consistant")

        # ------------ Group by type length
        if repeatclasses:
            outfile.write("\n-------- Number of SSRs by repeat class (experimental) ------\n")
            outfile.write("Perfect SSRs: (p)                : {:>8}\n".format(num_perfect_ssr))
            outfile.write("Imperfect SSRs (c)               : {:>8}\n".format(num_imperfect_ssr))
            outfile.write("Compound SSRs (c*)               : {:>8}\n".format(num_compound_ssr))
            subtotal = num_perfect_ssr + num_imperfect_ssr + num_compound_ssr
            outfile.write("                        Subtotal : {:>8}\n".format(subtotal))
            outfile.write("SSRs in compound formation (MISA): {:>8}\n".format(total_ssr_in_compound_misa))
            outfile.write("                           Total : {:>8}\n".format(subtotal + total_ssr_in_compound_misa))

        outfile.write("\n======== Longest motives =======\n")
        outfile.write("Mono  A: %d, T: %d, C: %d, G: %d\n" % (
            longest_motifs["A"], longest_motifs["T"], longest_motifs["C"], longest_motifs["G"]))
        # Di- to Deca
        for i in range(2, (max_repeat_unit_length + 1)):
            outfile.write("{0:<5} {1:>5}: {2}\n".format(
                repeat_len_name[i], longest_repeat[repeat_len_name[i]],
                list_as_string(longest_repeat_motifs[repeat_len_name[i]])))

        outfile.write("\n======= Total abundance of repeat types =======\n")
        out_format = "{0}: {1:5.1f}% ({2:>5})\n"

        outfile.write("---------  Mono-nucleotides\n")
        percent_repeat_type = percent_repeat(dict_type_mono)
        for repeat_type in ['A/T', 'C/G']:
            abundance = dict_type_mono[repeat_type]
            outfile.write(out_format.format(repeat_type, percent_repeat_type[repeat_type], abundance))

        outfile.write("---------  Di-nucleotides\n")
        percent_repeat_type = percent_repeat(dict_type_di)
        for repeat_type in ['AC/GT', 'AG/CT', 'AT/AT', 'CG/CG']:
            if repeat_type in dict_type_di:
                abundance = dict_type_di[repeat_type]
                outfile.write(out_format.format(repeat_type, percent_repeat_type[repeat_type], abundance))
        # ----- Tri- to Dec-nucleotides
        repeat_types = [dict_type_tri, dict_type_tetra, dict_type_penta, dict_type_hexa, dict_type_septa,
                        dict_type_octa, dict_type_nona, dict_type_deca]
        unit_length = 2
        for dictType in repeat_types:
            unit_length += 1
            if unit_length > max_repeat_unit_length:
                break
            else:
                type_length = len(list(dictType.keys())[0])
            if type_length == 7:
                outfile.write("---------  Tri-nucleotides\n")
            elif type_length == 9:
                outfile.write("---------  Tetra-nucleotides\n")
            elif type_length == 11:
                outfile.write("---------  Penta-nucleotides\n")
            elif type_length == 13:
                outfile.write("---------  Hexa-nucleotides\n")
            elif type_length == 15:
                outfile.write("---------  Septa-nucleotides\n")
            elif type_length == 17:
                outfile.write("---------  Octa-nucleotides\n")
            elif type_length == 19:
                outfile.write("---------  Nona-nucleotides\n")
            elif type_length == 21:
                outfile.write("---------  Deca-nucleotides\n")
            else:
                outfile.write("--------- \n")

            percent_repeat_type = percent_repeat(dictType)
            # sort first by length of k and then by k
            od = OrderedDict(sorted(dictType.items(), key=lambda t: (len(t[0]), t[0])))
            for repeat_type in od:
                abundance = od[repeat_type]
                percent = percent_repeat_type[repeat_type]
                if percent > 1.0:
                    outfile.write(out_format.format(repeat_type, percent, abundance))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("statisticfile", help="MISA statistics file")
    parser.add_argument("misafile", help="MISA file")
    parser.add_argument("-rpc", "--repeatclasses", help="include SSR repeat classes", action="store_true")
    args = parser.parse_args()
    if args.statisticfile and args.misafile:
        main(args.statisticfile, args.misafile, args.repeatclasses)
