# Directory Example

## Contents

The directory *Example* contains all input and output files of an analysis with MISA and the PySSRstat programs:

File with all sequences ("database") in FASTA format:

    sequences.fasta

Input and output files of MISA:

    misa.ini                           - configuration for MISA
    misa.pl                            - MISA script
    sequences.fasta.misa               - Output, list of repeats
    sequences.fasta.statistics         - Output, statstical data

Input and output files of PySSRstat:

    repeats_analysis.txt               - Output of statistics_misa.py
    longest-sequences-list.txt         - Output of statgetlongest.py
    filtered-repeats-sequence-list.txt - Output of filterrepeatsmisa.py
    border.txt                         - Output of getsequences.py
    getsequences-info.txt              - Output of getsequences.py
    repeats-sequences-border.fas       - Output of getsequences.py
    index.txt                          - Temporary file of getsequences.py
    
## Run analysis as batch

Besides the files mentioned above the directory Example contains also two MS-DOS batch files to run all the programs of MISA and PySSRstat as a batch and to clean-up:

    run.bat                             - Example batch to run analysis as batch
    clean.bat                           - clean-up all output files

Running clean.bat should leave you just with the starter files:
    
    sequences.fasta
    misa.pl
    misa.ini
    run.bat
    clean.bat
   
Launching run.bat runs first MISA and then the four main PySSRstat programs creating the output files described above in the section *Content*.
