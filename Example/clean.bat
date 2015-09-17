@ECHO OFF
REM ********** clean-up  *********
REM
REM ====  MISA output files  ======
REM
IF EXIST sequences.fasta.misa DEL sequences.fasta.misa
IF EXIST sequences.fasta.statistics DEL sequences.fasta.statistics
REM
REM ====== PySSRstat  ========
REM 
REM ---  statistics_misa
IF EXIST repeats_analysis.txt DEL repeats_analysis.txt
REM
REM ---  statgetlongest
IF EXIST longest-sequences-list.txt DEL longest-sequences-list.txt
REM
REM --- filterrepeatsmisa
IF EXIST filtered-repeats-sequence-list.txt DEL filtered-repeats-sequence-list.txt
REM
REM --- getsequences
IF EXIST index.txt DEL index.txt
IF EXIST getsequences-info.txt DEL getsequences-info.txt
IF EXIST repeats-sequences.fas DEL repeats-sequences.fas
IF EXIST repeats-sequences-border.fas DEL repeats-sequences-border.fas
IF EXIST border.txt DEL border.txt
