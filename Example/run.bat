@ECHO OFF
SET SEQFILE=sequences.fasta
REM
ECHO ===============  MISA  ===============
ECHO run misa.pl
perl misa.pl %SEQFILE%
REM
ECHO =============== PySSRstat =============
REM
REM --------- set path ---------------------
Echo Set path to PySSRstat
REM PYSSRSTAT_HOME or 
IF "%PYSSRSTAT_HOME%"=="" GOTO LOneup
IF NOT "%PYSSRSTAT_HOME%"=="" GOTO LHome
:LHome
SET PSS=%PYSSRSTAT_HOME%
GOTO Running
:LOneup
SET PSS=%~dp0..
GOTO Running
REM
REM
REM --------------  Run  -------------------
:Running
ECHO Run statistics_misa
call %PSS%\statistics_misa.bat %SEQFILE%.statistics %SEQFILE%.misa -rpc
REM 
ECHO Run statgetlongest
call %PSS%\statgetlongest.bat repeats_analysis.txt %SEQFILE%.misa
REM 
ECHO Run filterrepeatsmisa
call %PSS%\filterrepeatsmisa.bat %SEQFILE%.misa 30 200 repeat
REM 
ECHO Run getsequences
call %PSS%\getsequences.bat filtered-repeats-sequence-list.txt %SEQFILE% -b 200
