#!/bin/bash
#SBATCH --job-name=gp2x               # Job name
#SBATCH --partition=sixhour           # Partition Name (Required)
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b047m507@ku.edu   # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --time=0-06:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=gp2x_%j.log          # Standard output and error log

module list

FNAME=$1
NUMGEN=$2
BES1=$3
BES2=$4
GNAME=$5
KNAME=$6

echo 'File Name '$FNAME
echo 'Number of events to generate: '$NUMGEN
echo 'Electron Beam Energy Spread [GeV]: '$BES1
echo 'Positron Beam Energy Spread [GeV]: '$BES2
echo 'Guinea Pig Input File: '$GNAME
echo 'KKMC Input File: '$KNAME

echo 'Running job as '
echo 'User '$USER

pwd
hostname
date

echo 'PATH '
echo $PATH
 
echo 'LD_LIBRARY_PATH'
echo $LD_LIBRARY_PATH

echo 'Run GP2X'

MWORK=/panfs/pfs.local/work/wilson/b047m507/GP2KKMC/GP2KKMC_Run

echo 'Script defines'
echo 'FNAME:  '${FNAME}
echo 'MWORK:  '${MWORK}

cd ${MWORK}
root -l 'GP2X.C('${NUMGEN}','${BES1}','${BES2}',"'${FNAME}'","'${GNAME}'","'${KNAME}'")'

echo 'Finished running GP2X'

date

exit
