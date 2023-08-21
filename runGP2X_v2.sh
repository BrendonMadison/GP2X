#!/bin/bash
#SBATCH --job-name=gp2x               # Job name
#SBATCH --partition=bigjay           # Partition Name (Required)
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b047m507@ku.edu   # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=4gb                     # Job memory request
#SBATCH --time=0-12:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=gp2x_%j.log          # Standard output and error log

module list

WNAME=$1
FNAME=$2
NUMGEN=$3
BES1=$4
BES2=$5
GNAME=$6
KNAME=$7

echo 'Work directory name '$WNAME
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

echo 'Run GP2X_v2'

MWORK=/panfs/pfs.local/work/wilson/b047m507/GP2X_dev

mkdir ${MWORK}/${WNAME}

echo 'Script defines'
echo 'MWORK:  '${MWORK}

cd ${MWORK}
#Copy over to our work directory of this run
cp GP2X_v2.C ${WNAME}/GP2X_v2.C
cd ${WNAME}

root -l 'GP2X_v2.C('${NUMGEN}','${BES1}','${BES2}',"'${FNAME}'","'${GNAME}'","'${KNAME}'")'

echo 'Finished running GP2X_v2'

date

exit
