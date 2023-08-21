#!/bin/bash
#SBATCH --job-name=dcsvgval           # Job name
#SBATCH --partition=bigjay            # Partition Name (Required)
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b047m507@ku.edu   # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=0-04:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=dcsvgval_%j.log      # Standard output and error log

module load python
module list

FNAME=$1
RNAME=$2
WNAME=$3
VERS=$4
HNAME=$5
ONAME=$6
NBINS=$7

#datafile needs to be absolute path
#reference data needs to be absolute path
#version needs to be local (NOT ANYMORE IT IS NOW LOCAL)
#work directory name needs to be absolute path
#others should be local (no path needed)
echo 'Datafile name'$FNAME
echo 'Reference data name'$RNAME
echo 'Work directory name '$WNAME
echo 'Version name '$VERS
echo 'Histogram Name '$HNAME
echo 'Savgolified File Name '$ONAME
echo 'Deconvolution number of bins'$NBINS

echo 'Running job as '
echo 'User '$USER

pwd
hostname
date

echo 'PATH '
echo $PATH

echo 'LD_LIBRARY_PATH'
echo $LD_LIBRARY_PATH

MWORK=/panfs/pfs.local/work/wilson/b047m507/GP2X_dev

#Where the repo for SVG python is located                                                                                                                               
REPO=/panfs/pfs.local/work/wilson/b047m507/GP2X_dev/MUMUUNP

cd ${MWORK}

mkdir ${WNAME}

cd ${WNAME}

cp ${REPO}/SavGolHist3.py ${WNAME}/SavGolHistWork.py
cp ${REPO}/SavGolFilt3.py ${WNAME}/SavGolFiltWork.py
cp ${REPO}/SavGolSave3.py ${WNAME}/SavGolSaveWork.py
cp ${REPO}/FFTDetDeconv_Val.C ${WNAME}/FFTDetDeconv_Val.C

echo 'Script defines'
echo 'MWORK:  '${MWORK}

echo 'Running FFT Deconvolution of Detector level to Generator level'

root -l 'FFTDetDeconv_Val.C(180,121.0,125.5,'${NBINS}',0.0001,250.0,124.960,0.24,"'${RNAME}'","'${FNAME}'","'${VERS}'","0.5*(DiMup4->E() + DiMup4->Pz())","0.5*(SDiMup4->E() + SDiMup4->Pz())")'

echo 'Running python SVG filter algorithms'

date

python2 SavGolHistWork.py HistMC_${VERS}.root ${HNAME} OutHist.pkl

python SavGolFiltWork.py OutHist.pkl OutFilt.pkl 71 3 25

python2 SavGolSaveWork.py OutFilt.pkl ${ONAME} hs1 180 121.0 125.5 407222.0

echo 'Finished SVG filter algorithms'

date

echo 'Run DECONSVG_VAL'

cd ${MWORK}
#Copy over to our work directory of this run
cp DECONSVG_VALM.C ${WNAME}/DECONSVG_VALM.C
cd ${WNAME}

root -l 'DECONSVG_VALM.C("'${FNAME}'","'${ONAME}'")'

echo 'Finished running DECONSVG_VALM'

date
exit
