#!/bin/bash
#SBATCH --job-name=whiz               # Job name
#SBATCH --partition=bigjay            # Partition Name (Required)
#SBATCH --mail-type=ALL               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=b047m507@ku.edu   # Where to send mail	
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=2gb                     # Job memory request
#SBATCH --time=0-10:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=whiz_%j.log          # Standard output and error log

module list

VERSION=$1
WNAME=$2
ELO=$3
EHI=$4
NEVT=$5

echo 'Version '$VERSION
echo 'Lower sqrt(s) bound (GeV) '$ELO
echo 'Upper sqrt(s) bound (GeV) '$EHI
echo 'Number of events to generate '$NEVT
echo 'Name of Sindarin file '$WNAME

echo 'Running job as '
echo 'User '$USER

pwd
hostname
date

echo 'PATH '
echo $PATH
 
echo 'LD_LIBRARY_PATH'
echo $LD_LIBRARY_PATH

echo 'Run KKMCee dimuon demo script'

MWORK=/panfs/pfs.local/work/wilson/b047m507/WHIZARD3
MYBIN=${MWORK}/BatchWhiz
MYWDIR=${MWORK}/WhizResults/Run${VERSION}

echo 'Creating directory '${MYWDIR}
mkdir ${MYWDIR}

echo 'Script defines'
echo 'MWORK :  '${MWORK}
echo 'MYWDIR:  '${MYWDIR}
echo 'MYBIN :  '${MYBIN}

cd ${MYBIN}

pwd

for ENE in $(seq ${ELO} 0.05 ${EHI})
do 

    echo 'CoM Energy ' ${ENE}
    
    cd ${MYBIN}

    MYWHIZ=${MYBIN}/${WNAME}_${ENE}
    mkdir ${MYWHIZ}
    echo 'Constructing local directory at ' ${MYWHIZ}
    cp ${MYBIN}/${WNAME}.sin ${MYWHIZ}/${WNAME}.sin

    cd ${MYWHIZ}

    sed 's/VCOM/'${ENE}'/g' ${MYWHIZ}/${WNAME}.sin > ${MYWHIZ}/whiz.temp
    sed 's/VNEV/'${NEVT}'/' ${MYWHIZ}/whiz.temp > ${MYWHIZ}/whiz.tmp

    cp ${MYWHIZ}/whiz.tmp ${MYWHIZ}/${WNAME}.sin
    
    rm ${MYWHIZ}/whiz.tmp
    rm ${MYWHIZ}/whiz.temp

    echo 'Running whizard'

    whizard ${MYWHIZ}/${WNAME}.sin

    sleep 10

    date

    echo 'Finished running WHIZARD for ' ${ENE}

    echo 'Scraping the log file for XSEC, XSECERR'

    XSEC=$(grep -B2 "Time estimate for" whizard.log --text | head -n 1 | awk '{print $3}')
    XSECERR=$(grep -B2 "Time estimate for" whizard.log --text | head -n 1 | awk '{print $4}')

    echo ${XSEC}
    echo ${XSECERR}

    echo 'Running python LHE->ROOT converter...'

    python ${MYBIN}/lhereader_whiz_ww.py ${ENE} ${WNAME}.lhe ${WNAME}_${ENE}.root ${XSEC} ${XSECERR}

    cp ${WNAME}_${ENE}.root ${MYWDIR}/${WNAME}_${ENE}.root
    cp whizard.log ${MYWDIR}/whiz_${ENE}.log
    
    echo 'Finished running python converter for ' ${ENE}

    cd ..
    
    rm -r ${MYWHIZ}

    ls -lrt
done

echo 'Finished running batch WHIZARD for all energies!'

date

exit
