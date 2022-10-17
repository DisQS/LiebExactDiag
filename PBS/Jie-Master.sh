#!/bin/bash

# settings from input

size=${1:-10}
seed=${2:-1}
config=${3:-2}
keep=${4:-1}

echo "LED: making for M=" $size "with starting seed=" $seed "and" $config "samples"

# settings for files

binary=LEDdiag.IC

# settings for directories

currdir=`pwd`
#jobdir=$currdir

binarydir=$currdir/../EXE
[ -d ${binarydir} ] || mkdir ${binarydir}
#binarydir=/storage/disqs/LiebSparseDiag/EXE

cp $currdir/../src/$binary $binarydir



for CubeConstPoten in 0.0 2.0 #5.0 10.0 20.0 50.0 100.0
do

echo "--- CubeConstPoten=" $CubeConstPoten

for disorder in 1.0 2.0 3.0 #4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0 50.0 60.0 70.0 80.0
do

echo "--- hDis=" $disorder

jobname="LED-$size-hD$disorder-CubeP$CubeConstPoten"
echo $jobname

jobfile=`printf "$jobname.sh"`
logfile=`printf "$jobname.log"`
jobdir="LED-$size"
mkdir -p $jobdir

echo "binarydir=" $binarydir " jobdir=" $jobdir 

# settings for parallel submission

cd $jobdir

cat > ${jobfile} << EOD
#!/bin/bash
#PBS -l nodes=${nodes}:ppn=16
#PBS -l pmem=${memory}
#PBS -l walltime=04:00:00

#       The jobname
#PBS -N ${jobname}

##############################
## comment out when INTERACTIVE
## SLURM_NTASKS=8
## SLURM_JOBID=1
##############################

# construct the input file

for iseed in {1..$config..1}
do

myseed=\$(( $seed + \$iseed - 1))
echo "--- working on config" \$iseed "with seed" \$myseed

# create the input file
echo create the input file
inpfile=LEDdiag-$disorder-$CubeConstPoten-\$iseed.inp
touch \$inpfile

echo "ISeed         = \$myseed       ">  \$inpfile #
echo "NConfig       = 1        ">>  \$inpfile #
echo "Dim           = 3            ">>  \$inpfile #
echo "Nx            = 1            ">>  \$inpfile #
echo "IBCFlag       = 1             ">>  \$inpfile #
echo "IRNGFlag      = 0             ">>  \$inpfile #
echo "IKeepFlag     = $keep      ">>  \$inpfile #
echo "IWriteFlag    = 2       ">>  \$inpfile #
echo "IStateFlag    = -1       ">>  \$inpfile #
echo "Width0        = $size       ">>  \$inpfile #
echo "Width1        = $size       ">>  \$inpfile #
echo "dWidth        = 2          ">>  \$inpfile #
echo "HubDis0       = $disorder      ">>  \$inpfile #
echo "HubDis1       = $disorder           ">>  \$inpfile #
echo "dHubDis       = 1.0           ">>  \$inpfile #
#echo "RimDis0       = $disorder      ">>  \$inpfile #
echo "RimDis0       = 0.0            ">>  \$inpfile #
echo "CubeConstPoten = $CubeConstPoten  ">>  \$inpfile #
echo "LiebConstPoten = 0.0           ">>  \$inpfile #
 
cat \$inpfile

${binarydir}/${binary} <\$inpfile >& ${logfile}
#${binarydir}/${binary} <\$inpfile >& ${logfile}

done


wait
#exit 0

EOD

chmod 755 ${jobfile}
#(msub -q devel $jobdir/${jobfile}) # for queueing system
#(sbatch -q devel $jobdir/${jobfile}) # for queueing system
#sbatch ${jobfile} # for queueing system
#(source $jobdir/${jobfile} ) >& $jobdir/${logfile} & # for parallel shell execution
#source ${jobfile} #>& ${logfile} # for sequential shell execution
(source ${jobfile} ) &

sleep 1

cd ..

done

done
