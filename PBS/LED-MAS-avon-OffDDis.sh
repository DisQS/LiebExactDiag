#!/bin/bash

# settings from input

CPflag=${1:-0}
size=${2:-10}
seed=${3:-1}
config=${4:-2}
keep=${5:-1}

echo "LED: making for M=" $size "with starting seed=" $seed "and" $config "samples"

# settings for files

binary=LEDdiag.IC

# settings for directories

currdir=`pwd`
jobdir=$currdir

binarydir=$HOME/Projects/LiebExactDiag/EXE
#binarydir=/storage/disqs/LiebSparseDiag/EXE

for OffDDis in 0.1 1.0 2.0 5.0 6.0 #10.0 20.0 50.0 100.0
    #1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0 50.0 60.0 70.0 80.0
do

echo "--- OffDDis=" $OffDDis

jobname="LED-$CPflag-$size-oD$OffDDis"
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
#!/bin/bash
#SBATCH --nodes=1
##SBATCH --ntasks-per-node=24
#SBATCH --ntasks-per-node=4
##SBATCH --time=48:00:00
#SBATCH --time=00:00:20
##SBATCH --mem-per-cpu=31418
##SBATCH --partition=hmem

module purge
module load GCC/10.2.0 parallel 
module load intel

for iseed in {1..$config..1}
do

myseed=\$(( $seed + \$iseed - 1))
echo "--- working on config" \$iseed "with seed" \$myseed

# create the input file
echo create the input file
inpfile=LEDdiag-oD$OffDDis-\$iseed.inp
touch \$inpfile

echo "ISeed          = \$myseed       ">  \$inpfile #
echo "NConfig        = 1              ">>  \$inpfile #
echo "Dim            = 3              ">>  \$inpfile #
echo "Nx             = 1              ">>  \$inpfile #
echo "IBCFlag        = 1              ">>  \$inpfile #
echo "IRNGFlag       = $CPflag        ">>  \$inpfile #
echo "IKeepFlag      = $keep          ">>  \$inpfile #
echo "IWriteFlag     = 2              ">>  \$inpfile #
echo "IStateFlag     = -1             ">>  \$inpfile #
echo "IFluxFlag      = 5              ">>  \$inpfile #
echo "Width0         = $size          ">>  \$inpfile #
echo "Width1         = $size          ">>  \$inpfile #
echo "dWidth         = 2              ">>  \$inpfile #
echo "0 CubeConPot0  = 0.0            ">>  \$inpfile #
echo "0 CubeConPot1  = 0.1            ">>  \$inpfile #
echo "0 dCubeConPot  = 0.25           ">>  \$inpfile #
echo "1 CubeDis0     = 0.0            ">>  \$inpfile #
echo "1 CubeDis1     = 0.1            ">>  \$inpfile #
echo "1 dCubeDis     = 0.25           ">>  \$inpfile #
echo "2 LiebConPot   = 0.0            ">>  \$inpfile #
echo "3 LiebDis      = 0.0            ">>  \$inpfile #
echo "4 OffDShift0   = 0.0            ">>  \$inpfile #
echo "4 OffDShift1   = 0.1            ">>  \$inpfile #
echo "4 dOffDShift   = 0.5            ">>  \$inpfile #
echo "5 OffDDis0     = $OffDDis       ">>  \$inpfile #
echo "5 OffDDis1     = $OffDDis       ">>  \$inpfile #
echo "5 dOffDDis     = 0.5            ">>  \$inpfile #

cat \$inpfile

#$binarydir/$binary <\$inpfile

done

MY_PARALLEL_OPTS="-N 1 --delay .2 -j \$SLURM_NTASKS --joblog parallel-\${SLURM_JOBID}.log"
MY_SRUN_OPTS="-N 1 -n 1 --exclusive"
MY_EXEC="$binarydir/$binary <LEDdiag-oD$OffDDis-{}.inp"

echo parallel \$MY_PARALLEL_OPTS srun \$MY_SRUN_OPTS \$MY_EXEC ::: {1..$config}
parallel \$MY_PARALLEL_OPTS srun \$MY_SRUN_OPTS \$MY_EXEC ::: {1..$config}

#zip -mv LED-$CubeConPot.zip Evec*.raw
zip -um inp.zip LEDdiag-CP$CubeConPot-*.inp
zip -um sh.zip *.sh

exit 0

EOD

chmod 755 ${jobfile}
#(msub -q devel $jobdir/${jobfile}) # for queueing system
#(sbatch -q devel $jobdir/${jobfile}) # for queueing system
#sbatch --account=su007 ${jobfile} # for queueing system
sbatch ${jobfile} # for queueing system
#(source $jobdir/${jobfile} ) >& $jobdir/${logfile} & # for parallel shell execution
#source ${jobfile} #>& ${logfile} # for sequential shell execution

sleep 1

cd ..

done

