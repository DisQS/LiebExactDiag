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
jobdir=$currdir

binarydir=$HOME/Projects/LiebExactDiag/EXE
#binarydir=/storage/disqs/LiebSparseDiag/EXE

for CubeConPot in 0.0 10.0 # 1.0 2.0 5.0 10.0 20.0 50.0 100.0
    #1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 20.0 50.0 60.0 70.0 80.0
do

echo "--- CCPot=" $CubeConPot

jobname="LED-$size-CP$CubeConPot"
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
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=7700
#SBATCH --partition=hmem
#SBATCH --account=su007-rr

module purge
module load intel/2019b # sulis
module load GCCcore/10.3.0 parallel/20210622 # sulis

#module load GCC/10.2.0 parallel intel #avon

for iseed in {1..$config..1}
do

myseed=\$(( $seed + \$iseed - 1))
echo "--- working on config" \$iseed "with seed" \$myseed

# create the input file
echo create the input file
inpfile=LEDdiag-$CubeConPot-\$iseed.inp
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
echo "CubeConPot    = $CubeConPot      ">>  \$inpfile #
echo "CubeDis0      = 0.0        ">>  \$inpfile #
echo "CubeDis1      = 10.0            ">>  \$inpfile #
echo "dCubeDis      = 5.0           ">>  \$inpfile #
#echo "LiebDis0      = $CubeConPot      ">>  \$inpfile #
echo "LiebConPot    = 0.0            ">>  \$inpfile #
echo "LiebDis0      = 0.0            ">>  \$inpfile #

cat \$inpfile

#$binarydir/$binary <\$inpfile

done

MY_PARALLEL_OPTS="-N 1 --delay .2 -j \$SLURM_NTASKS --joblog parallel-\${SLURM_JOBID}.log"
MY_SRUN_OPTS="-N 1 -n 1 --exclusive"
MY_EXEC="$binarydir/$binary <LEDdiag-$CubeConPot-{}.inp"

parallel \$MY_PARALLEL_OPTS srun \$MY_SRUN_OPTS \$MY_EXEC ::: {1..$config}

#zip -mv LED-$CubeConPot.zip Evec*.raw
zip -m inp.zip LEDdiag-$CubeConPot-*.inp
zip -m sh.zip *.sh

exit 0

EOD

chmod 755 ${jobfile}
#(msub -q devel $jobdir/${jobfile}) # for queueing system
#(sbatch -q devel $jobdir/${jobfile}) # for queueing system
sbatch ${jobfile} # for queueing system
#(source $jobdir/${jobfile} ) >& $jobdir/${logfile} & # for parallel shell execution
#source ${jobfile} #>& ${logfile} # for sequential shell execution

sleep 1

cd ..

done

