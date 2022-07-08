#!/bin/bash

# settings from input

size=${1:-10}
seed=${2:-1}
config=${3:-2}
keep=${4:-1}

echo "LSD: making for M=" $size "with starting seed=" $seed "and" $config "samples"

# settings for files

binary=LSDdiag.IC

# settings for directories

currdir=`pwd`
jobdir=$currdir

binarydir=$HOME/Projects/LiebSparseDiag/EXE
#binarydir=/storage/disqs/LiebSparseDiag/EXE

for disorder in 1.0 2.0 #20.0 25.0 30.0 35.0 40.0 10.0 9.0 8.0 7.0 6.0 5.0 4.0 3.0 2.0 1.0 15.0 14.0 13.0 12.0 11.0 0.2 0.4 0.6 0.7 0.8 1.2 1.4 1.6 1.8 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5
#50.0 60.0 70.0 80.0 90.0 100.0 
do

energy=`echo "$disorder/2.0 + 4.05"| bc`

echo "--- hDis=" $disorder ", Min_Eng=" $energy

jobname="LSD-$size-hD$disorder"
echo $jobname

jobfile=`printf "$jobname.sh"`
logfile=`printf "$jobname.log"`
jobdir="LSD-$size"
mkdir -p $jobdir

echo "binarydir=" $binarydir " jobdir=" $jobdir 

# settings for parallel submission

cd $jobdir

cat > ${jobfile} << EOD
#!/bin/bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=3700

module purge
module load GCC/10.2.0 parallel 
module load intel

for iseed in {1..$config..1}
do

myseed=\$(( $seed + \$iseed - 1))
echo "--- working on config" \$iseed "with seed" \$myseed

# create the input file
echo create the input file
inpfile=LSDdiag-$disorder-\$iseed.inp
touch \$inpfile

echo "ISeed         = \$myseed       ">  \$inpfile #
echo "NConfig       = 1        ">>  \$inpfile #
echo "Dim           = 3            ">>  \$inpfile #
echo "Nx            = 1            ">>  \$inpfile #
echo "IBCFlag       = 1             ">>  \$inpfile #
echo "IRNGFlag      = 0             ">>  \$inpfile #
echo "IKeepFlag     = $keep      ">>  \$inpfile #
echo "IWriteFlag    = 2       ">>  \$inpfile #
echo "IStateFlag    = 0       ">>  $inpfile #
echo "Width0        = $size       ">>  \$inpfile #
echo "Width1        = $size       ">>  \$inpfile #
echo "dWidth        = 2          ">>  \$inpfile #
echo "HubDis0       = $disorder      ">>  \$inpfile #
echo "HubDis1       = $disorder           ">>  \$inpfile #
echo "dHubDis       = 1.0           ">>  \$inpfile #
#echo "RimDis0       = $disorder      ">>  \$inpfile #
echo "RimDis0       = 0.0            ">>  \$inpfile #
echo "Energy0       = -$energy       ">>  \$inpfile #
echo "Energy1       = $energy       ">>  \$inpfile #
echo "dEnergy       = 1.0       ">>  \$inpfile #
echo "NEvals        = 100           ">>  \$inpfile #
echo "Memory        = 100          ">>  \$inpfile #

cat \$inpfile

#$binarydir/$binary <\$inpfile

done

MY_PARALLEL_OPTS="-N 1 --delay .2 -j \$SLURM_NTASKS --joblog parallel-\${SLURM_JOBID}.log"
MY_SRUN_OPTS="-N 1 -n 1 --exclusive"
MY_EXEC="$binarydir/$binary <LSDdiag-$disorder-{}.inp"

parallel \$MY_PARALLEL_OPTS srun \$MY_SRUN_OPTS \$MY_EXEC ::: {1..$config}

exit 0

EOD

chmod 755 ${jobfile}
#(msub -q devel $jobdir/${jobfile}) # for queueing system
#(sbatch -q devel $jobdir/${jobfile}) # for queueing system
sbatch ${jobfile} # for queueing system
#(source $jobdir/${jobfile} ) >& $jobdir/${logfile} & # for parallel shell execution
#source ${jobfile} #>& ${logfile} # for sequential shell execution

#echo "<return>"
sleep 1

cd ..

done

