#!/bin/bash

# settings from input

energy=${1:-1.0}
size=${2:-10}
seed=${3:-1}
config=${4:-2}
keep=${5:-1}

echo "LSD: making for E=" $energy "with M=" $size "and starting seed=" $seed "and" $config "samples"

# settings for files

binary=LSDdiag.IC

# settings for directories

currdir=`pwd`
jobdir=$currdir

binarydir=$HOME/Projects/LiebSparseDiag/EXE
#binarydir=/storage/disqs/LiebSparseDiag/EXE

for disorder in 16.0 16.05 16.1 16.15 16.2 16.25 16.3 16.35 16.4 16.45 16.5 16.55 16.6 16.7 16.8 16.9
    #16.0 16.1 16.15 16.2 16.25 16.3 16.35 16.4 16.45 16.6 16.7 16.8 16.9
    #8.0 7.5 7.0 6.5 6.0 5.5 5.0 4.5 4.0 3.5 3.0 2.5 2.0 1.5 1.0 0.5
    #18.0 17.75 17.5 17.25 17.0 16.75 16.5 16.25 16.0 15.75 15.5 15.25 15.0
do

#energy=`echo "$disorder/2.0 + 4.05"| bc`
#energy=1.0

echo "--- hDis=" $disorder ", Eng=" $energy

jobdir="L31-$size-E$energy-hD$disorder"
mkdir -p $jobdir

jobname=$jobdir-$seed-$config
echo $jobname

jobfile=`printf "$jobname.sh"`
logfile=`printf "$jobname.log"`

#inpfile=LSDdiag-$size-$energy-$disorder.inp

echo "binarydir=" $binarydir " jobdir=" $jobdir 

# settings for parallel submission

cd $jobdir

cat > ${jobfile} << EOD
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=48:00:00
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
echo "IStateFlag    = 0       ">>  \$inpfile #
echo "Width0        = $size       ">>  \$inpfile #
echo "Width1        = $size       ">>  \$inpfile #
echo "dWidth        = 2          ">>  \$inpfile #
echo "HubDis0       = $disorder      ">>  \$inpfile #
echo "HubDis1       = $disorder           ">>  \$inpfile #
echo "dHubDis       = 1.0           ">>  \$inpfile #
#echo "RimDis0       = $disorder      ">>  \$inpfile #
echo "RimDis0       = 0.0            ">>  \$inpfile #
echo "Energy0       = $energy    ">>  \$inpfile #
echo "Energy1       = $energy       ">>  \$inpfile #
echo "dEnergy       = 0.1      ">>  \$inpfile #
echo "NEvals        = 100           ">>  \$inpfile #
echo "Memory        = 100          ">>  \$inpfile #

cat \$inpfile

#$binarydir/$binary <$inpfile

done

MY_PARALLEL_OPTS="-N 1 --delay .2 -j \$SLURM_NTASKS --joblog parallel-\${SLURM_JOBID}.log"
MY_SRUN_OPTS="-N 1 -n 1 --exclusive"
MY_EXEC="$binarydir/$binary <LSDdiag-$disorder-{}.inp"

parallel \$MY_PARALLEL_OPTS srun \$MY_SRUN_OPTS \$MY_EXEC ::: {1..$config}

zip -m inp.zip *.inp
zip -m sh.zip *.sh

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

