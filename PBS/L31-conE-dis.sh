#!/bin/bash

# settings from input

energy=${1:-1.0}
size=${2:-10}
disorder0=${3:-1.0}
seed=${4:-1}
config=${5:-2}
keep=${6:-1}

echo "LSD: making for E=" $energy "with M=" $size " DISORDER=" $disorder0 "and starting seed=" $seed "and" $config "samples"

# settings for files

binary=LSDdiag.IC

# settings for directories

currdir=`pwd`
jobdir=$currdir

#binarydir=$HOME/Projects/LiebSparseDiag/EXE
binarydir=/storage/disqs/LiebSparseDiag/EXE

for disorder in $disorder0
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

inpfile=LSDdiag-$size-$energy-$disorder.inp

echo "binarydir=" $binarydir " jobdir=" $jobdir 

# settings for parallel submission

cd $jobdir

cat > ${jobfile} << EOD
#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2012
#SBATCH --time=48:00:00

module purge
module load intel
#module load DMTCP/2.5.2

for iseed in {1..$config..1}
do

myseed=\$(( $seed + \$iseed - 1))
echo "--- working on config" \$iseed "with seed" \$myseed

# create the input file
echo create the input file
touch $inpfile

echo "ISeed         = \$myseed       ">  $inpfile #
echo "NConfig       = 1        ">>  $inpfile #
echo "Dim           = 3            ">>  $inpfile #
echo "Nx            = 1            ">>  $inpfile #
echo "IBCFlag       = 1             ">>  $inpfile #
echo "IRNGFlag      = 0             ">>  $inpfile #
echo "IKeepFlag     = $keep      ">>  $inpfile #
echo "IWriteFlag    = 2       ">>  $inpfile #
echo "IStateFlag    = 0       ">>  $inpfile #
echo "Width0        = $size       ">>  $inpfile #
echo "Width1        = $size       ">>  $inpfile #
echo "dWidth        = 2          ">>  $inpfile #
echo "HubDis0       = $disorder      ">>  $inpfile #
echo "HubDis1       = $disorder           ">>  $inpfile #
echo "dHubDis       = 1.0           ">>  $inpfile #
#echo "RimDis0       = $disorder      ">>  $inpfile #
echo "RimDis0       = 0.0            ">>  $inpfile #
echo "Energy0       = $energy    ">>  $inpfile #
echo "Energy1       = $energy       ">>  $inpfile #
echo "dEnergy       = 0.1      ">>  $inpfile #
echo "NEvals        = 100           ">>  $inpfile #
echo "Memory        = 100          ">>  $inpfile #

cat $inpfile

$binarydir/$binary <$inpfile
rm -f $inpfile

done

zip -m inp.zip *.inp
#zip -m sh.zip *.sh

#exit 0

EOD

chmod 755 ${jobfile}
##(msub -q devel $jobdir/${jobfile}) # for queueing system
##(sbatch -q devel $jobdir/${jobfile}) # for queueing system
#sbatch ${jobfile} # for queueing system
##(source $jobdir/${jobfile} ) >& $jobdir/${logfile} & # for parallel shell execution
source ${jobfile} >& ${logfile} # for sequential shell execution

#echo "<return>"
sleep 1

cd ..

done

