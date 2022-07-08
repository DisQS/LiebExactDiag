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

#binarydir=$HOME/Projects/LiebSparseDiag/EXE
binarydir=/storage/disqs/LiebSparseDiag/EXE

for disorder in 10.0 9.0 8.0 7.0 6.0 5.0 4.0 3.0 2.0 1.0 0.2 0.4 0.6 0.7 0.8 1.2 1.4 1.6 1.8 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 20.0 17.5 15.0 14.0 13.0 12.0 11.0 25.0 30.0 35.0 40.0 #50.0 60.0 70.0 80.0 90.0 100.0 
do

energy=`echo "$disorder/2.0 + 4.05"| bc`

echo "--- hDis=" $disorder ", Min_Eng=" $energy

jobname="LSD-$size-hD$disorder"
echo $jobname

jobfile=`printf "$jobname.sh"`
logfile=`printf "$jobname.log"`
jobdir="LSD-$size"
mkdir -p $jobdir

inpfile=LSDdiag-$disorder.inp

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
echo "Nx            = 4            ">>  $inpfile #
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
echo "Energy0       = 0.01       ">>  $inpfile #
echo "Energy1       = $energy       ">>  $inpfile #
#echo "Energy1       = 2.5           ">>  $inpfile #
echo "dEnergy       = 0.2      ">>  $inpfile #
#echo "dEnergy       = 0.05      ">>  $inpfile #
echo "NEvals        = 100           ">>  $inpfile #
echo "Memory        = 100          ">>  $inpfile #

cat $inpfile

$binarydir/$binary <$inpfile

done

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

