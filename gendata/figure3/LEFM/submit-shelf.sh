#!/bin/bash 

#loop over shelves
#declare -a arr=("PineIslandFast" "Thwaites" "Amery")
declare -a arr=("PineIslandFast")

STEP=8
for SHELF in "${arr[@]}"
do
   echo "$SHELF"
   #for each one, submit a job 

   JOBNO="CT_$SHELF-$STEP"
   OUTNAME="CT_$SHELF-$STEP.out"
   echo $JOBNO
   sbatch -J $JOBNO \
	  --output=$OUTNAME \
	  --export JOBNO=$JOBNO,SHELF=$SHELF,STEP=$STEP \
	  ./run-shelf.sh
done



