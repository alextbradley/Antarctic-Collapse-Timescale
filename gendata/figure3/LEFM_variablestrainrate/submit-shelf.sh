#!/bin/bash 

#loop over shelves

STEP=1
for SHELFNO in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42
do
   echo "$SHELFNO"
   #for each one, submit a job 

   JOBNO="CT_$SHELFNO-$STEP"
   OUTNAME="CT_$SHELFNO-$STEP.out"
   echo $JOBNO
   sbatch -J $JOBNO \
	  --output=$OUTNAME \
	  --export JOBNO=$JOBNO,SHELFNO=$SHELFNO,STEP=$STEP \
	  ./run-shelf.sh
done



