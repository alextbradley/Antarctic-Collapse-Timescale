#!/bin/bash 

%loop over shelves

%for each one, submit a job 
SHELF="PineIslandFast"
STEP=1

JOBNO="CT_$SHELF-$STEP"
OUTNAME="CT_$SHELF-$STEP.out"
echo $JOBNO
sbatch -J $JOBNO \
	--output=$OUTNAME \
	--export JOBNO=$JOBNO,SHELF=$SHELF,STEP=$STEP \
	./run-shelf.sh
