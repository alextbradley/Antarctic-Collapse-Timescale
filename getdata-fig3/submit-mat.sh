#!/bin/bash

# set parameters
NR=25
NROW=11

JOBNO="CT_$NR-$NROW"
OUTNAME="CT_$NR-$NROW-.out"
echo $JOBNO
sbatch -J $JOBNO \
	--output=$OUTNAME \
	--export JOBNO=$JOBNO,NR=$NR,NROW=$NROW \
	./run.sh



