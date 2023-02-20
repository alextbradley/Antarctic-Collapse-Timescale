NR=10000
NROW=1

FNAME="circum_Antarctic_collapse_time($NR,$NROW)"
echo $FNAME
matlab -nodisplay -nodesktop -r $FNAME
