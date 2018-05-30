#$ -S /bin/bash
#$ -V


# Actual work
echo "Running: ${0}"

#get the enviroment
HOST=`hostname -s`

#echo the parameters
echo "################################################################################"
echo "${0}"
echo "Run at: $HOST"
date
echo "###############################################################################"
echo ""


cd *rootPath*
Rscript *algorithmScript* $algorithmParams $datasetParams


echo "Finished" 
echo ""
echo ""
echo ""
