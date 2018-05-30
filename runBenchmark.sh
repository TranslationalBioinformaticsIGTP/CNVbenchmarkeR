#!/bin/bash
# Description: Run the benchmark analysis for those algorithms and dataset defined in the algorithms.yaml and datasets.yaml file

# Read params
source ./utils/parse_yaml.sh
eval $(parse_yaml algorithms.yaml "pars_")

# create logs/output folder if not exists
mkdir -p logs
mkdir -p output

# Execute algorithms over selected datasets
if [ "$pars_algorithms_panelcn" == "true" ]; then
    echo "[$(date)] Executing panelcn.MOPS"
	Rscript ./algorithms/panelcnmops/runPanelcnmops.r ./algorithms/panelcnmops/panelcnmopsParams.yaml datasets.yaml  > logs/panelcnmops.log 2>&1
fi

if [ "$pars_algorithms_decon" == "true" ]; then
    echo "[$(date)] Executing DECoN"
	Rscript ./algorithms/decon/runDecon.r ./algorithms/decon/deconParams.yaml datasets.yaml  > logs/decon.log 2>&1
fi

if [ "$pars_algorithms_exomedepth" == "true" ]; then
    echo "[$(date)] Executing ExomeDepth"
	Rscript ./algorithms/exomedepth/runExomedepth.r ./algorithms/exomedepth/exomedepthParams.yaml datasets.yaml  > logs/exomedepth.log 2>&1
fi

if [ "$pars_algorithms_codex2" == "true" ]; then
    echo "[$(date)] Executing CODEX2"
	Rscript ./algorithms/codex2/runCodex2.r ./algorithms/codex2/codex2Params.yaml datasets.yaml  > logs/codex2.log 2>&1
fi

if [ "$pars_algorithms_convading" == "true" ]; then
    echo "[$(date)] Executing CoNVaDING"
	Rscript ./algorithms/convading/runConvading.r ./algorithms/convading/convadingParams.yaml datasets.yaml  > logs/convading.log 2>&1
fi


# Generate summary file
echo "[$(date)] Generating summary file"
Rscript ./utils/summary.r algorithms.yaml datasets.yaml  > logs/summary.log 2>&1