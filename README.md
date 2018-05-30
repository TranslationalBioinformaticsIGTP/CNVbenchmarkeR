# CNVbenchmarkeR #

CNVbenchmarkeR is a framework to benchmark algorithms when detecting germline copy number variations (CNVs) against different NGS datasets. Current version supports DECoN, CoNVaDING, panelcn.MOPS, ExomeDepth and CODEX2 tools.


### Prerequisites ###

Algorithms have to be properly installed. Links for algorithms installation:

- https://github.com/bioinf-jku/panelcn.mops
- https://molgenis.gitbooks.io/convading/
- https://github.com/RahmanTeam/DECoN
- https://github.com/yuchaojiang/CODEX2
- https://cran.r-project.org/web/packages/ExomeDepth/index.html


### How to use
1. Get Code
```
git clone https://github.com/TranslationalBioinformaticsIGTP/CNVbenchmarkeR 
```

2. Configure algorithms.yaml to set which algortithms will be benchmarked, and configure datasets.yaml to set against which datasets the algorithms will be executed. 
In case of executing DECoN, modify algorithms/decon/deconParams.yaml by setting deconFolder to your DECoN folder installation. In case of executing CoNVaDING, modify algorithms/convading/convadingParams.yaml by setting convadingFolder param.


3. Launch CNVbenchmarker
```
cd CNVbenchmarkerR
./runBenchmark.sh
```


### Output ###

A summary file and a .csv results file will be generated at output/summary folder. Stats include sensitivity, specificity, no-call rate, precision (PPV), NPV, F1, MCC and kappa coefficient.

Stats are calculated per ROI, per gene and at whole strategy level (gene level including no-calls, i. e., low quality regions)

Logs files will be generated at logs folder. Output for each algorithm and dataset will be generated at output folder.



## Extra feature: optimizer ##

An optimizer is also attached in the framework. It executes a CNV calling algorithm against a dataset with many different values for each param.
Up to 22 values are evaluated for each param. It is implemented using a greedy algorithm which starts from each different param. The CNV algorithm will be executed a maximum of (n_params^2)\*22 times. 

It will be improve sensitivity allowing drops of specificity defined at optimizerParams.yaml.


### Prerequisites ###

An SGE cluster system has to be available.

### How to use

1. Configure optimizers/optimizerParams.yaml by defining optimizer params, dataset and algorithm to be optimized.
2. Execute optimizer:
```
cd optimizers
Rscript optimizer.r optimizerParams.yaml
```
