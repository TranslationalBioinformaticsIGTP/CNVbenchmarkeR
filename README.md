**NOTE**: New version supporting 12 tools (CNVbenchmarkeR2) can be found [here](https://github.com/jpuntomarcos/CNVbenchmarkeR2). 

# CNVbenchmarkeR #

CNVbenchmarkeR is a framework to benchmark algorithms when detecting germline copy number variations (CNVs) against different NGS datasets. Current version supports DECoN, CoNVaDING, panelcn.MOPS, ExomeDepth and CODEX2 tools.

It is part of our [publication](https://www.nature.com/articles/s41431-020-0675-z) in which we performed a benchmark of germline CNV calling tools for targeted gene-panel data. Citation:
Moreno-Cabrera, J.M., del Valle, J., Castellanos, E. et al. Evaluation of CNV detection tools for NGS panel data in genetic diagnostics. Eur J Hum Genet (2020). https://doi.org/10.1038/s41431-020-0675-z

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

2. **Configure algorithms.yaml** to set which algortithms will be benchmarked. In case of executing DECoN, modify algorithms/decon/deconParams.yaml by setting deconFolder to your DECoN folder installation. In case of executing CoNVaDING, modify algorithms/convading/convadingParams.yaml by setting convadingFolder param.

3. **Configure datasets.yaml** to define against which datasets the algorithms will be executed. Within this file, it is important to provide files with the exact expected format (**special attention** to `validated_results_file` and `bed_file` that are **tab-delimited** files). To do so, please **check the [examples](https://github.com/TranslationalBioinformaticsIGTP/CNVbenchmarkeR/tree/master/examples) folder**.


4. Launch CNVbenchmarker
```
cd CNVbenchmarkerR
./runBenchmark.sh
```


### Output ###

A summary file and a .csv results file will be generated at output/summary folder. Stats include sensitivity, specificity, no-call rate, precision (PPV), NPV, F1, MCC and kappa coefficient.

Stats are calculated per ROI, per gene and at whole strategy level (gene level including no-calls, i. e., low quality regions)

Logs files will be generated at logs folder. Output for each algorithm and dataset will be generated at output folder.


### Troubleshooting  ###

Two important checks to ensure that metrics are computed correctly:

- The **sample names in the `validated_results_file` should match the file names of your bam files** (excluding the .bam extension). For example, if the `validated_results_file` contains sample names like mySample2312, your bam files should have file names like mySample2312.bam .
- Provide and use chromosomes names with the same format, that is, do not use "chr5" and "5" in you bed and `validated_results_file` files, for example.


## Extra feature: optimizer ##

An optimizer is also attached in the framework. It executes a CNV calling algorithm against a dataset with many different values for each param.
Up to 22 values are evaluated for each param. It is implemented using a greedy algorithm which starts from each different param. The CNV algorithm will be executed a maximum of (n_params^2)\*22 times. 

It will be improve sensitivity allowing drops of specificity defined at optimizerParams.yaml.


### Prerequisites ###

An SGE cluster system has to be available.

### How to use

1. Configure optimizers/optimizerParams.yaml by defining optimizer params, dataset and algorithm to be optimized. Note: it is recommended to optimize over a random subset (training subset) of the original subset. Then, performance can be compared on the validation subset.
2. Execute optimizer:
```
cd optimizers
Rscript optimizer.r optimizerParams.yaml
```
