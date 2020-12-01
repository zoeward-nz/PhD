# A Pipeline developed by Zoe Ward 

To analyse short read RNA-Seq for (novel)lncRNA and circRNA analysis


##Prereqs

* The workflow depends on cpat files that have to be predownloaded from cpat 
sourceforge. 
* Circexplorer2 files are needed (see data) and can be created with circexplorer2 fetch_ucsc.py 
script.
* Adjust the config file where necessary.


## Usage
```python
snakemake --use-conda --jobs 10 --snakefile Snakefile --configfile config.yaml 
```

## Workflow bottlenecks
* 2-pass STAR mapping takes time
* 1st pass on average 10min / sample
* Salmon expression quantification takes around 15min / sample with --gcBias

