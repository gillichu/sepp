
Introduction 
===

UPP2 is a wrapper around UPP, and is subsequently also a method for performing alignments of ultra-large and fragmentary datasets. Please find more details about [UPP on its respective repository](https://github.com/smirarab/sepp/blob/master/README.UPP.md). We have two key strategies to note: 

* UPP-Hierarchical
* UPP-EarlyStop

For a particular query sequence and a particular HMM, we've made some changes to allow UPP to recalculate the bitscore by which it finds the "best" HMM match for a particular query sequence. This recalculation now takes into account the number of sequences summarized in the HMM. We call this the adjustedUPP strategy. 

According to this new weighting the, we also present several strategies to redesign how the ensembles of HMMs are organized, taking into account the phylogenetic information. That is, we explored several strategies that allow us to incorporate more sequences or fewer sequences, and generate overlapping HMMs that each query sequence can be tested against. 

We have also implemented a fastUPP strategy which does a fast search down the hierarchy of trees. That is, recursively from the root, we choose the child HMM which best fits our current query sequence, continuing until the leaf. 

### Inputs

The main pipeline for UPP2 can be found in the script `run_upp.py`. We recommend the use of a config file as followin
```
[commandline]
sequence_file=<query sequence file>
alignment=<backbone alignment>
tree=<backbone tree>
alignmentSize=<decomposition size>
molecule=<amino/dna/rna>
cpu=16
tempdir=./tmp/
outdir=./output/

[upp2]
decomp_only=True
bitscore_adjust=True
hier_upp=True
early_stop=<True/False>
```
The command to run UPP2 therefore follows the format: 
`python run_upp.py -c <Config file>`

### Outputs
The main output of UPP2 is `output_alignment.fasta` located in the output directory. UPP2 also outputs the time it took to perform its different stages.
