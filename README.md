# RNA-seq Analysis for NCG vs H37Rv (in re: Saelans, Sweeney, Viswanathan, et al. 2022)

This contains the scripts and other files needed to analyze the RNA sequencing files contained in (). The general pipeline is as follows:

Parse a freshly downloaded gene list from [Mycobrowser](https://mycobrowser.epfl.ch/releases) for H37Rv with cdnaParser.py. This script will take the downloaded gene list and remove the appended information in the FASTA header (the headers are formatted ">Rv#|gene name|type|loc|strand|long name" and for downstream analysis, we only need the second piece (because the Rv number is simply duplicated in the second slot if there is no generally recognized gene name. For instance:

```
>Rv1792|esxM|CDS|2030347-2030643|+|ESAT-6 like protein EsxM
```

This script will parse all the genes into a Python dictionary, replace esxM, esxJ, esxW, esxK, and esxP with their respective 5' UTR sequences (which facilitiates differentiation of these highly similar genes) and then write it out into a new fasta file, which is then used to generate the index file in Kallisto. It can be run piecewise in a Python terminal or accessed directly from the terminal by:

```
chmod +x cdnaParser.py
./cdnaParser.py parser [reference cDNA file location] [replacement esx sequence location]
```

Small modifications could be make to accept any gene list, but this was not required for the current project. Optimally, the reference and replacement sequence files are all in the same directory as the script itself as this is where the output will go.

Kallisto can be run manually from the command line or from the kallistoRun.py script. This script serves as a wrapper to automate both populating needed information in Kallisto and then looping over all of the files. The reads for this experiment were single-end, but the script should properly support paired-end reads as well. The output will be a series of folders with abundance.tsv files in them that are then ready for processing in R.

For their own reasons, Ensembl has discontinued the production of biomaRt entities for bacterial genomes due to the immense number of bacterial genome sequences available. This introduces a new challenge for the analysis of bacterial RNA-seq data. Specifically for mycobacteria, we can use the TSVs provided with each release to construct a simplified mart-like object that can be use for Sleuth. This is an imperfect option (as it is somewhat more manual and organism-specific) but is more than adequately servicable for the time being. 

Sleuth is largely just used to process and output the desired genes into a matrix, which is then plotted with pheatmap. 
