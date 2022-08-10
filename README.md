# RNA-seq Analysis for NCG vs H37Rv and Image Analysis for BLaER1 Cells (in re: Saelans, Sweeney, Viswanathan, et al. 2022)

## RNA-seq Pipeline

This contains the scripts and other files needed to analyze the RNA sequencing files contained in (). The general pipeline is as follows:

Parse a freshly downloaded gene list from [Mycobrowser](https://mycobrowser.epfl.ch/releases) for H37Rv with cdnaParser.py. This script will take the downloaded gene list and remove the appended information in the FASTA header (the headers are formatted ">Rv#|gene name|type|loc|strand|long name" and for downstream analysis, we only need the second piece (because the Rv number is simply duplicated in the second slot if there is no generally recognized gene name. For instance:

```
>Rv1792|esxM|CDS|2030347-2030643|+|ESAT-6 like protein EsxM
```

This script will parse all the genes into a Python dictionary, replace esxM, esxJ, esxW, esxK, and esxP with their respective 5' UTR sequences (which facilitiates differentiation of these highly similar genes) and then write it out into a new fasta file, which is then used to generate the index file in Kallisto. It can be run piecewise in a Python terminal by running the contents and running with:

```
parser(ref , esx_utr) # Providing the paths as strings
```

or accessed directly from the terminal by:

```
chmod +x cdnaParser.py
./cdnaParser.py parser [reference cDNA file location] [replacement esx sequence location]
```

Small modifications could be make to accept any gene list, but this was not required for the current project. Optimally, the reference and replacement sequence files are all in the same directory as the script itself as this is where the output will go.

[Kallisto](https://github.com/pachterlab/kallisto) 0.48 was built on a 2017 5K iMac running macOS Monterey 12.5 with a modified CMakeList.txt requiring inclusion of HDF5 support. This is the latest version of Kallisto available at time of writing and standard versions of Kallisto are unable to perform bootstrapping for Sleuth quantification due to removal of this HDF5 dependency. The underlying rationale for this change is unclear, but the provided binary should work for others. This Kallisto binary is provided under the same license terms as the source (BSD 2-Clause "Simplified" License).

Kallisto can be run manually from the command line or from the kallistoRun.py script. This script serves as a wrapper to automate both populating needed information in Kallisto and then looping over all of the files. The reads for this experiment were single-end, but the script should properly support paired-end reads as well. The output will be a series of folders with abundance.tsv files in them that are then ready for processing in R.

For their own reasons, Ensembl has discontinued the production of biomaRt entities for bacterial genomes due to the immense number of bacterial genome sequences available. This introduces a new challenge for the analysis of bacterial RNA-seq data. Specifically for mycobacteria, we can use the TSVs provided with each release on Mycobrowser to construct a simplified mart-like object that can be use for Sleuth. This is an imperfect option (as it is somewhat more manual and organism-specific) but is more than adequately servicable for the time being.

Sleuth is largely just used to process and output the desired genes into a matrix, which is then plotted with pheatmap. Sleuth also includes useful interactive features for *ad hoc* visualization of the sequencing data. The plot itself was generated by loading the esx_5_tpm.csv into R, reading it into a matrix, and then plotting.

## Image Analysis in FIJI/ImageJ

This repository also contains the necessary files and information for analyzing the physical protein distribution within a stellate cell, in this case for ARPC2. The basic question revolves around whether this protein is more centrally or distally located in particular cells under particular conditions and these scripts (./surface_plot) are designed to facilitate the analysis of this sort of pixel-intensity based distribution through a combination of automated processing and manual adjustments within the 3D Surface Plot plug-in in FIJI/ImageJ. The raw measurement files are provided in ./measurements. A fuller description and the most up-to-date versions of these scripts will be able to be found at [the origin repo](https://github.com/jaredbrewer/image-analysis).

Methodological details for these images (acquisition, processing, and analysis) are provided in the source publication. Briefly, macrophage differentiated BLaER1 cells were infected with *Mycobacterium tuberculosis* mc<sup>2</sup>6020 expressing either full-length EsxM or the mCerulean fluorescent protein. These were then imaged singly at 100x. Using either __autoIsolator.py__ or __manIsolator.py__, these are then maximum intensity projected and a cellular outline is calculated from a channel that is well-distributed across the cell body. This outline is then used to mask the cell and clear any external background. Using __launchSP.py__, these cells are then visualized in the 3D Surface Plot plugin in FIJI/ImageJ, where the longest X axis is identified and the cells arranged in the XZ orientation along that axis. Smoothing is applied (8.5 in this analysis, but any is acceptable if consistent) and all axes are removed and the background is changed to black before exporting ("Save Plot" - it is advised to then save these to disk (a __plotSaver.py__ convenience script is provided), but this is not required if using __manSurfacePlotMeasure.py__ immediately afterwards) to a static 2D image, which is then used for quantitation (with __autoSurfacePlotMeasure.py__ or __manSurfacePlotMeasure.py__) by scanning for pixels with intensity values >0 on a binarized image and writing those coordinates to a .csv, which can then be further analyzed in R. The __auto__ scripts work well together while the __man__ scripts cooperate well. However, it is principally possible to use them in different arrangements. The __man__ scripts all work from the currently selected image while the __auto__ scripts work from provided directories of compatible files.

The analysis in R consists of identifying the maximum Y coordinate at each pixel and using that to calculate the area under the curve for each cell across a defined interval of the cell (in percent X distance). Here, our analysis compared the first 25% of the AUC for each condition, normalized to the AUC for that cell from 25-75%. This allows comparison of the distal/central fluorescence intensity ratio for each cell, which is then plotted in ggplot2. The script allows for arbitrary numbers of cells to be analyzed and compared at once, provided each cell is given a discrete identifier.
