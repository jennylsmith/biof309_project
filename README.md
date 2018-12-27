# biof309_project

## Final Project

This program will merge two ChIP-seq or ATAC-seq peak files to identify shared peaks between the two sets, and then merge this list with RNA-seq expression data.

The purpose of this program is to easily identify genes which are present in both ChIP-seq and ATAC-seq datasets and then determine if these genes have associated differential expression based on an RNA-seq dataset.  
Peaks_DEG command line tool will integrate ChIP-sequencing, ATAC-sequencing and RNA-sequencing datasets into a single CSV file as the output and then  create a Venn diagram illustrating the overlap between them. 

Peaks_DEG is compatible with all peaklists and RNA-seq datasets that are in CSV file format (comma separated values).
However, this program was initially designed for use on PAPST peaklists. For more information on PAPST, please see https://github.com/paulbible/papst. 
 
The user is prompted first to choose the column headers to merge the peaklists, and then then is prompted to choose the RNA-seq column headers to merge on.  Please always ensure to select the same number of columns to merge on for the peaklists, e.g. if you merge peaklist1 on genbank accession and official gene symbol columns, you must select both the corresponding columns in peaklist2. It will not affect the merging if the headers do not have the same names, just that  both columns contain the same information. 

