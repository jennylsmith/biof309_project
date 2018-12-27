# -*- coding: utf-8 -*-
"""
@author: Jenny Smith

This program will merge two ChIP-seq or ATAC-seq peak files to identify shared peaks between the two sets, and then merge this list with RNA-seq expression data. 
The purpose of this program is to easily identify genes which are present in both ChIP-seq and ATAC-seq datasets and then determine if these genes have
associated differential expression based on an RNA-seq dataset.  


Peaks_DEG command line tool will integrate ChIP-sequencing, ATAC-sequencing and RNA-sequencing datasets into a single CSV file as the output and then 
create a Venn diagram illustrating the overlap between them. 

Peaks_DEG is compatible with all peaklists and RNA-seq datasets that are in CSV file format (comma separated values).
However, this program was initially designed for use on PAPST peaklists. For more information on PAPST, please see https://github.com/paulbible/papst. 
 

The user is prompted first to choose the column headers to merge the peaklists, and then then is prompted to choose the RNA-seq column headers to 
merge on.  Please always ensure to select the same number of columns to merge on for the peaklists, e.g. if you merge peaklist1 on genbank accession
and official gene symbol columns, you must select both the corresponding columns in peaklist2. It will not affect the merging if the headers do not 
have the same names, just that  both columns contain the same information. 


"""

import sys

try:
    import pandas as pd
except NameError:
    print("You must install missing PANDAS package. Please see http://pandas.pydata.org/pandas-docs/stable/install.html for more information.")
    raise SystemExit
try:
    import matplotlib_venn
except NameError:
    print("You must install missing matplot-lib venn package. Please see https://pypi.python.org/pypi/matplotlib-venn for more information.")
    raise SystemExit
    
from matplotlib_venn import venn3 

try:
    import matplotlib.pyplot as plt 
except NameError:
    print("You must install missing matplotlib package pyplot. Please see http://matplotlib.org/users/installing.html for more information.")
    raise SystemExit
    
    
def getColumnsToMergeOn(header):
    #print the columns with an associated interger index
    print ("\n Choose the indices to merge the table on. The indices chosen must be present in all three datasets: peakslist1, peaklist2, and the rna-seq file. \n\
    \n For multiple columns, separate using a comma: \n")
    for col in enumerate(header):
        print("%d\t%s" % (col[0], col[1]))  

    #read the input from user
    column_choices = input('Indices: ')    
    print("")
    
    #split into a list of integers
    column_choices = list(map(int, column_choices.split(",")))
    
    #convert into column names
    result = []
    for i in column_choices:
        result.append(header[i])
    return result


if len(sys.argv) != 4:
    print ("Usage: %s peakfile1.csv peakfile2.csv rna-seqfile.csv" %(sys.argv[0]))
    exit()

else:
    #read in the files
    peaklist1= ""
    peaklist2= ""
    rna_seqfile= ""
    try:
         peaklist1 = pd.read_csv(sys.argv[1], header='infer', delimiter=',')
         peaklist2 = pd.read_csv(sys.argv[2], header='infer', delimiter=',')
         rna_seqfile = pd.read_csv(sys.argv[3], header='infer', delimiter=',')
    except:
        print("An error occured while reading the files. Please check that you have input a csv file.")
        raise SystemExit

    #get header from peaklist files
    header1 = list(peaklist1.columns)
    header2 = list(peaklist2.columns)
    
    try:
        columns1 = getColumnsToMergeOn(header1)
    except ValueError:
        print("ValueError: Input only the interger values that correspond to the column(s) you wish to merge on. Please choose again below.")
        columns1 = getColumnsToMergeOn(header1)   #I want the program to simply re-promt the user to start over again at columns 1, rather than crash    
    
    columns2 = getColumnsToMergeOn(header2)    
        
    if len(columns1) != len(columns2):
        print("Merge Error: you must merge peaklists on the same number of columns. Please choose again below.")
        columns2 = getColumnsToMergeOn(header2)   
    
   #merge the peaklists first
    try:    
        shared = pd.merge(peaklist1, peaklist2, left_on=columns1, right_on=columns2, how='inner')
        print("Shared peaks identified.")
        #I want a pause or wait command to allow used to read this 
    except KeyError:
        print("Please check header of peak files; both files must have indices in common to merge.")
        raise SystemExit
    
    index_shared = shared.set_index(columns1, drop=False, inplace=False) #no copywarning if I dont have this but just to make sure. 
         
    
    #merge the shared peaks with the rna_seq expression data
    header3 = list(rna_seqfile.columns) 
    
    columns3 = getColumnsToMergeOn(header3)
   
    if len(columns1) != len(columns3):
        print("Merge Error: you must merge the RNA-seq file on the same number of columns as your peaklists . Please choose again below.")
        columns3 = getColumnsToMergeOn(header3) 
        
    peaks_DEG = pd.DataFrame()
    
    try:
        peaks_DEG = pd.merge(left=index_shared, right=rna_seqfile, left_on=columns1, right_on=columns3, how='inner')
        print("Peak files merged with RNA-seq data.")
    except ValueError:
        print("Please ensure that you select the same number of columns to merge on for the shared peaks and rna-seq file.")
        raise SystemExit

    peaks_DEG_index = peaks_DEG.set_index(columns1, drop=False, inplace=False)                   
    
    
    #remove duplicates from peaklists and merged files - in preparation for Venn Diagram to be based on one differentally expressed gene, for one peak. 
    peaklist1_rmdup = peaklist1.drop_duplicates(subset=columns1, inplace=False)
    peaklist2_rmdup = peaklist2.drop_duplicates(subset=columns2, inplace=False)    
    index_shared_rmdup = index_shared.drop_duplicates(subset=columns1, inplace=False)
    peaks_DEG_rmdup = peaks_DEG.drop_duplicates(subset=columns1, inplace=False)    
     
     
    #define variables for the venn diagram 
    num_peaks_DEG = len(peaks_DEG_rmdup) 
        
    num_peaks_shared = len(index_shared_rmdup)  
        
    peaks_notDEG = num_peaks_shared - num_peaks_DEG
    
    #merge peaklist1-rna_seqfile and remove duplicates
    peaklist1_DEG = pd.merge(left=peaklist1, right=rna_seqfile, left_on=columns1, right_on=columns3, how='inner')  
    peaklist1_DEG_rmdup = peaklist1_DEG.drop_duplicates(subset=columns1, inplace=False)    
    excl_peaklist1_DEG = len(peaklist1_DEG_rmdup) - num_peaks_DEG
                               
    exclusively_peaklist1 = len(peaklist1_rmdup) - num_peaks_shared - excl_peaklist1_DEG  
    
    
    #merge peaklist2-rna_seqfile and remove duplicates
    peaklist2_DEG = pd.merge(left=peaklist2, right=rna_seqfile, left_on=columns2, right_on=columns3, how='inner')  
    peaklist2_DEG_rmdup = peaklist2_DEG.drop_duplicates(subset=columns1, inplace=False)      
    excl_peaklist2_DEG = len(peaklist2_DEG_rmdup) - num_peaks_DEG 
    
    exclusively_peaklist2 = len(peaklist2_rmdup) - num_peaks_shared - excl_peaklist2_DEG     
       
    exclusively_DEGs = len(rna_seqfile) - len(peaklist1_DEG_rmdup) - excl_peaklist2_DEG
      
    #"7-element list of subset sizes (Abc, aBc, ABc, abC, AbC, aBC, ABC), and draw a three-circle area-weighted venn diagram" from python docs
    venn = venn3(subsets = (exclusively_peaklist1, exclusively_DEGs, excl_peaklist1_DEG, exclusively_peaklist2, peaks_notDEG, excl_peaklist2_DEG, num_peaks_DEG), set_labels = (sys.argv[1],sys.argv[3],sys.argv[2] ))    
    
    
    
    outfile1 = input("Please specify the output filename [peaks_DEG]: ")
    if (outfile1==""):
        outfile1="peaks_DEG"    
    outfile1_csv = outfile1 + ".csv"
    peaks_DEG_index.to_csv(outfile1_csv, index=False)  
    print("your csv will be saved as", outfile1_csv)
    
    outfile2 = input("please specify name of venn diagram [peaks_DEG]: ")    
    if (outfile2==""):
        outfile2="peaks_DEG"
    outfile2_png = outfile2 + ".png"
    print("your venn diagram will be saved as", outfile2_png)    
    plt.savefig(outfile2_png)
    
    
    print("Saving Complete. Exiting.")
 
 
