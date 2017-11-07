This code deconvolutes FASTQ files from Illumina sequencing for screens that use the Big Papi vector.<br/>
Requires:<br/>
Python-2.7<br/>
pandas <br/>


Files required:<br/>
• Gzipped FASTQ file with construct reads<br/>
• Gzipped FASTQ file with sample barcodes<br/>
• Reference file with U6 guides. This file should be in .csv format and should not have any headers. The first column should have the barcodes and the second column should have the construct IDs.<br/>
• Reference file with H1 guides. This file should be in .csv format and should not have any headers. The first column should have the barcodes and the second column should have the construct IDs.<br/>
• File with sample conditions. This file should be in .csv format and should not have any headers. The first column should have the barcodes and the second column should have the sample condition.<br/>


To run this code type the following on the terminal:<br/>
<b>python big_papi_deconvolution.py --read-file \<FASTQ construct reads file\> --barcode-file \<FASTQ sample barcodes file\> --ref-u6 \<File with U6 guides\> --ref-h1 \<File with H1 guides\> --cond-file \<File with sample conditions\> --outputfile \<.txt output file name\> <br/></b>

This code prints summary statistics to standard output at the end of the run. It also outputs a tab-delimited file with read counts for all combinations of U6-H1 guides in all samples. 


