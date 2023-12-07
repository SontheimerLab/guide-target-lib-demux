###System requirements				
##Software dependencies and operating systems (including version numbers)			
##Versions the software has been tested on			
Red Hat Enterprise Linux 8		
python3 v.3.11.4		
biopython v.1.78 		
python-levenshtein v.0.12.2		
##Any required non-standard hardware			
No		
				
###Installation guide				
	"install python3, biopython, python-levenshtein into a conda environment"			
	"Typical install time on a """"normal"""" desktop computer"			
	<1 hour			
				
Demo and instructions for use				
	Instructions to run on data			
		Copy demux.py and demuxRef.csv into the directory with a sample FASTQ files for read 1 and read 2. Run the following command		
		python3 demux.py <sample_R1_001.fastq.gz> <sample_R2_001.fastq.gz>		
				
	Expected output			
		"demultiplexed FASTQ files labeled by library member ID, <Nme2_A001_matched.fastq>, <Nme2_A002_matched.fastq>"		
		"CSV file containing the IDs of matched spacer and barcode sequnces, <output_list.csv>"		 
		"CSV file containing the nucleotide sequences of unmatched spacer and barcode sequences, <no_match_list.csv>"		
				
	Expected run time for demo on a normal desktop computer			
	4 hours			
