# QASDRA
#Quality Assessment of Sequencing Data via Range Analysis
Ali Fotouhi
Research Assistant - Istanbul Technical University  fotouhi@itu.edu.tr
  
GENERAL INFO:

QASDRA_cp1 is a tool for analyzing fastq files using their phred quality scores. this tool gives you general and overall insight about the file you have in hand. It is well known that long intervals having fewer errors improve the performances of the post-processing tools in the down-stream analysis of the DNA sequencing data. We investigate another dimension for the quality assessment motivated with the fact that reads including long intervals having fewer errors improve the performances of the post-processing tools in the down-stream analysis.

this tool has been implemented using Inverse Range Query algorithm, which we are trying to find the longest read segments with at most k base under user defined quality value v. at the end of the program this tool will provide graphical and analytical results in a single page pdf file.
this tool has been written in python 2.7 under Linux Ubuntu 16.04.1 LTS and before running this program make sure your system has requiremnts to run this tool.

Requirments are as below:

In order to provide graphical results, this program uses matplotlib (http://matplotlib.org/1.3.1/index.html) make sure this library is already installed on your system. A complete installation guide can be found at the link below. but most of the platforms have this tool by default   
install documentation: http://matplotlib.org/1.3.1/users/installing.html

In order to generate an official report, this tool uses report lab (http://www.reportlab.com/) which most of the time is already included in python.
install documentation: http://www.reportlab.com/software/downloads/

and finally, this tool uses python imaging library. (PIL)
 
HOW IT WORKS:

In order to use this tool efficiently, it is better to get a pre-insight about your file by running "QASDRAinfo.py" on your data. this tool will tell you about the maximum, minimum, and average lengths of reads in the file and also maximum and minimum Phred quality score. so you will choose 'v' value wisely by knowing about the ranges of Phred quality scores. and by knowing read lengths you can decide about the best value for k.
in order to run this tool you should write this in command line on mac or Linux: python QAVARAinfo.py filename.fastq
which file name is the name of your file that you want to be processed.

After finding out about k and v values now you can run the QASDRA_cp1 tool in the command line:
python QASDRA_cp1.py filename.fastq k v [-r rvalue -fl flvalue -fs fsvalue -fa favalue -b bvalue -p]
an example can be like this:
python QASDRA_cp1.py filename.fastq 2 20 -r 50 -fl 20 -fs 10 -fa 15 -b 33 -p
which we will explain them separately

before using this tool you can also get help by this command line:
python QASDRA_cp1.py -h
this command briefly will tell you about the values you should enter.

Here we will explain each of these values

Mandatory values:

1. tool name: QASDRA_cp1.py should be written in order to run this tool
2. file name: fastq file's name that you want to be processed should be written after tool name with its format.fastq
3. k value: this value indicates the number of low-quality values that user wants to be in read segments which is a number.
4. v value: this is the lowest quality score that user wants to use as the threshold.
for example k=2 and v=20 indicates that there can be at most just 2 values less than or equal to 20 in longest ranges that this tool will detect using inverse range query algorithm.

we call these longest ranges "Maximal Ranges" which means that these ranges are the only possible ranges that can be long enough with k kvalues less than or equal to v value.
further information about this algorithm can be found in this paper: "https://www.researchgate.net/profile/M_Kuelekci/publications?sorting=recentlyAdded"


Optional values:

-r rvalue: random sampling, this program can be slow for largest files and in order to get a quicker response user can apply random sampling. for this, user will use "-r rvalue" command. which rvalue is a number (integer) between 0 and 100, which approximately is the percentage of the sampling. for instance rvalue=50 means 50% of reads will be processed or rvalue=30 means that this tool will process approximately 30% of reads.

-fl flvalue: longest length filtration, by including this feature user aims to filter some unwanted reads. flvalue is the minimum longest segment length (Maximal Range Length) that user wants values bigger than this value to be included during process, this means that segments smaller than lfvalue will be filtered. again fl is a number (an integer) which can vary according to length of reads.

-fs fsvalue: shortest length filtration, fsvalue is the minimum shortest segment length that user wants values bigger than this value to be included during the process, this means that segments smaller than lsvalue will be filtered. again fs is a number (an integer) which can vary according to the length of reads.

-fa favalue: average length filtration, favalue is the minimum average segment length that user wants values bigger than this value to be included during the process, this means that segments smaller than lavalue will be filtered. again fa is a number (an integer) which can vary according to length of reads.

-b bvalue: ASCII-BASE Phred quality scores, as you know most of the fastq files are 33 ASCII-BASE but still some of the centers produces 64 ASCII-BASE files. this tool works on both of the. if the user knows the file's ASCII-BASE, can enter using -b command. if no, the tool itself will find the exact ASCII-BASE. for instance for 33 ASCII-BASE one can enter common like this: -b 33.

-p: this command which gets no value is for plotting results. by default this tool uses percentage base plotting, but if one wants to know about the distribution of these segments lengths by including this command instead of percentage base plots can have distribution results. percentage base also depicts distributions of lengths.
 
"a positive aspect of this tool is filtration of longest, shortest and the average length of maximal ranges would not affect the overall analysis. it will just affect the analysis of related metrics. it means longest length filtration will just affect the analysis of longest maximal ranges, shortest length filtration will just affect the analysis of shortest maximal ranges and average length filtration will just affect the analysis of average maximal ranges, for instance, changes will ba applied to associated plots and figures."

finally, this tool will generate a pdf file encompassing results of the quality assessment of sequencing data. furthermore, terminal based results will be printed and graphical results will be produced.

Cite this paper as:
Fotouhi A., Majidi M., Külekci M.O. (2018) Quality Assessment of High-Throughput DNA Sequencing Data via Range Analysis. In: Rojas I., Ortuño F. (eds) Bioinformatics and Biomedical Engineering. IWBBIO 2018. Lecture Notes in Computer Science, vol 10813. Springer, Cham DOI https://doi.org/10.1007/978-3-319-78723-7_37

in the case of any ambiguity do not hesitate to mail me at fotouhi@itu.edu.tr. I will be grateful to know about this tool's flaws and will try to fix them.

#Ali Fotouhi
#Research Assistant - Istanbul Technical University
#office: DAMGA Lab - UYBHM Building - ITU Ayazaga Campus - 34469 Istanbul-Turkey
#E-mail: fotouhi@itu.edu.tr
