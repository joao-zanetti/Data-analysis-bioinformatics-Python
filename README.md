# Data-analysis-bioinformatics-Python
Scientific Research that analyzes genome data (text files that describe and show the genome) 

## main.py:
Call all the functions and static methods.

## reagff.py:
Read the archive GFF,separate into different files: genome.csv,the genes.csv,exon.csv and introns.csv.<br/>
Each one with their respectives positions,size,start,end,region,sinal of tape  in genome fasta (archive.fasta)

## intronclass.py:
Read the intron.csv with the informations of introns (created by readgff.py).<br/>
Create the intronfasta.csv with their respectives sequences (from .fasta archive).<br/>
Find the not minimums and not minimums introns. (minimum  is the introns that appears 25% times compared with modal intron) <br/>
Creates other archives with the edges of introns (minimum, not minimum and  in intron) and their respectives sequences<br/>
In intron is the 20 intron tip positions<br/>
Calculate the proportion of each bases and the quantity of information (Shannon Entropy) for each position in all introns (minimum, not minimum and in intron)<br/>



## To compile the code:
gcc mulmatrix.c -o mulmatrix -pthread<br/>

## To run the code:
./mulmatrix<br/>

## Details:
Inside the code, edit the numbers of defines at the beginning of the file, with numbers whatever you want.<br/>

### example:<br/>
```
#define nt 4    //is the number of threads is created in execution.
#define nt 1000 //N:is the order of the matrix square
#define seed 100 //is the  seed of rand(), or seed ==100 ->  numbers between 0 - 100 in matrix
```

### Performance:<br/>
For better load balancing, the numbers of threads (nt) selected in define must be higher than number of processors:<br/> nt > number of processors.<br/>
This improve lowest idleness of processors<br/>
