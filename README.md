# Data-analysis-bioinformatics-Python
Scientific Research that analyzes genome data (text files that describe and show the genome) 
Text archives that were analyzed:
schmidtea_mediterranea.gff3 
schmidtea_mediterranea.fasta
https://drive.google.com/drive/folders/1VSzoSqiezh12fcjMK9AYbkOHmtNbdLXO?usp=sharing

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
