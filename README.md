# Data-analysis-bioinformatics-Python
Scientific Research that analyzes genome data (text files that describe and show the genome) <br/>
From these archives, extract,search and analyzes informations.<br/>
Find the intron, exon,gene, primes and intron edges sequences and informations.<br/>
Calculates the proportion of each position of nucleotides, respective quantity of informations (Shannon Entropy) and plot graphics<br/>

## Text archives that were used to analyze:
#### schmidtea_mediterranea.gff3:
3821230 lines that describe  the start and end positions of regions, gene, exons and primes  of the other archive schmidtea_mediterranea.fasta<br/>
### schmidtea_mediterranea.fasta:
15091651 lines wich show the nucleotide sequence separated for region of genome<br/> 
#### the files link:
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
