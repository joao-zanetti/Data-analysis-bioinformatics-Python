##Developed by João Zanetti
##https://github.com/joao-zanetti

import re

#This function read the archive GFF
#Separate in the genome,the genes,exon and introns
#Their respectives positions,size,start,end,region,sinal of tape  in genome fasta (archive.fasta)
def readgff(limite):#limite == how many lines of GFF it will be read
	archivegff=open('schmidtea_mediterranea.gff3','r')
	lc=0 #loop controller
	nexons=0 #number of exon/gene
	nintron=0
	primerfive='' #5' primer controller
	primerthree=''#3' primer controller
	genepc=0 #gene protein coding controller
	introns= []
	ngenes=0
	v=','
	genedata='' #gene data
	genearchive=open('gene.csv','w') #open write gene archive
	exonarchive=open('exon.csv','w') #open write exon archive
	intronarchive=open('intron.csv','w') #open write intron archive

	for line in archivegff: 
		lineclear= re.sub('\n','',line) #remove \n for ''
		linesplit= re.split('\t|;| ',lineclear) #remove '\t' ';' ','' espace from read line
		result= re.match("#",lineclear) #check '#'  
		if(result!=None): #if is '#'
			print('#')
		elif((linesplit[1]=="dust") or (linesplit[1]=='tandem') or (linesplit[1]=='RepeatMasker')): #check dust or tandem or repeatmasker 
				print("D T R M")
		elif(linesplit[2]=='exon' and genepc==1): #if is a exon and gene PC
			nexons= nexons+1
			if(nexons>1): #if is a second or high exon, readed in gene
				idexon= re.sub('ID=exon:','',linesplit[8])
				parentexon= re.sub('Parent=transcript:','',linesplit[9])
				size= int(linesplit[4])-int(linesplit[3])
				size= str(size)
				exondata=idexon+v+parentexon+v+size+v+linesplit[3]+v+linesplit[4]+v+linesplit[6]+'\n'
				exonarchive.writelines(exondata) #WRITE EXON
				print ('EXON'+linesplit[6])
				if(linesplit[6]=='-'): #if exon in a negative tape AND CALCULATE/WRITE INTRON
					finalintron=int(exonanterior[3]) #finalpos of intron = initial pos of previus exon
					initintron=int(linesplit[4]) #initial pos of intron = final pos of current exon
					tamintron=finalintron-initintron
					initintron=initintron+1
					finalintron=finalintron-1
					tamintron=str(tamintron)
					finalintron=str(finalintron)
					initintron=str(initintron)
					introndata=sequence_region+v+parentexon+v+tamintron+v+initintron+v+finalintron+v+linesplit[6]
					introns.append(introndata)
					print ('INTRON-')
					exonanterior=linesplit
				else: #if exon in a positive tape
					nintron=nexons-1
					nintron=str(nintron)
					initintron=int(exonanterior[4]) #initial pos of intron = final pos of previous exon
					finalintron=int(linesplit[3]) #finalpos of intron = initial pos of current exon 
					tamintron=finalintron-initintron
					initintron=initintron+1
					finalintron=finalintron-1
					tamintron=str(tamintron)
					finalintron=str(finalintron)
					initintron=str(initintron)
					introndata=sequence_region+v+parentexon+v+tamintron+v+initintron+v+finalintron+v+linesplit[6]+v+nintron+v+'oi'+'\n'
					intronarchive.writelines(introndata) #WRITE INTRON+
					print ('INTRON+')
					exonanterior=linesplit
			else: #if its first exon read from gene
				idexon= re.sub('ID=exon:','',linesplit[8])
				parentexon= re.sub('Parent=transcript:','',linesplit[9])
				size= int(linesplit[4])-int(linesplit[3])
				size= str(size)
				exonanterior=linesplit
				exondata=idexon+v+parentexon+v+size+v+linesplit[3]+v+linesplit[4]+v+linesplit[6]+'\n'
				exonarchive.writelines(exondata) #WRITE FIRST EXON OF GENE
				print('EXON'+linesplit[6])
		elif(linesplit[2]=='gene'): # if is a gene
					if(linesplit[10]=='biotype=protein_coding'): #if is gene PC
						sequence_region=linesplit[0]
						print(sequence_region) 
						if(nexons>1 and exonanterior[6]=='-'):
							isize=len(introns)
							print(isize) 
							count=isize
							counti=0
							while(counti<isize):
								nintronw=str(count)
								print(nintronw)
								intronwrite=introns[counti]+v+nintronw+v+'oi'+'\n'
								intronarchive.writelines(intronwrite) #WRITE INTRON-
								count=int(count)
								count-=1
								counti+=1
						introns= []
						genepc=1
						nexons=str(nexons)
						genedata=genedata+nexons+primerfive+primerthree+'\n'
						if(ngenes>0):
							genearchive.writelines(genedata) #Write a gene when the second gene is being reading
						ngenes+=1
						primerfive=',N' #start primeflag of a current gene read
						primerthree=',N' #start primeflag of a current gene read
						nexons='0' #start nexons of a current gene read
						nexons=int(nexons)
						nintron=0
						idgene= re.sub('ID=gene:','',linesplit[8])
						namegene= re.sub('Name=','',linesplit[9])
						size= int(linesplit[4])-int(linesplit[3])
						size= str(size)
						genedata=idgene+v+namegene+v+size+v+linesplit[3]+v+linesplit[4]+v+linesplit[6]+v
						print('GENE PC') 
					else: #IF IS NOT PC
						genepc=0
						print('GENE NOT PC')		
		elif(linesplit[2]=='five_prime_UTR'): #if is a five_prime_UTR
			primerfive=',Y' # of a current gene read
			print('PRIMERFIVE') 
		elif(linesplit[2]=='three_prime_UTR'): #if is a five_prime_UTR
			primerthree=',Y' # of a current gene read
			print('PRIMERTHREE')
					

		lc=lc+1
		if(lc==limite):
			if(nexons>1 and exonanterior[6]=='-'):
				isize=len(introns)
				count=isize-1
				counti=0
				while(counti<isize):
					nintronw=str(count)
					intronwrite=introns[counti]+v+nintronw+'\n'
					intronarchive.writelines(intronwrite) #WRITE INTRON-
					count=int(count)
					count-=1
					counti+=1
			nexons=str(nexons)
			genedata=genedata+nexons+primerfive+primerthree
			genearchive.writelines(genedata) #write the last gene
			archivegff.close()
			genearchive.close()
			exonarchive.close()
			intronarchive.close()
			break