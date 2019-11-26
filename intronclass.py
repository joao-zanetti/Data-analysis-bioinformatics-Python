##Developed by João Zanetti
##https://github.com/joao-zanetti

import collections
import matplotlib.pyplot as plt
import numpy as np
import operator
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import re
import math
 
#Object for manipulation of intron data
class Intron(object): 
	def __init__(self,sequencereg,gene,size,start,end,sense,nintron,seq):
		self.sequencereg=sequencereg
		self.gene=gene
		self.size=size
		self.start=start
		self.end=end
		self.sense=sense
		self.nintron=nintron
		self.seq=seq


	#############
	#method that create and return one vector introns
	@staticmethod
	def introns_create(): 
		archive=open('intron.csv','r')
		introns=[]
		for line in archive:
			if(line!='\n'):
				intronread=line.split(',')
				print(intronread[7])
				intronread[2]=int(intronread[2])
				intronread[3]=int(intronread[3])
				intronread[4]=int(intronread[4])
				introns.append(Intron(*intronread))
		archive.close()
		return introns


	#############
	#Method that create and return one vector with all minimun introns
	@staticmethod
	def introns_createfasta(): 
		archive=open('intronfastaNotMinimum.csv','r')
		introns=[]
		for line in archive:
			if(line!='\n'):
				intronread=line.split(',')
				intronread[2]=int(intronread[2])
				intronread[3]=int(intronread[3])
				intronread[4]=int(intronread[4])
				introns.append(Intron(*intronread))
		archive.close()
		return introns


	#############
	#Method that creates one dictionary wich all sizes of introns with your respective apparition
	#AND plot the line graphic
	@staticmethod
	def introns_counter(): 
		introns= Intron.introns_create()
		dic={}
		for intro in introns:
			flag=0
			for elemdic in dic:
				if(elemdic==intro.size):
					flag=1
			if(flag==1):
				dic[intro.size]=dic[intro.size]+1
			else:
				dic[intro.size]=1

		dic=collections.Counter(dic)
		#print dic.most_common(10)# number of elements that most appeared
		nintronmodal=max(dic.iterkeys(), key=lambda k: dic[k])
		minintron= dic[nintronmodal]/4

		maxsize=max(dic.iteritems(), key=operator.itemgetter(0))[0]
		minsize=min(dic.iteritems(), key=operator.itemgetter(0))[0]
		maxnumb=max(dic.iteritems(), key=operator.itemgetter(1))[1]
		minnumb=min(dic.iteritems(), key=operator.itemgetter(1))[1]
		maxsize+=100
		maxnumb+=5

		lists = sorted(dic.items())#sorted by key, return a list of tuples
		#print lists
		y,x= zip(*lists) # unpack a list of pairs into two tuples

		plt.ylim([0,maxsize]) #NUMBER x SIZE
		plt.xlim([minnumb,100]) #NUMBER x SIZE
		plt.title('Introns from schmidtea mediterranea')
		plt.xlabel('Size(bp)')
		plt.ylabel('Numbers')
		plt.grid(True)
		plt.fill_between(y, 0, x, alpha=0.2)
		plt.plot(y,x,linewidth=1.0)
		#plt.show()

		d= {k:v for (k,v) in dic.items() if v > minintron}

		return d
		
	
	#############
	#Plot the histogram of introns
	#Wich all sizes of introns with your respective apparition
	@staticmethod 
	def introns_hist(): 
		introns= Intron.introns_create()
		tamanhos=[]
		
		for intro in introns:
			tamanhos.append(intro.size)

		n, bins, patches = plt.hist(tamanhos,25000, alpha=0.7, rwidth=0.85)
		plt.title('Introns from schmidtea mediterranea')
		plt.xlabel('Size(bp)')
		plt.ylabel('Numbers')

		plt.xlim((0, 100)) 
		plt.grid(True)
		plt.show()

		return


	#############
	#EXTRACT INTRONS EDGE Minimum 
	@staticmethod
	def introns_edgeM():
		intronearchive=open('intronEdgeMinimum.csv','w')
		introns= Intron.introns_createfasta()
		for intro in introns:
			intro.seq=re.sub('\n','',intro.seq)
			if(intro.sense=='+'):
				intronfwrite1=intro.seq[0:15]
				intronfwrite2=intro.seq[intro.size-6:]
				intronfwrite=intronfwrite1+','+intronfwrite2+'\n'
				intronearchive.writelines(intronfwrite)
			else:
				intronfwrite1=intro.seq[0:15]
				intronfwrite2=intro.seq[intro.size-6:]
				intronfwrite=intronfwrite1+','+intronfwrite2+'\n'
				intronearchive.writelines(intronfwrite)
		return


	#############
	#EXTRACT INTRONS EDGE not minimum
	@staticmethod
	def introns_edgeNM():
		intronearchive=open('intronEdgeNotMinimum.csv','w')
		introns= Intron.introns_createfasta()
		for intro in introns:
			intro.seq=re.sub('\n','',intro.seq)
			if(intro.sense=='+'):
				intronfwrite1=intro.seq[0:15]
				intronfwrite2=intro.seq[intro.size-6:]
				intronfwrite=intronfwrite1+','+intronfwrite2+'\n'
				intronearchive.writelines(intronfwrite)
			else:
				intronfwrite1=intro.seq[0:15]
				intronfwrite2=intro.seq[intro.size-6:]
				intronfwrite=intronfwrite1+','+intronfwrite2+'\n'
				intronearchive.writelines(intronfwrite)
		return


	#############
	#EXTRACT INTRONS EDGE in the intron
	@staticmethod
	def introns_INedge():
		intronearchive=open('intronInEdge.csv','w')
		introns= Intron.introns_createfasta()
		for intro in introns:
			intro.seq=re.sub('\n','',intro.seq)
			if(intro.sense=='+'):
				intronfwrite1=intro.seq[5:25]
				intronfwrite2=intro.seq[intro.size-16:intro.size+4]
				intronfwrite=intronfwrite1+','+intronfwrite2+'\n'
				intronearchive.writelines(intronfwrite)
			else:
				intronfwrite1=intro.seq[5:25]
				intronfwrite2=intro.seq[intro.size-16:intro.size+4]
				intronfwrite=intronfwrite1+','+intronfwrite2+'\n'
				intronearchive.writelines(intronfwrite)
		return


	#############
	#Create one archive with the intron fasta sequence
	@staticmethod
	def introns_fastaAll():
		v=','
		intronfarchive=open('intronfastaAll.csv','w')
		record_dict=SeqIO.index("schmidtea_mediterranea.PRJNA12585.WBPS11.genomic.fasta", "fasta")
		introns= Intron.introns_create()
		for intro in introns:
			if(intro.sense=='+'):
				intronrecord=record_dict[intro.sequencereg].seq[intro.start-6:intro.end+5]
				introsize=str(intro.size)
				introstart=str(intro.start)
				introend=str(intro.end)
				intronfwrite=intro.sequencereg+v+intro.gene+v+introsize+v+introstart+v+introend+v+intro.sense+v+intro.nintron+v+intronrecord+'\n'
				intronfarchive.writelines(intronfwrite)
			else:
				intronrecord=record_dict[intro.sequencereg].seq[intro.start-6:intro.end+5]
				intronrecord=intronrecord.reverse_complement()
				introsize=str(intro.size)
				introstart=str(intro.start)
				introend=str(intro.end)
				intronfwrite=intro.sequencereg+v+intro.gene+v+introsize+v+introstart+v+introend+v+intro.sense+v+intro.nintron+v+intronrecord+'\n'
				intronfarchive.writelines(intronfwrite)
		return


	#############
	#Create one archive with information of sequence of the minimum intron 
	@staticmethod
	def introns_fastaM(dicminmodal):
		v=','
		intronfarchive=open('intronfastaMinimum.csv','w')
		record_dict=SeqIO.index("schmidtea_mediterranea.PRJNA12585.WBPS11.genomic.fasta", "fasta")
		introns= Intron.introns_create()
		for intro in introns:
			if(intro.size in dicminmodal):
				if(intro.sense=='+'):
					intronrecord=record_dict[intro.sequencereg].seq[intro.start-6:intro.end+5]
					introsize=str(intro.size)
					introstart=str(intro.start)
					introend=str(intro.end)
					intronfwrite=intro.sequencereg+v+intro.gene+v+introsize+v+introstart+v+introend+v+intro.sense+v+intro.nintron+v+intronrecord+'\n'
					intronfarchive.writelines(intronfwrite)
				else:
					intronrecord=record_dict[intro.sequencereg].seq[intro.start-6:intro.end+5]
					intronrecord=intronrecord.reverse_complement()
					introsize=str(intro.size)
					introstart=str(intro.start)
					introend=str(intro.end)
					intronfwrite=intro.sequencereg+v+intro.gene+v+introsize+v+introstart+v+introend+v+intro.sense+v+intro.nintron+v+intronrecord+'\n'
					intronfarchive.writelines(intronfwrite)
		return


	#############
	#Create one archive with the information of sequence of not minimum intron 
	@staticmethod
	def introns_fastaNM(dicminmodal):
		v=','
		intronfarchive=open('intronfastaNotMinimum.csv','w')
		record_dict=SeqIO.index("schmidtea_mediterranea.PRJNA12585.WBPS11.genomic.fasta", "fasta")
		introns= Intron.introns_create()
		for intro in introns:
			if(intro.size in dicminmodal):
				print('not in')
			else:
				if(intro.sense=='+'):
					intronrecord=record_dict[intro.sequencereg].seq[intro.start-6:intro.end+5]
					introsize=str(intro.size)
					introstart=str(intro.start)
					introend=str(intro.end)
					intronfwrite=intro.sequencereg+v+intro.gene+v+introsize+v+introstart+v+introend+v+intro.sense+v+intro.nintron+v+intronrecord+'\n'
					intronfarchive.writelines(intronfwrite)
				else:
					intronrecord=record_dict[intro.sequencereg].seq[intro.start-6:intro.end+5]
					intronrecord=intronrecord.reverse_complement()
					introsize=str(intro.size)
					introstart=str(intro.start)
					introend=str(intro.end)
					intronfwrite=intro.sequencereg+v+intro.gene+v+introsize+v+introstart+v+introend+v+intro.sense+v+intro.nintron+v+intronrecord+'\n'
					intronfarchive.writelines(intronfwrite)
				
		return


	#############
	#create one archive that shows the proportion base of each position of intron not minimum
	@staticmethod
	def proportionSbases():
		intronearchive=open('intronEdgeNotMinimum.csv','r')
		proportionEdgeBases=open('proportionEdgeBasesNotMinimum.csv','w')
		ddic=[]
		ddicp=[]
		introns=[]
		size=30
		for introneread in intronearchive:
			intronsplit=re.sub('\n|,','',introneread)
			introns.append(intronsplit)
		nbases=len(introns)
		nbases=float(nbases)
		dicmain= {'A':0,'T':0,'C':0,'G':0,'N':0}
		for i in range(0,size):
			dic= {'A':0,'T':0,'C':0,'G':0,'N':0}
			for introneread in introns:
				dic[introneread[i]]+=1
				dicmain[introneread[i]]+=1
			ddic.append(dic.copy())

		dic= {'A':0.0,'T':0.0,'C':0.0,'G':0.0,'N':0.0}
		for i in range(0,size):
			dic['A']=float(ddic[i]['A'])/nbases
			dic['T']=float(ddic[i]['T'])/nbases
			dic['C']=float(ddic[i]['C'])/nbases
			dic['G']=float(ddic[i]['G'])/nbases
			ddicp.append(dic.copy())

		#for i in range(0,30):
			#print ddic[i]

		for i in range(0,size):
			j=str(i)
			prop= j+','+str(ddicp[i])+'\n'
			proportionEdgeBases.writelines(prop)
			print('POSIIITIOONN :'+j) 
			print(ddicp[i]) 
			print('\n\n\n') 
		#print dicmain

		archiveqofi= open('qofiEdgeBasesNotMinimum','w')
		Intron.qofinformation(ddicp,size,archiveqofi)


	#############
	#create one archive that shows the proportion base of each position of edge in intron
	@staticmethod
	def proportionINbases():
		intronearchive=open('intronInEdge.csv','r')
		proportionINBases=open('proportionInBases.csv','w')
		ddic=[]
		ddicp=[]
		introns=[]
		size=40
		for introneread in intronearchive:
			intronsplit=re.sub('\n|,','',introneread)
			introns.append(intronsplit)
		nbases=len(introns)
		nbases=float(nbases)
		dicmain= {'A':0,'T':0,'C':0,'G':0,'N':0}
		for i in range(0,size):
			dic= {'A':0,'T':0,'C':0,'G':0,'N':0}
			for introneread in introns:
				dic[introneread[i]]+=1
				dicmain[introneread[i]]+=1
			ddic.append(dic.copy())

		dic= {'A':0.0,'T':0.0,'C':0.0,'G':0.0,'N':0.0}
		for i in range(0,size):
			dic['A']=float(ddic[i]['A'])/nbases
			dic['T']=float(ddic[i]['T'])/nbases
			dic['C']=float(ddic[i]['C'])/nbases
			dic['G']=float(ddic[i]['G'])/nbases
			ddicp.append(dic.copy())

		#for i in range(0,30):
			#print ddic[i]
		
		for i in range(0,size):
			j=str(i)
			prop= j+','+str(ddicp[i])+'\n'
			proportionINBases.writelines(prop)
			print('POSIIITIOONN :'+j)
			print(ddicp[i]['A'])
			print(ddicp[i]['T'])
			print(ddicp[i]['C'])
			print(ddicp[i]['G'])
			print('\n\n\n')
		#print dicmain

		archiveqofi= open('qofiInBases','w')
		Intron.qofinformation(ddicp,size,archiveqofi)
	
	
	#############
	#Calcule the quantity of information of intron
	#From the archive that was fill by proportion of bases
	@staticmethod##fu
	def qofinformation(ddicp,size,archiveqofi):
		qi= []
		for i in range(0,size):
			#b=str(i)
			#print 'POSITIIIOOONNN'+b
			hi=0
			hi=ddicp[i]['A']*(math.log(ddicp[i]['A'],2))
			hi+=ddicp[i]['C']*(math.log(ddicp[i]['C'],2))
			hi+=ddicp[i]['T']*(math.log(ddicp[i]['T'],2))
			hi+=ddicp[i]['G']*(math.log(ddicp[i]['G'],2))
			hi=float(hi)
			#print hi
			en=(1.0/0.6931472)*((4.0-1.0)/(2.0*22177.0))
			qix=2.0+hi-en
			qi.append(qix)
			#print '\n'
		print(qi)
		j=0
		j= str(j)
		for elem in qi:
			aux= str(elem)
			qofi= j+','+aux+'\n'
			archiveqofi.writelines(qofi)
			j=int(j)
			j+=1
			j=str(j)

	
	#############
	#Calcule the quantity of information of in intron
	@staticmethod
	def edgein_tofasta():
		edgein= open('intronInEdge.csv','r')
		edgeinfasta= open('intronInEdge.fasta','w')
		for i in edgein:
			edgesread= re.sub(',','',i)
			linefasta='>\n'+edgesread
			edgeinfasta.writelines(linefasta)

	
	#############
	#Create one archive with the fasta intron not minimum fasta sequence
	@staticmethod
	def edge_tofasta():
		edgein= open('intronEdgeMinimum.csv','r')
		edgeinfasta= open('intronEdgeMinimum.fasta','w')
		for i in edgein:
			edgesread= re.sub(',','',i)
			linefasta='>\n'+edgesread
			edgeinfasta.writelines(linefasta)


	#############
	# method that plot the graphic of quantity of information for each intron and your edge
	@staticmethod
	def plotQoiEdge(): # AND plot the graphic
		archivem= open('qofiEdgeBasesMinimum','r')
		archivenm= open('qofiEdgeBasesNotMinimum','r')
		dicm={}
		dicnm={}
		for line in archivem:
			linesplit=line.split(',')
			#print linesplit[1]
			linesplit[1]=re.sub('\n','',linesplit[1])
			linesplit[0]=int(linesplit[0])+1
			linesplit[1]=float(linesplit[1])
			dicm[linesplit[0]]=linesplit[1]
		print(dicm)

		for line in archivenm:
			linesplit=line.split(',')
			#print linesplit[1]
			linesplit[1]=re.sub('\n','',linesplit[1])
			linesplit[0]=int(linesplit[0])+1
			linesplit[1]=float(linesplit[1])
			dicnm[linesplit[0]]=linesplit[1]
		print(dicnm)


		dicm=collections.Counter(dicm)
		dicnm=collections.Counter(dicnm)

		listsm = sorted(dicm.items())#sorted by key, return a list of tuples
		#print lists
		y,x= zip(*listsm)

		listsnm = sorted(dicnm.items())#sorted by key, return a list of tuples
		#print lists
		y2,x2= zip(*listsnm)

		plt.ylim([0,2]) #NUMBER x SIZE
		plt.xlim([1,15]) #NUMBER x SIZE
		plt.title('Schmidtea mediterranea 5')
		plt.xlabel('Position')
		plt.ylabel('Information Entropy(bits)')
		plt.grid(True)
		plt.xticks(np.arange(1,16,step=1))
		plt.plot(y,x,linewidth=2.0,label='Minimal')
		plt.plot(y2,x2,linewidth=2.0,label='Not Minimal')
		plt.legend()

		#plt.annotate('size=47', xy=(47,2050), xytext=(3, 1.5),arrowprops=dict(facecolor='black', shrink=0.01),)
		#plt.plot(y,x,color='blue',lw=3)
		plt.show()


		plt.ylim([0,2]) #NUMBER x SIZE
		plt.xlim([16,30]) #NUMBER x SIZE
		#plt.xlim([minnumb,maxnumb])#NUMBER x SIZE
		plt.title('Schmidtea mediterranea 3')
		plt.xlabel('Position')
		plt.ylabel('Information Entropy(bits)')
		plt.grid(True)
		plt.xticks(np.arange(16,31,step=1))
		#plt.savefig(images_path + file_name.replace('.txt','-commiter-contribution.eps'), format='eps')
		#plt.fill_between(y, 0, x, alpha=0.2)
		plt.plot(y,x,linewidth=2.0,label='Minimal')
		plt.plot(y2,x2,linewidth=2.0,label='Not Minimal')
		plt.legend()

		#plt.annotate('size=47', xy=(47,2050), xytext=(3, 1.5),arrowprops=dict(facecolor='black', shrink=0.01),)
		#plt.plot(y,x,color='blue',lw=3)
		plt.show()


	#############
	#Method that creates four dictionarys
	#Wich all proportion of each position in edge of introns
	@staticmethod
	def plotPropInEdge(): # AND plot the graphic
		archive= open('proportionInBases.csv','r')
		dicA={}
		dicT={}
		dicC={}
		dicG={}
		for line in archive:
			linesplit=re.sub('\n|{|}','',line)
			linesplit=re.sub(':',',',linesplit)
			linesplit=linesplit.split(',')
			linesplit[0]=int(linesplit[0])+1
			dicA[linesplit[0]]=float(linesplit[2])
			dicT[linesplit[0]]=float(linesplit[6])
			dicC[linesplit[0]]=float(linesplit[4])
			dicG[linesplit[0]]=float(linesplit[8])

		
		dicA=collections.Counter(dicA)
		dicT=collections.Counter(dicT)
		dicC=collections.Counter(dicC)
		dicG=collections.Counter(dicG)

		listsa = sorted(dicA.items())#sorted by key, return a list of tuples
		ya,xa= zip(*listsa)
		listst = sorted(dicT.items())#sorted by key, return a list of tuples
		yt,xt= zip(*listst)
		listsc = sorted(dicC.items())#sorted by key, return a list of tuples
		yc,xc= zip(*listsc)
		listsg = sorted(dicG.items())#sorted by key, return a list of tuples
		yg,xg= zip(*listsg)

		plt.ylim([0,1]) #NUMBER x SIZE
		plt.xlim([1,20]) #NUMBER x SIZE
		plt.title('Schmidtea mediterranea 5')
		plt.xlabel('Position')
		plt.ylabel('Proportion of Bases')
		plt.grid(True)
		plt.xticks(np.arange(1,21,step=1))
		plt.plot(ya,xa,linewidth=2.0,label='A')
		plt.plot(yt,xt,linewidth=2.0,label='T')
		plt.plot(yc,xc,linewidth=2.0,label='C')
		plt.plot(yg,xg,linewidth=2.0,label='G')
		plt.legend()

		#plt.annotate('size=47', xy=(47,2050), xytext=(3, 1.5),arrowprops=dict(facecolor='black', shrink=0.01),)
		#plt.plot(y,x,color='blue',lw=3)
		plt.show()

		plt.ylim([0,1]) #NUMBER x SIZE
		plt.xlim([21,40]) #NUMBER x SIZE
		plt.title('Schmidtea mediterranea 3')
		plt.xlabel('Position')
		plt.ylabel('Proportion of Bases')
		plt.grid(True)
		plt.xticks(np.arange(21,41,step=1))
		#plt.savefig(images_path + file_name.replace('.txt','-commiter-contribution.eps'), format='eps')
		#plt.fill_between(y, 0, x, alpha=0.2)
		plt.plot(ya,xa,linewidth=2.0,label='A')
		plt.plot(yt,xt,linewidth=2.0,label='T')
		plt.plot(yc,xc,linewidth=2.0,label='C')
		plt.plot(yg,xg,linewidth=2.0,label='G')
		plt.legend()

		#plt.annotate('size=47', xy=(47,2050), xytext=(3, 1.5),arrowprops=dict(facecolor='black', shrink=0.01),)
		#plt.plot(y,x,color='blue',lw=3)
		plt.show()