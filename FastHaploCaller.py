#!/apps/python/2.7.11/gcc-5.3.0/bin/python
from timeit import default_timer as timer
import sys, subprocess 
import time

from collections import Counter
from argparse import ArgumentParser
from collections import defaultdict
from itertools import combinations
import cProfile

from fastBamRead import Read

def parseArgs():
	parser = ArgumentParser(description='Parse a bam file to make a depthDict from the reads for a specific contig.')
	parser.add_argument("-b",help="Input Bam File", nargs = '*', action="store", dest="bamFile", required=True)
	parser.add_argument("-v",help="Input VCF File", action="store", dest="vcfFile", required=True)
	parser.add_argument("-c",help="Name of the Contig of Interest",action="store", dest="contigOfInterest", required=True)
	args = parser.parse_args()
	return args

def getReads(bamFile, contigOfInterest, positions):
	bashCommand = "samtools view " + bamFile + " " + contigOfInterest
	try:
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	except:
		print "load samtools module"
		sys.exit()
	output, error = process.communicate()
	
	lines = output.strip().split('\n')
	file = open('read.sam', 'w')
	file.write(output)
	reads = []
	for line in lines:
		ll = line.split('\t')
		if ll[5] != '*' and ll[4] != '0': # no empty CIGAR string and no mapq 0
			newRead = Read(ll)
			reads.append(newRead)
	return reads

def getPositionList(vcf, startPos, stopPos):
	# 21    48112747    .    T    TCTC    103.69    PASS    AC=1;AN=2    GT:AD:DP:GQ:PL    0/1:18,5:23:99:153,0,963
	# 21    48112773    .    T    C    223.48    PASS    AC=1;AN=2    GT:AD:DP:GQ:PL    0/1:16,3:19:78:78,0,711
	posDict = {}
	if stopPos == 0:
		stopPos = float('inf')
	with open(vcf,'r') as f:
	        for line in f:
	                if line[0] != '#':
	                        ll = line.strip().split()
	                        pos = int(ll[1])
				if '/' in ll[9]:
	                        	genotype = ll[9].split(':')[0].split('/')
				else:
	                        	genotype = ll[9].split(':')[0].split('|')
				if genotype[0]  == genotype[1]: # no Homozygous
					continue
				if startPos <= pos <= stopPos:
					possible = [ll[3]] + ll[4].split(',')
					first = possible[int(genotype[0])]
					second = possible[int(genotype[1])]

	                		posDict[pos -1] = (first, second)
				elif pos > stopPos:
					break

					#if abs(len(first) - len(second)) > 2:
	                		#	posDict[pos -1] = (first, second)
					#elif len(first) == 1 and len(second) == 1: #SNP
	                		#	posDict[pos -1] = (first, second)
					#elif len(first) > 1 and len(second) == 1: # deletion
	                		#	posDict[pos -1] = (first, second)
					#elif len(first) == 1 and len(second) > 1: # insertion
	                		#	posDict[pos -1] = (first, second)
	#print posDict
	return posDict

#def isInformative(read, positions):
#        posList = read.rangesDict.keys()
#        count = 0
#        for i in posList:
#		if count >= 2:
#			return True
#		if i in positions:
#			count += 1
#        return False

def createHaploStrings(reads, keys, posDict):
	haploDict = defaultdict(dict) 
	subsetReads = []
	for read in reads:
		numVars = 0
		for position in keys:
			if position in read.rangesDict:
				tup = posDict.get(position)
				first = tup[0]
				second = tup[1]
				
				readStartInd = read.rangesDict.get(position)
				readVal1 = read._seq[readStartInd:readStartInd + len(first)]
				readVal2 = read._seq[readStartInd:readStartInd + len(second)]
				len1 = len(first)
				len2 = len(second)

				if len1 == 1 and len2 == 1: #SNP
					if readVal1 == first:
						haploDict[read.qname][position] = '0'
						numVars += 1
					elif readVal2 == second:
						haploDict[read.qname][position] = '1'
						numVars += 1
				else:
					if len(readVal1) == len1 and len(readVal2) == len2:
						if len1 > len2:
							if readVal1 == first:
								haploDict[read.qname][position] = '0'
								numVars += 1
							elif readVal2 == second:
								haploDict[read.qname][position] = '1'
								numVars += 1
						else:
							if readVal2 == second:
								haploDict[read.qname][position] = '1'
								numVars += 1
							elif readVal1 == first:
								haploDict[read.qname][position] = '0'
								numVars += 1
		read._seq = None
	for r in reads:
		#print r.qname
		#if r.qname == "HSQ1004:134:C0D8DACXX:2:2101:2267:43445":
			#print 'yes'
		if len(haploDict[r.qname]) >= 2:
			subsetReads.append(r)
	return haploDict, subsetReads

def pickBest(myDict):
	combos = list(combinations(myDict, 2))
	max = 0
	finish = "NONE"
	for tup in combos:
		geno1 = tup[0]
		geno2 = tup[1]
		sumGeno = myDict[geno1] + myDict[geno2]
		if int(geno1) ^ int(geno2) == int('1' * len(geno1)):
			if max < sumGeno:
				max = sumGeno
				finish = tup
	return finish

def mergeWithPrev(h1, h2, tuple):
	if tuple == "NONE":
		return '',''
	if len(h1) == 0:
		return tuple[0], tuple[1]
	
	one = tuple[0]
	two = tuple[1]
	if h1[-1] == one[0] and h2[-1] == two[0]:
		h1 += one[1:]
		h2 += two[1:]
		return h1,h2
	elif h2[-1] == one[0] and h1[-1] == two[0]:
		h2 += one[1:]
		h1 += two[1:]
		return h1,h2
	else:
		return '',''
		
		

if __name__ == "__main__":		

	start = timer()

	args = parseArgs()
#	pr = cProfile.Profile()
#	pr.enable()
	
	bamFiles = args.bamFile
	contigOfInterest = args.contigOfInterest
	vcfFile = args.vcfFile
	if ':' in contigOfInterest:
		positions = contigOfInterest.split(':')[1].split('-')
		startPosition = int(positions[0]) - 1
		stopPosition = int(positions[1])
	else: 
		startPosition = 0
		stopPosition = 0
	

	# get list of variant positions within contigOfInterest
	posDict = getPositionList(vcfFile, startPosition + 1, stopPosition)
	positions = sorted(posDict.keys())
	print "read vcf file: "
	end = timer()
	print(end - start)
		
	# samtools grabs all reads in bamfiles within contigOfInterest
	reads = []
	for file in bamFiles:
		reads += getReads(file, contigOfInterest, positions)
	print "got reads    : "
	end = timer()
	print(end - start)

	# print summary of input reads
	print "contig: {}\n\tlength={}\n\t{} reads\n\t{} variants".format(contigOfInterest, stopPosition - startPosition, len(reads), len(posDict))

	haploDict, subsetReads = createHaploStrings(reads, positions, posDict)
	print "haploDicts   : "
	end = timer()
	print(end - start)

	print "sort position: "
	end = timer()
	print(end - start)
#	#{ position : ('A','AT'), ... }
	# zero-based position

	haplotype1 = ''
	haplotype2 = ''
	start = positions[0]
	del reads
	subsetReads.sort(key=lambda x: x.pos)
	print "\tsubset of reads =", len(subsetReads)
	print "sort reads   : "
	end = timer()
	print(end - start)



	count = 0
	for num in range(len(positions) - 1):
		occurenceList = []
		pos1 = positions[num]
		pos2 = positions[num + 1]
		for read in subsetReads:
			if read.first > pos2 + 117: # Read is beyond the position
				break
			if read.last < pos1 + 1: # Read is too far before the position
				continue
			count += 1
			readDict = haploDict[read.qname]
			if pos1 in readDict and pos2 in readDict:
				pos1Val = readDict[pos1]
				pos2Val = readDict[pos2]
				occurenceList.append(pos1Val + pos2Val)
				
		countDict = Counter(occurenceList)
		haploTup = pickBest(countDict)
		if haploTup == "NONE":
			stop = pos1
			print start,'-',stop
			print haplotype1
			print haplotype2
			start = pos2
		haplotype1, haplotype2 = mergeWithPrev(haplotype1,haplotype2, haploTup)


	stop = positions[num]
	print start,'-',stop
	print haplotype1
	print haplotype2
		
	print count
	print "haplotypes   : "
	end = timer()
	print(end - start)
#	pr.disable()
#	pr.print_stats(sort = 'time')

	
