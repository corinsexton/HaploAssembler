#!/apps/python/2.7.11/gcc-5.3.0/bin/python

class Read(object):

	"""A Read entry from a bam file accessed through `samtools view`

	Attributes: (correspond to the samtools headers)
		qname: (str)
		flag: (int)
		rname: (str)
		pos: (int)
		mapq: (int)
		cigar: (str)
		rnext: (str)
		pnext: (int)
		tlen: (int)
		_seq: (str)
		qual: (str)

		QAranges: (list<int>) query indexes (based on refSeq and CIGAR) where sequences aligned
		RAranges: (list<int>) reference indexes (based on refSeq and CIGAR) where sequences aligned
		rangesDict: (dict) corresponding matching indices { refIndex : queryIndex }

		rev: (bool) determines whether the 0x10 flag is present = reversed seq matched
		queryAlignment: (str) this read's string with spaces included for insertions and deletions
		refSeqAlignment: (str) refSeq string that corresponds to this read
		fullSeq: (str) full size of the refseq
		refSeq = converted to a global of the ref seq in question
	"""
	
	def __init__(self, ll):
	
		self.qname = ll[0]
		self.flag = int(ll[1])
		self.rname = ll[2]
		self.pos = int(ll[3])
		self.mapq = int(ll[4])
		self.cigar = ll[5]
		self.rnext = ll[6]
		self.pnext = int(ll[7])
		self.tlen = int(ll[8])
		self._seq = ll[9]
		self._seq2 = None
		self.qual = ll[10]

		self.length = len(self._seq)
		self.rangesDict1 = {} #{ refIndex : queryIndex }
		self.rangesDict2 = {} #{ refIndex : queryIndex }

		# populate RA and QA ranges
		self.getCigarAlignment(self.pos, self.cigar, self.rangesDict1)

	@classmethod
	def ReadPair(self,ll, ll1):
		read = Read(ll)
		read.getCigarAlignment(int(ll1[3]), ll1[5], read.rangesDict2)
		read._seq2 = ll[9]
		return read

	
	def findIfRev(self):
		if self.flag & 16:
			return True
		else:
			return False

	def toString(self):
		string = ''
		string += self.qname
		string += '\n'
		string +=  "\t{}\n\t{}\n\n".format(
				self.pos, self.cigar)
		return string
	
	def getCigarAlignment(self, position, cigar, rangesDict):
		cigarList = list(cigar)
		typeList = []
		refStart = position - 1 #samtools = 1-based position
		queStart = 0

		for c in cigarList: #convert to ints and strings in CIGAR
			try:
				typeList.append(int(c))
			except ValueError:
				typeList.append(c)

		count = 0
		lent = len(typeList)
		while count < lent:
			c = typeList[count]
			if c != '*':
				num =''
				while isinstance(c, int):
					num += str(c)
					count += 1
					c = typeList[count]
				num = int(num)
				if c in ["M", "=", "X"]:
					for i in range(num):
						rangesDict[refStart + i] = queStart + i
						
					queStart = queStart + num
					refStart = refStart + num
				elif c in ["I","S"]:
					queStart = queStart + num
				elif c in ["D", "N"]:
					refStart = refStart + num
			count += 1

