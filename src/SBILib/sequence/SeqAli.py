'''
IMPORTANT! Differientatie between sequence position (counts gaps and initial index) and alingment position

Iterables of SeqAli refer to the alignment

jbonet @ oliva's lab 2013
'''
import re
from collections import Counter

from SBILib.error    import SeqAliError as SAE
from SBILib.sequence import Sequence    as Sequence

class SeqAli(object):

	def __init__(self, sequences, sequenceInits, identities = None, positives = None, gaps = None):

		if not isinstance(sequences, list): 	 raise AttributeError('Sequences must be added in a list\n')
		if not isinstance(sequenceInits, list):	 raise AttributeError('Sequence inits must be added in a list\n')
		if len(sequences) != len(sequenceInits): raise AttributeError('One Sequence Init is required for each sequence\n')

		self._seq     = []
		for aliseq in sequences:
			self._seq.append(Sequence(sequence = aliseq))
		self._num_seq = len(self._seq)
		self._segment = []

		self._idx	  = self._set_index(sequenceInits)

		self._segmentation_ok = self._search_segments()

		if not self._segmentation_ok:
			raise SAE(code = 1)

		self._identities = identities
		self._positives  = positives
		self._gaps       = gaps
		self._alipatt    = None
		self._aliptmeth  = None

		self._aligned_aa = self._get_aligned_aa()

	'''ATTRIBUTES'''
	@property
	def sequences(self):			return self._seq
	@property
	def identities(self):			return int(self._identities)
	@property
	def identities_pec(self):		return 100 * float(self._identities)/self._aligned_aa
	@property
	def positives(self):			return int(self._positives)
	@property
	def positives_pec(self):		return 100 * float(self._positives)/self._aligned_aa
	@property
	def gaps(self):					return int(self._gaps)
	@property
	def alignment_pattern(self):	return self._alipatt
	@property
	def is_multiple(self):			return self._num_seq > 2
	@property
	def aligned_aminoacids(self):	return self._aligned_aa
	@property
	def number_of_sequences(self):	return self._num_seq
	@property
	def are_segments_ok(self):		return self._segmentation_ok

	'''INTERNAL GETTERS/SETTERS'''
	def _sequence_position_id(self, refseq, pos): #Only for ungapped (basically to use in _search_segments())
		return self._idx[refseq] + pos - 1

	def _aligned_sequence_length(self, refseq):
		return len(re.sub('-','',self.sequences[refseq].sequence))

	def _sequence_position_from_alignment_position(self, refseq, pos): 
		#Considering over alignment with everything
		#if the alignemnt position is gap, we move up to the next position
		data = self._seq[refseq][:pos]
		c    = 1
		if bool(re.match(r'\S+\-\w\Z',data)): cr = 1
		else:							cr = 0	
		while (data.endswith('-') and (data <= len(self._seq[refseq]))):
			data = self._seq[refseq][:pos + c]
			c += 1
		c = Counter(data)['-']
		if bool(c):
			return self._idx[refseq] + len(data) - Counter(data)['-'] - cr
		else:
			return self._idx[refseq] + pos - 1

	def _alignment_position_from_sequence_position(self, refseq, pos):
		if not self._seq[refseq].is_gapped:
			try: return pos - self._idx[refseq] + 1
			except ValueError as e: raise IndexError(e.message)
		else:
			c = 0
			pregaps = True
			for x in range(len(self._seq[refseq])):
				if pregaps and self._seq[refseq][x] == '-': continue
				pregaps = False
				if pos == self._idx[refseq] + c:
					return  x  + 1
				if x+1 < len(self._seq[refseq]) and self._seq[refseq][x+1] != '-':
					c += 1
			raise IndexError
 
	def _sequence_position_to_sequence_position(self, refseq, destseq, pos):
 		pos = self._alignment_position_from_sequence_position(refseq, pos)
 		return self._sequence_position_from_alignment_position(destseq, pos)

	def _first_segment_value(self, refseq):
		for x in range(len(self._segment[refseq])):
			if self._segment[refseq][x] != '-':
				return self._segment[refseq][x]

	def _last_segment_value(self, refseq):
		for x in reversed(list(range(len(self._segment[refseq])))):
			if self._segment[refseq][x] != '-':
				return self._segment[refseq][x]

	'''FUNCTIONS'''
	def overlap(self, qrefseq, mappedrefseq, otherAli):
		# measured over the SHORTEST sequence
		# ranges 0<->1
		length = min(self._aligned_sequence_length(qrefseq), otherAli._aligned_sequence_length(qrefseq))

		me = self.get_section_from_sequence_position(qrefseq,     
													 self._first_segment_value(qrefseq),     
													 self._last_segment_value(qrefseq))
		he = otherAli.get_section_from_sequence_position(qrefseq, 
														 otherAli._first_segment_value(qrefseq), 
														 otherAli._last_segment_value(qrefseq))

		if me._first_segment_value(mappedrefseq) > he._last_segment_value(mappedrefseq) or\
		   me._last_segment_value(mappedrefseq)  < he._first_segment_value(mappedrefseq):
			return 0

		if me._first_segment_value(mappedrefseq) == he._first_segment_value(mappedrefseq) and\
		   me._last_segment_value(mappedrefseq)  == he._last_segment_value(mappedrefseq):
			return 1

		upper_ini = max(me._first_segment_value(mappedrefseq), he._first_segment_value(mappedrefseq))
		lower_end = min(me._last_segment_value(mappedrefseq),  he._last_segment_value(mappedrefseq))
		cover = lower_end - upper_ini + 1

		return float(cover)/float(length)

	def get_sequence_section_from_sequence_position(self, refseq, ini, end, gapped = True):

		qseq = self.get_section_from_sequence_position(refseq, ini, end)._seq[refseq]

		if gapped: 	return qseq
		else:		
			qseq.do_ungap()
			return qseq

	def get_segment_section_from_sequence_position(self, refseq, ini, end):

		return self.get_section_from_sequence_position(refseq, ini, end)._segment[refseq]

	def increment_sequence_index(self, refseq, new_index):
		self._idx[refseq] += (new_index - 1)
		for x in range(len(self._segment[refseq])):
			self._segment[refseq][x] = self._segment[refseq][x] - 1 + self._idx[refseq]

	def get_respective_coordinate(self, refseq, destseq, pos):
		return self._sequence_position_to_sequence_position(refseq, destseq, pos)

	def get_coverage_of_sequence(self, refseq,  original_seqlength, ini = None, end = None):
		#if a range is given (sequence positions), it only calculates the coverage
		#according to that section (coverage of the section)
		#else, gives the coverage over the full seq length (original_seqlength)
		is_section = (ini is not None or end is not None)
		
		if ini is None: 
			ini = self._first_segment_value(refseq)
		if end is None: 
			end = self._last_segment_value(refseq)

		if ini is None and end is None: # this sequence is all gaps
			return 0

		section = self.get_section_from_sequence_position(refseq, ini, end)

		if is_section:
			original_seqlength = section._aligned_sequence_length(refseq)

		coverage = 0
		tk_seqs = section._tokenize() 
		for x in range(len(tk_seqs[refseq])):
			if bool(tk_seqs[refseq][x]) and sum([int(n[x]) for n in tk_seqs]) > 1:
				coverage += 1

		return float(coverage)/int(original_seqlength)

	def get_section_from_sequence_position(self, refseq, ini, end):
		align_ini = self._alignment_position_from_sequence_position(refseq, ini)
		align_end = self._alignment_position_from_sequence_position(refseq, end)
		if align_ini > len(self) or align_end > len(self): raise IndexError
		if align_ini <= 0 or align_end <= 0:               raise IndexError
		if align_ini == align_end: 
			if align_ini == len(self):
				align_end = None
			else:
				align_end += 1
		return self.__getitem__(slice(align_ini,align_end,1))

	def add_alignment_pattern(self, alipatt, method = 'blast'):
		self._alipatt   = alipatt
		self._aliptmeth = method
		if self._aliptmeth == 'blast' and not self.is_multiple:
			if self._positives is None and self._identities is None and self._gaps is None:
				self._positives  = 0
				self._identities = 0
				self._gaps       = 0
				for x in range(len(self._alipatt)):
					position = self._alipatt[x]
					if position == '+':
						self._positives  += 1
					elif position == ' ':
						if bool(re.search('-',"".join(self[x + 1]))):
							self._gaps += 1
					else:
						self._identities += 1
						self._positives  += 1
		else:
			return NotImplemented

	"""
		FORMULA DEFINITIONS FOR ROST TWILIGHT ZONE:
		Rost, B. (1999). Twilight zone of protein sequence alignments. Protein Engineering, 12(2), 85 94.
		MODIFIED FROM J.Planas-Iglesias IMPLEMENTATION

		INCLUDES 1 public and 2 private functions.
		Rost's Twilight Zone can only be applyed to alignments between TWO (2) SEQUENCES.
		This means is not allowed for MULTIPLE SEQUENCE ALIGNMENT
	"""
	def evaluate_Rost_twilight_zone(self, equation = "ID", parameter = 0):
		if self.is_multiple: raise SAE(code = 2, value = self._num_seq)

		if equation == "ID" :
			isHomolog = self.identities_pec >= self._get_Rost_ID_threshold(param=parameter)
		elif equation == "SIM": 
			isHomolog = self.positives_pec  >= self._get_Rost_SIM_threshold(param=parameter)
		return isHomolog

	def _get_Rost_ID_threshold(self, param=0):
		import decimal, math
		L = self._aligned_aa
		return decimal.Decimal(param)+(480*pow(L,decimal.Decimal('-0.32')*(1+pow(decimal.Decimal(repr(math.e)),decimal.Decimal(repr(float(-L)/1000))))))

	def _get_Rost_SIM_threshold(self, param=0):
		import decimal, math
		L = self._aligned_aa
		return decimal.Decimal(param)+(420*pow(L,decimal.Decimal('-0.335')*(1+pow(decimal.Decimal(repr(math.e)),decimal.Decimal(repr(float(-L)/2000))))))

	'''PRINTING FORMAT'''
	def formatPositions(self, human_readable = False):
		"""
		METHOD: formatPositions()
			 __STRING_GETTER__
			=>Returns:
				self.qpos
				self.hpos
			as a single string in which fragment definitions will be represented as follows:
				qpos[0]:hpos[0],qpos[1]:hpos[1];qpos[2]:hpos[2],qpos[3]:hpos[3]; [...]
				^       ^       ^       ^       ^       ^       ^       ^
				|       |       |       |       |       |       |       |Hit segment 2 end
				|       |       |       |       |       |       |Query segment 2 end
				|       |       |       |       |       |Hit segment 2 start
				|       |       |       |       |Query segment 2 start
				|       |       |       |Hit segment 1 end
				|       |       |Query segment 1 end
				|       |Hit segment 1 start
				|Query segment 1 start

			IN MULTIPLE ALIGNMENTS, a section which is not aligned while other are will display the same ini-start for a
			given segment.
		"""
		outdata = []
		if not human_readable:
			for x in range(len(self._segment[0])):
				text = []
				for i in range(len(self._segment)):
					text.append("{0}".format(self._segment[i][x]))
				outdata.append(":".join(text))
				if x % 2:
					if x != len(self._segment[0]) - 1:
						outdata.append(";")
				else:
					outdata.append(",")
			return "".join(outdata)
		else:
			for i in range(len(self._segment)):
				text = []
				for x in range(len(self._segment[i])):
					text.append("{0:>4}".format(self._segment[i][x]))
					if not x % 2:
						text.append(" :")
					elif x != len(self._segment[0]) - 1:
						text.append(" |")
				text = "".join(text)
				outdata.append(text)
			return "\n".join(outdata)

	'''FUNCTIONS CALLED ON LOAD'''
	def _set_index(self, seqini):
		return seqini

	def _tokenize(self):
		'''
		returns a binary representation of the alignment in an array
		'''
		token_seqs = []
		for seq in self._seq:
			token_seqs.append(seq.tokenize('binary'))
		return token_seqs

	def _get_aligned_aa(self):
		this_len = 0
		tk_seqs = self._tokenize()

		for x in range(len(self._seq[0])):
			profile_str = "".join([n[x] for n in tk_seqs])
			profile_vle = sum([int(i) for i in profile_str])
			if profile_vle > 1:
				this_len += 1
		return this_len

	def _check_segments(self):
		static_length = len(self._segment[0])
		if static_length %2 != 0: return False
		for i in range(len(self._segment)):
			if len(self._segment[i]) != static_length:
				return False
		return True

	def _search_segments(self):
		counter = []
		previous_profile = ''
		for i in range(self._num_seq):
			counter.append(1)
			previous_profile += '0'
			self._segment.append([])

		tk_seqs = self._tokenize()
		for x in range(len(self._seq[0])):
			profile_str = "".join([n[x] for n in tk_seqs])

			for i in range(self._num_seq):
				if x == 0:
					if bool(int(profile_str[i])): self._segment[i].append(self._sequence_position_id(i,counter[i]))
					else:						  self._segment[i].append('-')	
				else:
					if profile_str != previous_profile:
						if bool(int(previous_profile[i])): self._segment[i].append(self._sequence_position_id(i, (counter[i] - 1)))
						else:							   self._segment[i].append('-')
						if bool(int(profile_str[i])): 	   self._segment[i].append(self._sequence_position_id(i, counter[i]))
						else:						  	   self._segment[i].append('-')
					if x == (len(self._seq[0]) - 1):
						if bool(int(profile_str[i])): self._segment[i].append(self._sequence_position_id(i, counter[i]))
						else:						  self._segment[i].append('-')
				if bool(int(profile_str[i])): counter[i] += 1

			previous_profile = profile_str

		x = 0
		while (x < len(self._segment[0])):
			counter = 0
			for y in range(self._num_seq):
				if self._segment[y][x] != '-': counter += 1
			if counter == 1:
				for y in range(self._num_seq):
					self._segment[y].pop(x)
					self._segment[y].pop(x)
				x = x - 2
			x += 2
		if len(self) == 1:
			for y in range(self._num_seq):
				self._segment[y].append(self._segment[y][0])
		return self._check_segments()

	'''OVEWRITTE OBJECT FUNCTIONS'''
	def __len__(self):
		return len(self._seq[0])

	def __getitem__(self, key):
		# We correct the count as one would try to ask for the first position as 1
		try: 		
			int(key)
			if int(key) > len(self): raise IndexError
			return [n[int(key) - 1] for n in self._seq]
		except:
			if not isinstance(key, slice): 	
				raise TypeError()
			else:
				index = key.indices(len(self) + 1)
				if key.stop is not None:
					index = (index[0],index[1]+1,index[2]) #Remember! we go according to align position, not array!!
				else:
					index = (index[0],index[1],index[2])
				sequences = ['' for x in range(self.number_of_sequences)]
				seqInits  = ['' for x in range(self.number_of_sequences)]
				for i in range(*index):
					data = self[i]
					for x in range(len(data)):
						sequences[x] += data[x]
				for x in range(self.number_of_sequences):
					seqInits[x] = self._sequence_position_from_alignment_position(x,index[0])
				newali = self.__class__(sequences, seqInits) 
				if not self._alipatt is None:
					newali.add_alignment_pattern(self._alipatt[index[0]-1:index[1]-1:index[2]],self._aliptmeth)
				return newali
				
	def __iter__(self):
		for x in range(1, len(self) + 1):
			yield self[x]

	def __repr__(self):
		data = []
		for x in range(self._num_seq):
			data.append('{0}\t{1}\t{2}'.format(self._first_segment_value(x), self._seq[x].sequence, self._last_segment_value(x)))
		if not self._alipatt is None:
			data.append('\t{0}'.format(self._alipatt))
		if not self._identities is None and not self._positives is None and not self._gaps is None:
			data.append('I: {0._identities:<4} P: {0._positives:<4} G: {0._gaps:<4}'.format(self))
		data.append(self.formatPositions(True))
		return "\n".join(data)
