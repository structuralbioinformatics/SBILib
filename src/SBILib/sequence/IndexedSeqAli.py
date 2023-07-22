'''
jbonet @ oliva's lab 2013
'''
import re

from SBILib.error    import SeqAliError as SAE
from SBILib.sequence import SeqAli      as SeqAli
from SBILib.beans.IndexedNum    import IndexedNum

class IndexedSeqAli(SeqAli):

	def __init__(self, sequences, sequenceInits, identities = None, positives = None, gaps = None):

		super(IndexedSeqAli, self).__init__(sequences, sequenceInits, identities, positives, gaps)

		self._has_complex_index = self._check_index_kind(sequenceInits)

	'''FUNCTIONS CALLED ON LOAD'''
	def _check_index_kind(self, index):
		for idx in index:
			if isinstance(idx, list): return True
		return False

	'''INTERNAL GETTERS/SETTERS'''
	def _sequence_position_id(self, refseq, pos): #Only for ungapped (basically to use in _search_segments())
		if not isinstance(self._idx[refseq], list):
			return super(IndexedSeqAli, self)._sequence_position_id(refseq, pos)
		else:
			return self._idx[refseq][pos - 1]

	def _sequence_position_from_alignment_position(self, refseq, pos):
		#Considering over alignment with everything
		#if the alignemnt position is gap, we move up to the next position
		if not isinstance(self._idx[refseq], list):
			return super(IndexedSeqAli, self)._sequence_position_from_alignment_position(refseq, pos)
		else:

			data = self._seq[refseq][:pos]
			c    = 1
			while (data.endswith('-') and (data <= len(self._seq[refseq]))):
				data = self._seq[refseq][:pos + c]
				c += 1
			return self._idx[refseq][len(data) - 1]

	def _alignment_position_from_sequence_position(self, refseq, pos):
		if not isinstance(self._idx[refseq], list):
			return super(IndexedSeqAli, self)._alignment_position_from_sequence_position(refseq, pos)

		pos = IndexedNum(pos)
		if not self._seq[refseq].is_gapped:
			try: return self._idx[refseq].index(pos) + 1
			except ValueError as e: raise IndexError(e.message)
		else:
			c = 0
			pregaps = True
			for x in range(len(self._seq[refseq])):
				if pregaps and self._seq[refseq][x] == '-': continue
				pregaps = False
				if pos == self._idx[refseq][c]:
					return x + 1
				if x+1 < len(self._seq[refseq]) and self._seq[refseq][x+1] != '-':
					c += 1
			raise IndexError

	'''FUNCTIONS'''
	def add_complex_index(self, refseq, index, pregaps = 0):
		if isinstance(self._idx[refseq], list):
			raise AttributeError('A complex index cannot be modified')
		if not isinstance(index, list):
			raise AttributeError('A complex index is a list of IndexedNum or integers')
		if len(index) < len(re.sub('-','',self._seq[refseq].sequence)):
			raise AttributeError('The given index is not long enought to cover que assigned sequence')

		init_value  = self._idx[refseq]
		part_to_add = index[init_value - 1 - pregaps:]
		self._idx[refseq] = []
		count             = 0

		for x in range(len(self._seq[refseq])):
			if self._seq[refseq][x] != '-':
				self._idx[refseq].append(IndexedNum(part_to_add[count]))
				count += 1

		for x in range(len(self._segment[refseq])):
			if self._segment[refseq][x] != '-':
				self._segment[refseq][x] = self._sequence_position_id(refseq, self._segment[refseq][x] - init_value + 1)

		self._has_complex_index = True

	'''OVEWRITTE OBJECT FUNCTIONS (FROM PARENT)'''
	def __getitem__(self, key):
		# We correct the count as one would try to ask for the first position as 1
		# Runs by alignment position
		# if not self._has_complex_index:
		# 	return super(IndexedSeqAli, self).__getitem__(key)
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
				complex_indexes = {}
				complex_pregaps = {}
				for i in range(*index):
					data = self[i]
					for x in range(len(data)):
						sequences[x] += data[x]
					if i == len(self): break

				for x in range(self.number_of_sequences):
					if not isinstance(self._idx[x], list):
						seqInits[x] = self._sequence_position_from_alignment_position(x,index[0])
					else:
						seqInits[x] = self._idx[x].index(self._sequence_position_from_alignment_position(x,index[0])) + 1
						complex_indexes[x] = self._idx[x]
						complex_pregaps[x] = Counter(str(self.sequences[x]).split(str(sequences[x]))[0])['-']
				newali = self.__class__(sequences, seqInits)
				for refident in complex_indexes:
					newali.add_complex_index(refident, complex_indexes[refident], complex_pregaps[refident])
				if not self._alipatt is None:
					newali.add_alignment_pattern(self._alipatt[index[0]-1:index[1]-1:index[2]],self._aliptmeth)
				return newali
