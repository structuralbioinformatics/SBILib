"""
Sequence

author: jbonet
date:   05/2013


@oliva's lab
"""

import re, copy
from collections import Counter

global AVAILABLE_FORMATS
AVAILABLE_FORMATS = set(['TAB', 'FASTA'])

global gapdefinition
gapdefinition = '[-x]'

class Sequence(object):

	def __init__(self, seqID = '', sequence = ''):
		global gapdefinition

		if bool(re.search(r'\s',seqID)):
			self._seqID = seqID.split()[0]
			self._info  = " ".join(seqID.split()[1:])
		else:
			self._seqID = seqID
			self._info  = ''

		self._sequence = sequence
		self._gapped   = bool(re.search(gapdefinition,self._sequence))

	'''ATTRIBUTES'''
	@property
	def id(self): 			return self._seqID
	@id.setter
	def id(self,value): 	  
		if bool(re.search(r'\s',value)):
			self._seqID = value.split()[0]
			self._info  = " ".join(value.split()[1:])
		else:
			self._seqID = value

	@property
	def sequence(self):	return self._sequence
	@sequence.setter
	def sequence(self, value): 
		self._sequence = value
		self._gapped   = bool(re.search('-',self._sequence))

	@property
	def info(self):		return self._info
	@info.setter
	def info(self, value):
		self._info = value

	@property
	def is_gapped(self):return self._gapped

	'''FUNCTIONS'''
	def contains(self, sequence):
		if isinstance(sequence, str):
			return bool(re.search(sequence, self.sequence))
		elif isinstance(sequence, Sequence):
			return bool(re.search(sequence.sequence, self.sequence))
		else:
			return NotImplemented

	def contained(self, sequence):
		if isinstance(sequence, str):
			return bool(re.search(self.sequence, sequence))
		elif isinstance(sequence, Sequence):
			return bool(re.search(self.sequence, sequence.sequence))
		else:
			return NotImplemented

	def format(self, format = 'TAB'):
		global AVAILABLE_FORMATS
		if format.upper() not in AVAILABLE_FORMATS:
			raise AttributeError('format option not available')

		if format.upper() == 'TAB':
			return '{0.id}\t{0.sequence}'.format(self)
		elif format.upper() == 'FASTA':
			return '>{0.id}\n{0.sequence}'.format(self)

	def tokenize(self, token_coding = None):
		if token_coding is None:
			return copy.deepcopy(self._sequence)
		if token_coding == 'binary':
			tc = ['\D','1','-','0']
			return re.sub(tc[0],tc[1],re.sub(tc[2],tc[3],self._sequence))

	def aa_frequency(self):
		global gapdefinition
		return dict([(x,float(y)/len(self)) for x,y in Counter(re.sub(gapdefinition,'',self._sequence)).items()])

	def duplicate(self):
		return copy.deepcopy(self)

	def do_ungap(self):
		global gapdefinition
		if self._gapped:
			self._sequence = re.sub(gapdefinition,'',self._sequence)

	'''DEFAULT FUNCTIONS OVERWRITTE'''
	def append(self, sequence):
		if isinstance(sequence, str):
			self._sequence += sequence
			if not self._gapped:
				self._gapped = bool(re.search('-',self._sequence))
		elif isinstance(sequence, (tuple, list)):
			self._sequence += "".join(sequence)
			if not self._gapped:
				self._gapped = bool(re.search('-',self._sequence))
		else:
			raise AttributeError()

	def __len__(self):
		return len(self._sequence)

	def __eq__(self, other):
		if isinstance(other, Sequence):
			return self.sequence == other.sequence
		return NotImplemented

	def __ne__(self, other):
		result = self.__eq__(other)
		if result is  NotImplemented:
			return result
		return not result

	def __lt__(self, other):
		if isinstance(other, Sequence):
			return len(self) < len(other)
		return NotImplemented

	def __gt__(self, other):
		if isinstance(other, Sequence):
			return len(self) > len(other)
		return NotImplemented

	def __getitem__(self, key):
		try: 		
			int(key)
			return self._sequence[int(key)]
		except:
			if not isinstance(key, slice): 	
				raise TypeError
			else:
				return self._sequence[key]

	def __iter__(self):
		for s in self._sequence:
			yield s

	def __repr__(self):
		return '<{0.__class__.__name__}: [{0.id}, {0.sequence}]>'.format(self)

	def __str__(self):
		return '{0.id}\t{0.sequence}'.format(self)
