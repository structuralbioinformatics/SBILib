"""
BlastHit

The BlastHit object stores all the data related to a blast XML output hit.
    
That means the data inside each <Hit> tag of the blast XML file which corresponds to a single blast match to the query.

WARNING!: It has to be a blastpgp from a single sequence; multi-fasta blasts are not compatible
          with xml.dom.minidom parsing capabilities.

                                   ############################################ 
                                        jbonet @ boliva's lab        2011
"""
import re

from SBILib.sequence import IndexedSeqAli as IndexedSeqAli
from collections      import Counter
class BlastHit(IndexedSeqAli):    
    """
    The BlastHit object stores all the data related to a blast XML output hit.
    
    That means the data inside each <Hit> tag of the blast XML file which corresponds to a single blast match to the query.
    
    WARNING!: It has to be a blast from a single sequence; multi-fasta blasts are not compatible
              with xml.dom.minidom parsing capabilities.
    """

    global VALID_ALIGNMENT_SEQ_TYPES
    VALID_ALIGNMENT_SEQ_TYPES = {'query': 0, 0:0, 'hit': 1, 1:1}

    def __init__(self, name, length, iteration = 0,    e_value   = None, align_length      = None,
                 identities = None,  positives = None, gaps      = None, qseq = None, hseq = None, 
                 qpos       = None,  hpos      = None, score_seq = None):
        """
        BlastResult Constructor [PARAMETERS]: [ALL MANDATORY]
            (str)    name           -> name of the hit sequence              -> self.sequenceID
            (int)    length         -> length of the hit sequence            -> self.length
            (int)    iteration      -> iteration to which the hit belongs    -> self.iteration
            (str)    e_value        -> alignment e-value                     -> self.e_value
            (int)    align_length   -> alignment length                      -> self.align_length
            (int)    identities     -> alignment identities                  -> self.identities
            (int)    positives      -> alignment positives                   -> self.positives
            (int)    gaps           -> alignment gaps                        -> self.gaps
            (str)    qseq           -> query sequence aligned                -> self.query_seq
            (str)    hseq           -> hit sequence aligned                  -> self.hit_seq
            (arr)    qpos           -> positions of the aligned query        -> self.qpos[]
            (arr)    hpos           -> positions of the aligned hit          -> self.hpos[]
            
            
            BlastResult [INTERNAL VARIABLES]
            (str)    score_seq      -> score sequence summary
        """
        super(BlastHit, self).__init__(sequences     = [str(qseq),str(hseq)], 
                                       sequenceInits = [qpos,hpos], 
                                       identities    = identities,
                                       positives     = positives,
                                       gaps          = gaps)

        self.add_alignment_pattern(alipatt = score_seq, method = 'blast')

        self._sequenceID   = str(name)                # Hit Name
        self._length       = int(length)              # Hit Length
        self._iteration    = int(iteration)           # Iteration to which the hit belongs
        
        self._e_value      = e_value                  # Alignment e-value
        self._align_length = int(align_length)        # Alignment Length

        self._used         = False
    
    '''ATTRIBUTES (RENAMED FROM PARENT)'''
    @property
    def query_seq(self):            return str(self._seq[0].sequence)
    @property
    def ungapped_query_seq(self):   return str(re.sub('-','',self.query_seq))
    @property
    def hit_seq(self):              return str(self._seq[1].sequence)
    @property
    def ungapped_hit_seq(self):     return str(re.sub('-','',self.hit_seq))
    @property
    def query_pos(self):            return self._segment[0]
    @property
    def hit_pos(self):              return self._segment[1]

    '''ATTRIBUTES'''
    @property
    def sequenceID(self):           return self._sequenceID
    @property
    def length(self):               return self._length
    @property
    def iteration(self):            return self._iteration
    @property
    def e_value(self):    
        try:          
            return float(self._e_value)
        except:
            return None
    @property
    def align_length(self):         return self._align_length
    @property
    def is_used(self):              return self._used
    @is_used.setter
    def is_used(self, value):       self._used = value

    '''MIX ATTRIBUTES'''
    @property
    def query_threshold_coord(self):        return [self._first_segment_value(0), self._last_segment_value(0)]
    @property
    def query_align_segment_length(self):   return self._aligned_sequence_length(0)
    @property
    def hit_threshold_coord(self):          return [self._first_segment_value(1), self._last_segment_value(1)]
    @property
    def hit_align_segment_length(self):     return self._aligned_sequence_length(1)

    def _get_valid_alignment_seq_types(self):
        global VALID_ALIGNMENT_SEQ_TYPES
        return VALID_ALIGNMENT_SEQ_TYPES

    '''OVERWRITTE PARENT GETTERS TO ADAPT TO refseq=["query","hit"]'''
    def get_respective_coordinate(self, refseq, coordinate):
        refseq = self._check_valid_alignment_seq_type(refseq)
        return super(BlastHit, self).get_respective_coordinate(refseq, int(not refseq), coordinate)

    def get_coverage_of_sequence(self, refseq,  fulllength = None, ini = None, end = None):
        """
        METHOD: get_coverage_of_sequence(refseq =["query","hit", fulllength = None, ini = None, end = None))

            __DEVEL_INFO__
            => Overwrittes from parent (SeqAli)
        
            __EVALUATOR__
            =>Returns coverage ratio (0<->1) of the quey/hit over the full sequence (for query fulllength is required)
            
            =>If a *range* is given, returns __the sequence coverage for the given range__
             (adjusted to actually aligned regions)

            =>[TIP] To obtain __the sequence coverage for the aligned region over itself__:
                ini = self.query_threshold_coord[0] or self.hit_threshold_coord[0]
                AND
                end = self.query_threshold_coord[1] or self.hit_threshold_coord[1]
        """
        refseq = self._check_valid_alignment_seq_type(refseq)
        if bool(refseq): #hit
            fulllength = self.length 
        else:            #query
            if fulllength is None:
                raise AttributeError('To check the query coverage the query length must be given')

        return super(BlastHit, self).get_coverage_of_sequence(VALID_ALIGNMENT_SEQ_TYPES[refseq],fulllength, ini, end)
    
    def get_min_coverage_of_sequence(self, fulllength):
        return min(self.get_coverage_of_sequence('hit'),self.get_coverage_of_sequence('query',fulllength))

    def get_coverage_of_full_sequence_segment(self, refseq):
        refseq = self._check_valid_alignment_seq_type(refseq)
        if bool(refseq): #hit
            ini = self.hit_threshold_coord[0]
            end = self.hit_threshold_coord[1]
            length = len(self.hit_seq.replace('-',''))
        else: #query
            ini = self.query_threshold_coord[0]
            end = self.query_threshold_coord[1]
            length = len(self.query_seq.replace('-',''))

        return super(BlastHit, self).get_coverage_of_sequence(VALID_ALIGNMENT_SEQ_TYPES[refseq],length)

    def get_min_coverage_of_full_sequence_segment(self):
        return min(self.get_coverage_of_full_sequence_segment('hit'),self.get_coverage_of_full_sequence_segment('query'))

    def get_section_from_sequence_position(self, refseq, ini, end):
        refseq = self._check_valid_alignment_seq_type(refseq)
        return super(BlastHit, self).get_section_from_sequence_position(refseq, ini, end)

    def correct_hit_count(self, new_index):
        if isinstance(self._idx[1], list):
            raise AttributeError('Complex index of hit cannot be modified')

        if not isinstance(new_index, list):
            self.increment_sequence_index(1, new_index)
        else:
            self.add_complex_index(1, new_index)

    def correct_query_count(self, new_index):
        if isinstance(self._idx[0], list):
            raise AttributeError('Complex index of querycannot be modified')

        if not isinstance(new_index, list):
            self.increment_sequence_index(0, new_index)
        else:
            self.add_complex_index(0, new_index)
    
    def _check_valid_alignment_seq_type(self, seq_type): 
        """
        METHOD: _check_valid_alignment_seq_type()
        
            __STATUS_CHECKER__
            =>Returns BOOLEAN:
                                True if the seq_type is accepted
                                False otherwise
        """
        if not seq_type in self._get_valid_alignment_seq_types():
            raise ValueError( "'ref_seq' must be in %s. Recieved 'ref_seq' is %s" %(repr(self._get_valid_alignment_seq_types()), seq_type) )
        return self._get_valid_alignment_seq_types()[seq_type]
        
    def overlap(self, blastHit):
        #It returns a range between 0<->1 of how much the actual blastHit and the one passed
        #overlap in the query sequence.
        #The overlap is measured OVER THE SHORTEST lenght
        return super(BlastHit, self).overlap(1, 0, blastHit)
    
    """
        toString
    """
    def __str__(self):
        """
        METHOD: __str__() IS DEFINED
        """

        return "\t".join((self.sequenceID, str(self.length), str(self.identities), str(self.positives), str(self.gaps),
                          "{:.3}".format(self.e_value), self.query_seq, self.hit_seq, self.formatPositions()))

        '''OVEWRITTE OBJECT FUNCTIONS (FROM PARENT)'''
    def __getitem__(self, key):
        # We correct the count as one would try to ask for the first position as 1
        # Runs by alignment position
        # if not self._has_complex_index:
        #     return super(IndexedSeqAli, self).__getitem__(key)
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
                newali = self.__class__(name = self._sequenceID, length = self._length, iteration = self._iteration,
                                        qseq = sequences[0], hseq = sequences[1], hpos = seqInits[1], qpos = seqInits[0],
                                        align_length = len(sequences[0]), score_seq = self._alipatt[index[0]-1:index[1]-1:index[2]])
                for refident in complex_indexes:
                    newali.add_complex_index(refident, complex_indexes[refident],complex_pregaps[refident])
                if not self._alipatt is None:
                    newali.add_alignment_pattern(self._alipatt[index[0]-1:index[1]-1:index[2]],self._aliptmeth)
                return newali

            
