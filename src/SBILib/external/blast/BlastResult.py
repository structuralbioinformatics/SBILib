"""
BlastResult

The BlastResult object contains the data from a single blastpgp XML file.

WARNING!: It has to be a blastpgp from a __single sequence__; multi-fasta blasts are not compatible
          with xml.dom.minidom parsing capabilities.

                                   ############################################
                                        jbonet @ boliva's lab        2011
"""

'''
    Supporting Modules
'''
import os, sys, pydoc, warnings, gzip
from SBILib.error import BlastError as BE
from SBILib.beans.StorableObject import StorableObject
from SBILib           import SBIglobals

try:
    import pickle as pickle
except:
    import pickle

class BlastResult(StorableObject):

    def __init__(self, queryname   = None, querylength     = 0,    blastversion = None,
                       blastmatrix = None, gap_open        = 0,    gap_extend   = 0,
                       blastdb     = 0,    BlastResultFile = None, queryseq     = None):
        """
        The BlastResult object contains the data from a single blastpgp XML file.

        @type  queryname:       String
        @param queryname:       Name of the query sequence.
                                B{MANDATORY} unless I{BlastResultFile} is given.

        @type  querylength:     Integer
        @param querylength:     Length of the query sequence.
                                B{MANDATORY} unless I{BlastResultFile} is given.

        @type  blastversion:    String
        @param blastversion:    Version of blastpgp that created the L{BlastResult}.

        @type  blastmatrix:     String
        @param blastmatrix:     Weight matrix used in the blastpgp run.

        @type  gap_open:        Integer
        @param gap_open:        Score used in the blastpgp run.

        @type  gap_extend:      Integer
        @param gap_extend:      Score used in the blastpgp run.

        @type  blastdb:         String
        @param blastdb:         Sequence database over which the search was performed.

        @type  BlastResultFile: String
        @param BlastResultFile: File containing a L{BlastResult} object to upload.

        @type  queryseq:        String
        @param queryseq:        Query sequence.

        @raise AtributeError: IF mandatory paremeters are not met and BlastResultFile is not given.
        """

        if BlastResultFile != None:
            self.__dict__ = BlastResult.load(object_file = BlastResultFile).__dict__

        else:
            self._query          = queryname         # Query Name
            self._query_length   = int(querylength)  # Query Sequence Length
            self._query_sequence = queryseq          # Query Sequence
            self._version        = blastversion      # Blast version that has been executed for the query
            self._matrix         = blastmatrix       # Similarity matrix used
            self._gap_open       = int(gap_open)     # Gap opening penalty
            self._gap_extend     = int(gap_extend)   # Gap extension penalty
            self._db             = blastdb           # Database used in the search
            self._lastiteration  = 0                 # Keep the last iteration
            self._hits           = []                # List of BlastHit objects

            self._correctedHits  = False
            self._correctFile    = None

    '''ATTRIBUTES'''
    @property
    def query(self):                return self._query
    @property
    def query_length(self):         return int(self._query_length)
    @property
    def query_sequence(self):       return self._query_sequence
    @property
    def version(self):              return self._version
    @property
    def matrix(self):               return self._matrix
    @property
    def gap_open(self):             return int(self._gap_open)
    @property
    def gap_extend(self):           return int(self._gap_extend)
    @property
    def database(self):             return self._db
    @property
    def lastiteration(self):        return int(self._lastiteration)
    @property
    def are_hits_corrected(self):    return self._correctedHits


    def get_hits(self, iteration = None, evalue = None, tz_type = None, tz_parameter = 0, overlap = 1):
        """
        Browse through the BlastHit of the BlastResult

        @type  iteration:    Integer
        @param iteration:    Iteration of interest. B{DEFAULT = last iteration}.
                             IF iteration is not a number returns as with DEFAULT.
                             IF iteration is out of range returns as with DEFAULT.
                             IF iteration == 0 returns the FIRST iteration

        @type  evalue:       Float
        @param evalue:       Maximum e-value threshold. B{DEFAULT = No e-value threshold}.

        @type  tz_type:      String
        @param tz_type:      Rost Evaluation Curve. B{DEFAULT = NO-ACTIVE}. I{OPTIONS = "ID" for identity curve
                                                                                        "SIM" for similarity curve}.

        @type  tz_parameter: Integer
        @param tz_parameter: Rost Evaluation Curve "Strength" Parameter. B{DEFAULT = 0}.

        @type  overlap:      Float
        @param overlap:      Amount of allowed overlap (0<->1) between single results. B{DEFAULT = 1}.

        @rtype: Array of L{BlastHist}

        B{TIP:} Rost Evaluation Curve is based on Rost's twilight zone:
                I{Rost, B. (1999). Twilight zone of protein sequence alignments. Protein Engineering, 12(2), 85 94.}
        """

        hits_of_interest = []

        '''
            First we check that the iteration asked is correct.
            Any error will redirect to either the first or last iteration (depending on the error)
        '''
        try: iteration = int(iteration)
        except: iteration = self.lastiteration

        if iteration == None:
            iteration = self.lastiteration
        elif iteration > self.lastiteration:
            iteration = self.lastiteration
        elif iteration == 0:
            iteration = 1

        '''
            We check that the asked options for the Rost curve are OK
            When the type is set to None, the curve is not applyed, even if tz_parameter is defined
        '''
        if tz_type != None:
            if not (tz_type == "ID" or tz_type == "SIM"):
                tz_type = None

        for hit in self._hits:
            '''
                We check out only hits of the iteration of interest
            '''
            if hit.iteration == iteration:
                '''
                    If applyes, we check the evalue
                '''
                if (evalue == None or hit.e_value <= float(evalue)):
                    '''
                        If applyes, we check Rost's curves
                    '''
                    if(tz_type == None or hit.evaluate_Rost_twilight_zone(equation = tz_type, parameter = tz_parameter)):
                        hits_of_interest.append(hit)
            '''
                As the hits are ordered, once we are over the iteration of interest, we can skip
            '''
            if hit.iteration > iteration:
                break

        if overlap == 1 or len(hits_of_interest) == 0:
            return hits_of_interest
        else:
            return self._get_nonoverlaping_hits(hitlist = hits_of_interest, overlap = overlap)

    def get_used_hits(self):
        """
        @rtype: Array of L{BlastHist}
        """
        hits_of_interest = []
        for hit in self._hits:
            if hit.is_used:
                hits_of_interest.append(hit)
        return hits_of_interest

    def _get_nonoverlaping_hits(self, hitlist, overlap = 1):
        """
        Due to the way this needs to be calculated, it will be a last step of the get_hits() function
        Returns all hits whose overlap with the already accepted ones is LOWER OR EQUAL than the given parameter
        It asumes that the hits are ordered by reliability, thus gives priority to the initial ones

        @type  hitlist: Array of L{BlastHist}
        @param hitlist: List of BlastHit to filter

        @type  overlap: Float
        @param overlap: Amount of allowed overlap (0<->1) between single results. B{DEFAULT = 1}.

        @rtype: Array of L{BlastHist}

        """
        final_hits = []
        final_hits.append(hitlist[0])
        for x in range(1,len(hitlist)):
            accepted = True
            for y in range(0,len(final_hits)):
                if hitlist[x].overlap(blastHit = final_hits[y]) > overlap:
                    accepted = False
                    break
            if accepted:
                final_hits.append(hitlist[x])

        return final_hits

    '''
        FUNCTIONS
    '''
    def correct_hit_count(self, countFile = None, countqueryFile = None, return_correction_dict = False):
        """
        Corrects the starting point of the hits.

        B{RATIONALE:}
        When blasting vs. PDB (for example), sometimes the hit positions given by blast are wrong,
        as the blast always consider the first position of the hit sequence as 1 and PDB does not.
        Given a tabulated file with:
        C{HIT_ID    FIRST_POSITION_REAL_NUMBER}
        (obviously the HIT_ID should match the format of the blasted database)

        @type  countFile:              String
        @param countFile:              Name of the "correction file"

        @type  return_correction_dict: Boolean
        @param return_correction_dict: When active the correction is not applied over the BlastHit objects and
                                       correctedHits is not set to True. The correction dictionary is returned

        @rtype:  Dictionary {String}[Integer]
        @return: Only when it is asked for through the return_correction_dict parameters

        @raise IOError:        IF the correction index file does not exist.
        @raise AttributeError: IF the BlastResult does not contain any BlastHit.
        @raise BlastError:     IF L{correctHitCount} has been called already.
        """

        if not os.path.isfile(countFile):
            raise IOError("The file %s does not exist and can not be used to correct the BlastHit positions.\n" %countFile)

        if not self.has_hits():
            warnings.warn("BlastResult of {0} has no hits to correct\n".format(self._query))

        if self.are_hits_corrected:
            raise BE(code = 3, value = "HitCounts has already been corrected\n")

        codes_of_interest = set([hit.sequenceID for hit in self._hits])
        if countqueryFile == countFile:
            codes_of_interest.add(self._query)
        start_index_dic = {}
        file_fd = open(countFile,'r')
        for line in file_fd:
            if len(line.strip()) > 0:
                if line.split('\t')[0].lstrip('>') in codes_of_interest:
                    start_index_dic[line.split('\t')[0].lstrip('>')] = line.split('\t')[1].strip().split(';')
        file_fd.close()

        if countqueryFile is not None and countqueryFile != countFile:
            file_fd = open(countqueryFile,'r')
            for line in file_fd:
                if len(line.strip()) > 0:
                    if line.split('\t')[0] == self._query:
                        start_index_dic[line.split('\t')[0]] = line.split('\t')[1].strip().split(';')
            file_fd.close()

        if return_correction_dict: return start_index_dic


        for hit in self._hits:
            '''
            This tests between the options PDB/PDB_ID or PDB_ID in case the TAB file has different codification
            '''
            try:
                hit.correct_hit_count(new_index = start_index_dic[hit.sequenceID])
                if countqueryFile is not None:
                    hit.correct_query_count(new_index = start_index_dic[self._query])
            except:
                hit_ID = hit.sequenceID.split("/")[-1]
                hit.correct_hit_count(new_index = start_index_dic[hit_ID])
                if countqueryFile is not None:
                    hit.correct_query_count(new_index = start_index_dic[self._query])
        self._correctedHits = True

    def set_last_iteration(self):
        """
        Modifies the value of self.lastiteration with the value of iteration of the last BlastHit found

        IF self.hits == [] : NO BLAST HITS HAVE BEEN FOUND then self.lastiteration remains at 0
        """
        if self._hits == []:
            self._lastiteration = 0
        else:
            self._lastiteration = self._hits[-1].iteration

    def has_hits(self):
        """
        @rtype: Boolean
        """
        if self.lastiteration == 0: return False
        else:                       return True

    def is_empty(self):
        """
        @rtype: Boolean
        """
        return not self.has_hits()

    '''
        BlastResult comparisson
    '''
    def same_blast_version(self, blastObject = None):
        """
        @type  blastObject: L{BlastResult}
        @param blastObject: BlastResult to which compare

        @rtype: Boolean

        @raise AttributeError: IF no blastObject is given
        """
        if blastObject == None:
            raise AttributeError("A second object needs to be given for the comparison\n")

        if self.version == blastObject.version: return True
        return False

    def same_blast_database(self, blastObject = None):
        """
        @type  blastObject: L{BlastResult}
        @param blastObject: BlastResult to which compare

        @rtype: Boolean

        @raise AttributeError: IF no blastObject is given
        """
        if blastObject == None:
            raise AttributeError("A second object needs to be given for the comparison\n")

        if self.database == blastObject.database: return True
        return False

    def same_blast_matrix(self, blastObject = None):
        """
        @type  blastObject: L{BlastResult}
        @param blastObject: BlastResult to which compare

        @rtype: Boolean

        @raise AttributeError: IF no blastObject is given
        """
        if blastObject == None:
            raise AttributeError("A second object needs to be given for the comparison\n")

        if self.matrix == blastObject.matrix: return True
        return False

    def same_blast_iterations(self, blastObject = None):
        """
        @type  blastObject: L{BlastResult}
        @param blastObject: BlastResult to which compare

        @rtype: Boolean

        @raise AttributeError: IF no blastObject is given
        """
        if blastObject == None:
            raise AttributeError("A second object needs to be given for the comparison\n")

        if self.lastiteration == blastObject.lastiteration: return True
        return False

    def same_blast_gap_properties(self, blastObject = None):
        """
        @type  blastObject: L{BlastResult}
        @param blastObject: BlastResult to which compare

        @rtype: Boolean

        @raise AttributeError: IF no blastObject is given
        """
        if blastObject == None:
            raise AttributeError("A second object needs to be given for the comparison\n")

        if (self.gap_open == blastObject.gap_open) and (self.gap_extend == blastObject.gap_extend): return True
        return False

    def comparable(self, blastObject = None):
        """
        @type  blastObject: L{BlastResult}
        @param blastObject: BlastResult to which compare

        @rtype: Boolean

        @raise AttributeError: IF no blastObject is given
        """
        if blastObject == None:
            raise AttributeError("A second object needs to be given for the comparison\n")

        if self.same_blast_database(blastObject   = blastObject)  and self.same_blast_version(blastObject        = blastObject) and \
           self.same_blast_matrix(blastObject     = blastObject)  and self.same_blast_gap_properties(blastObject = blastObject) and \
           self.same_blast_iterations(blastObject = blastObject): return True
        return False

    '''
        TO_STRING
    '''
    def str_compacted_blast(self, iteration = None, evalue = None, tz_type = None, tz_parameter = 0, overlap = 1):
        """
        @type  iteration:    Integer
        @param iteration:    Iteration of interest. B{DEFAULT = last iteration}.
                             IF iteration is not a number returns as with DEFAULT.
                             IF iteration is out of range returns as with DEFAULT.
                             IF iteration == 0 returns the FIRST iteration

        @type  evalue:       Float
        @param evalue:       Maximum e-value threshold. B{DEFAULT = No e-value threshold}.

        @type  tz_type:      String
        @param tz_type:      Rost Evaluation Curve. B{DEFAULT = NO-ACTIVE}. I{OPTIONS = "ID" for identity curve
                                                                                        "SIM" for similarity curve}.

        @type  tz_parameter: Integer
        @param tz_parameter: Rost Evaluation Curve "Strength" Parameter. B{DEFAULT = 0}.

        @type  overlap:      Float
        @param overlap:      Amount of allowed overlap (0<->1) between single results. B{DEFAULT = 1}.

        @rtype: String

        B{TIP:} Rost Evaluation Curve is based on Rost's twilight zone:
                I{Rost, B. (1999). Twilight zone of protein sequence alignments. Protein Engineering, 12(2), 85 94.}
        """

        lines = []

        for hit in self.get_hits(iteration=iteration, evalue=evalue, tz_type=tz_type, tz_parameter=tz_parameter, overlap=overlap):
            lines.append("\t".join((self.query,str(self.query_length),str(hit))))

        return "\n".join(lines)

    def str_PIR(self, iteration = None, result = 0, evalue = None, tz_type = None, tz_parameter = 0, overlap = 1):
        """
        @type  iteration:    Integer
        @param iteration:    Iteration of interest. B{DEFAULT = last iteration}.
                             IF iteration is not a number returns as with DEFAULT.
                             IF iteration is out of range returns as with DEFAULT.
                             IF iteration == 0 returns the FIRST iteration

        @type  result:       Integer
        @param result:       Result of interest (in order of selected results)

        @type  evalue:       Float
        @param evalue:       Maximum e-value threshold. B{DEFAULT = No e-value threshold}.

        @type  tz_type:      String
        @param tz_type:      Rost Evaluation Curve. B{DEFAULT = NO-ACTIVE}. I{OPTIONS = "ID" for identity curve
                                                                                        "SIM" for similarity curve}.

        @type  tz_parameter: Integer
        @param tz_parameter: Rost Evaluation Curve "Strength" Parameter. B{DEFAULT = 0}.

        @type  overlap:      Float
        @param overlap:      Amount of allowed overlap (0<->1) between single results. B{DEFAULT = 1}.

        @rtype: String

        B{TIP:} Rost Evaluation Curve is based on Rost's twilight zone:
                I{Rost, B. (1999). Twilight zone of protein sequence alignments. Protein Engineering, 12(2), 85 94.}
        """
        lines = []
        hit = self.get_hits(iteration=iteration, evalue=evalue, tz_type=tz_type, tz_parameter=tz_parameter, overlap=overlap)[result]

        #Counting in the possibility that the ID of the PDB hit is:
        #    PDB/2ERD_A
        hit_name = hit.sequenceID.replace("/","|").upper()
        hit_chain = hit_name.split("_")[-1]
        hit_id = hit_name.split("|")[-1].split("_")[0]
        hit_sq = hit.hit_seq + "*"
        q_sq = hit.query_seq + "*"

        header1 = ">P1;%s\nsequence:%s:%d:.:%d:.:.:.:.:." % (self.query,self.query,int(hit.query_pos[0]),int(hit.query_pos[-1]))
        header2 = ">P1;%s\nstructureX:%s:%d:%s:%d:%s:.:.:.:." % (hit_id,hit_id,int(hit.hit_pos[0]),hit_chain,int(hit.hit_pos[-1]),
                                                                 hit_chain)

        lines.append(header1)
        lines.extend([q_sq[i:i+60] for i in range(0, len(q_sq), 60)])
        lines.append(header2)
        lines.extend([hit_sq[i:i+60] for i in range(0, len(hit_sq), 60)])

        return "\n".join(lines)

    def print_query_fasta(self):
        """
        @rtype: String. B{Through STDOUT}.

        @raise IOError: IF query_sequence is not known.
        """

        if self.query_sequence == "":
            raise IOError("The sequence for %s is not known and cannot be printed\n" %self.query)

        print(">%s\n%s" %(self.query, self.query_sequence))

    def print_compacted_blast(self, iteration = None, evalue = None, tz_type = None, tz_parameter = 0, overlap = 1, file = None):
        """
        @type  iteration:    Integer
        @param iteration:    Iteration of interest. B{DEFAULT = last iteration}.
                             IF iteration is not a number returns as with DEFAULT.
                             IF iteration is out of range returns as with DEFAULT.
                             IF iteration == 0 returns the FIRST iteration

        @type  evalue:       Float
        @param evalue:       Maximum e-value threshold. B{DEFAULT = No e-value threshold}.

        @type  tz_type:      String
        @param tz_type:      Rost Evaluation Curve. B{DEFAULT = NO-ACTIVE}. I{OPTIONS = "ID" for identity curve
                                                                                        "SIM" for similarity curve}.

        @type  tz_parameter: Integer
        @param tz_parameter: Rost Evaluation Curve "Strength" Parameter. B{DEFAULT = 0}.

        @type  overlap:      Float
        @param overlap:      Amount of allowed overlap (0<->1) between single results. B{DEFAULT = 1}.

        @type  file:         String
        @param file:         Name of the output file. Redirects function output to file

        @rtype: String. B{Through STDOUT} IF file is None.
                String. B{Through FILE} otherwise.

        B{TIP:} Rost Evaluation Curve is based on Rost's twilight zone:
                I{Rost, B. (1999). Twilight zone of protein sequence alignments. Protein Engineering, 12(2), 85 94.}
        """
        if file:
            try:
                file_out = open(file, 'w')
            except:
                raise IOError("%s cannot be opened to write\n" %file)
            file_out.write("%s\n" % self.str_compacted_blast(iteration=iteration, evalue=evalue, tz_type=tz_type, tz_parameter=tz_parameter, overlap=overlap))
            file_out.close()
        else:
            print(self.str_compacted_blast(iteration=iteration, evalue=evalue, tz_type=tz_type, tz_parameter=tz_parameter, overlap=overlap))

    def str_blast_details(self):
        """
        @rtype: String
        """
        lines = []

        lines.append("Blast Version: %s" %self.version)
        lines.append("Search on matrix: %s" %self.matrix)
        lines.append("Gap open penalty: %s\nGap extension penalty: %s" %(str(self.gap_open),str(self.gap_extend)))
        lines.append("Database searched: %s\n" %self.database)

        return "\n".join(lines)

    def __str__(self):
        """
        Strings the String version of ALL the BlastHit objects contained (with all iterations!!)

        @rtype: String
        """
        lines = []

        for hit in self._hits:
            lines.append("\t".join((self.query,str(self.query_length),str(hit))))

        return "\n".join(lines)
