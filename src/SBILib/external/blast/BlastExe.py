'''
Check BLAST user's guide at: http://www.ncbi.nlm.nih.gov/books/NBK1763/

jbonet @ boliva's lab 2013
'''
import os, sys, time, copy, re, configparser

from SBILib.beans.Executable    import Executable
from SBILib.beans.File     import File
from SBILib.sequence  import Fasta      as Fasta
from SBILib.error     import BlastError as BE
from SBILib           import SBIglobals

from SBILib.external.blast import BlastResult as BR
from SBILib.external.blast import BlastHit    as BH

global default_configuration_file
default_configuration_file = os.path.join(os.path.normpath(os.path.join(os.path.dirname(__file__),'..')),'configSBI.txt')

class BlastExe(object):

    def __init__(self, database, search_type = 'prot'):

        #Search Type Check
        if search_type not in set(['prot','nucl']):
            raise BE(-10)
        self._search_type = search_type

        #Blast executable configuration
        self._configurator = configparser.RawConfigParser(allow_no_value=True)
        self._configurator.read(os.getenv('SBI_CONFIG_FILE',default_configuration_file))
        self._exe    = Executable(executable    = self._configurator.get('blast','executable'),
                                  path          = self._configurator.get('blast','path'),
                                  variable_path = self._configurator.get('blast','variable_path'))

        #Database Configuration
        self._database = self._check_database(os.path.abspath(database))
        if os.path.isfile(self._database.file.full + ".idx"):
            self._idx = File(file_name = self._database.file.full + ".idx", action = 'r')
        else:
            self._idx = None

        #Adding fixed blast parameters
        self._exe.add_attribute(self._database.file.full, '-db')
        self._exe.add_attribute('5', '-outfmt')
        self._exe.add_parameter('-lcase_masking')

        SBIglobals.alert('debug', self, 'New Blast Executable created.\nBlast executable at {0}\n'.format(self._exe.full_executable))

        self._selfHit     = False
        self._hitIDformat = 'single'
        self._overwritte  = False
        self._clean_files = True

    '''ATTRIBUTES'''
    @property
    def database(self): return self._database

    @property
    def selfHit(self):        return self._selfHit
    @selfHit.setter
    def selfHit(self, vaule): self._selfHit = value

    @property
    def hitIDformat(self): return self._hitIDformat
    @hitIDformat.setter
    def hitIDformat(self, value):
        if value not in set(['single','double','all']):
            self._hitIDformat = 'single'
        else:
            self._hitIDformat = value

    @property
    def overwritte(self):       return self._overwritte
    @overwritte.setter
    def overwritte(self, value): self._overwritte = value

    @property
    def clean_files(self):        return self._clean_files
    @clean_files.setter
    def clean_files(self, value): self._clean_files = value

    '''FUNCTIONS'''
    def add_attribute(self, attribute_value, attribute_id):
        if attribute_id in set(['-db','-outfmt']):
            raise AttributeError('The parameters in {0} cannot be altered'.format(set(['-db','-outfmt'])))

        self._exe.add_attribute(str(attribute_value), attribute_id)

    def execute_query_seq(self, sequenceID = None, sequence = None, blast_input_file = None, blast_output_file = None):
        if sequenceID is None and sequence is None:
            raise AttributeError('Either a sequence or ID is needed to perform the blast')

        if isinstance(sequenceID, (list,set,tuple)):
            raise AttributeError('Blasts can only be executed one at a time due to XML output restrictions')

        if sequenceID is None:
            sequenceID = 'QuerySequence'
        #Given only a code implies that the protein of interest is in the database itself
        if sequence is None:
            grabbedSequence = self._database.retrieve(sequenceID)
            sequenceID = grabbedSequence.id
            sequence   = grabbedSequence.sequence

        if len(re.sub(r'[Xx]','',sequence)) == 0: #All the sequence is unknown, it will crash blast
            return BR.BlastResult(queryname=sequenceID, querylength=len(sequence))

        file_prefixes    = ".".join([str(os.getpid()), str(int(time.process_time()*100000))])
        if blast_input_file is None:
            temp_input_name  = os.path.join(os.getcwd(), file_prefixes + ".tmp.fa")
        else:
            temp_input_name  = blast_input_file
        if blast_output_file is None:
            temp_output_name = os.path.join(os.getcwd(), file_prefixes + ".blast.xml.out")
        else:
            temp_output_name = blast_output_file

        QueryFasta = Fasta.build(file_name = temp_input_name, sequenceID = sequenceID, sequence = sequence, force = True)

        self._execute(input_file = QueryFasta, output_file = temp_output_name)

        BlastResult = self._parse_blast(sequence, temp_output_name)

        if self.clean_files:
            self._clean([temp_input_name,temp_output_name])

        return BlastResult

    def execute_query(self, query_file = None, blast_output_file = None):
        if isinstance(query_file, str) or isinstance(query_file, File):
            newFasta = Fasta(fasta_file = query_file)
        elif isinstance(query_file, Fasta):
            newFasta = query_file

        if newFasta.is_multifasta:
            raise BE(code = -4, value = newFasta.file.full)

        newFasta.load()
        if len(re.sub(r'[Xx]','',newFasta.sequences[0].sequence)) == 0: #All the sequence is unknown, it will crash blast
            return BR.BlastResult(queryname=newFasta.sequences[0].id, querylength=len(newFasta.sequences[0].sequence))

        if blast_output_file is None:
            temp_output_name = os.path.join(os.getcwd(), newFasta.file.prefix + "." + str(os.getpid()) + ".blast.xml.out")
        else:
            temp_output_name = blast_output_file

        self._execute(input_file = newFasta, output_file = temp_output_name)

        BlastResult = self._parse_blast(newFasta.sequences[0].sequence, temp_output_name)

        if self.clean_files:
            self._clean([temp_output_name])

        return BlastResult

    '''PRIVATE FUNCTIONS'''
    def _execute(self, input_file, output_file):
        if not os.path.isfile(output_file) or self.overwritte:
            final_executable = copy.deepcopy(self._exe)
            final_executable.add_attribute(input_file.file.full, '-query')
            final_executable.add_attribute(output_file, '-out')

            try:
                final_executable.execute()
            except SystemError as e:
                psiblast_default_warning = 'Warning: Composition-based score adjustment conditioned on sequence properties and unconditional composition-based score adjustment is not supported with PSSMs, resetting to default value of standard composition-based statistics'
                selenocysteine_warning   = 'Selenocysteine \(U\) at position'
                if not bool(re.search(psiblast_default_warning,str(e))) and not bool(re.search(selenocysteine_warning,str(e))):
                    raise BE(code = -1, value = str(e))

    def _clean(self, files):
        if not SBIglobals.debug:
            for temp_file in files:
                os.unlink(temp_file)

    def _parse_blast(self, query_sequence, blast_output_file):
        from .blast_parser import parse_blast
        return parse_blast(query_sequence, blast_output_file, self._selfHit, self._hitIDformat)

    def _check_database(self, database):
        #Database file does not exist
        if not os.path.isfile(database):
            raise BE(code = -5, value = database)

        #Database is not formated (if dbformatexe is added in the configuration path it will be autoformated)
        formatdb_files = []
        if self._search_type == 'nucl':     formatdb_sufix = ['.nhr','.nin','.nsq']
        elif self._search_type == 'prot':   formatdb_sufix = ['.phr','.pin','.psq']
        for sufix in formatdb_sufix:
            if not os.path.isfile(database + sufix):
                formatdb_files.append(database + sufix)

        if len(formatdb_files) > 0:
            try:
                dbexe  = Executable(executable    = self._configurator.get('blast','dbformatexe'),
                                    path          = self._configurator.get('blast','path'),
                                    variable_path = self._configurator.get('blast','variable_path'))

                SBIglobals.alert('debug', self, 'Trying to format de DB {0} to perform a blast search.\n'.format(database))

                dbexe.add_attribute(database, '-in')
                dbexe.add_attribute(self._search_type, '-dbtype')

                SBIglobals.alert('debug', self, 'Executing command {0}\n'.format(dbexe))

                dbexe.execute()

            except configparser.NoOptionError as e:
                raise BE(code = -6, value = formatdb_files)
            except SystemError as e:
                raise BE(code = -11, value = e)

        return Fasta(fasta_file = database)


