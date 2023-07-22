'''
jbonet @ boliva's lab 2013
'''

import os, sys, time, copy, re, configparser

from SBILib.beans.Executable import Executable
from SBILib.beans.File       import File
from SBILib.external.CDhit   import CDhitList  as CHL

global default_configuration_file
default_configuration_file = os.path.join(os.path.normpath(os.path.join(os.path.dirname(__file__),'..')),'configSBI.txt')

class CDhitExe(object):

    def __init__(self, fasta, threshold, output_dir = None, execute = True):
        fasta = os.path.abspath(fasta)
        if output_dir is None:
            output_dir = os.path.split(fasta)[0]

        if   threshold >= 0.7: word = 5
        elif threshold >= 0.6: word = 4
        elif threshold >= 0.5: word = 3
        elif threshold >= 0.4: word = 2
        else:                  word = 1

        #CDhit executable configuration
        self._configurator = configparser.RawConfigParser(allow_no_value=True)
        self._configurator.read(os.getenv('SBI_CONFIG_FILE',default_configuration_file))
        self._exe    = Executable(executable    = self._configurator.get('cd-hit','executable'),
                                  path          = self._configurator.get('cd-hit','path'))

        self._input     = fasta
        output_file     = os.path.split(fasta)[1] + '.' + str(threshold).replace('.','_')
        self._output    = os.path.join(output_dir, output_file)
        self._threshold = str(threshold)
        self._word      = word

        if execute:
            self._execute()

    @property
    def output_file(self):
         return self._output + '.clstr'

    def _execute(self):
        self._exe.add_attribute(self._input,     '-i')
        self._exe.add_attribute(self._output,    '-o')
        self._exe.add_attribute(self._threshold, '-c')
        self._exe.add_attribute('1',             '-g')
        self._exe.add_attribute(self._word,      '-n')

        try:
            self._exe.execute(silent = True)
        except SystemError as e:
            sys.stderr.write('Some error occurred while executing cd-hit\n{0}\n'.format(e))

    def parse(self):
        return CHL.CDhitList(self._output + '.clstr')

    @staticmethod
    def parse_cdhit_output(self, cdhitfile):
        return CHL.CDhitList(cdhitfile)
