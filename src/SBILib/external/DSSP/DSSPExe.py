'''
jbonet @ boliva's lab 2013
'''

import os
import sys
import configparser

from SBILib.beans.Executable import Executable
from .DSSP                 import DSSP       as DSSP


class DSSPExe(object):
    
    default_config_file = os.path.join(
        os.path.normpath(
            os.path.join(
                os.path.dirname(__file__), '..')),
        'configSBI.txt')

    #DSSP executable configuration
    _config = configparser.RawConfigParser(allow_no_value=True)
    _config.read(os.getenv('SBI_CONFIG_FILE', default_config_file))
    _exe = None

    def __init__(self, pdbfile, dsspfile, cleanpdb=False, cleandssp=False):

        self._pdbfile  = pdbfile
        self._dsspfile = dsspfile
        self._dsspdata = []
        self._gapped   = False

        if DSSPExe._exe is None:
            DSSPExe._exe = Executable(
                executable=DSSPExe._config.get('dssp', 'executable'),
                path=DSSPExe._config.get('dssp', 'path'))

        self._execute()
        self._parse()
        self._clean(cleanpdb, cleandssp)

    @property
    def dsspdata(self):
        return self._dsspdata

    @property
    def gapped(self):
        return self._gapped

    @property
    def empty_dssp(self):
        return DSSP(secondary_structure = '-',
                    accessibility       = 1,
                    AAtype              = 'X')

    def _execute(self):
        print((self._pdbfile,self._dsspfile))
        self._exe.add_parameter(self._pdbfile)
        self._exe.add_parameter(self._dsspfile)
        try:
            self._exe.execute(silent=True)
        except SystemError as e:
            msg = 'Some error occurred while executing dssp\n{0}\n'.format(e)
            sys.stderr.write(msg)

        self._exe.clean_command()

    @staticmethod
    def dynamic(executable, path):
        DSSPExe._exe = Executable(executable=executable, path=path)

    def _parse(self):
        file_fd    = open(self._dsspfile)
        read       = False
        continuity = -1000
        readline   = 0
        for line in file_fd:
            if line.startswith("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC"):
                read = True
                continue
            if read:
                if line[13:14] != '!':
                    res_num = int(line[6:10].strip())
                    ss      = line[16:17]
                    buried  = int(line[35:38].strip())
                    aa      = line[13:15].strip()
                    if ss == " ":
                        ss = "-"
                    self._dsspdata.append(DSSP(secondary_structure = ss,
                                               accessibility       = buried,
                                               AAtype              = aa))
                    self._dsspdata[-1].add_hydrogen_links(line[39:50],
                                                          line[50:61],
                                                          line[61:72],
                                                          line[72:84])
                    if readline > 0:
                        if res_num != continuity + 1:
                            self._gapped = True
                        continuity = res_num
                    readline += 1
                else:
                    msg = "truncated chain!{0}\n".format(self._dsspfile)
                    sys.stderr.write(msg)
                    self._gapped = True
        file_fd.close()

    def _clean(self, cleanpdb, cleandssp):
        if cleanpdb:
            os.unlink(self._pdbfile)
        if cleandssp:
            os.unlink(self._dsspfile)
