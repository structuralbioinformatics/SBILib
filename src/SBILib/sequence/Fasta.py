"""
Fasta

author: jbonet
date:   05/2013


@oliva's lab
"""

import subprocess, sys

from SBILib.beans.File    import File
from SBILib.sequence import Sequence
from SBILib              import SBIglobals

class Fasta(object):

    def __init__(self, fasta_file):

        if isinstance(fasta_file, str):
            self._file = File(file_name = fasta_file, action = 'r')
        elif isinstance(fasta_file, File):
            self._file = fasta_file
            self._file.action = 'r'
        else:
            raise AttributeError('Check the input of the Fasta object')

        self._is_multifasta = self._check_multifasta()

        self._sequences = []
        self._seqfinder = {}

    '''ATTRIBUTES'''
    @property
    def file(self):          return self._file

    @property
    def is_multifasta(self): return self._is_multifasta

    @property
    def sequences(self):     return self._sequences

    '''FUNCTIONS'''
    def load(self):
        for line in self.file.descriptor:
            if line.startswith('>'):
                self._sequences.append(Sequence(seqID = line.lstrip('>').strip()))
                self._seqfinder[self._sequences[-1].id] = len(self._sequences) - 1
            elif len(line.strip()) > 0:
                self._sequences[-1].append(line.strip())
        self.file.close()

    def live_show(self):
        n = 0
        for line in self.file.descriptor:
            if line.startswith('>'):
                if n > 0:
                    yield s
                n += 1
                s = Sequence(seqID = line.lstrip('>').strip())
            elif len(line.strip()) > 0:
                s.append(line.strip())
        self.file.close()
        yield s

    def retrieve(self, seqID, allbut = False, prefix_size = None):
        SBIglobals.alert('debug', self, 'Getting sequence for {0}'.format(seqID))

        if isinstance(seqID, str):
            if len(self) == 0:
                sequence = Sequence()
                read = False
                for line in self.file.descriptor:
                    if line.startswith('>'):
                        if prefix_size is not None: sid = line.lstrip('>').split()[0].strip()[:prefix_size]
                        else:                       sid = line.lstrip('>').split()[0].strip()
                        if sid == seqID:
                            sequence.id = line.lstrip('>').split()[0].strip()
                            read = True
                        elif read: break
                    elif read and len(line.strip()) > 0:
                        sequence.append(line.strip())
                self.file.close()
            else:
                if not seqID in self._seqfinder:
                    raise KeyError(seqiID)
                return self._sequences[self._seqfinder[seqID]]

        if isinstance(seqID, list): seqID = set(seqID)
        if isinstance(seqID, set):
            sequence = []
            if len(self) == 0:
                read  = False
                for line in self.file.descriptor:
                    if line.startswith('>'):
                        if prefix_size is not None: sid = line.lstrip('>').split()[0].strip()[:prefix_size]
                        else:                       sid = line.lstrip('>').split()[0].strip()
                        if (not allbut and sid in seqID) or \
                           (allbut and sid not in seqID):
                            newSeq = Sequence(seqID = line.lstrip('>').split()[0].strip())
                            sequence.append(newSeq)
                            if not allbut and prefix_size is None:
                                seqID.remove(newSeq.id)
                            read = True
                        elif read: read = False
                        elif len(seqID) == 0: break
                    elif read and len(line.strip()) > 0:
                        sequence[-1].append(line.strip())
                self.file.close()
            else:
                for queryID in seqID:
                    if not queryID in self._seqfinder:
                        raise KeyError(queryID)
                    sequence.append(self._sequences[self._seqfinder[seqID]])

        return sequence

    def get_all_ids(self):
        identifiers = ()
        if len(self) == 0:
            for line in self.file.descriptor:
                if line.startswith('>'):
                    identifiers.append(line.lstrip('>').split()[0].strip())
            self.file.close()
        else:
            for seq in self:
                identifiers.append(seq.id)
        return identifiers

    @staticmethod
    def build(file_name, sequenceID, sequence, force = False):
        newFasta = File(file_name,'w', overwrite = force)
        newSeq   = Sequence(seqID = sequenceID, sequence = sequence)
        file_dsc = newFasta.descriptor
        file_dsc.write(newSeq.format('FASTA'))
        newFasta.close()
        return Fasta(fasta_file = newFasta.full)

    @staticmethod
    def build_multifasta(file_name, sequenceList, force = False):
        newFasta = File(file_name,'w', overwrite = force)
        file_dsc = newFasta.descriptor
        for sequence in sequenceList:
            file_dsc.write(sequence.format('FASTA')+"\n")
        newFasta.close()
        return Fasta(fasta_file = newFasta.full)

    '''PRIVATE FUNCTIONS'''
    def _check_multifasta(self):

        grep_search  = ['grep','-c','>',self.file.full]
        grep_process = subprocess.Popen(grep_search, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        g_out, g_err = grep_process.communicate()

        if int(g_out.strip()) > 1: return True
        else:                      return False

    '''DEFAULT FUNCTIONS OVERWRITTE'''
    def __len__(self):
        return len(self._sequences)

    def __getitem__(self, key):
        try:    int(key)
        except: raise TypeError
        return self._sequences[int(key)]

    def __iter__(self):
        for s in self._sequences:
            yield s

    def __repr__(self):
        return '<{0.__class__.__name__}: [{0.file}]>'.format(self)

    def __str__(self):
        if len(self) == 0:
            return self.file.full
        else:
            text = []
            for seq in self._sequences:
                text.append(seq.format('FASTA'))
            return "\n".join(text)
