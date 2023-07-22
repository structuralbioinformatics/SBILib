# -*-
#
# @author: jaumebonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.github.io
#
# @date:   2016-02-15 11:21:10
#
# @Last modified by:   bonet
# @Last modified time: 12-Oct-2018
#
# -*-
from .Element import PeriodicTable
e = PeriodicTable()


class WrongAtomCodingName(Exception):
    pass


class AtomVariants(object):
    """Translator between how atoms are called in different formats"""
    def __init__(self, names):
        self.bmrb  = names[0]
        self.pdb   = names[1]
        self.ucfs  = names[2]
        self.msi   = names[3]
        self.xplor = names[4]
        self.sybyl = names[5]
        self.midas = names[6]
        self.diana = names[7]

    @staticmethod
    def known_format(self, format):
        if not format.lower() in self.__dict__:
            raise WrongAtomCodingName

    def get_name(self, name):
        return self.__getattribute__(name.lower())

    def get_atom_type(self):
        return e[self.bmrb[0]]


C    = AtomVariants(['C',    'C',    '',     'C',    'C',    'C',    'C',   'C'])
CA   = AtomVariants(['CA',   'CA',   '',     'CA',   'CA',   'CA',   'CA',  'CA'])
CB   = AtomVariants(['CB',   'CB',   '',     'CB',   'CB',   'CB',   'CB',  'CB'])
CD   = AtomVariants(['CD',   'CD',   '',     'CD',   'CD',   'CD',   'CD',  'CD'])
CD1  = AtomVariants(['CD1',  'CD1',  '',     'CD1',  'CD1',  'CD1',  'CD1', 'CD1'])
CD2  = AtomVariants(['CD2',  'CD2',  '',     'CD2',  'CD2',  'CD2',  'CD2', 'CD2'])
CE   = AtomVariants(['CE',   'CE',   '',     'CE',   'CE',   'CE',   'CE',  'CE'])
CE1  = AtomVariants(['CE1',  'CE1',  '',     'CE1',  'CE1',  'CE1',  'CE1', 'CE1'])
CE2  = AtomVariants(['CE2',  'CE2',  '',     'CE2',  'CE2',  'CE2',  'CE2', 'CE2'])
CG   = AtomVariants(['CG',   'CG',   '',     'CG',   'CG',   'CG',   'CG',  'CG'])
CG2  = AtomVariants(['CG2',  'CG2',  '',     'CG2',  'CG2',  'CG2',  'CG2', 'CG2'])
CZ   = AtomVariants(['CZ',   'CZ',   '',     'CZ',   'CZ',   'CZ',   'CZ',  'CZ'])
H    = AtomVariants(['H',    'H',    'HN',   'HN',   'HN',   'H',    '',    'HN'])
HA   = AtomVariants(['HA',   'HA',   'HA',   'HA',   'HA',   'HA',   '',    'HA'])
HB   = AtomVariants(['HB',   'HB',   'HB',   'HB',   'HB',   'HB',   '',    'HB'])
HB1  = AtomVariants(['HB1',  '1HB',  'HB1',  'HB1',  'HB1',  'HB1',  '',    'HB1'])
HB2  = AtomVariants(['HB2',  '1HB',  'HB1',  'HB1',  'HB2',  'HB2',  '',    'HB2'])
HB2b = AtomVariants(['HB2',  '2HB',  'HB2',  'HB2',  'HB2',  'HB2',  '',    'HB2'])
HB3  = AtomVariants(['HB3',  '2HB',  'HB2',  'HB2',  'HB1',  'HB1',  '',    'HB3'])
HB3b = AtomVariants(['HB3',  '3HB',  'HB3',  'HB3',  'HB3',  'HB3',  '',    'HB3'])
HD1  = AtomVariants(['HD1',  'HD1',  'HD1',  'HD1',  'HD1',  'HD1',  '',    'HD1'])
HD11 = AtomVariants(['HD11', '1HD1', 'HD11', 'HD11', 'HD11', 'HD11', '',    'HD11'])
HD12 = AtomVariants(['HD12', '2HD1', 'HD12', 'HD12', 'HD12', 'HD12', '',    'HD12'])
HD13 = AtomVariants(['HD13', '3HD1', 'HD13', 'HD13', 'HD13', 'HD13', '',    'HD13'])
HD2  = AtomVariants(['HD2',  'HD2',  'HD2',  'HD2',  'HD2',  'HD2',  '',    'HD2'])
HD2b = AtomVariants(['HD2',  '1HD',  'HD1',  'HD1',  'HD2',  'HD2',  '',    'HD2'])
HD3  = AtomVariants(['HD3',  '2HD',  'HD2',  'HD2',  'HD1',  'HD1',  '',    'HD3'])
HE1  = AtomVariants(['HE1',  'HE1',  'HE1',  'HE1',  'HE1',  'HE1',  '',    'HE1'])
HE2  = AtomVariants(['HE2',  'HE2',  'HE2',  'HE2',  'HE2',  'HE2',  '',    'HE2'])
HE2b = AtomVariants(['HE2',  '',     '',     'HE2',   '',    '',     '',    'HE2']),
HG   = AtomVariants(['HG',   'HG',   'HSG',  'HG',   'HG',   'HG',   '',    'HG'])
HG2  = AtomVariants(['HG2',  '1HG',  'HG1',  'HG1',  'HG2',  'HG2',  '',    'HG2'])
HG3  = AtomVariants(['HG3',  '2HG',  'HG2',  'HG2',  'HG1',  'HG1',  '',    'HG3'])
HG21 = AtomVariants(['HG21', '1HG2', 'HG21', 'HG21', 'HG21', 'HG21', '',    'HG21'])
HG22 = AtomVariants(['HG22', '2HG2', 'HG22', 'HG22', 'HG22', 'HG22', '',    'HG22'])
HG23 = AtomVariants(['HG23', '3HG2', 'HG23', 'HG23', 'HG23', 'HG23', '',    'HG23'])
N    = AtomVariants(['N',    'N',    '',     'N',    'N',    'N',    'N',   'N'])
NE2  = AtomVariants(['NE2',  'NE2',  '',     'NE2',  'NE2',  'NE2',  'NE2', 'NE2'])
O    = AtomVariants(['O',    'O',    '',     'O',    'O',    'O',    'O',   'O'])
OD1  = AtomVariants(['OD1',  'OD1',  '',     'OD1',  'OD1',  'OD1',  'OD1', 'OD1'])
OE1  = AtomVariants(['OE1',  'OE1',  '',     'OE1',  'OE1',  'OE1',  'OE1', 'OE1'])
OE2  = AtomVariants(['OE2',  'OE2',  '',     'OE2',  'OE2',  'OE2',  'OE2', 'OE2'])
SG   = AtomVariants(['SG',   'SG',   '',     'SG',   'SG',   'SG',   'SG',  'SG'])
