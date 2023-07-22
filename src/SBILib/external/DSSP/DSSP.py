"""
DSSP

author: jbonet
date:   03/2013

@oliva's lab
"""

global exposed_threshold
exposed_threshold = 2.5 #Why? DSSP exposition refers only to secondary chain (aprox.) see PMID: 15768403

from SBILib.data import aminoacids3to1, aminoacids_surface

class DSSP(object):
    """
    The {DSSP} stores the dssp prediction
    """
    def __init__(self, secondary_structure, accessibility, AAtype):

        self._ss         = secondary_structure
        self._access     = int(accessibility)
        self._type       = AAtype
        if AAtype in aminoacids_surface:
            self._access10   = int((10 * self._access) / aminoacids_surface[AAtype])
            self._accesscode = self._codifyaccess((float(self._access)/ aminoacids_surface[AAtype]) * 100)
        else:
            self._access10   = 1
            self._accesscode = 1
        self._exposed    = self._access10 > exposed_threshold
        self._rnhoa      = None
        self._enhoa      = None
        self._rohna      = None
        self._eohna      = None
        self._rnhob      = None
        self._enhob      = None
        self._rohnb      = None
        self._eohnb      = None

    #
    # ATTRIBUTES
    #
    @property
    def secondary_structure(self): return self._ss

    @property
    def aminoacid(self):            return self._type

    @property
    def accessibility(self):       return self._access

    @property
    def accessibility10(self):     return self._access10

    @property
    def accesscode(self):          return str(self._accesscode)

    @staticmethod
    @property
    def exposition_threshold():
        global exposed_threshold
        return exposed_threshold

    @property
    def exposed(self):             return self._exposed

    @property
    def nhoa(self):
        return self._rnhoa, self._enhoa

    @property
    def ohna(self):
        return self._rohna, self._eohna

    @property
    def nhob(self):
        return self._rnhob, self._enhob

    @property
    def ohnb(self):
        return self._rohnb, self._eohnb

    #
    # FUNCTIONS
    #
    def add_hydrogen_links(self, nhoa, ohna, nhob, ohnb):
        self._rnhoa = int(  nhoa.split(',')[0].strip())
        self._enhoa = float(nhoa.split(',')[1].strip())
        self._rohna = int(  ohna.split(',')[0].strip())
        self._eohna = float(ohna.split(',')[1].strip())
        self._rnhob = int(  nhob.split(',')[0].strip())
        self._enhob = float(nhob.split(',')[1].strip())
        self._rohnb = int(  ohnb.split(',')[0].strip())
        self._eohnb = float(ohnb.split(',')[1].strip())

    def _codifyaccess(self, value):
        value = float(value)
        if value == 0:
            return '*'
        elif value > 100:
            return '?'
        else:
            if   value > 0  and value <= 10:  return '1'
            elif value > 10 and value <= 20:  return '2'
            elif value > 20 and value <= 30:  return '3'
            elif value > 30 and value <= 40:  return '4'
            elif value > 40 and value <= 50:  return '5'
            elif value > 50 and value <= 60:  return '6'
            elif value > 60 and value <= 70:  return '7'
            elif value > 70 and value <= 80:  return '8'
            elif value > 80 and value <= 90:  return '9'
            elif value > 90 and value <= 100: return '#'

