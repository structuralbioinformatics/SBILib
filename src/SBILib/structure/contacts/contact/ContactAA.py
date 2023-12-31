"""
ContactAA

author: jbonet
date:   02/2014

@oliva's lab
"""

from .   import Contact
from SBILib import SBIglobals

class ContactAA(Contact):

    available_distance_types = set(["min","ca","cb","geometric","backbone"])
    atomic_types             = set(["min","ca","cb","backbone"])
    description              = 'AminoAmino'

    def __init__(self, aminoacid1, aminoacid2, threshold_type = "cb", threshold_distance = 12):
        
        self._check_type(requested_type = threshold_type)

        self._threshold_type     = threshold_type
        self._threshold_distance = threshold_distance

        super(ContactAA, self).__init__(residue1 = aminoacid1, residue2 = aminoacid2)

    #
    # ATTRIBUTES
    #
    @property
    def aminoacid1(self):         return self._residue1

    @property
    def aminoacid2(self):         return self._residue2

    @property
    def min_distance(self):       return self._distance["min"][2]
    @property
    def min_atoms(self):          return self._distance["min"][0:2]

    @property
    def ca_distance(self):        return self._distance["ca"][2]
    @property
    def ca_atoms(self):           return self._distance["ca"][0:2]

    @property
    def cb_distance(self):        return self._distance["cb"][2]
    @property
    def cb_atoms(self):           return self._distance["cb"][0:2]

    @property
    def geometric_distance(self): return self._distance["geometric"][2]

    @property
    def backbone_distance(self):  return self._distance["backbone"][2]
    @property
    def backbone_atoms(self):     return self._distance["backbone"][0:2]

    #
    # OVERWRITE PARENT FUNCTION
    #
    def _build(self):
        SBIglobals.alert('deepdebug', self, 'Analyzing AA Contact {0.type}:{0.number} - {1.type}:{1.number}'.format(self.aminoacid1, self.aminoacid2))

        for dist_type in self.available_distance_types:
            self._distance.setdefault(dist_type, None)

        self._distance[self._threshold_type] = self.aminoacid1.distance(self.aminoacid2, dist_type = self._threshold_type)
        SBIglobals.alert('deepdebug', self, '\tEvaluating distance {0:.3f} of {1}'.format(self._distance[self._threshold_type][2], self._threshold_type))
        if float(self._distance[self._threshold_type][2]) <= self._threshold_distance and \
           float(self._distance[self._threshold_type][2]) >= 0:
            SBIglobals.alert('deepdebug', self, '\tDistance under threshold.')
            self._underthreshold = True
            for dist_type in self._distance:
                if dist_type != self._threshold_type:
                    SBIglobals.alert('deepdebug', self, '\tGathering {0} distance'.format(dist_type))
                    self._distance[dist_type] = self.aminoacid1.distance(self.aminoacid2, dist_type = dist_type)
                    SBIglobals.alert('deepdebug', self, '\t\tDistance {0:.3f} of {1}'.format(self._distance[dist_type][2], dist_type))      

    #
    # TOSTRING
    #
    def toString(self, requested_type = "cb"):
        if requested_type != "all":
            self._check_type(requested_type = requested_type)
            ddata = self._distance[requested_type]
            if requested_type in self.atomic_types:
                return '{0.aminoacid1.type}:{0.aminoacid1.number}:{1[0].name}\t{0.aminoacid2.type}:{0.aminoacid2.number}:{1[1].name}\t{1[2]:.3f}'.format(self, ddata)
            else:
                return '{0.aminoacid1.type}:{0.aminoacid1.number}\t{0.aminoacid2.type}:{0.aminoacid2.number}\t{1[2]:.3f}'.format(self, ddata)
        else:
            return '{0.aminoacid1.type}:{0.aminoacid1.number}\t{0.aminoacid2.type}:{0.aminoacid2.number}\t{0.min_distance:.3f}\t{0.ca_distance:.3f}\t{0.cb_distance:.3f}\t{0.geometric_distance:.3f}\t{0.backbone_distance:.3f}'.format(self)
        
