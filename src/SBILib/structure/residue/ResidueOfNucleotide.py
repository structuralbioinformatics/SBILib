"""
Nucleotide

author: jbonet
date:   02/2013

@oliva's lab
"""

from .          import Residue
from SBILib.data   import nucleic3to1, nucleic1to3
import numpy    as np

class ResidueOfNucleotide(Residue):
    """
    A {Nucleotide} collects a series of {NucleotideAtom}s
    """
    def __init__(self, number, version, Rtype, mode):
        """
        @type  number: Integer
        @param number: Residue number

        @type  version: Char
        @param version: Optional char used on pdbs to change count

        @type  type: String
        @param type: Residue type

        @type  mode: String
        @param mode: Residue mode: ATOM or HETATM
        """
        super(ResidueOfNucleotide, self).__init__(number = number, version = version, Rtype = Rtype, mode = mode)
        
        self._p        = None
        self._o3       = None

    #
    # ATTRIBUTES
    #
    @property
    def p(self):
        """
        Returns the Phosphate {NucleotideAtom}
        @rtype: {NucleotideAtom}
        """
        return self._p

    @property
    def o3(self):
        """
        Returns the Sugar O3' {NucleotideAtom}
        @rtype: {NucleotideAtom}
        """
        return self._o3

    @property
    def backbone(self):
        """
        Returns a list of the backbone {NucleotideAtom}s
        @rtype: List of {NucleotideAtom}s
        """
        return self._backbone_atoms

    @property
    def single_letter(self):
        """
        Returns the AminoAcid identifier as a single letter code
        @rtype: String
        """
        return nucleic3to1[self.type]

    @property
    def standard_type(self):
        """
        For some Nucleotide are HETATM, we can return the standard Nucleotide
        @rtype: String
        """
        return nucleic1to3[self.single_letter]

    @property
    def nitrogenous_base(self):
        """
        Returns the Nucleotide base as a single letter code
        @rtype: String
        """
        return nitrogenous_bases[self.single_letter]

    #
    # BOOLEANS
    #
    @property
    def has_p(self):
        """
        Checks if the Nucleotide has Phosphate
        @rtype: Boolean
        """
        return self._p is not None

    @property
    def has_o3(self):
        """
        Checks if the Nucleotide has Sugar O3'
        @rtype: Boolean
        """
        return self._o3 is not None

    #
    # METHODS
    #
    def normalize(self):
        self._type = self.standard_type
        self._mode = 'ATOM'
        
    def remove_side_chain(self):
        # self._atoms = self._backbone_atoms
        self._sidechain_atoms = []
        self._sidechain_coordinates = None

    def _distance_p(self, nucleotide):
        """
        Calculates the p-p distance between two residues

        @type  nucleotide: {Nucleotide}
        @param nucleotide: Nucleotide to which we want to find the distance

        @rtype: {Contact}
        @rtype: None if one of them does not have p
        """
        if not self.has_p or not nucleotide.has_p:
            return (None, None, -1)

        return (self.p, nucleotide.p, self.p.distance(atom = nucleotide.p))

    #
    # OVERRIDE PARENT'S METHODS
    #
    def add_atom(self, atom): 
        """
        Adds a new {NucleotideAtom} to the {Nucleotide}
        This includes:
            filling the all_coordinates attribute (parent)
            assign P, BB and backbone
        @type  atom: {AminoAtom}
        @param atom: New {AminoAtom} added to the {AminoAcid}
        """
        super(ResidueOfNucleotide, self).add_atom(atom = atom)

        if atom.is_Phosphate:
            self._p = atom

        if atom.is_SugarOxygen3:
            self._o3 = atom

    def distance(self, residue, dist_type):
        """
        Calculates the distance between two residues

        @type  residue: {Nucleotide}
        @param residue: Nucleotide to which we want to find the distance

        @type  dist_type: String
        @param dist_type: Type of distance to evaluate
        @limit dist_type: Accepted: distance.types (p, min, geometric)

        @rtype: (AT1, AT2, DISTANCE)
        """
        if not isinstance(residue, ResidueOfNucleotide):
            if isinstance(residue, ResidueOfAminoAcid):
                return super(ResidueOfNucleotide, self).distance(residue = residue, dist_type = dist_type)
            else:
                return super(ResidueOfAminoAcid, self).distance(residue = residue, dist_type = dist_type)
        if dist_type.strip().lower() == 'p':
            return self._distance_p(nucleotide = residue)
        else:
            return super(ResidueOfNucleotide, self).distance(residue = residue, dist_type = dist_type)


    def follows(self, residue):
        return residue.is_followed(self)

    def is_followed(self, residue):
        p0, o30 = self.p,    self.o3
        p1, o31 = residue.p, residue.o3
        a, b, distance = self.distance(residue, "geometric")

        if p0 is not None and p1 is not None and o30 is not None and o31 is not None:
            return o30.distance(p1) < o31.distance(p0) and distance <= 7.0 # measured distance P-P is around 7.0
        else:
            return self._backbone_atoms[-1].distance(residue._backbone_atoms[0]) < residue._backbone_atoms[-1].distance(self._backbone_atoms[0]) and distance <= 7.0 # measured distance P-P is around 7.0
