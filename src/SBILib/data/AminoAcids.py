# -*-
#
# @author: jaumebonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.github.io
#
# @date:   2015-06-18 10:01:28
# @lab:    LPDI/EPFL
#
# @last modified by:   jaumebonet
# @last modified time: 2016-02-15 11:37:34
#
# @source: atom_nom.tbl of Eldon Ulrich (elu@nmrfam.wisc.edu)
#
# -*-
from pynion import Singleton
from AtomVariants import *
from .SetDict import SetDict


class AminoAcidsTable(object):
    """It knows about Amino Acid properties"""
    __metaclass__ = Singleton

    def __init__(self):
        self._atom_composition = {
            'A': [ C, H, CA, O, N, CB, HA, HB1, HB2b, HB3b],
            'C': [ C, H, CA, O, N, CB, HA, HB3, HB2, HG, SG],
            'E': [ C, H, CA, O, N, CB, HA, HG2, CD, CG, OE1, HB3, HB2, HG3, HE2b, OE2],
            'D': [ C, H, CA, O, N, CB, HA, CG,  OD1, HB3, HB2,
                    AtomVariants(['HD2', '', '', 'HD2', '', '', '', 'HD2']),
                    AtomVariants(['OD2', 'OD2', '', 'OD2', 'OD2', 'OD2', 'OD2', 'OD2'])],
            'G': [ C, H, CA, O, N,
                    AtomVariants(['HA2', '1HA', 'HA1', 'HA1', 'HA2', 'HA2', '', 'HA1']),
                    AtomVariants(['HA3', '2HA', 'HA2', 'HA2', 'HA1', 'HA1', '', 'HA2'])],
            'F': [ C, H, CA, O, N, CB, HA, HE2, HD2, HD1, CG, HE1, CE1, CZ, CE2, CD1, CD2,
                    HB3, HB2,
                    AtomVariants(['HZ', 'HZ', 'HZ', 'HZ', 'HZ', 'HZ', '', 'HZ'])],
            'I': [ C, H, CA, O, N, CB, HA, HG22, CG2, HG21, HD13, CD1, HB, HD12, HG23, HD11,
                    AtomVariants(['CG1', 'CG1', '', 'CG1', 'CG1', 'CG1', 'CG1', 'CG1']),
                    AtomVariants(['HG12', '1HG1', 'HG11', 'HG11', 'HG12', 'HG12', '', 'HG12']),
                    AtomVariants(['HG13', '2HG1', 'HG12', 'HG12', 'HG11', 'HG11', '', 'HG13'])],
            'H': [ C, H, CA, O, N, CB, HA, HD2, CG, HE1, CE1, CD2, NE2, HB3, HB2,
                    AtomVariants(['HE2', '', 'HNE2', 'HE2', '', '', '', 'HE2']),
                    AtomVariants(['HD1', 'HD1', 'HND1', 'HD1', 'HD1', 'HD1', '', 'HD1']),
                    AtomVariants(['ND1', 'ND1', '', 'ND1', 'ND1', 'ND1', 'ND1', 'ND1'])],
            'K': [ C, H, CA, O, N, CB, HA, HG2, HG3, CG, CE, CD, HD3, HD2b, HB3, HB2,
                    AtomVariants(['HE2', '1HE', 'HE1', 'HE1', 'HE2', 'HE2', '', 'HE2']),
                    AtomVariants(['HE3', '2HE', 'HE2', 'HE2', 'HE1', 'HE1', '', 'HE3']),
                    AtomVariants(['NZ', 'NZ', '', 'NZ', 'NZ', 'NZ', 'NZ', 'NZ']),
                    AtomVariants(['HZ1', '1HZ', 'HNZ1', 'HZ1', 'HZ1', 'HZ1', '', 'HZ1']),
                    AtomVariants(['HZ3', '3HZ', 'HNZ3', 'HZ3', 'HZ3', 'HZ3', '', 'HZ3']),
                    AtomVariants(['HZ2', '2HZ', 'HNZ2', 'HZ2', 'HZ2', 'HZ2', '', 'HZ2'])],
            'M': [ C, H, CA, O, N, CB, HA, HG2, HG3, CE, CG, HB3, HB2,
                    AtomVariants(['HE1', '1HE', 'HE1', 'HE1', 'HE1', 'HE1', '', 'HE1']),
                    AtomVariants(['HE2', '2HE', 'HE2', 'HE2', 'HE2', 'HE2', '', 'HE2']),
                    AtomVariants(['HE3', '3HE', 'HE3', 'HE3', 'HE3', 'HE3', '', 'HE3']),
                    AtomVariants(['SD', 'SD', '', 'SD', 'SD', 'SD', 'SD', 'SD'])],
            'L': [ C, H, CA, O, N, CB, HA, CG, CD1, CD2, HD13, HD12, HD11, HB3, HB2,
                    AtomVariants(['HD22', '2HD2', 'HD22', 'HD22', 'HD22', 'HD22', '', 'HD22']),
                    AtomVariants(['HD23', '3HD2', 'HD23', 'HD23', 'HD23', 'HD23', '', 'HD23']),
                    AtomVariants(['HD21', '1HD2', 'HD21', 'HD21', 'HD21', 'HD21', '', 'HD21']),
                    AtomVariants(['HG', 'HG', 'HG', 'HG', 'HG', 'HG', '', 'HG'])],
            'N': [ C, H, CA, O, N, CB, HA, CG, OD1, HB3, HB2,
                    AtomVariants(['ND2', 'ND2', '', 'ND2', 'ND2', 'ND2', 'ND2', 'ND2']),
                    AtomVariants(['HD22', '2HD2', 'HN22', 'HD22', 'HD22', 'HD22', '', 'HD22']),
                    AtomVariants(['HD21', '1HD2', 'HN21', 'HD21', 'HD21', 'HD21', '', 'HD21'])],
            'Q': [ C, H, CA, O, N, CB, HA, HG2, HG3, CD, CG, NE2, OE1, HB3, HB2,
                    AtomVariants(['HE22', '2HE2', 'HN22', 'HE22', 'HE22', 'HE22', '', 'HE22']),
                    AtomVariants(['HE21', '1HE2', 'HN21', 'HE21', 'HE21', 'HE21', '', 'HE21'])],
            'P': [ C, CA, O, N, CB, HA, HG2, CD, CG, HB3, HB2,
                    AtomVariants(['HD3', '2HD', 'HD1', 'HD2', 'HD1', 'HD1', '', 'HD3']),
                    AtomVariants(['HD2', '1HD', 'HD2', 'HD1', 'HD2', 'HD2', '', 'HD2']),
                    AtomVariants(['H2', 'H2', '', 'HN2', 'HT2', '', '', '']),
                    AtomVariants(['H3', 'H1', '', 'HN1', 'HT1', '', '', ''])],
            'S': [ C, H, CA, O, N, CB, HA, HB3, HB2,
                    AtomVariants(['OG', 'OG', '', 'OG', 'OG', 'OG', 'OG', 'OG']),
                    AtomVariants(['HG', 'HG', 'HOG', 'HG', 'HG', 'HG', '', 'HG'])],
            'R': [ C, H, CA, O, N, CB, HA, HG2, HG3, CG, CD, CZ, HD3, HD2b, HB3, HB2,
                    AtomVariants(['NE', 'NE', '', 'NE', 'NE', 'NE', 'NE', 'NE']),
                    AtomVariants(['HE', 'HE', 'HNE', 'HE', 'HE', 'HE', '', 'HE']),
                    AtomVariants(['HH22', '2HH2', 'HN22', 'HH21', 'HH22', 'HH22', '', 'HH22']),
                    AtomVariants(['HH21', '1HH2', 'HN21', 'HH22', 'HH21', 'HH21', '', 'HH21']),
                    AtomVariants(['NH1', 'NH1', '', 'NH1', 'NH1', 'NH1', 'NH1', 'NH1']),
                    AtomVariants(['NH2', 'NH2', '', 'NH2', 'NH2', 'NH2', 'NH2', 'NH2']),
                    AtomVariants(['HH12', '2HH1', 'HN12', 'HH12', 'HH12', 'HH12', '', 'HH12']),
                    AtomVariants(['HH11', '1HH1', 'HN11', 'HH11', 'HH11', 'HH11', '', 'HH11'])],
            'T': [ C, H, CA, O, N, CB, HA, HG21, CG2, HB, HG23, HG22,
                    AtomVariants(['OG1', 'OG1', '', 'OG1', 'OG1', 'OG1', 'OG1', 'OG1']),
                    AtomVariants(['HG1', 'HG1', 'HOG1', 'HG1', 'HG1', 'HG1', '', 'HG1'])],
            'W': [ C, H, CA, O, N, CB, HA, CD1, CD2, CG, CE2, HD1, HB3, HB2,
                    AtomVariants(['HH2', 'HH2', 'HH2', 'HH2', 'HH2', 'HH2', '', 'HH2']),
                    AtomVariants(['CZ2', 'CZ2', '', 'CZ2', 'CZ2', 'CZ2', 'CZ2', 'CZ2']),
                    AtomVariants(['CZ3', 'CZ3', '', 'CZ3', 'CZ3', 'CZ3', 'CZ3', 'CZ3']),
                    AtomVariants(['HE1', 'HE1', 'HNE1', 'HE1', 'HE1', 'HE1', '', 'HE1']),
                    AtomVariants(['HE3', 'HE3', 'HE3', 'HE3', 'HE3', 'HE3', '', 'HE3']),
                    AtomVariants(['CH2', 'CH2', '', 'CH2', 'CH2', 'CH2', 'CH2', 'CH2']),
                    AtomVariants(['HZ3', 'HZ3', 'HZ3', 'HZ3', 'HZ3', 'HZ3', '', 'HZ3']),
                    AtomVariants(['HZ2', 'HZ2', 'HZ2', 'HZ2', 'HZ2', 'HZ2', '', 'HZ2']),
                    AtomVariants(['CE3', 'CE3', '', 'CE3', 'CE3', 'CE3', 'CE3', 'CE3']),
                    AtomVariants(['NE1', 'NE1', '', 'NE1', 'NE1', 'NE1', 'NE1', 'NE1'])],
            'V': [ C, H, CA, O, N, CB, HA, HB,
                    AtomVariants(['HG22', '2HG2', 'HG22', 'HG12', 'HG22', 'HG22', '', 'HG22']),
                    AtomVariants(['HG21', '1HG2', 'HG21', 'HG11', 'HG21', 'HG21', '', 'HG21']),
                    AtomVariants(['CG1', 'CG1', '', 'CG2', 'CG1', 'CG1', 'CG1', 'CG1']),
                    AtomVariants(['HG11', '1HG1', 'HG11', 'HG21', 'HG11', 'HG11', '', 'HG11']),
                    AtomVariants(['HG12', '2HG1', 'HG12', 'HG22', 'HG12', 'HG12', '', 'HG12']),
                    AtomVariants(['HG13', '3HG1', 'HG13', 'HG23', 'HG13', 'HG13', '', 'HG13']),
                    AtomVariants(['CG2', 'CG2', '', 'CG1', 'CG2', 'CG2', 'CG2', 'CG2']),
                    AtomVariants(['HG23', '3HG2', 'HG23', 'HG13', 'HG23', 'HG23', '', 'HG23'])],
            'Y': [ C, H, CA, O, N, CB, HA, HE2, HD2, HD1, CG, HE1, CE1, CZ, CD1, CD2, CE2, HB3, HB2,
                    AtomVariants(['OH', 'OH', '', 'OH', 'OH', 'OH', 'OH', 'OH']),
                    AtomVariants(['HH', 'HH', 'HOH', 'HH', 'HH', 'HH', '', 'HH'])],
            'X': [
                AtomVariants(['H2', '2H', 'HN2', 'HN2', 'HT2', '', '', '']),
                AtomVariants(['H3', '3H', 'HN3', 'HN3', 'HT3', '', '', '']),
                AtomVariants(['H1', '1H', 'HN1', 'HN1', 'HT1', 'HNCAP', '', '']),
                AtomVariants(['O\'', 'O', '', 'O', 'OT1', 'O', 'O', '']),
                AtomVariants(['H\'\'', '', '', '', '', 'HOCAP', '', '']),
                AtomVariants(['O\'\'', 'OXT', '', 'OXT', 'OT2', 'OXT', 'OXT', ''])]}

    def rename_atom(self, residue_name, atom_name, current_format, new_format):
        AtomVariants.known_format(current_format)
        AtomVariants.known_format(new_format)
        for av in self._atom_composition[residue_name]:
            if av.get_name(current_format) == atom_name:
                return av.get_name(new_format)
        return atom_name

AminoAcids = AminoAcidsTable()
