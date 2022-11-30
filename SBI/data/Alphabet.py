# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>
"""
# Standard Libraries
import copy

# External Libraries
import pandas as pd

# This Library
import SBI.core as core

__all__ = ['alphabet']


class Alphabet( pd.DataFrame ):

    _aminoacids_main3    = []
    _aminoacids_extr3    = []
    _aminoacids_main1    = []
    _aminoacids_3to1     = {}
    _aminoacids_1to3     = {}
    _protein_backbone    = ['CA', 'N', 'O', 'C']

    _nucleotide_main     = []
    _nucleotide_backbone = ['P',
                            'O1P', 'O2P', 'O3P', "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                            'OP1', 'OP2', 'OP3', "O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*"]

    def __init__( self, *args, **kwargs ):
        super(Alphabet, self).__init__(*args, **kwargs)
        if self.shape[0] > 0:
            self._init_values()

    @property
    def protein_backbone( self ):
        return self._protein_backbone

    @property
    def aminoacids_main3( self ):
        return copy.copy(self._aminoacids_main3)

    @property
    def aminoacids_main1( self ):
        return copy.copy(self._aminoacids_main1)

    @property
    def nucleotide_main( self ):
        return copy.copy(self._nucleotide_main)

    @property
    def nucleotide_backbone( self ):
        return self._nucleotide_backbone

    def is_aminoacid( self, code3 ):
        if code3 in self._aminoacids_main3 or code3 in self._aminoacids_extr3:
            return True
        data = self[(self['three_letter_code'] == code3) &
                    (self['type'].str.contains('PEPTIDE LINKING'))
                    ]
        if data.shape[0] == 0:
            return False
        elif data.shape[0] == 1:
            data = data.iloc[0]
        if data['one_letter_code'] != ' ' or data['mon_nstd_parent_comp_id'] != ' ':
            self._aminoacids_extr3.append(code3)
            return True
        else:
            return False

    def is_nucleotide( self, code3 ):
        # @TODO: Improve nucleotide detection of non-main types
        # @BODY: Do a similar approach as in the aminoacid.
        if code3 in self._nucleotide_main:
            return True
        return False

    def aminoacids3to1( self, code3 ):
        try:
            return self._aminoacids_3to1[code3]
        except KeyError:
            data = self[(self['three_letter_code'] == code3) &
                        (self['type'].str.contains('PEPTIDE LINKING'))
                        ]
            if data.shape[0] == 0:
                if core.get_option('data', 'strict'):
                    raise Unknown3LetterCodeError('Unknown code {}'.format(code3))
                else:
                    self._aminoacids_3to1[code3] = 'X'
            elif data.shape[0] > 1:
                raise Ambiguous3LetterCodeError('Multiple hits for {}'.format(code3))
            else:
                data = data.squeeze()
                if data['one_letter_code'] != ' ':
                    if len(data['one_letter_code']) == 1:
                        self._aminoacids_3to1[code3] = data['one_letter_code']
                    else:
                        self._aminoacids_3to1[code3] = '[' + data['one_letter_code'] + ']'
                elif data['mon_nstd_parent_comp_id'] != ' ':
                    self._aminoacids_3to1[code3] = self.aminoacids3to1(data['mon_nstd_parent_comp_id'])
                else:
                    if core.get_option('data', 'strict'):
                        raise Unknown3LetterCodeError('Unknown code {}'.format(code3))
                    else:
                        self._aminoacids_3to1[code3] = 'X'
            return self._aminoacids_3to1[code3]

    def aminoacids1to3( self, code1 ):
        try:
            return self._aminoacids_1to3[code1]
        except KeyError:
            data = self[(self['one_letter_code'] == code1) &
                        (self['type'].str.contains('PEPTIDE LINKING'))
                        ]
            if data.shape[0] == 0:
                if core.get_option('data', 'strict'):
                    raise Unknown1LetterCodeError('Unknown code {}'.format(code1))
                else:
                    self._aminoacids_1to3[code1] = 'UNK'
            else:
                self._aminoacids_1to3[code1] = data.iloc[0]['three_letter_code']
            return self._aminoacids_1to3[code1]

    def _init_values( self ):
        aas = self[((self['type'].str.contains('L-PEPTIDE LINKING'))  | (self['type'] == 'PEPTIDE LINKING')) &
                   (self['pdbx_initial_date'] == '1999-07-08') &
                   (self['mon_nstd_parent_comp_id'] == ' ') &
                   (~self['one_letter_code'].isin(['X', ' ']))
                   ]
        self._aminoacids_main1 = list(aas['one_letter_code'])
        self._aminoacids_main3 = list(aas['three_letter_code'])
        self._aminoacids_3to1 = pd.Series(aas.one_letter_code.values, index=aas.three_letter_code).to_dict()
        self._aminoacids_1to3 = pd.Series(aas.three_letter_code.values, index=aas.one_letter_code).to_dict()
        self._aminoacids_1to3['X'] = 'UNK'

        nns = self[(self['type'].str.contains('NA LINKING')) & (self['mon_nstd_parent_comp_id'] == ' ') &
                   (self['pdbx_type'] == 'ATOMN') & (self['name'].str.isupper()) &
                   ((self['three_letter_code'].str.len() == 1) | ((self['three_letter_code'].str.len() == 2) &
                    (self['three_letter_code'].str.startswith('D')))) ]

        self._nucleotide_main = list(nns['three_letter_code'])


class Unknown3LetterCodeError( KeyError ):
    """Raised when a new 3to1 code translation is requested on an unknown 3 letter code.
    """


class Unknown1LetterCodeError( KeyError ):
    """Raised when a new 1to3 code translation is requested on an unknown 1 letter code.
    """


class Ambiguous3LetterCodeError( KeyError ):
    """Raised when a 3 letter code generates more than one hit.
    """


alphabet = Alphabet(pd.read_csv(core.get_option('data', 'alphabet')))
"""
:data:`.alphabet` is the main container to hold static data necessary to interpret PDB codes for compound.

The data contained in this variable is loaded at runtime and generated from the
`PDBeChem <ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/>`_ database.
"""
