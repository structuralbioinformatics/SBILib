# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>

.. class:: ChainFrame
.. class:: ProteinChainFrame
.. class:: NucleotideChainFrame
"""

# Standard Libraries
import warnings

# External Libraries
import pandas as pd

# This Library
from SBILib.data import aminoacids_main3, nucleic_main
from ..Frame3D import Frame3D
from ..residue import ResidueFrame

__all__ = ['ChainFrame', 'ProteinChainFrame', 'NucleotideChainFrame']


class ChainFrame( Frame3D ):
    """A DataFrame container of single chain data. A protein fragment
    would also fall into this category.
    """
    _subtyp = 'chain_frame'

    def __init__( self, *args, **kw ):
        super(ChainFrame, self).__init__(*args, **kw)

    # Attributes
    @property
    def globalID( self ):
        """
        Full identifier of the Chain.
        By default, it would be ``<pdb-name>_<chain-id>``

        :returns: :py:class:`str`
        """
        return "{0}_{1}".format(self.id, self.chain)

    @property
    def pdb( self ):
        warnings.warn(".pdb will be deprecated. use .id instead", FutureWarning, stacklevel=2)
        return self.id

    @property
    def chain( self ):
        """
        Chain identifier.

        :returns: :py:class:`str`
        """
        return str(self._current_model[self._current_asym].unique()[0])

    # Magic Methods
    def __len__(self):
        """Number of residues in the Chain (of the first -main- entity of the Chain).
        """
        entity = self._current_model.iloc[0]["label_entity_id"]
        entity = self._current_model[self._current_model["label_entity_id"] == entity]
        return len(entity.groupby([self._current_seq, "pdbx_PDB_ins_code"]))

    # Private Pandas Methods
    def _inheritance( self, other ):
        self.id = other.id
        return self

    # Pandas Magic Methods
    @property
    def _constructor( self ):
        def f(*args, **kwargs):
            df = pd.DataFrame(*args, **kwargs)
            if len(df[self._current_seq].unique()) == 1:
                return ResidueFrame(*args, **kwargs).__finalize__(self, method='inherit')
            else:
                return ChainFrame(*args, **kwargs).__finalize__(self, method='inherit')
        return f

    def __finalize__(self, other, method=None, **kwargs):
        if method == 'inherit':
            if (self[self._current_comp].isin(aminoacids_main3) == True).any():
                return ProteinChainFrame(self).__finalize__(other, 'inherit')
            if (self[self._current_comp].isin(nucleic_main) == True).any():
                return NucleotideChainFrame(self).__finalize__(other, 'inherit')
            return self._inheritance(other)
        if method == 'concat':
            if len(self._current_model[self._current_asym].unique()) > 1:
                from SBILib.structure import PDBFrame
                return PDBFrame(self)
        return self


class ProteinChainFrame( ChainFrame ):
    _subtyp = 'protein_chain_frame'

    def __init__( self, *args, **kw ):
        super(ProteinChainFrame, self).__init__(*args, **kw)

    # Attributes
    @property
    def aminoacids( self ):
        """Iterate over the :class:`.ProteinResidueFrame` of the :class:`.ChainFrame`.

        :return: Iterator[:class:`.ProteinResidueFrame`]
        """
        pass

    # Pandas Magic Methods
    @property
    def _constructor( self ):
        def f(*args, **kwargs):
            df = pd.DataFrame(*args, **kwargs)
            if len(df[self._current_seq].unique()) == 1:
                return ResidueFrame(*args, **kwargs).__finalize__(self, method='inherit')
            else:
                return ProteinChainFrame(*args, **kwargs).__finalize__(self, method='inherit')
        return f

    def __finalize__(self, other, method=None, **kwargs):
        if method == 'inherit':
            return self._inheritance(other)
        else:
            return super(ProteinChainFrame, self).__finalize__(other, method, **kwargs)
        return self


class NucleotideChainFrame( ChainFrame ):
    _subtyp = 'nucleotide_chain_frame'

    def __init__( self, *args, **kw ):
        super(NucleotideChainFrame, self).__init__(*args, **kw)

    # Pandas Magic Methods
    @property
    def _constructor( self ):
        def f(*args, **kwargs):
            df = pd.DataFrame(*args, **kwargs)
            if len(df[self._current_seq].unique()) == 1:
                return ResidueFrame(*args, **kwargs).__finalize__(self, method='inherit')
            else:
                return NucleotideChainFrame(*args, **kwargs).__finalize__(self, method='inherit')
        return f

    def __finalize__(self, other, method=None, **kwargs):
        if method == 'inherit':
            return self._inheritance(other)
        else:
            return super(NucleotideChainFrame, self).__finalize__(other, method, **kwargs)
        return self
