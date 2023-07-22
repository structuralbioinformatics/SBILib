# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>

.. class:: AtomSerie
.. class:: VirtualAtom
"""

# Standard Libraries

# External Libraries
import numpy as np
from numpy.linalg import norm
import pandas as pd

# This Library
from SBILib.structure.io import mmNone
from SBILib.structure.Frame3D import mandatory_fields, frame3D_classify

__all__ = ['AtomSeries', 'VirtualAtom']


class AtomSeries( pd.Series ):
    """A single Atom data.
    """
    _subtyp = 'atom_series'
    _atmtyp = 'ATOM'

    def __init__( self, *args, **kwargs ):
        fields = kwargs.pop('empty_fields', mandatory_fields)
        super(AtomSeries, self).__init__(*args, **kwargs)
        if self.empty:
            data = {}
            for f in fields:
                data.setdefault(f, np.nan)
            super(AtomSeries, self).__init__(data)

        if not set(mandatory_fields).issubset(self.index):
            raise AtomSeriesCreationError('AtomSeries requires a given set of index to be present: {}'.format(','.join(mandatory_fields)))

    # Attributes
    @property
    def coordinates( self ):
        """Spatial position of the :class:`.AtomSeries` as coordinates ``X``, ``Y``, ``Z``.

        :return: :class:`~numpy.ndarray`
        """
        return self[['Cartn_x', 'Cartn_y', 'Cartn_z']].values

    # Boolean Checks
    @property
    def is_empty( self ):
        return (self.isna()).all()

    # Methods
    def distance( self, atom=None ):
        """
        Euclidean distance between two atoms
        """
        if atom is None:
            return norm(self.coordinates - np.zeros(3))
        else:
            return norm(self.coordinates - atom.coordinates)

    def simple_atom( self, atom_id, coordinates ):
        """
        """
        data = mmNone().read()['empty']['_atom_site']

        for i, xyz in enumerate(['x', 'y', 'z']):
            data['Cartn_{}'.format(xyz)] = coordinates[i]
        data['group_PDB'] = self._atmtyp
        for x in ['id', 'label_entity_id', 'label_seq_id', 'auth_seq_id', 'pdbx_PDB_model_num']:
            data[x] = 1
        for x in ['auth_comp_id', 'label_comp_id']:
            data[x] = 'UNK'
        for x in ['auth_asym_id', 'label_asym_id']:
            data[x] = 'A'
        for x in ['auth_atom_id', 'label_atom_id']:
            data[x] = atom_id
        for x in ['occupancy', 'B_iso_or_equiv']:
            data[x] = 1.0
        data['type_symbol'] = data['auth_atom_id'][0] if data['auth_atom_id'] != data['auth_comp_id'] else data['auth_atom_id']

        for k in data:
            if isinstance(data[k], list):
                data[k] = ' '

        return AtomSeries(data)

    # Pandas Magic Methods
    @property
    def _constructor( self ):
        # return AtomSeries
        def f(*args, **kwargs):
            try:
                return frame3D_classify(AtomSeries(*args, **kwargs))
            except AtomSeriesCreationError:
                return frame3D_classify(pd.Series(*args, **kwargs))
        return f

    @property
    def _constructor_expanddim( self ):
        from SBILib.structure import Frame3D

        def f(*args, **kwargs):
            return frame3D_classify(Frame3D(*args, **kwargs))
        return f


class VirtualAtom( AtomSeries ):
    """A single Atom data.
    """
    _subtyp = 'virtual_atom_series'
    _atmtyp = 'VIRATOM'

    def __init__( self, *args, **kw ):
        super(VirtualAtom, self).__init__(*args, **kw)

    @property
    def _constructor( self ):
        return VirtualAtom


class AtomSeriesCreationError( Exception ):
    """Raises an error when a :class:`.AtomSeries` cannot be created.
    """