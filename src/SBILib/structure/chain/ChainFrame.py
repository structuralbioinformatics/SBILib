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
import os
import re
from typing import Optional, Union, TypeVar
from tempfile import gettempdir
from itertools import groupby

# External Libraries
import pandas as pd
import numpy as np

# This Library
from ..Frame3D import Frame3D
from SBILib.core import core
from SBILib.data import alphabet

__all__ = ['ChainFrame', 'ProteinChainFrame', 'NucleotideChainFrame']

CF = TypeVar('CF', bound='ChainFrame')


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

        :returns: :class:`str`
        """
        warnings.warn(".globalID will be deprecated. use .id instead", FutureWarning, stacklevel=2)
        return self.id

    @property
    def id( self ):
        """Structure identifier.

        :returns: :class:`str`
        """
        return "{0}_{1}".format(self.pdb, self.chain)

    @property
    def pdb( self ):
        """
        PDB identifier.

        :returns: :class:`str`
        """
        return self._id

    @property
    def chain( self ):
        """Chain identifier.

        :returns: :class:`str`
        """
        return str(self._current_model[self._current_asym].unique()[0])

    @property
    def residue_count( self ):
        """Number of residues in the :class:`.ChainFrame` (of the first *-main-* entity).

        .. warning::
            This substitutes :func:`len` in pre-pandas versions of the code,
            as overwritting :func:`len` is problematic.

        :return: :class:`int`
        """
        try:
            entity = self._current_model.iloc[0]['label_entity_id']
            entity = self._current_model[self._current_model['label_entity_id'] == entity]
            return len(entity.groupby(self._residue_columns, sort=False))
        except IndexError:
            # Produced for no .iloc[0], which means there is no data here.
            return 0

    @property
    def compounds( self ):
        """List the compounds in the :class:`.ChainFrame`.

        :yields: :class:`.ResidueFrame`
        """
        for _ in self.groupby(self._residue_columns + [self._current_comp, self._current_asym], sort=False):
            yield _[1]

    @property
    def first_compound( self ):
        """Obtain the first compound of the :class:`.ChainFrame`.

        :return: :class:`.ResidueFrame`
        """
        return list(self.compounds)[0]

    @property
    def last_compound( self ):
        """Obtain the first compound of the :class:`.ChainFrame`.

        :return: :class:`.ResidueFrame`
        """
        return list(self.compounds)[-1]

    # Boolean
    @property
    def is_ordered( self ):
        """Checks if the residue numbering of the chain is ordered; considering insertion codes too.

        Ignores waters and ligands.

        .. note::
            This functionality depends on the :ref:`global configuration options <configuration>` ``structure.source``.

        :return: :class:`bool`
        """
        if core.get_option('structure', 'source') == 'label':
            return True  # label NEVER is unoerdered
        else:
            if not self.has_insertion_codes:
                df = self.dehydrate(False).remove_heteroatoms(False)
                return pd.Index(df[df._current_seq].loc[df[df._current_seq].shift() != df[df._current_seq]]).is_monotonic
            return False  # Insertion codes will always mean some type of "disorder" in the sequence.

    # Methods
    def renumber( self, start=1, inplace=True ):
        """Renumber the residues of :class:'.ChainFrame'.

        .. note::
            This functionality depends on the :ref:`global configuration options <configuration>` ``structure.source``.
        """
        src = core.get_option('structure', 'source')

        df = self

        gf = [list(group) for key, group in groupby(df[df._residue_columns].values.tolist())]
        max_value = len(gf)
        atoms = np.asarray([(len(x)) for x in gf])

        numbers = [[x, ] * atoms[i] for i, x in enumerate(range(start, start + max_value))]

        df[df._current_seq] = np.asarray(numbers).flatten()
        if src == 'auth' and 'pdbx_PDB_ins_code' in df.columns:
            df['pdbx_PDB_ins_code'] = np.asarray(['', ] * df.shape[0])

        df['id'] = list(range(1, df.shape[0] + 1))

        return self._inplace(df, inplace, False)

    def B_factor( self, values, inplace=False ):
        """Substitue the B-factor values for the provided ones.

        This is useful as coloring by B-factor is allowed in multiple protein representation
        softwares and allows plotting through other data of interest.

        :param values: Values to replace the current B-factor. If a single value is provided,
            all residues are given the same B-factor. If a number of values equal to the number
            of residues is provided, each residue gets the same value for all atoms. If a number
            of values equal to the atoms in the residue are provided, each atom gets its own value.
        :type values: Union[:class:`int`, :func:`list`]

        :return: :class:`.ResidueFrame`

        :raise:
            :AttributeError: If the number of values given do not match any of the accepted
                conditions.
        """
        df = self
        if isinstance(values, (int, float)):
            df['B_iso_or_equiv'] = [values, ] * self.shape[0]
        elif isinstance(values, (list, np.ndarray)):
            if len(values) == 1:
                df['B_iso_or_equiv'] = [values[0], ] * self.shape[0]
            elif len(values) == self.shape[0]:
                df['B_iso_or_equiv'] = values
            elif len(values) == self.compounds:
                df2 = []
                for i, aa in enumerate(df.compounds):
                    df2.append(aa.B_factor(values[i], inplace=True))
                df = pd.concat(df2)
            else:
                raise AttributeError('Number of values cannot be fitted to any definition')
        return self._inplace(df, inplace)

    def pop_out( self, key: Optional[Union[int, str]] = 0 ) -> CF:
        """Delete the requested position by order from the :class:`.ChainFrame`.

        .. note::
            This functionality depends on the :ref:`global configuration option <configuration>` ``structure.model``.
            This functionality depends on the :ref:`global configuration option <configuration>` ``structure.source``.

        :param key: Array like access of the target position.

        :return: The new :class:`.ChainFrame` without the requested position.

        :raise:
            :IndexError: If the position requested cannot be found (bigger than ``count(self)``).
        """
        df = self._current_model
        cl = self._residue_columns
        dd = [x for x, _ in df.groupby(cl, sort=False)][key]
        if len(cl) == 1:
            idx = df[df[cl[0]] == dd[0]].index
        else:
            idx = df[(df[cl[0]] == dd[0]) & (df[cl[1]] == dd[1])].index
        return df.drop(idx)

    # Private Pandas Methods
    def _inheritance( self, other ):
        self._id = other._id
        return self


class ProteinChainFrame( ChainFrame ):
    _subtyp = 'protein_chain_frame'

    def __init__( self, *args, **kw ):
        super(ProteinChainFrame, self).__init__(*args, **kw)
        self._aalimit = None

    # Attributes
    @property
    def aminoacids( self ):
        """Iterate over the compounds in the :class:`.ChainFrame`.

        :yields: :class:`.ProteinResidueFrame`
        """
        df = self.remove_heteroatoms(False).dehydrate(False)
        df.set_index(df._residue_columns, drop=False, inplace=True)
        for _ in df.groupby(level=[0, 1], axis=0, sort=False):
            yield _[1]

    @property
    def first_aminoacid( self ):
        """Obtain the first amino acid.

        :return: :class:`.ProteinResidueFrame`
        """
        return self.first_compound

    @property
    def last_aminoacid(self):
        """Obtain the first amino acid.

        :return: :class:`.ProteinResidueFrame`
        """
        for residue in reversed(list(self.aminoacids)):
            return residue

    @property
    def protein_sequence(self):
        """Returns the sequence of the crystalized residues.

        It considers all the residues with density in the 3D structure and has no
        consideration for whether or not there might be gaps in the structure.

        :return: :class:`str`

        .. seealso::
            :meth:`.ProteinChainFrame.gapped_protein_sequence`
            :meth:`.ProteinChainFrame.full_protein_sequence`
        """
        df = self.remove_heteroatoms(False).dehydrate(False)
        df = df.set_index(df._residue_columns, drop=False)
        lvl = list(range(len(df._residue_columns)))
        seq = df.groupby(level=lvl, axis=0, sort=False).first()[df._current_comp].values
        return "".join([alphabet.aminoacids3to1(code3) for code3 in seq])

    @property
    def gapped_protein_sequence(self):
        """Returns the sequence of the crystalized residues filling with 'X' missed densities.

        Missing densities are calculated by evaluating for each two residues ir they follow each other or not.

        :return: :class:`str`

        .. seealso::
            :meth:`.ProteinChainFrame.protein_sequence`
            :meth:`.ProteinChainFrame.full_protein_sequence`
        """
        from ..protein import Sequencer
        return Sequencer._sequencer(self, 'seq')

    @property
    def full_protein_sequence( self ):
        """Returns the sequence of the crystalized residues filling with 'X' missed densities and starting
        non-crystalized residues.

        Takes into account the starting number of the chain and considers that all previous residues are missing.

        :return: :class:`str`

        .. seealso::
            :meth:`.ProteinChainFrame.protein_sequence`
            :meth:`.ProteinChainFrame.gapped_protein_sequence`
        """
        if self.first_aminoacid.number <= 1:
            return self.gapped_protein_sequence

        seq = ''
        for x in range(self.first_aminoacid.number - 1):
            seq += 'x'
        return seq + self.gapped_protein_sequence

    # Methods
    def calculate_dssp( self, tmppdb=None, tmpdssp=None, cleanfiles=True, minimize=False, simplify='BGITS' ):
        """Generate the [DSSP](https://github.com/cmbi/xssp) assignation for the :class:`.ProteinChainFrame`.

        .. warning::
            If the object already contains DSSP data,

        :param str tmppdb: Name for the PDB file. If :data:`None` provided, a temporary file is created.
        :param str tmpdssp: Name for the DSSP file. If :data:`None` provided, a temporary file is created.
        :param bool cleanfiles: If :data:`True`, PDB and DSSP files are deleted after execution.
        :param bool minimize: If :data:`True`, ignore data from DSSP output that is not secondary structure
            assignation or accessibility.
        :param str simplify: Secondary structure codes to simplify. ``H``, ``E`` and ``L`` are ignored if provided,
            as they cannot be simplified. The default represents the maximum possible simplification.

        :return: :class:`.DSSPFrame`
        """
        from SBILib.external import DSSPExe

        tmppdb = tmppdb if tmppdb is not None else os.path.join(gettempdir(), '{}.pdb'.format(os.getpid()))
        tmpdssp = tmpdssp if tmpdssp is not None else os.path.join(gettempdir(), '{}.dssp'.format(os.getpid()))
        if not os.path.isfile(tmppdb):
            self.dehydrate(False).remove_heteroatoms(False).write(tmppdb, format='pdb', force=False, clean=False)
        dssp = DSSPExe(pdb=tmppdb, dssp=tmpdssp, cleanpdb=cleanfiles, cleandssp=cleanfiles, minimize=minimize).dsspdata.simplify(simplify)

        return dssp.copy()

    def split_discontinuous( self ):
        """Find the discontinuous pieces of :class:`.ProteinChain` and separate them.

        :return: :func:`list` of :class:`.Frame3D`
        """
        gps = re.split(r'x+', self.gapped_protein_sequence)
        ranges = [[1, len(gps[0])]]
        base = 0
        for i in range(1, len(gps)):
            base += len(gps[i - 1])
            ranges.append([base + 1, base + len(gps[i])])
        with core.on_option_value('structure', 'source', 'label'):
            for i, r in enumerate(ranges):
                ranges[i] = self['Residue:{0}-{1}'.format(r[0], r[1])]
        return ranges


class NucleotideChainFrame( ChainFrame ):
    _subtyp = 'nucleotide_chain_frame'

    def __init__( self, *args, **kw ):
        super(NucleotideChainFrame, self).__init__(*args, **kw)
