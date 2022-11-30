# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>

.. class:: ResidueFrame
.. class:: ProteinResidueFrame
.. class:: NucleotideResidueFrame
"""

# Standard Libraries

# External Libraries
import numpy as np
import pandas as pd

# This Library
from SBI.core import core
from SBI.data import alphabet
from SBI.data import nucleic_main
from SBI.structure.Frame3D import Frame3D
from SBI.structure.atom import AtomSeries

__all__ = ['ResidueFrame', 'ProteinResidueFrame', 'NucleotideResidueFrame']


class ResidueFrame( Frame3D ):
    """A DataFrame container of single residue data.
    """
    _subtyp = 'residue_frame'

    def __init__( self, *args, **kw ):
        super(ResidueFrame, self).__init__(*args, **kw)

    # Attributes
    @property
    def number(self):
        """Residue count identifier.

        :return: :class:`int`
        """
        return self[self._current_seq].unique()[0]

    @number.setter
    def number( self, value ):
        """Residue count identifier.
        """
        self[self._current_seq] = [value] * len(self)

    @property
    def version(self):
        """Residue insertion code.

        :return: :class:`str`
        """
        return self.pdbx_PDB_ins_code.unique()[0]

    @property
    def type( self ):
        """Residue compound type.

        :return: :class:`str`
        """
        return self[self._current_comp].unique()[0]

    @type.setter
    def type( self, value ):
        """
        """
        self[self._current_comp] = [value] * len(self)

    @property
    def chain( self ):
        """Chain identifier.

        :returns: :class:`str`
        """
        return str(self._current_model[self._current_asym].unique()[0])

    @property
    def atom_count( self ):
        """Returns the number of atoms in the :class:`.ResidueFrame`.

        This is equivalent to use :func:`len`, but is included for consistency.

        :return: :class:`int`
        """
        return len(self)

    # @property
    # def secondary_structure( self ):
    #     if self._dssp is None:
    #         raise AttributeError("To call secondary structure DSSP needs to be calculated")
    #
    #     return self._dssp.secondary_structure
    #
    # @property
    # def accessibility( self ):
    #     if self._dssp is None:
    #         raise AttributeError("To call accessibility DSSP needs to be calculated")
    #
    #     return self._dssp.accessibility
    #
    # @property
    # def accessibility10( self ):
    #     if self._dssp is None:
    #         raise AttributeError("To call accessibility DSSP needs to be calculated")
    #
    #     return self._dssp.accessibility10
    #
    # @property
    # def accessibilitycoded( self ):
    #     if self._dssp is None:
    #         raise AttributeError("To call accessibility DSSP needs to be calculated")
    #
    #     return self._dssp.accesscode
    #
    # @property
    # def exposed( self ):
    #     if self._dssp is None:
    #         raise AttributeError("To call exposition DSSP needs to be calculated")
    #
    #     return self._dssp.exposed
    #
    # @property
    # def exposed_text( self ):
    #     """Encodes exposition as either ''
    #     """
    #     if self.exposed:
    #         return 'E'
    #     else:
    #         return 'B'

    # Methods
    def types( self, value ):
        """
        """
        self[['auth_comp_id', 'label_comp_id']] = [[value, value]] * len(self)
        return self

    def numbers( self, value ):
        """
        """
        self[['auth_seq_id', 'label_seq_id']] = [[value, value]] * len(self)
        return self

    def simple_residue( self, name, atoms ):
        """
        """
        residue = []
        for at in atoms:
            residue.append(AtomSeries().simple_atom(at, atoms[at]))
        a = ResidueFrame(pd.concat(residue, ignore_index=True, axis=1).T)
        a.type = name
        return a.__finalize__(a, 'inherit')

    def follows( self, residue ):
        """Evaluate if the current :class:`.ResidueFrame` follows another
        in labeled sequence order.

        :param residue: Compound with which to evaluate if it follows.
        :type residue: :class:`.ResidueFrame`

        :return: :class:`bool` - Are the two residues consecutive according to
            sequence naming
        """
        number0  = self.number
        version0 = self.version if self.version != ' ' else '@'
        number1  = residue.number
        version1 = residue.version if residue.version != ' ' else '@'
        if number0 == number1 + 1:
            return True
        if number0 == number1:
            if ord(version0) == ord(version1) + 1:
                return True
        return False

    def is_followed( self, residue ):
        """Evaluate if the current :class:`.ResidueFrame` is followed by another
        in labeled sequence order.

        Notice that this is *not* the inverse of the :meth:`~ResidueFrame.follows`
        method.

        :param residue: Compound with which to evaluate if it follows.
        :type residue: :class:`.ResidueFrame`

        :return: :class:`bool` - Are the two residues consecutive according to
            sequence naming
        """
        number0  = self.number
        version0 = self.version if self.version != ' ' else '@'
        number1  = residue.number
        version1 = residue.version if residue.version != ' ' else '@'
        if number0 == number1 - 1:
            return True
        if number0 == number1:
            if ord(version0) == ord(version1) - 1:
                return True
        return False

    def identifier_distance( self, residue ):
        """Evaluate the identifier distance with another :class:`.ResidueFrame`.

        :return: :class:`int`
        """
        number0  = self.number
        version0 = self.version if self.version != ' ' else '@'
        number1  = residue.number
        version1 = residue.version if residue.version != ' ' else '@'

        if number0 != number1:
            return number1 - number0
        else:
            return abs(ord(version1) - ord(version0))

    def B_factor( self, values, inplace=False ):
        """Substitue the B-factor values for the provided ones.

        This is useful as coloring by B-factor is allowed in multiple protein representation
        softwares and allows plotting through other data of interest.

        :param values: Values to replace the current B-factor. If a single value is provided,
            all atoms are given the same B-factor. If a number of values equal to the atoms in
            the residue are provided, each atom gets its own value.
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
            else:
                raise AttributeError('Number of values cannot be fitted to any definition')
        return self._inplace(df, inplace)

    # Private Methods
    def _inheritance( self, other ):
        self.id = other.id
        if len(self['occupancy'].unique()) != 1:
            if core.get_option('structure', 'occupancy') == 'max':
                idx = self.groupby(self._current_atom)['occupancy'].transform(max) == self.occupancy
            else:
                idx = self.groupby(self._current_atom)['occupancy'].transform(min) == self.occupancy
            return self[idx]
        return self


class ProteinResidueFrame( ResidueFrame ):
    _subtyp = 'protein_residue_frame'

    def __init__( self, *args, **kw ):
        super(ProteinResidueFrame, self).__init__(*args, **kw)

    # Attributes
    @property
    def ca( self ):
        """Returns the C-alpha of the :class:`.ProteinResidueFrame`

        :return: :class:`.AtomSeries`
        """
        return self[self[self._current_atom] == 'CA'].squeeze()

    @property
    def cb( self ):
        """Returns the C-beta of the :class:`.ProteinResidueFrame`

        :return: :class:`.AtomSeries`
        """
        return self[self[self._current_atom] == 'CB'].squeeze()

    @property
    def n( self ):
        """Returns the backbone N of the :class:`.ProteinResidueFrame`

        :return: :class:`.AtomSeries`
        """
        return self[self[self._current_atom] == 'N'].squeeze()

    @property
    def c( self ):
        """Returns the backbone C of the :class:`.ProteinResidueFrame`

        :return: :class:`.AtomSeries`
        """
        return self[self[self._current_atom] == 'C'].squeeze()

    @property
    def o( self ):
        """Returns the backbone O of the :class:`.ProteinResidueFrame`

        :return: :class:`.AtomSeries`
        """
        return self[self[self._current_atom] == 'O'].squeeze()

    @property
    def single_letter( self ):
        """Amino acid identifier in 1 letter code.

        :return: :class:`str`
        """
        return alphabet.aminoacids3to1(self.type)

    # Methods
    def follows( self, residue ):
        """Evaluate if a residue follows another according to their fisical distances.

        :param residue: Residue with which to evaluate if it follows.
        :type residue: :class:`.ProteinResidueFrame`

        :return: :class:`bool`
        """
        return residue.is_followed(self)

    def is_followed( self, residue ):
        """Evaluate if a residue is followed by another according to their fisical distances.

        :param residue: Residue with which to evaluate if it follows.
        :type residue: :class:`.ProteinResidueFrame`

        :return: :class:`bool`
        """
        from ..atom import AtomSeries

        c0, ca0, n0 = self.c, self.ca, self.n
        c1, ca1, n1 = residue.c, residue.ca, residue.n

        if isinstance(c0, AtomSeries) and isinstance(n1, AtomSeries):
            return c0.distance(n1) <= 1.5  # measured distance is 1.3
        if isinstance(ca0, AtomSeries) and isinstance(ca1, AtomSeries):
            return ca0.distance(ca1) <= 4  # measured distance around 3.8
        if isinstance(c0, AtomSeries) and isinstance(c1, AtomSeries):
            return c0.distance(c1) <= 4
        if isinstance(n0, AtomSeries) and isinstance(n1, AtomSeries):
            return n0.distance(n1) <= 4
        return False


class NucleotideResidueFrame( ResidueFrame ):
    _subtyp = 'nucleotide_residue_frame'

    def __init__( self, *args, **kw ):
        super(NucleotideResidueFrame, self).__init__(*args, **kw)
