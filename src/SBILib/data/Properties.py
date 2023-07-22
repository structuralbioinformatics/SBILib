# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>
"""
# Standard Libraries

# External Libraries
import pandas as pd

# This Library
import SBI.core as core

__all__ = ['properties']


class Properties( pd.DataFrame ):
    def get_surface( self, residue_type ):
        """Get the expected **surface** for a query amino acid.

        .. note::
            Depends on global configuration option ``data.surface`` to determine
            the source from which the surface data is obtained. There are 4 main
            data sources for surface information:

            * _tien2013t_ and _tien2013e_: Experimental (e) and theoretical (t) surface \
                data provided by Tien, M.Z. _et al._ (2013). [Maximum allowed solvent accessibilites \
                of residues in proteins.](https://doi.org/10.1371%2Fjournal.pone.0080635) **PLoS ONE**.
            * _miller87_: Surface data obtained from Miller, S. _et al._ (1987). [Interior and surface \
                of monomeric proteins.](https://doi.org/10.1016%2F0022-2836%2887%2990038-6) **JMB**.
            * _rose85_: Surface data obtained from Rose, G.D. _et al._ (1985). [Hydrophobicity of amino \
                acid residues in globular proteins.](https://doi.org/10.1126%2Fscience.4023714) **Science**.

        :param str residue_type: Amino acid type (1 or 3 letter code).

        :returns: :class:`.float`

        :raises:
            :IndexError: If no surface can be found for the specified amino acid code
                (i.e. if the code is not found)
        """
        surf_source = core.get_option('data', 'surface')
        residue_type = residue_type.strip()
        col = 'one_letter_code' if len(residue_type) == 1 else 'three_letter_code'
        return self[self[('residue', col)] == residue_type][('surface', surf_source)].values[0]


colnames = pd.MultiIndex.from_tuples([
    ('residue', 'one_letter_code'), ('residue', 'three_letter_code'),
    ('surface', 'tien2013t'), ('surface', 'tien2013e'), ('surface', 'miller87'), ('surface', 'rose85')
])
properties = Properties([
    ['A', 'ALA', 129.0, 121.0, 113.0, 118.1],
    ['R', 'ARG', 274.0, 265.0, 241.0, 256.0],
    ['N', 'ASN', 195.0, 187.0, 158.0, 165.5],
    ['D', 'ASP', 193.0, 187.0, 151.0, 158.7],
    ['C', 'CYS', 167.0, 148.0, 140.0, 146.1],
    ['E', 'GLU', 223.0, 214.0, 183.0, 186.2],
    ['Q', 'GLN', 225.0, 214.0, 189.0, 193.2],
    ['G', 'GLY', 104.0,  97.0,  85.0,  88.1],
    ['H', 'HIS', 224.0, 216.0, 194.0, 202.5],
    ['I', 'ILE', 197.0, 195.0, 182.0, 181.0],
    ['L', 'LEU', 201.0, 191.0, 180.0, 193.1],
    ['K', 'LYS', 236.0, 230.0, 211.0, 225.8],
    ['M', 'MET', 224.0, 203.0, 204.0, 203.4],
    ['F', 'PHE', 240.0, 228.0, 218.0, 222.8],
    ['P', 'PRO', 159.0, 154.0, 143.0, 146.8],
    ['S', 'SER', 155.0, 143.0, 122.0, 129.8],
    ['T', 'THR', 172.0, 163.0, 146.0, 152.5],
    ['W', 'TRP', 285.0, 264.0, 259.0, 266.3],
    ['Y', 'TYR', 263.0, 255.0, 229.0, 236.8],
    ['V', 'VAL', 174.0, 165.0, 160.0, 164.5]
], columns=colnames)
"""
:data:`.properties` is the main container to hold static data necessary to interpret PDB codes for compound.

The data contained in this variable is loaded at runtime and generated from the
`PDBeChem <ftp://ftp.ebi.ac.uk/pub/databases/msd/pdbechem/>`_ database.
"""