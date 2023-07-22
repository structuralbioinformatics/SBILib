"""
PDB

author: jbonet
date:   02/2013

@oliva's lab
"""

import numpy
import math

from .chain import ChainOfProtein, ChainOfNucleotide

from SBILib.beans.StorableObject import StorableObject
from SBILib.beans.File           import File
from SBILib.structure.geometry import RMSD


class PDB(StorableObject):

    """
    A {PDB} is a collection of {Chain}
    """
    def __init__(self, pdb_file=None, cif_file=None, dehydrate=False, header=False,
                 onlyheader=False, biomolecule=False):
        """
        @type  pdb_file: String
        @param pdb_file: PDB formated file to read

        @raise IOError if pdb_file does not exist and it is not an empty object
        """
        if biomolecule or onlyheader:
            header = True

        self._pdb_file      = pdb_file
        # Additional cif variable added -Patrick Gohl
        self._cif_file      = cif_file
        self._chains        = []
        self._NMR           = False
        self._NMR_chains    = []
        self._chain_id      = set()

        self._biomol_id     = -1    # -1 -> original
                                    #  0 -> symmetry
                                    # >0 -> biomolecule

        self._header        = None

        self._has_prot      = False
        self._has_nucl      = False

        self._COMPND        = None

        if self.pdb_file is not None:
            self._pdb_file  = File(file_name=self._pdb_file, action='r')
            self._read_PDB_file(header=header,
                                onlyheader=onlyheader,
                                biomolecule=biomolecule)

        """
        Adding the option to read a mmCIF file
        """
        if self._cif_file is not None:
            self._cif_file  = File(file_name=self._cif_file, action='r')
            self._read_mmCIF_file(header=header,
                                onlyheader=onlyheader,
                                biomolecule=biomolecule)
        """
        -Patrick Gohl
        """


        if dehydrate:
            self.dehydrate()

    #
    # ATTRIBUTES
    #
    @property
    def pdb_file(self):
        """
        PDB file name
        @rtype: String
        """
        return self._pdb_file

    @pdb_file.setter
    def pdb_file(self, value):
        """
        Sets a PDB file if none has been given
        @raise UsedAttributeError
        """
        if self._pdb_file is not None:
            raise AttributeError(
                "The PDB object is loaded from file {0}. To load the new file {1} create a new PDB object".format(self._pdb_file.full, value))

        if isinstance(value, File):
            self._pdb_file = value
        else:
            self._pdb_file = File(file_name=value, type='r')

    @property
    def chain_identifiers(self):
        return self._chain_id

    @property
    def id(self):
        return self._chains[0].pdb

    @property
    def chains(self):
        """
        List of {Chain} contained in the PDB w/out NMR replicas
        @rtype: List of {Chain}
        """
        return self._chains

    @property
    def proteins(self):
        """
        List of {ProteinChain} contained in the PDB w/out NMR replicas
        @rtype: List of {ProteinChain} (iterator)
        """
        for chain in self.chains:
            if isinstance(chain, ChainOfProtein):
                yield chain

    @property
    def nucleotides(self):
        """
        List of {NucleotideChain} contained in the PDB w/out NMR replicas
        @rtype: List of {NucleotideChain} (iterator)
        """
        for chain in self.chains:
            if isinstance(chain, ChainOfNucleotide):
                yield chain

    @property
    def non_standard_chains(self):
        """
        List of non {NucleotideChain}/ non {ProteinChain} contained in the PDB w/out NMR replicas
        @rtype: List of non {NucleotideChain}/ non {ProteinChain} (iterator)
        """
        for chain in self.chains:
            if not isinstance(chain, ChainOfNucleotide) and not isinstance(chain, ChainOfProtein):
                yield chain

    @property
    def all_models(self):
        """
        List of {Chain} contained in the PDB w/ NMR replicas
        @rtype: List of {Chain}
        """
        return self._chains + self._NMR_chains

    @property
    def header(self):
        if self._header is None:
            return ''
        else:
            return self._header

    @property
    def biomolecule_identifier(self):
        return self._biomol_id

    #
    # COMPLEX GETTERS & SETTERS
    #
    def get_chain_by_id(self, id):
        """
        Returns a chain according to its id or None if no chain with that id is found
        @rtype: {Chain}
        """
        for chain in self._chains:
            if chain.chain == id:
                return chain
        return None

    def add_chain(self, chain, NMR=False):
        """
        Adds a new chain to the PDB
        """
        if not NMR:
            self._chains.append(chain)
        elif NMR and self._NMR:
            self._NMR_chains.append(chain)

        self._chain_id.add(chain.chain)

    def add_chains(self, chains, NMR=False):
        """
        Adds a new chains to the PDB
        """
        for chain in chains:
            self.add_chain(chain=chain, NMR=NMR)

    def _get_chain_position_by_id(self, id):
        """
        Returns the position in the chain array where the chain is
        @rtype: Integer
        """
        for x in range(len(self._chains)):
            if self._chains[x].chain == id:
                return x
        return None

    #
    # BOOLEANS
    #
    @property
    def is_NMR(self):
        """
        Identifies if the PDB contains NMRs
        @rtype: Boolean
        """
        return self._NMR

    def chain_exists(self, chain):
        """
        Confirms if a given chain exists in the PDB
        @rtype: Boolean
        """
        return chain in self._chain_id

    @property
    def has_protein(self):
        """
        Checks if the PDB contains a protein (not only)
        @rtype: Boolean
        """
        return self._has_prot

    @property
    def has_nucleotide(self):
        """
        Checks if the PDB contains a nucleotide chain (not only)
        @rtype: Boolean
        """
        return self._has_nucl

    @property
    def repeated_chain_ids(self):
        """
        Checks if more than one {Chain} has the same assigned ID
        @rtype: Boolean
        """
        return len(self._chain_id) < len(self._chains)

    @property
    def is_all_ca(self):
        for p in self.proteins:
            if p.is_only_ca:
                return True
        return False

    #
    # METHODS
    #
    def dehydrate(self):
        recheck_chains = False
        for c in self.chains:
            c.dehydrate()
            if c.is_empty:
                recheck_chains = True
        if recheck_chains:
            c = []
            for ch in self.chains:
                if not ch.is_empty:
                    c.append(ch)
                else:
                    self._chain_id.remove(ch.chain)
            self._chains = c

    def duplicate(self, hetero=True, water=False, NMR=False):
        """
        Returns a {PDB} identical to the original but as a new object
        @rtype: {PDB}
        """
        new_PDB = PDB()
        new_PDB.pdb_file = self.pdb_file

        for chain in self.chains:
            new_PDB.add_chain(
                chain=chain.duplicate(hetero=hetero, water=water))

        if NMR:
            for chain in self._NMR_chains:
                new_PDB.add_chain(chain=chain.duplicate(
                    hetero=hetero, water=water), NMR=True)

        new_PDB._NMR = self._NMR
        new_PDB._has_prot = self._has_prot
        new_PDB._has_nucl = self._has_nucl

        return new_PDB

    def apply_symmetry_matrices(self):
        """
        Only works if the PDB file is an original PDB file
        or the matrices have been added in the correct PDB format
        @rtype: {PDB}
        """
        if self._header is None:
            self._read_PDB_file(header=True, onlyheader=True)
        return self._apply_matrix(matrix=self.header.symmetry_matrix)

    def apply_biomolecule_matrices(self, keepchains=False, water=True):
        """
        Only works if the PDB file is an original PDB file or
        the matrices have been added in the correct PDB format
        @rtype: {PDB}
        """
        if self._header is None:
            self._read_PDB_file(header=True, onlyheader=True)
        PDB_list = []
        for matrix in self.header.biomolecules:
            PDB_list.append(self._apply_matrix(matrix=matrix,
                                               keepchains=keepchains,
                                               water=water))
        return PDB_list

    def _apply_matrix(self, matrix, keepchains=False, water=True):
        new_PDB            = PDB()
        new_PDB._biomol_id = matrix.identifier

        for chain in self.chains:
            if chain.chain in matrix.chains:
                for mat in matrix.matrices:
                    new_chain = chain.duplicate(water=water)
                    new_chain.reposition(matrix=mat.matrix, vector=mat.vector)
                    if len(new_chain) >= 1:
                        new_PDB.add_chain(chain=new_chain)
        if not keepchains:
                new_PDB.tmpclean(cluster_by_alternative_id=True)
        return new_PDB

    def clean(self):
        first_atom = 1
        for c in self.chains:
            c.clean(initatom=first_atom)
            first_atom = c.last_residue.last_atom_number + 1

    def tmpclean(self, cluster_by_alternative_id=False):
        """
        Makes a clean version of the PDB, rechaining in order and renumerating atoms.
        Renumbering residues is optional
        """
        pchainsIDs = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
        chainsIDs = ""
        chainsNIDs = ""
        chainID = 0
        atom_count = 1

        for x in range(len(pchainsIDs)):
            if not self.chain_exists(chain=pchainsIDs[x]):
                chainsIDs += pchainsIDs[x]
            else:
                chainsNIDs += pchainsIDs[x]

        chain_change = len(self) <= len(chainsIDs)

        for chain in self.chains:
            if (not chain.chain in chainsNIDs) and chain_change:
                self._chain_id.add(chain.chain)
                chain.chain = chainsIDs[chainID]
                chainID += 1
                self._chain_id.add(chain.chain)
                if cluster_by_alternative_id:
                    if self._COMPND is None:
                        self._COMPND = {}
                    if chain.alternative_id not in self._COMPND:
                        self._COMPND.setdefault(
                            chain.alternative_id, []).append(chain.alternative_id)
                    self._COMPND[chain.alternative_id].append(chain.chain)
            else:
                chainsNIDs = chainsNIDs.replace(chain.chain, '')

            chain.renumerate_atoms(init=atom_count)
            atom_count += (chain.atom_length)

    def fuse_chains(self, chains_ids):
        """
        Fuses several chains into the first one.
        It will not allow to fuse different structural chains.
        It does not alter the {PDB}, but provides a new one
        @rtype: {Chain}

        @raise AttributeError if:
            a) A given chain ID is not present
            b) Try to fuse different structural chains
        """
        if len(self._chain_id.intersection(set(chains_ids))) < len(chains_ids):
            raise AttributeError(
                "Some of the given chains to fues do not exist")

        error_counter = 0
        error_control = [False, False]
        new_PDB = PDB()
        for c in chains_ids:
            chain = self.get_chain_by_id(id=c)
            new_PDB.add_chain(chain=chain.duplicate())
            if isinstance(chain, ChainOfProtein) and not error_control[0]:
                error_counter += 1
                error_control[0] = True
            elif isinstance(chain, ChainOfNucleotide) and not error_control[1]:
                error_counter += 1
                error_control[1] = True
            if error_counter == 2:
                raise AttributeError(
                    "Fuse different kinds of structural chain is not possible\n")

        init_chain_num = new_PDB.chains[0].last_residue.number
        for x in range(1, len(new_PDB.chains)):
            new_PDB.chains[x].renumerate_residues(init=init_chain_num + 1)
            init_chain_num = new_PDB.chains[0].last_residue.number
            new_PDB.chains[0].fuse(chain=new_PDB.chains[x])

        return_PDB = PDB()
        return_PDB.add_chain(chain=new_PDB.chains[0])
        return return_PDB

    # def calculate_dssp(self, out_dir = None, store = True):
    #     """
    #     Executes DSSP and assigns the prediction to each chain

    #     @param  out_dir: directory to save the output
    #     @defaut out_dir: None

    #     @param store: Save the dssp output(?)
    #     """

    #     for chain in self.proteins:
    #         if out_dir is None:
    #             pdb_file  = chain.globalID + ".pdb2dssp"
    #             dssp_file = chain.globalID + ".dssp"
    #         else:
    #             Path.mkdir(newdir = out_dir)
    #             pdb_file  = os.path.join(os.path.abspath(out_dir), chain.globalID + ".pdb2dssp")
    #             dssp_file = os.path.join(os.path.abspath(out_dir), chain.globalID + ".dssp")

    #         pdb_fd = open(pdb_file, 'w')
    #         pdb_fd.write(chain.PDB_format())
    #         pdb_fd.close()

    #         dssp_calc = DSSPexec(pdb_file = pdb_file, dssp_file = dssp_file,
    #                              chain    = chain,    store     = store)

    def rotate(self, matrix=None):
        """
        Rotates each {Chain} according to a given matrix

        @type matrix: numpy.matrix
        """
        if matrix is None:
            matrix = numpy.identity(3, float)
        for chain in self.all_models:
            chain.rotate(matrix=matrix)

    def translate(self, vector=None):
        """
        Translates each {Chain} according to a translational vector

        @type vector: numpy.array
        """
        if vector is None:
            vector = numpy.zeros(3, float)
        for chain in self.all_models:
            chain.translate(vector=vector)

    def reposition(self, matrix=None, vector=None):
        """
        Rotates and Translates each {Chain} according to a matrix and a translational vector

        @type matrix: numpy.matrix

        @type vector: numpy.array
        """
        if matrix is None:
            matrix = numpy.identity(3, float)
        if vector is None:
            vector = numpy.zeros(3, float)
        for chain in self.all_models:
            chain.reposition(matrix=matrix, vector=vector)

    # def calculate_protein_heteroatom_contacts(self, distance = 6):
    #     """
    #     Returns a {HeteroatomContacts} list with the contacts between a protein and its heteroatoms
    #     at a maximum given distance
    #     @type distance: Integer
    #     @rtype: list of {HeteroatomContacts}
    #     """
    #     data = []
    #     for protein in self.proteins:
    #         data.append(HeteroatomContacts(chain = protein, max_distance = distance))
    #     return data

    #
    # OVERRIDE PARENT'S FUNCTIONS
    #
    @staticmethod
    def read(input_file, format='PDB'):
        """
        Reads a file of data in a specific format and returns the object

        @type  input_file: String
        @param input_file: File to read

        @type  format: String
        @param format: Format of the file to read
        """
        if format == 'PDB':
            pdb = PDB(pdb_file=input_file)
            return pdb

    def write(self, output_file=None, format='PDB', force=False, clean=False):
        """
        Writes the object in a specific format

        @type  output_file: String
        @param output_file: File to write

        @type  format: String
        @param format: Format of the file to print
        """
        outfile = File(file_name=output_file, action='w', overwrite=force)
        if format == 'PDB':
            self._write_PDB_file(pdb_file=outfile, clean=clean)

    #
    # IO
    #
    def _read_PDB_file(self, header=False, onlyheader=False, biomolecule=False):
        """
        Process and load crystal data from a PDB formated file
        """
        from .parse_pdb import read_PDB_file, read_PDB_header
        if header:
            read_PDB_header(self)
            self._pdb_file.close()
        if not onlyheader:
            # read_PDB_file(self, biomolecule=biomolecule)
            read_PDB_file(self)
        self._pdb_file.close()

    """
    Reading a mmCIF file
    """
    def _read_mmCIF_file(self, header=False, onlyheader=False, biomolecule=False):
        """
        Process and load crystal data from a PDB formated file
        """
        from . import parse_mmCIF
        MM = parse_mmCIF.mmCIF(self)
        MM.read(self)
        self._cif_file.close()
    """
    -Patrick Gohl
    """

    # def _represent_COMPND(self):
    #     if self._COMPND is None: return ''

    #     data = []
    #     mol_counter = 1
    #     for chain in self._COMPND:
    #         data.append("COMPND    MOL_ID: %d;" %mol_counter)
    #         data.append("COMPND   2 CHAIN: " + ",".join(self._COMPND[chain]) + ";")
    #         if len(self._biomolecA) > 0:
    #             matrices = []
    #             for mat in self._biomolecA:
    #                 if mat[1] == chain: matrices.append(str(mat[0]))
    #             data.append("COMPND   3 MATRICES: " + ",".join(sorted(matrices)))
    #         mol_counter += 1
    #     return "\n".join(data) + "\n"

    def _write_PDB_file(self, pdb_file, clean=False):
        """
        Print a crystal into a PDB formated file
        """
        out_fd = pdb_file.descriptor
        # out_fd.write(self._represent_COMPND())
        out_fd.write(self.PDB_format(clean=clean) + "\n")
        pdb_file.close()

    def PDB_format(self, clean=False, terminal=True):
        """
        Strings a {PDB} in PDB format
        @rtype: String
        """
        lines = []
        if clean:
            self.clean()
        for chain in self._chains:
            lines.append(chain.PDB_format(terminal=terminal))
        lines.append("END")

        return "\n".join(lines)

    def FASTA_format(self, gapped=True, protein=True, nucleotide=False):
        lines = []
        for c in self.chains:
            if isinstance(c, ChainOfProtein) and protein:
                lines.append(
                    ">{0}\t{1}".format(c.globalID, c.aminoacids[0].identifier))
                if gapped:
                    lines.append("{0}".format(c.gapped_protein_sequence))
                else:
                    lines.append("{0}".format(c.protein_sequence))
            if isinstance(c, ChainOfNucleotide) and nucleotide:
                lines.append(
                    ">{0}\t{1}".format(c.globalID, c.nucleotides[0].identifier))
                if gapped:
                    lines.append("{0}".format(c.gapped_nucleotide_sequence()))
                else:
                    lines.append("{0}".format(c.nucleotide_sequence()))
        if len(lines) == 0:
            return ""
        else:
            return "\n".join(lines) + "\n"

    def IDX_format(self, protein=True, nucleotide=False):
        lines = []
        for c in self.chains:
            if isinstance(c, ChainOfProtein) and protein:
                lines.append(">{0}\t{1}".format(c.globalID, c.protein_idx))
            if isinstance(c, ChainOfNucleotide) and nucleotide:
                lines.append(
                    ">{0}\t{1}".format(c.globalID, c.nucleotide_idx()))
        if len(lines) == 0:
            return ""
        else:
            return "\n".join(lines) + "\n"

    def FASTA_IDX(self, protein=True, nucleotide=False):
        data = {}
        data.setdefault('FASTA', [])
        data.setdefault('IDX', [])
        for c in self.chains:
            if isinstance(c, ChainOfProtein) and protein:
                data['FASTA'].append(
                    ">{0}\n{1}".format(c.globalID, c.gapped_protein_sequence))
                data['IDX'].append(
                    ">{0}\t{1}".format(c.globalID, c.protein_idx))
            if isinstance(c, ChainOfNucleotide) and nucleotide:
                data['FASTA'].append(
                    ">{0}\n{1}".format(c.globalID, c.gapped_nucleotide_sequence()))
                data['IDX'].append(
                    ">{0}\t{1}".format(c.globalID, c.nucleotide_idx()))

        return data

    """
    Gathering the loop similarities between two proteins
    """
    def compare_loops(self, PDB2, rho_thresh=10, delta_thresh=5, theta_thresh=5, distance_thresh=0.5):
        Similarities = []
        for target in PDB2.proteins:
            target.clean()
            target.calculate_dssp()
            target.calculate_archs()

        for query in self.proteins:
            query.clean()
            query.calculate_dssp()
            query.calculate_archs()
            for query_arch in query.archs:
                for target in PDB2.proteins:
                    for target_arch in target.archs:
                        Rho = math.sqrt((query_arch.rho - target_arch.rho)**2)
                        Delta = math.sqrt((query_arch.delta - target_arch.delta)**2)
                        Theta = math.sqrt((query_arch.theta - target_arch.theta)**2)
                        Distance = math.sqrt((query_arch.cartesian_distance - target_arch.cartesian_distance)**2)
                        if Distance <= distance_thresh  and Rho <= rho_thresh and Delta <= delta_thresh and Theta <= theta_thresh:
                            Similarities.append("QueryProt:"+query.chain+","+str(int(query_arch.initial_position)+query.first_aminoacid.number -1 )+","+str(int(query_arch.end_position)+query.first_aminoacid.number -1)+":TargetProt:"+target.chain+","+str(int(target_arch.initial_position)+target.first_aminoacid.number -1)+","+str(int(target_arch.end_position)+ target.first_aminoacid.number-1)+":"+ query_arch.aminoacid_sequence + "," + target_arch.aminoacid_sequence )

        if len(Similarities) == 0:
            print("No similar loops at the given thresholds were found")

        return Similarities
    """
    Patrick Gohl
    """

    """ 
    Grafting a selected loop pair 
    """
    def graft(self, PDB2, selection):
        # gathering loop info from selection
        query_chain_id = selection.split(":")[1].split(",")[0]
        query_start = selection.split(":")[1].split(",")[1]
        query_end = selection.split(":")[1].split(",")[2]
        target_chain_id = selection.split(":")[3].split(",")[0]
        target_start = selection.split(":")[3].split(",")[1]
        target_end = selection.split(":")[3].split(",")[2]

        # grabbing the chains with selected loops
        query_chain = self.get_chain_by_id(query_chain_id)
        target_chain = PDB2.get_chain_by_id(target_chain_id)

        #grabbing the guiding loop 
        nguide = query_chain.extract(int(query_start)-1 , int(query_end) - 1 , backbone=True)
        # retriveing vector and center of guide
        nguide_vector = []
        for r in nguide.aminoacids:
            for a in r.atoms:
                nguide_vector.append(numpy.copy(a.coordinates))

        nguide_center = numpy.mean(nguide_vector, axis = 0, dtype = numpy.float64)

        # Doing the same for the target loop
        nquery = target_chain.extract(int(target_start) -1 , int(target_end)-1 , backbone=True)
        nquery_vector = []
        for r in nquery.aminoacids:
            for a in r.atoms:
                nquery_vector.append(numpy.copy(a.coordinates))

        nquery_center = numpy.mean(nquery_vector, axis = 0, dtype = numpy.float64)

        # Calculating rotation and translation and Super positioning 
        M, nrmsd = RMSD.rigid_transform_3D(numpy.copy(nquery_vector), numpy.copy(nguide_vector))
        target_chain.translate(-nquery_center)
        target_chain.rotate(M)
        target_chain.translate(nguide_center)

        #Creating Grafted Protein
        NewChain = query_chain.extract(query_chain.first_aminoacid.number,query_chain.first_aminoacid.number+int(query_start)-2, backbone=True)
        extension = target_chain.extract(target_chain.first_aminoacid.number + int(target_start)-1 , target_chain.first_aminoacid.number + int(target_end) -1 , backbone=True)
        NewChain.join(extension)
        extension = query_chain.extract(query_chain.first_aminoacid.number+int(query_end) , -1 , backbone=True)
        NewChain.join(extension)
        NewChain.chain = query_chain_id
        GraftedProtein = PDB()
        for chain in self.chain_identifiers:
            if chain != query_chain_id:
                GraftedProtein.add_chain(self.get_chain_by_id(chain))
                if (self.get_chain_by_id(chain).chaintype == "P"):
                    GraftedProtein._has_prot = "TRUE"
                elif (self.get_chain_by_id(chain).chaintype == "N"):
                    GraftedProtein._has_nucl = "TRUE"
            else:
                GraftedProtein.add_chain(NewChain)

        return GraftedProtein



    """
    Patrick Gohl
    """





    #
    # OVERRIDE DEFAULT METHODS
    #
    def __len__(self):
        return len(self._chains)
