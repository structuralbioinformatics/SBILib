from collections import OrderedDict
from SBILib.data       import aminoacids3to1, nucleic3to1
from .chain         import Chain, ChainOfProtein, ChainOfNucleotide
import gzip
import re

columns_order = [
    "group_PDB",
    "id",
    "type_symbol",
    "label_atom_id",
    "label_alt_id",
    "label_comp_id",
    "label_asym_id",
    "label_entity_id",
    "label_seq_id",
    "pdbx_PDB_ins_code",
    "Cartn_x",
    "Cartn_y",
    "Cartn_z",
    "occupancy",
    "B_iso_or_equiv",
    "pdbx_formal_charge",
    "auth_seq_id",
    "auth_comp_id",
    "auth_asym_id",
    "auth_atom_id",
    "pdbx_PDB_model_num"
]


def specialSplit( content, monoliners ):
    output = [["", False]]
    quote  = False
    length = len(content)

    if re.search(r'^(?:"\w*)\'\S{1}\'^(?:"\w*)', content) and not content.startswith('_'):
        output[-1][0] == content
        return output

    if monoliners:
        if content.startswith(';') and content.endswith(';'):
            output[-1][0] = content.strip(';')
            return output

    for c in range(length):
        isWS   = content[c] == " " or content[c] == "\t"
        wasWS  = c == 0 or content[c - 1] == " " or content[c - 1] == "\t"
        willWS = c == length - 1 or content[c + 1] == " " or content[c + 1] == "\t"
        braket = frozenset(["'", '"'])
        if (content[c] in braket) and (wasWS or willWS):
            quote = not quote
        elif not quote and isWS and output[-1][0] != "":
            output.append(["", False])
        elif not quote and content[c] == "#":
            # Beware of identifiers with # sign!
            # Problem PBD ID: 2XQB
            if c == 0 or content[c - 1] == ' ':
                break
        elif not isWS or quote:
            output[-1][0] += content[c]
            output[-1][1] = quote
    if output[-1][0] == "":
        output.pop()
    return output


def typefy( content, only_empties=False ):
    if content == "?":
        return " "
    if content == ".":
        return " "

    if only_empties:
        return content

    if "." in content:
        try:
            return float(content)
        except Exception:
            pass
    else:
        try:
            return int(content)
        except Exception:
            pass
    return content


class mmCIF( object ):
    """Read mmCIF format as defined in http://mmcif.wwpdb.org/ and transform it
    into a json object.
    This parser is based on the one created by Gert-Jan Bekker
    <http://github.com/gjbekker/cif-parsers>, but highly simplified.
    """
    def __init__( self, pdbobject, monoliners=False ):
        self.multiline  = False
        self.buffer     = ""
        self.data       = OrderedDict()
        self.loopstatus = False
        self.keyname    = None
        self.lastkey    = None
        self.lastsubkey = None
        self.headers    = frozenset(['loop_', 'save_', 'global_', 'data_'])
        self.monoliners = monoliners
        self.WriteCIF = True

    def read( self, pdbobject ):
        """Parse the provided file.
        :param str file_name: Name of the file to parse.
        :return: json containing the file's data
        """
        file_name = pdbobject._cif_file.name
        fd = gzip.open(file_name, 'rt') if file_name.endswith("gz") else open(file_name)
        return self._from_filehandle(fd, pdbobject)

    def _from_filehandle( self, fhandle, pdbobject ):
        for line in [l.rstrip() for l in fhandle]:
            if len(line) == 0:
                continue
            if self._check_multiline(line):
                self._assign(specialSplit(self.buffer, self.monoliners), pdbobject)
                self.buffer = ""
        info = self.data
        self._set_default(True)
        return info

    def _is_header( self, data ):
        for x in self.headers:
            if data[0].startswith(x) and not data[1]:
                return True
        return False

    def _manage_header( self, data ):
        if data[0][0] == 'loop_':
            self.loopstatus = True
        elif data[0][0].startswith('data_'):
            self.keyname = data[0][0][5:].strip()
            self.data.setdefault(self.keyname, {})
        elif data[0][0].startswith('save_'):
            pass  # TODO
        elif data[0][0].startswith('global_'):
            pass  # TODO

    def _set_default( self, strong=False ):
        self.loopstatus = False
        self.lastkey    = None
        if strong:
            self.multiline  = False
            self.buffer     = ""
            self.data       = OrderedDict()
            self.keyname    = None

    def _assign( self, data, pdbobject ):
        if len(data) == 0:
            self._set_default()
        elif self._is_header(data[0]):
            self._manage_header(data)
        elif not self.loopstatus:
            self._no_loop_status(data, pdbobject)
        else:
            self._loop_status(data, pdbobject)

    def _loop_status( self, data, pdbobject ):
        if data[0][0].startswith('_'):
            k = data[0][0].split('.')
            if self.lastkey is None:
                self.lastkey = {}
            self.lastkey.setdefault(k[0], []).append(k[1])
            self.data[self.keyname].setdefault(k[0], OrderedDict()).setdefault(k[1], [])
            self.lastsubkey = 0
        else:
            k = list(self.lastkey.keys())[0]
            for x in range(len(data)):
                #try:
                """
                subtracted 1 from [self.lastsubkey + x] for proper indexing
                """
                self.data[self.keyname][k][self.lastkey[k][self.lastsubkey + x - 1]].append(typefy(data[x][0], True))
                #except:
                     #print(len(data))
                     #print(x)
                     #exit(0)
                """
                -Patrick Gohl
                """
            self.lastsubkey = self.lastsubkey + x + 1
            
            if k == '_pdbx_poly_seq_scheme' and self.WriteCIF:
                self._write_to_pdb_object(pdbobject)
                self.WriteCIF = False
            if self.lastsubkey >= len(self.lastkey[k]):
                self.lastsubkey = 0

    def _no_loop_status( self, data, pdbobject ):
        if data[0][0].startswith('_'):
            self.lastkey = data[0][0].split('.')
            self.data[self.keyname].setdefault(self.lastkey[0], OrderedDict()).setdefault(self.lastkey[1], [])
            if len(data) == 2:
                d = typefy(data[1][0], True)
                self.data[self.keyname][self.lastkey[0]][self.lastkey[1]].append(d)
                self.lastkey = None
        elif self.lastkey is not None:
            if len(data) == 1:
                d = typefy(dataprotein_key[0][0], True)
                self.data[self.keyname][self.lastkey[0]][self.lastkey[1]].append(d)
            else:
                d = typefy(" ".join([x[0] for x in data]))
                self.data[self.keyname][self.lastkey[0]][self.lastkey[1]].append(d)
            self.lastkey = None

    # def _check_multiline( self, line ):
    #     if line[0] == ';':
    #         self.multiline = not self.multiline
    #     if self.multiline or line.strip() == ';':
    #         self.buffer += line  # .lstrip(';')
    #         return False
    #     else:
    #         if len(self.buffer) == 0:
    #             self.buffer = line
    #         return True

    def _check_multiline( self, line ):
        if line[0] == ';':
            self.buffer += ';'
            self.multiline = not self.multiline
        if self.multiline:
            if not self.monoliners:
                self.buffer += line.lstrip(';')
            else:
                self.buffer += line
            return False
        else:
            if len(self.buffer) == 0:
                self.buffer = line
            return True

    def _write_to_pdb_object(self, pdbobject):
            for key in list(self.data.keys()):
                old_chain = "Old_Chain"
                
                for index in range(len(self.data[key]['_atom_site']["label_entity_id"])):
                    """
                    The ordering is off in the column keys, everything was shifted one position to the right 
                    """
                    pos1 =  self.data[key]['_atom_site']['pdbx_PDB_model_num'][index]  
                    pos2 = self.data[key]['_atom_site']['group_PDB'][index]
                    pos3 = self.data[key]['_atom_site']['type_symbol'][index]
                    pos4 = self.data[key]['_atom_site']['label_alt_id'][index]
                    pos5 = self.data[key]['_atom_site']['auth_comp_id'][index]
                    pos6 = self.data[key]['_atom_site']['pdbx_formal_charge'][index]
                    pos7 = self.data[key]['_atom_site']['pdbx_PDB_ins_code'][index]
                    pos8 = self.data[key]['_atom_site']['Cartn_x'][index]
                    pos9 = self.data[key]['_atom_site']['Cartn_y'][index]
                    pos10 = self.data[key]['_atom_site']['Cartn_z'][index]
                    pos11 = self.data[key]['_atom_site']['occupancy'][index]
                    pos12 = self.data[key]['_atom_site']['id'][index]
                    pos13 = self.data[key]['_atom_site']['B_iso_or_equiv'][index]
                    pdbline = pos1.ljust(6) + pos2.rjust(5) + "  " + pos3.ljust(4) + pos4.rjust(3) + pos5.rjust(2) + pos6.rjust(4) + "    " + pos7.rjust(8) + pos8.rjust(8) + pos9.rjust(8) + pos10.rjust(6) + pos11.rjust(6) + pos12.rjust(12) + pos13.rjust(2)
  
                    chain = self.data[key]['_atom_site']['auth_comp_id'][index]
                    
                    if self.data[key]['_atom_site']['auth_atom_id'][index] != "1":
                        old_chain = "Old_Chain"
                        pdbobject._NMR = True
                    if chain != old_chain:
                        chain_type = self.data[key]['_atom_site']['label_alt_id'][index]
                        if chain_type in aminoacids3to1:
                            obj_chain = ChainOfProtein(pdb = key ,chain = chain)
                            pdbobject._has_prot = True
                        elif chain_type in nucleic3to1:
                            obj_chain = ChainOfNucleotide(pdb = key ,chain = chain)
                            pdbobject._has_nucl = True
                        else:
                            obj_chain = Chain(pdb = key ,chain = chain)
                            pdbobject.add_chain(obj_chain, NMR = pdbobject.is_NMR) 
                    else:
                        putative_old_chain = pdbobject.get_chain_by_id(id = chain)
                        obj_chain = putative_old_chain
                    
                    old_chain = chain

                    chain_type = self.data[key]['_atom_site']['label_alt_id'][index]
                    if not isinstance(obj_chain, ChainOfNucleotide) and not isinstance(obj_chain, ChainOfProtein):
                        if chain_type in nucleic3to1:
                            newobj_chain = ChainOfNucleotide(pdb = key ,chain = obj_chain.chain)
                            pdbobject._has_nucl = True
                            newobj_chain.fuse(chain = obj_chain, lapl = True)
                            del(obj_chain)
                            obj_chain = newobj_chain
                            pdbobject._chains[pdbobject._get_chain_position_by_id(id = obj_chain.chain)] = obj_chain
                        if chain_type in aminoacids3to1:
                            newobj_chain = ChainOfProtein(pdb = key ,chain = obj_chain.chain)
                            pdbobject._has_prot = True
                            newobj_chain.fuse(chain = obj_chain, lapl = True)
                            del(obj_chain)
                            obj_chain = newobj_chain
                            pdbobject._chains[pdbobject._get_chain_position_by_id(id = obj_chain.chain)] = obj_chain
                    read_translated_line(obj_chain, pdbline)
                    if not pdbobject.chain_exists(chain):
                        pdbobject.add_chain(obj_chain)
                    
            pdbobject._chain_id = set([c.chain for c in pdbobject.chains if not c.is_empty])
            pdbobject._chains   = [c for c in pdbobject.chains if not c.is_empty]

def read_translated_line(thischain, line, keep_version = "A"):
    """
    Given a PDB-formated line, creates an atom to add to a new or a pre-existent residue

    @type  line: String
    @param line: PDB formated line

    @type  keep_version: String
    @param keep_version: Some residues have two versions, codified in front of the residue_type (line[16:17])
                         By default we keep the A version of doubles Aa, but it can be changed through parameters
    """
    isOK = set([" ", keep_version, thischain._residue_version])

    #This skips several versions of the same residue. It always will grab " " and keep_version
    #but in chains that are all B (1EA4_Z) this will help grab them
    if line[16:17] not in isOK:
        if len(thischain) == 0:
            thischain._residue_version = line[16:17]
        else:
            return

    residue_num = int(line[22:26].strip())
    residue_ver = line[26:27]

    #Check the need of creating a new residue
    if thischain.is_empty or (len(thischain._last_appended_list) > 0 and thischain._last_appended_list[-1].identifier != str(residue_num) + residue_ver):
        residue = thischain._new_Residue(number = residue_num, version = residue_ver, Rtype = line[17:20].strip(), mode = line[:6].strip())
        thischain.add_residue(residue)


    #Create the new atom
    x, y, z = [float(line[30+8*i:38+8*i]) for i in range(3)]
    try:
        occupancy = float(line[54:60])
    except:
        occupancy = ""
    try:
        tempFactor = float(line[60:66])
    except:
        tempFactor = ""
    try:
        element = line[76:78].strip()
    except:
        element = ""
    try:
        charge = line[78:80].strip()
    except:
        charge = ""
    atom = thischain._new_Atom(number = line[6:12].strip(), name = line[12:16].strip(), x = x, y = y, z = z, occupancy = occupancy, tempFactor = tempFactor, element = element, charge = charge)

    #Add the atom to the last residue (it applies the same if the residue has been created now or before)
    thischain._last_appended_list[-1].add_atom(atom)
    












