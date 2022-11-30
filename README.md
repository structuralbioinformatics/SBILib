# Welcome to the StructuralBioInformatics pylib

The **SBI** python library is a continuous work in progress that clusters functionalities to work with data derived from [PDB][1], from creating a local PDB repository to do _real_ Blast searches over PDB sequences (i.e. sequences of the really crystallized protein sections) to perform a simply split chains (separate the sections of a PDB file).

## Requirements:


You can start to play by reading a PDB using the following command:

```python
from SBI.structure import PDB
newPDBobject = PDB('pdbfilename')
```

The library **is not exempt of bugs** (probably), and can not be considered a final version. This means that one must use it at its own risk. And, while new functionalities are bound to appear, others might disappear at any given moment. Any comments, complains or bug reports can be addressed in the corresponding sections.

[1]: http://www.pdb.org/




# SBI Library Manual 

This manual contains a short turorial split into 5 common Tasks.
Information on how to access more detailed help pages within python is provided at the end.
We have made available a list of slightly more complex scenarios on which the SBI library may be applied. In Scenarios.md

### Dependencies
python3 <br/>
scipy <br/>
numpy <br/>
DSSP <br/>
BLAST <br/>
CD-Hit <br/>


## Installation

Download the SBI library <br/>
Navigate to SBI/external/configSBI.txt <br/>
Alter the path variable for each external dependency to point to the location the executable is stored in on you system. <br/>
Once done the SBI library is ready to use either interactivly or through a script. <br/>


## TASK 1

### 1) loading in the appropriate modules

The following code demonstrates how a PDB file can be loaded: <br/>

The first step is to load in the modules necessary. <br/>
PDB link will allow us to grab a protein from the protein data bank. While the PDB module will allow us to do some basic tasks on the protein. <br/>
```{python}

from SBI.databases import PDBlink
from SBI.structure import PDB

```
### 2) Fetching PDBs

In the next step we create an instance of the PDBlink class and use it to get a pdb from the databank with the coresponding accession code. <br/>
Here we have selected protein "1a3q". <br/>
Note that the 'path/to/pdb' will be saved in the variable 'path'. This can then be read in as a new PDB object. <br/>
If the pdb file is stored locally the address can be provided directly without the need to fetch it from the databank. <br/>

```{python}

new = PDBlink()
path = new.get_PDB("1a3q")
protein1 = PDB(path)

```
### 3) Reading out chains in a protein

First I will load in a second protein this time providing the address directly, and then have a look at the chains contained within the protein. <br/>
We can look at exacty what chains are contained within the protein. <br/>
Here our protein contains 2 chains identified as A and B <br/>

```{python}

protein2 = PDB(path/to/pdbfile.pdb)
protein2.chain_identifiers
# set(['A', 'B'])

```

### 4) Duplicating protein and fusing chains

We can duplicate the protein as a new object with the duplicate command. <br/>
Then if we want to we could fuse chains "A" and "B" into a new chain "A" <br/>
Note that you will not be able to fuse two different kinds of structural chains. <br/>

```{python}

protein3 = protein2.duplicate()
protein3 = protein3.fuse_chains(["A","B"])
protein3.chain_identifiers
# set(["A"])

```


### 5) Adding a chain to a protein

Now let us imagine a scenario where we wanted to take chain A from protein2 and add it to protein1. <br/>
First we have to pull out the chain that we want (A) from protein2. <br/>
Then by checking the chain Identifiers in protein1 we can see that it already has a chain "A". <br/>
This means that we will have to change the identifier of the chain we wish to add. Lets make it "E". <br/>
First we will have to duplicate the chain, creating a new object (ChainE), so that we dont alter the chain Ids of protein2. <br/>
Now we can add the chain to protein1 <br/>
If we check the identifiers again we will see that a new chain "E" has been added to the protein. <br/>


```{python}

ChainA = protein2.get_chain_by_id("A")
protein1.chain_identifiers
# set(['A', 'C', 'B', 'D'])
ChainE = ChainA.duplicate()
ChainE.chain = "E"
protein1.add_chain(ChainE)
protein1.chain_identifiers
# set(['A', 'C', 'B', 'D', 'E'])

```

### 6) Removing chain(s) from a protein

Now we would like to delete this chain again: <br/>
First we set up a temporary empty PDB object and call it protein3. <br/>
Then we will define a list of the chain identifiers we wish to remove (in this case only "E") <br/>
First we loop through all of the chains ids in protein1. <br/>
If that protein is not in the deletion list we provided... <br/>
We add it to the newly created protein3. <br/>
We can then redifine protein1 to leave it without the unwanted chains. <br/>

```{python}
 
protein3 = PDB()
deletion = ["E"]
for chain in protein1.chain_identifiers:
     if chain not in deletion:
             protein3.add_chain(protein1.get_chain_by_id(chain))
             if (protein1.get_chain_by_id(chain).chaintype == "P"):
             	protein3._has_prot = "TRUE"
             elif (protein1.get_chain_by_id(chain).chaintype == "N"):
             	protein3._has_nucl = "TRUE"

protein1 = protein3

```

### 7) Adding multiple chains to a protein

Addition of multiple chains can also be accomplished. <br/>
First we have to grab the second chain from protein2 and change its id like we did to ChainA. <br/>
Then we change the identifiers of both chains as we have done previously <br/>
Then we create a list of the chains that we wish to add <br/>
We then add this list of chains to protein1 <br/>

```{python}

ChainB = protein2.get_chain_by_id("B")
ChainF = ChainB.duplicate()
ChainF.chain = "F"
chainlist = [ChainE,ChainF]
protein1.add_chains(chainlist)
protein1.chain_identifiers
# set(['A','B','C','D','E','F'])

```

You can revert back to the original protein by repeating step 6 and adding "F" to the deletion list. <br/>

### 8) Reading out protein sequence (gapped and exact)

We can now have a look at the protein sequence <br/>
First we can have a look at the exact sequence of the protein without gaps indicated (not containing gaps of non-crystalized positions).  <br/>
Then we look at the gapped version (default). Xs will be placed where the crystal sequence contains gaps. <br/>
We can also read out the sequence of each chain individually as a gapped sequence (as shown with ChainA) <br/>
We can also read out the non gapped sequence of each chain as shown <br/>


```{python}
ChainA = protein1.get_chain_by_id("A")
protein1.FASTA_format(gapped=False)
# "protein sequence in FASTA format"
protein1.FASTA_format()
# "gapped protein sequence in FASTA format"
ChainA.gapped_protein_sequence
# "gapped protein sequence for Chain A only"
ChainA.protein_sequence
# "The exact sequence in the crysal without non crystalized positions indicated"

```

### 9) Getting DNA sequence

In addition to protein sequences the PDB may contain DNA <br/>
We can check if this is the case for protein1 <br/>
In the case where this is true we can read out the sequences in Fasta format <br/>


```{python}
protein1.has_nucleotide
# True
for chain in protein1.nucleotides:
	print(">" + chain.globalID + "\n" + chain.nucleotide_sequence())
# >"globalID"
# "Nucleotide Sequence"
...

``` 

### 10) Cleaning protein and writing PDBfile

We can clean the protein with the following command <br/>
We can then write the PDB of the protein out to a file name that we provide <br/>

```{python}

protein1.clean()
protein1.write("protein1.pdb")

```

## TASK 2

### 1) Dependencies

The following application of this SBI library requires DSSP to be run <br/>

### 2) Configuration

In the configSBI.txt file are stored the paths to the executable necessary. <br/> 
This file is located in SBI/external and looks something like this. <br/>
Here only dssp is shown but the format holds for other dependencies (BLAST, HMMER etc.) <br/>
```text

[dssp]
executable     = mkdssp
path           = /path/to/the/executable/folder

```

If a config file is already available it can be used as well. <br/>
This file can be added to the scripts with the following commands <br/>

```{python}

import os 
os.environ['SBI_CONFIG_FILE'] = "path/to/new/configSBI.txt"

```


### 3) Calculating a chains secondary structure

First we load in a new protein and grab one of the chains as in the previous chapter <br/> 
We can then calculate the secondary structure of the chain <br/>
After this is done the secondary structure can be printed out. <br/>


```{python}

from SBI.structure import PDB
protein = PDB("pdb1a3q.ent.gz")
ChainA = protein.get_chain_by_id("A")
ChainA.calculate_dssp()
ChainA.gapped_protein_secondary_structure
# "gapped secondary structure"

```
### 4) printing out the entire proteins amino acid sequence and secondary structure

```{python}

for chain in protein.proteins:
	chain.calculate_dssp()
	print(">" + chain.globalID + "\n" + chain.gapped_protein_sequence + "\n" + chain.gapped_protein_secondary_structure)

# >"globalID"
# "gapped amino acid sequence"
# "gapped secondary structure"


```

## TASK 3

### 1) Set up a blast  


In order to blast you must have a database <br/>
This is a multifasta file containing protein or DNA sequences depending on the blast type. <br/>
Reformating of the multifasta file will be done automatically once it is provided. <br/>
The blast object will now know that the blast is to be conducted on this database. <br/>

```{python}

from SBI.external.blast import BlastExe
blast = BlastExe(database = "path/to/database/file.fa")

```

### 2) execute blast
Once the blast is set up we need a query sequence to pass. <br/>
For this example we will blast chainA of protein1 against our database. <br/>
With these variables we can now execute the blast and save the to "blastresults". <br/>
The result is an instance of type "BlastResult". <br/>

```{python}

querySequence = ChainA.protein_sequence
queryID = ChainA.globalID
blastresults = blast.execute_query_seq(sequenceID = queryID, sequence = querySequence)
type(blastresults)
# <class 'SBI.external.blast.BlastResult.BlastResult'>

```


### 3) View results
There are several ways that we can view the results of the blast. <br/>
Firstly we can print out a compacted version of all of the results. <br/>
We could also focus to select points of interest (eg. hit_id and evalue). <br/>

```{python}

blastresults.print_compacted_blast()
# query and target info, e-value, coverage etc. 

for hit in blastresults.get_hits():
	print(hit.sequenceID, " ", hit.e_value)
# list of all hits and their associated e-values.


```

### 4) Filter results
If we want to evaluate whether the proteins fall within the Rost Curve twilight zone, <br/>
we can parse through the list and print it off similarly to how it was done above. <br/>
If we want to only consider hits passing the rost evaluation we can pass that as an argument to the "get_hits" command <br/>



```{python}

for hit in blastresults.get_hits():
	print(hit.sequenceID, " ", hit.evaluate_Rost_twilight_zone())
# indication of wich proteins would pass the test

for hit in blastresults.get_hits(tz_type = "ID"):
	print(hit.sequenceID, " ", hit.evaluate_Rost_twilight_zone())
# only hits that pass the test will be shown


```

### 4) PIR 

PIR formatted alignments are available for every hit. <br/>
We can select hit number 13 with the result parameter <br/>
We could add in the filters to grab the PIR of alll of the best results <br/>

```{python}

blastresults.str_PIR(result = 13)
# PIR alignment
for hit in  blastresults.get_hits( evalue = 0, tz_type = "ID"):
	blastresults.str_PIR(result = hit._num_seq)
# All PIRs with hits that match the criteria given 

```


## Task 4

### 1) Complex
For this exercise we will use the protein DNA complex that we used in the previous tasks. <br/>
In order to extract the interfaces from a protein complex we will first have to create an object of the complex class. <br/> 
In createing a new complex class we need to provide the protein complex. <br/>
Optional arguments are PPI,PNI,PHI threshold distances (in Angstrom). At default these are 12, 8 and 6 respectively. <br/>
We can grab the individual chains of the complex straight from this class using the .pdb method if we wanted. <br/>

```{python}

from SBI.structure import Complex
from SBI.structure import PDB
protein = PDB("pdb1a3q.ent.gz")
complex = Complex(protein)
complex.pdb
# PDB object


```

### 2) PPI
The interfaces are already calculated and present within the Complex object. <br/>
The results of the Protein Protein Interfaces are stored in a list and accessible with the PPInterfaces method <br/>
We can count the number of interfaces produced by looking at the length of the list (here just 1) <br/>
To access the interfaces we can simply index the list <br/>
The results can quickly be summarized by converting to a string. <br/>
Alternatively we can enter the interface for more detailed info. <br/>
The interface contains information on all of the contacts in the form of objects of the contact class. We can find the amino acid identifiers of all of the contact residues by iterating through the interface object (here we simply selected the first element of the list because there are no more interfaces). <br/>
We can then print each amino acids position within its chain and the single letter code identifier as well as the coordinates for the alpha carbon of that amino acid. <br/>
(aminoacid1 is the first chain provided and aminoacid2 is the second chain/interactor) <br/>
We can then print the geometric distance between the two amino acids. <br/> 


```{python}
print(complex.PPInterfaces[0].toString())
#Chain1, Chain2, threshold_type(default cb), threshold distance. 

len(complex.PPInterfaces)
# 1

complex.PPInterfaces[0]
# Interface object

for contact in complex.PPInterfaces[0].contacts:
	print(contact.aminoacid1.identifier, contact.aminoacid1.single_letter, contact.aminoacid1.ca.coordinates)
	print(contact.aminoacid2.identifier, contact.aminoacid2.single_letter, contact.aminoacid2.ca.coordinates)
	print(contact.geometric_distance)
# complete list of the contacts, coordinates and distances in this interface

```

### 2) PNI

We can calculate Protein Nucleotide interfaces as well. <br/>
The interfaces are accessible with the PNInterfaces method. <br/> 
When we take a look at the interface list we see that this time we have multiple results. <br/>
These can be accessed by indexing or iterated through. <br/>
We can then gather the position and identifiers for the amino acid and nucleotides. As well as the alpha carbon coordinate in the case of the amino acid and the phosphate coordinate in the case of the nucleotide. <br/>


```{python}
len(complex.PNInterfaces)
# 4 

for contact in complex.PNInterfaces[0].contacts:
	print(contact.nucleotide.identifier, contact.nucleotide.single_letter, contact.nucleotide.p.coordinates)
	print(contact.aminoacid.identifier, contact.aminoacid.single_letter, contact.aminoacid.ca.coordinates)
	print(contact.geometric_distance)
# complete list of the contacts, coordinates and distances in this interface

```

### 3) PHI

In order to find the Protein Heteroatoms interface you must first check whether your chains contains heteroatoms. <br/>
This can be checked with the has_heteroatoms method from the chain class. <br/>
I will load in a pdb file that I know has heteroatoms <br/>
The PHInterface will mix the residue types between the two chains. <br/>
For this reason we will find one that is of type aminoacid so that we can check the type within the contact object. <br/>
Now that we know which is which we can control for the classes in the iterations. <br/>
In the case of an amino acid residue we can print the same information from the previous examples. <br/>
In the case that the residue is a heteroatom we will print the identifier, all of the residues atoms (which contain the coordinats as well), we can also print the coordinates seperately with the all_coordinates method. <br/>



```{python}
for chain in protein.chains:
	print(chain.has_heteroatoms)
# False
# False
# False
# False


protein = PDB("1a2w.pdb")
for chain in protein.chains:
	print(chain.has_heteroatoms)
# True
# True

complex = Complex(protein)
len(complex.PHInterfaces)
# 1

complex.PHInterfaces[0].contacts
# list of contacts

type(list(complex.PHInterfaces[0].contacts)[0].aminoacid)
# type of the first residue

type(list(complex.PHInterfaces[0].contacts)[0].heteroatom)
# type of the second Residue which will be opposite of the first (heteroatom residue vs aminoacid)


for contact in complex.PHInterfaces[0].contacts:
	if (contact.aminoacid.mode == "ATOM"):
		print(contact.aminoacid.identifier, contact.aminoacid.single_letter, contact.aminoacid.ca.coordinates)
		print(contact.heteroatom.identifier, contact.heteroatom.atoms, contact.heteroatom.all_coordinates)
		print(contact.geometric_distance)
	else:
		print(contact.aminoacid.identifier, contact.aminoacid.atoms, contact.aminoacid.all_coordinates)
		print(contact.heteroatom.identifier, contact.heteroatom.single_letter, contact.heteroatom.ca.coordinates)
		print(contact.geometric_distance)
# complete list of the contacts, coordinates and distances in this interface


```


## Task 5

### Protein loop geometry

We will load in the pdb we have worked with before and have a look at the loop geometry in its chain "A".

```{python}

protein = PDB("pdb1a3q.ent.gz")
ChainA = protein.get_chain_by_id("A")

```

In order to visualize the loop geometry we need to first have an idea of the proteins secondary structure. <br/>
Like before this is done with DSSP and can be called within the SBI library. <br/>
Once the secondary structures are resolved we can calculate the loop geometry. <br/>
In the SBI library loops are refered to as archs. <br/>

```{python}

ChainA.calculate_dssp()
ChainA.calculate_archs()
ChainA.archs
ChainA.superarchs

```

## Details on classes
If you wish to check the functions available for any class, they are available with pytons dir method. <br/>
first create an instance of the class you wish to inspect. <br/>
then call the dir method on that instance. What follows is a list of all of the functions of the PDB class <br/>
to get more detailed information call the help methdod. <br/>

```{python}

protein = PDB("pdb1a3q.ent.gz")
dir(protein)

#['FASTA_IDX', 'FASTA_format', 'IDX_format', 'PDB_format', '_COMPND', '_NMR', '_NMR_chains', '__abstractmethods__', '__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__len__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_abc_impl', '_apply_matrix', '_biomol_id', '_chain_id', '_chains', '_cif_file', '_get_chain_position_by_id', '_has_nucl', '_has_prot', '_header', '_pdb_file', '_read_PDB_file', '_read_mmCIF_file', '_write_PDB_file', 'add_chain', 'add_chains', 'all_models', 'apply_biomolecule_matrices', 'apply_symmetry_matrices', 'biomolecule_identifier', 'chain_exists', 'chain_identifiers', 'chains', 'clean', 'dehydrate', 'dump', 'duplicate', 'fuse_chains', 'get_chain_by_id', 'has_nucleotide', 'has_protein', 'header', 'id', 'is_NMR', 'is_all_ca', 'load', 'non_standard_chains', 'nucleotides', 'pdb_file', 'proteins', 'read', 'repeated_chain_ids', 'reposition', 'rotate', 'tmpclean', 'translate', 'write']

help(protein)

```

This is the result:

### PDB
class PDB(SBI.beans.StorableObject.StorableObject) <br/>
	A {PDB} is a collection of {Chain} <br/>

Methods defined here: <br/>

 |  FASTA_IDX(self, protein=True, nucleotide=False) <br/>
 |  <br/>
 |  FASTA_format(self, gapped=True, protein=True, nucleotide=False) <br/>
 |  <br/>
 |  IDX_format(self, protein=True, nucleotide=False) <br/>
 |  <br/>
 |  PDB_format(self, clean=False, terminal=True) <br/>
 |      Strings a {PDB} in PDB format <br/>
 |      @rtype: String <br/>
 |  <br/>
 |  __init__(self, pdb_file=None, dehydrate=False, header=False, onlyheader=False, biomolecule=False) <br/>
 |      @type  pdb_file: String <br/>
 |      @param pdb_file: PDB formated file to read <br/>
 |      <br/>
 |      @raise IOError if pdb_file does not exist and it is not an empty object <br/>
 |  __len__(self) <br/>
 |      # OVERRIDE DEFAULT METHODS <br/>
 |  <br/>
 |  add_chain(self, chain, NMR=False) <br/>
 |      Adds a new chain to the PDB <br/>
 |  <br/>
 |  add_chains(self, chains, NMR=False) <br/>
 |      Adds a new chains to the PDB <br/>
 |  <br/>
 |  apply_biomolecule_matrices(self, keepchains=False, water=True) <br/>
 |      Only works if the PDB file is an original PDB file or <br/>
 |      the matrices have been added in the correct PDB format <br/>
 |      @rtype: {PDB} <br/>
 |  <br/>
 |  apply_symmetry_matrices(self) <br/>
 |      Only works if the PDB file is an original PDB file <br/>
 |      or the matrices have been added in the correct PDB format <br/>
 |      @rtype: {PDB} <br/>
 |  <br/>
 |  chain_exists(self, chain) <br/>
 |      Confirms if a given chain exists in the PDB <br/>
 |      @rtype: Boolean <br/>
 |  <br/>
 |  clean(self) <br/>
 |  <br/>
 |  dehydrate(self) <br/>
 |      # METHODS <br/>
 |  <br/>
 |  duplicate(self, hetero=True, water=False, NMR=False) <br/>
 |      Returns a {PDB} identical to the original but as a new object <br/>
 |      @rtype: {PDB} <br/>
 |  fuse_chains(self, chains_ids) <br/>
 |      Fuses several chains into the first one. <br/>
 |      It will not allow to fuse different structural chains. <br/>
 |      It does not alter the {PDB}, but provides a new one <br/>
 |      @rtype: {Chain} <br/>
 |      <br/>
 |      @raise AttributeError if: <br/>
 |          a) A given chain ID is not present <br/>
 |          b) Try to fuse different structural chains <br/>
 |  <br/>
 |  get_chain_by_id(self, id) <br/>
 |      Returns a chain according to its id or None if no chain with that id is found <br/>
 |      @rtype: {Chain} <br/>
 |  <br/>
 |  reposition(self, matrix=None, vector=None) <br/>
 |      Rotates and Translates each {Chain} according to a matrix and a translational vector <br/>
 |      <br/>
 |      @type matrix: numpy.matrix <br/>
 |      <br/>
 |      @type vector: numpy.array <br/>
 |  <br/>
 |  rotate(self, matrix=None) <br/>
 |      Rotates each {Chain} according to a given matrix <br/>
 |      <br/>
 |      @type matrix: numpy.matrix <br/>
 |  <br/>
 |  tmpclean(self, cluster_by_alternative_id=False) <br/>
 |      Makes a clean version of the PDB, rechaining in order and renumerating atoms. <br/>
 |      Renumbering residues is optional <br/>
 |  <br/>
 |  translate(self, vector=None) <br/>
 |      Translates each {Chain} according to a translational vector <br/>
 |      <br/>
 |      @type vector: numpy.array <br/>
 |  <br/>
 |  write(self, output_file=None, format='PDB', force=False, clean=False) <br/>
 |      Writes the object in a specific format <br/>
 |      <br/>
 |      @type  output_file: String <br/>
 |      @param output_file: File to write <br/>
 |      <br/>
 |      @type  format: String <br/>
 |      @param format: Format of the file to print <br/>
 |  <br/>
 |  ---------------------------------------------------------------------- <br/>
 |  Static methods defined here: <br/>
 |  <br/>
 |  read(input_file, format='PDB') <br/>
 |      Reads a file of data in a specific format and returns the object <br/>
 |      <br/>
 |      @type  input_file: String <br/>
 |      @param input_file: File to read <br/>
 |      <br/>
 |      @type  format: String <br/>
 |      @param format: Format of the file to read <br/>
 |  <br/>
 |  ---------------------------------------------------------------------- <br/>









