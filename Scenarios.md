Scenarios to run the SBI library on


```{python}
from SBI.external.blast import BlastExe
from SBI.structure import PDB 
from SBI.databases import PDBlink
from SBI.structure import Complex
import re
import math
```


1) Blast and print alignment in pir format for modeller
```{python}
protein = PDB("pdb1a3q.ent.gz")
blast = BlastExe(database = "/home/pgohl/Work/pgohl/ModCRE/scripts/databases/sequences.fa")
ChainA = protein.get_chain_by_id("A")
sequence = ChainA.protein_sequence
blastresults = blast.execute_query_seq(sequenceID = "1A3Q_A", sequence = sequence)
##### The alignment that we want to match
PIR = blastresults.str_PIR(result=2)
print(PIR)
```
![Alt text](images/PIR.png?raw=true "Title")



2) Change the residue numbering of a protein to match that of an alignment.
```{python}
PIRlist = PIR.split(">P1;")
p1a3q = "".join(PIRlist[1].split("\n")[2:])
blasthit = "".join(PIRlist[2].split("\n")[2:])
hit_ID = "".join(PIRlist[2].split("\n")[0])
new = PDBlink()
path = new.get_PDB(hit_ID)
proteinmatch = PDB(path)
MatchChainStart = PIRlist[2].split(":")[2]
index = int(MatchChainStart) - 1
for x in range(0 , len(ChainA.aminoacids)):
	while p1a3q[index] == "-" :
		index = index + 1
	ChainA.aminoacids[x].number = index + 1
	index = index + 1
```





3) Calculate a proteins loops in a chain
```{python}
ChainA.calculate_dssp()
ChainA.calculate_archs()
ChainA.archs
ChainA.superarchs
for chain in proteinmatch.chain_identifiers:
	ChainB = proteinmatch.get_chain_by_id(chain)
	if str(type(ChainB)) != "<class 'SBI.structure.chain.ChainOfProtein.ChainOfProtein'>":
		continue
	ChainB.calculate_dssp()
	ChainB.calculate_archs()
	##### Here I will compare loops between ChainA and ChainB using FragRus´s most conservative thresholds
	for loop in ChainA.archs:
		for x in ChainB.archs:
			Rho = math.sqrt((loop.rho - x.rho)**2) <= 10 
			Delta = math.sqrt((loop.delta - x.delta)**2) <= 5
			Theta = math.sqrt((loop.theta - x.theta)**2) <= 5
			Distance = math.sqrt((loop.cartesian_distance - x.cartesian_distance)**2) <= 0.5 
			if Distance and Rho and Delta and Theta:
				print("\n{0}\n{1}\n".format(loop.line_format(),x.line_format()))
```

![Alt text](images/LoopMatches.png?raw=true "Title")

##### The line format is chain identifier \t loop type \t loop length \t loop distance \t theta \t rho \t delta \t internal sec. Str.
##### The order of loops comparisons may be different as the chains are stored in a set.

4) Getting the loop geometry values of all 1A3Q_A_19 BLAST hits 
```{python}
ChainA = protein.get_chain_by_id("A")
tempMatch = "<<<<<<<<" + ChainA.protein_sequence + ">>>>>>>>"
search = r'.{4}' + ChainA.archs[1].aminoacid_sequence + r'.{4}'
blastSearch = re.findall(search, tempMatch)[0].strip("<>")
tempblastresults = blast.execute_query_seq(sequenceID = ChainA.archs[1].identifier, sequence = blastSearch)
Rholist = []
Deltalist = []
Thetalist = []
identities = []
for hit in tempblastresults.get_hits():
	id,chain = hit.sequenceID.split("_")
	path = new.get_PDB(id)
	proteinmatch = PDB(path)
	ChainB = proteinmatch.get_chain_by_id(chain)
	ChainB.calculate_dssp()
	ChainB.calculate_archs()
	for x in ChainB.archs:
		if x.aminoacid_sequence in hit.hit_seq:
			Rho = math.sqrt((ChainA.archs[1].rho - x.rho)**2)
			Delta = math.sqrt((ChainA.archs[1].delta - x.delta)**2)
			Theta = math.sqrt((ChainA.archs[1].theta - x.theta)**2)
			Distance = math.sqrt((ChainA.archs[1].cartesian_distance - x.cartesian_distance)**2)
			Rholist.append(Rho)
			Deltalist.append(Delta)
			Thetalist.append(Theta)
			identities.append(hit.identities_pec)
```

In this example we ran a beta-link: <br/>
![Alt text](images/1A3Q_loop_geom.png?raw=true "Title")




With a repeat using an alpha helix - helix 310 from 5KRB_B: <br/>
![Alt text](images/5KRB_loop_geom.png?raw=true "Title")




5) Find the residues in an interface
```{python}
protein = PDB("pdb1a3q.ent.gz")
complex = Complex(protein)
##### Here I am getting the amount of interfaces that were detected
len(complex.PPInterfaces)
##### This is the interface on the first chain in the interaction
for residue in complex.PPInterfaces[0]._view_interface_from(1):
	print(str(residue.number) + " " + residue._type)
##### This is the interface on the second chain in the interaction
for residue in complex.PPInterfaces[0]._view_interface_from(2):
	print(str(residue.number) + " " + residue._type)
##### The above prints out the name and number of each residue if you want to get the actual residue object itself you simply call residue. 
```

6) Interface summary statistics
```{python}
## Finding the percentage of contact residues that are hydrophobic

# Define hydrophobic Amino Acids
hydrophobic = ['PHE', 'LEU', 'ILE', 'TYR', 'TRP', 'VAL', 'MET', 'PRO']
i = 0
x = 0
for contact in complex.PPInterfaces[0].contacts:
	if contact.aminoacid1.standard_type in hydrophobic:
		i+=1
	if contact.aminoacid2.standard_type in hydrophobic:
		x+=1
print("Contact residue hydrophobicity is:\nChain1:{0}% Chain2:{1}% ".format(i/len(complex.PPInterfaces[0].contacts),x/len(complex.PPInterfaces[0].contacts)))


## Percentage of contact residues that are exposed

# First the dssp for both chains in the complex interface must be calculated
complex.pdb.get_chain_by_id("A").calculate_dssp()
complex.pdb.get_chain_by_id("B").calculate_dssp()
i = 0
x = 0
for contact in complex.PPInterfaces[0].contacts:
	if contact.aminoacid1.accessibility:
		i+=1
	if contact.aminoacid2.accessibility:
		x+=1
print("Contact residue exposure is:\nChain1:{0}% Chain2:{1}% ".format(i/len(complex.PPInterfaces[0].contacts),x/len(complex.PPInterfaces[0].contacts)))

```

7) Read in a mmcif file and blast getting results that pass the rost curve
```{python}
protein = PDB(cif_file = "1a3q.cif.gz")
ChainA = protein.get_chain_by_id("A")
sequence = ChainA.protein_sequence
blastresults = blast.execute_query_seq(sequenceID = "1A3Q_A", sequence = sequence)
for hit in blastresults.get_hits(tz_type = "ID"):
	print(hit.sequenceID, “ “, hit.evaluate_Rost_twilight_zone() )
```
![Alt text](images/Rost_output.png?raw=true "Title")
















