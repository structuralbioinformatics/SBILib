import urllib.request

def download_pdb(uniprot_accession):
	url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_accession}-F1-model_v4.pdb"
	urllib.request.urlretrieve(url, f"AlphaPrediction_{uniprot_accession}.pdb")
	return(f"AlphaPrediction_{uniprot_accession}.pdb")

