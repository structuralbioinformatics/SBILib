import unittest, subprocess, gzip, os, shutil

from SBILib.structure import PDB, Complex 
from SBILib.databases import PDBlink
from SBILib.external.blast import BlastExe





class Test_structure(unittest.TestCase):
    def setUp(self):
        self.protein = PDB("Test_pdb1a3q.ent.gz")

    def test_download(self):
        new = PDBlink()
        path = new.get_PDB("1a3q")
        self.assertIsNot(path,False)

    def test_mmcif(self):
        protein = PDB(cif_file = "Test_1a3q.cif.gz")
        self.assertEqual(protein.chain_identifiers, {'A', 'C', 'B', 'D'})
        self.assertEqual(protein.FASTA_format(), '>1A3Q_A\t37 \nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRxxxxxxSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP\n>1A3Q_B\t37 \nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRxxxxxxSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP\n')

    def test_manipulatingProtein(self):
        self.assertEqual(self.protein.chain_identifiers, {'A', 'C', 'B', 'D'})
        protein2 = self.protein.duplicate()
        protein2 = protein2.fuse_chains(["A","B"])
        self.assertEqual(protein2.chain_identifiers , {"A"})
        ChainA = protein2.get_chain_by_id("A")
        ChainE = ChainA.duplicate()
        ChainE.chain = "E"
        protein2 = self.protein.duplicate()
        protein2.add_chain(ChainE)
        self.assertEqual(protein2.chain_identifiers, {'A', 'C', 'B', 'D', 'E'})
        protein3 = PDB()
        deletion = ["E"]
        for chain in protein2.chain_identifiers:
            if chain not in deletion:
                protein3.add_chain(protein2.get_chain_by_id(chain))
                if (protein2.get_chain_by_id(chain).chaintype == "P"):
                    protein3._has_prot = "TRUE"
                elif (protein2.get_chain_by_id(chain).chaintype == "N"):
                    protein3._has_nucl = "TRUE"
        
        self.assertEqual(protein3.chain_identifiers, self.protein.chain_identifiers)


    def test_reading_sequences(self):
        self.maxDiff = None
        ChainA = self.protein.get_chain_by_id("A")
        self.assertEqual(self.protein.FASTA_format(gapped=False), '>TEST_PDB1A3Q_A\t37 \nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP\n>TEST_PDB1A3Q_B\t37 \nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP\n' )
        self.assertEqual(self.protein.FASTA_format(), '>TEST_PDB1A3Q_A\t37 \nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRxxxxxxSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP\n>TEST_PDB1A3Q_B\t37 \nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRxxxxxxSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP\n' )
        self.assertEqual(ChainA.gapped_protein_sequence, 'GPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRxxxxxxSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP')
        self.assertEqual(ChainA.protein_sequence, 'GPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRSLPLKPVISQPIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP' )
        self.assertTrue(self.protein.has_nucleotide)
        first = next(self.protein.nucleotides)
        self.assertEqual(first.globalID, 'TEST_PDB1A3Q_C' )
        self.assertEqual(first.nucleotide_sequence(), 'GGGGAATCCCC' )

class Test_Blast(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.blastresults = None

    def setUp(self):
        url = "ftp://ftp.ebi.ac.uk/pub/databases/uniprot/knowledgebase/uniprot_sprot.fasta.gz"
        output = subprocess.run(["wget", url], capture_output=True, text=True)
        with gzip.open('uniprot_sprot.fasta.gz', 'rb') as f:
            with open('uniprot_sprot.fasta','wb') as f_out:
                shutil.copyfileobj(f,f_out)
        os.remove("uniprot_sprot.fasta.gz")
        self.protein = PDB("Test_pdb1a3q.ent.gz")
        self.chainA = self.protein.get_chain_by_id("A")

    def test_Compile_Blast(self):
        querySequence = self.chainA.protein_sequence
        queryID = self.chainA.globalID
        try:
            blast = BlastExe(database = "uniprot_sprot.fasta")
        except Exception as e:
            self.fail("BlastExe() raised ExceptionType unexpectedly!")
        try:
            self.__class__.blastresults = blast.execute_query_seq(sequenceID = queryID, sequence = querySequence)
        except Exception as e:
            self.fail("blast.execute_query_seq() raised ExceptionType unexpectedly!")
        first = self.__class__.blastresults.get_hits()[0]
        self.assertEqual(first.sequenceID, 'sp|Q00653|NFKB2_HUMAN')
        self.assertTrue(first.evaluate_Rost_twilight_zone())
    
    def test_PIR(self):
        self.maxDiff = None
        self.assertEqual(self.__class__.blastresults.str_PIR(result = 0), '>P1;TEST_PDB1A3Q_A\nsequence:TEST_PDB1A3Q_A:1:.:285:.:.:.:.:.\nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLV\nTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQR\nQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLR------SLPLKPVISQP\nIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDF\nSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP*\n>P1;NFKB2\nstructureX:NFKB2:37:HUMAN:327:HUMAN:.:.:.:.\nGPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLV\nTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQR\nQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRASDGSFSLPLKPVISQP\nIHDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDF\nSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP*'  )
        self.assertEqual(self.__class__.blastresults.get_hits( evalue = 0, tz_type = "ID")[0]._num_seq, 2)
        

class Test_Complex(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.Complex = None
    def setUp(self):
        self.protein = PDB("Test_pdb1a3q.ent.gz")
        self.chainA = self.protein.get_chain_by_id("A")
        

    def test_Complex_calc(self):
        try:
            self.__class__.Complex = Complex(self.protein)
        except Exception as e:
            self.fail("Complex() raised ExceptionType unexpectedly!")
        
    def test_Interactions(self):
        self.assertEqual(len(self.__class__.Complex.PPInterfaces[0]), 100 )
        self.assertEqual( len(self.__class__.Complex.PNInterfaces), 4)
        for chain in self.protein.chains:
            self.assertFalse(chain.has_heteroatoms)
        Hprotein = PDB("Test_pdb1a2w.ent.gz")
        try:
            HComplex = Complex(Hprotein)
        except Exception as e:
            self.fail("Complex() raised ExceptionType unexpectedly!")
        for chain in Hprotein.chains:
            self.assertTrue(chain.has_heteroatoms)



class Test_Loops(unittest.TestCase):
    def setUp(self):
        self.protein = PDB("Test_pdb1a3q.ent.gz")
        self.ChainA = self.protein.get_chain_by_id("A")

    def test_Calcs(self):
        try:
            self.ChainA.calculate_dssp()
        except Exception as e:
            self.fail("calculate_dssp() raised ExceptionType unexpectedly!")
        try:
            self.ChainA.calculate_archs()
        except Exception as e:
            self.fail("calculate_archs() raised ExceptionType unexpectedly!")

        self.assertEqual(self.ChainA.archs[0].aminoacid_sequence, 'YLVIVEQPKQRGFRFRY' )
        self.assertEqual(self.ChainA.archs[0].structure_sequence, 'EEEEEE-B-SSSB--EE' )



class Test_Grafting(unittest.TestCase):
    def setUp(self):
        self.protein = PDB("Test_pdb1a3q.ent.gz")
        self.ChainA = self.protein.get_chain_by_id("A")
        self.protein2 = PDB("Test_pdb2ram.ent.gz")

    def test_compare_loops(self):

        try:
            report = self.protein.compare_loops(self.protein2)
        except Exception as e:
            self.fail("compare_loops() raised ExceptionType unexpectedly!")
        
        self.assertEqual(report[0], 'QueryProt:A,176,182:TargetProt:A,160,166:VISQPIH,VLSHPIF')

    def test_graft(self):
        self.maxDiff = None
        graft = self.protein.graft(self.protein2, 'QueryProt:A,176,182:TargetProt:A,160,166:VISQPIH,VLSHPIF')
        self.assertEqual(graft.get_chain_by_id("A").protein_sequence, 'GPYLVIVEQPKQRGFRFRYGCEGPSHGGLPGASSEKGRKTYPTVKICNYEGPAKIEVDLVTHSDPPRAHAHSLVGKQCSELGICAVSVGPKDMTAQFNNLGVLHVTKKNMMGTMIQKLQRQRLRSRPQGLTEAEQRELEQEAKELKKVMDLSIVRLRFSAFLRSLPLKPVLSHPIFDSKSPGASNLKISRMDKTAGSVRGGDEVYLLCDKVQKDDIEVRFYEDDENGWQAFGDFSPTDVHKQYAIVFRTPPYHKMKIERPVTVFLQLKRKRGGDVSDSKQFTYYP')

        




if __name__ == '__main__':
    unittest.main()

