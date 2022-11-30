"""
#
# PROTEIN
#
"""

"""
CODING TRANSFORMATION
http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/
"""

aminoacids3to1 = {
    'ALA':"A", '02K': "A", '02O': "A", '02Y': "A", '0CS': "A", '0FL': "A", '0NC': "A", '175': "A", '1AC': "A",
               '1L1': "A", '1PI': "A", '23P': "A", '2AG': "A", '2RA': "A", '3GA': "A", '4U7': "A", '5AB': "A",
               '5OH': "A", 'AA3': "A", 'AA4': "A", 'ABA': "A", 'AHO': "A", 'AHP': "A", 'AIB': "A", 'ALC': "A",
               'ALM': "A", 'ALN': "A", 'ALS': "A", 'ALT': "A", 'ALV': "A", 'AN8': "A", 'AYA': "A", 'AYG': "A",
               'AZH': "A", 'B2A': "A", 'B3A': "A", 'BP5': "A", 'C22': "A", 'CLV': "A", 'CRW': "A", 'CRX': "A",
               'CWD': "A", 'DAB': "A", 'DBZ': "A", 'DDZ': "A", 'DNP': "A", 'DPP': "A", 'FB5': "A", 'FB6': "A",
               'FLA': "A", 'HAC': "A", 'HHK': "A", 'HQA': "A", 'HV5': "A", 'IAM': "A", 'KNB': "A", 'LAL': "A",
               'MA': "A", 'MAA': "A", 'MDO': "A", 'NA8': "A", 'NAL': "A", 'NAM': "A", 'NCB': "A", 'NPI': "A",
               'NWD': "A", 'ONH': "A", 'ORN': "A", 'PIA': "A", 'PYA': "A", 'S2P': "A", 'SE7': "A", 'SEG': "A",
               'TIH': "A", 'UM1': "A", 'UM2': "A", 'UMA': "A", 'X9Q': "A", 'XW1': "A", 'Z01': "A", 'ZZJ': "A",
    'ARG':"R", '0AR': "R", '0X9': "R", '2MR': "R", '4AR': "R", '9NR': "R", 'AAR': "R", 'ACL': "R", 'AGM': "R",
               'AHL': "R", 'AR2': "R", 'AR7': "R", 'ARM': "R", 'ARO': "R", 'BOR': "R", 'CIR': "R", 'DA2': "R",
               'DIR': "R", 'DP1': "R", 'FIO': "R", 'HAR': "R", 'HMR': "R", 'HR7': "R", 'HRG': "R", 'IAR': "R",
               'MAI': "R", 'MGG': "R", 'MMO': "R", 'MYN': "R", 'NMM': "R", 'NNH': "R", 'OAR': "R", 'ORQ': "R",
               'RGL': "R", 'THZ': "R", 'VR0': "R",
    'ASN':"N", '02L': "N", '9DN': "N", 'A5N': "N", 'AFA': "N", 'AHB': "N", 'AS7': "N", 'B3X': "N", 'DMH': "N",
               'MEN': "N", 'NYG': "N", 'SD4': "N", 'SNN': "N", 'XSN': "N",
    'ASP':"D", '0A0': "D", '0AK': "D", '0TD': "D", '3MD': "D", 'AEI': "D", 'AKL': "D", 'AKZ': "D", 'ASA': "D",
               'ASB': "D", 'ASI': "D", 'ASK': "D", 'ASL': "D", 'ASQ': "D", 'B3D': "D", 'BFD': "D", 'BH2': "D",
               'BHD': "D", 'DMK': "D", 'DOH': "D", 'DYG': "D", 'LAA': "D", 'OHS': "D", 'OXX': "D", 'PAS': "D",
               'PHD': "D", 'SIC': "D", 'SUI': "D", 'TAV': "D", 'XYG': "D",
    'ASX':"B",
    'CYS':"C", '03Y': "C", '07O': "C", '0A8': "C", '0QL': "C", '143': "C", '2CO': "C", '2XA': "C", '5CS': "C",
               'AGT': "C", 'BB6': "C", 'BB7': "C", 'BB9': "C", 'BBC': "C", 'BCS': "C", 'BCX': "C", 'BPE': "C",
               'BTC': "C", 'BUC': "C", 'C1T': "C", 'C3Y': "C", 'C4R': "C", 'C5C': "C", 'C6C': "C", 'CAF': "C",
               'CAS': "C", 'CAY': "C", 'CCS': "C", 'CCY': "C", 'CEA': "C", 'CFY': "C", 'CME': "C", 'CMH': "C",
               'CML': "C", 'CMT': "C", 'CS0': "C", 'CS1': "C", 'CS3': "C", 'CS4': "C", 'CSA': "C", 'CSB': "C",
               'CSD': "C", 'CSE': "C", 'CSJ': "C", 'CSO': "C", 'CSP': "C", 'CSR': "C", 'CSS': "C", 'CSU': "C",
               'CSW': "C", 'CSX': "C", 'CSZ': "C", 'CY0': "C", 'CY1': "C", 'CY3': "C", 'CY4': "C", 'CYA': "C",
               'CYD': "C", 'CYF': "C", 'CYG': "C", 'CYM': "C", 'CYQ': "C", 'CYR': "C", 'CYW': "C", 'CZ2': "C",
               'CZZ': "C", 'DC2': "C", 'DYS': "C", 'ECX': "C", 'EFC': "C", 'FFM': "C", 'FOE': "C", 'GT9': "C",
               'GYC': "C", 'HCM': "C", 'HNC': "C", 'HTI': "C", 'ICY': "C", 'JJJ': "C", 'JJK': "C", 'JJL': "C",
               'K1R': "C", 'M0H': "C", 'MCS': "C", 'MD3': "C", 'MD5': "C", 'NPH': "C", 'NYB': "C", 'NYS': "C",
               'OCS': "C", 'OCY': "C", 'P1L': "C", 'PBB': "C", 'PEC': "C", 'PR3': "C", 'PYX': "C", 'QCS': "C",
               'QPA': "C", 'R1A': "C", 'S2C': "C", 'SAH': "C", 'SCH': "C", 'SCS': "C", 'SCY': "C", 'SHC': "C",
               'SIB': "C", 'SIC': "C", 'SMC': "C", 'SNC': "C", 'SOC': "C", 'SYS': "C", 'TNB': "C", 'TQZ': "C",
               'TSY': "C", 'XCN': "C", 'YCM': "C", 'ZZD': "C",
    'GLN':"Q", 'CRQ': "Q", 'ECC': "Q", 'EPQ': "Q", 'GHG': "Q", 'GLH': "Q", 'GNC': "Q", 'LMQ': "Q", 'MEQ': "Q",
               'MGN': "Q", 'NLQ': "Q", 'QFG': "Q", 'QLG': "Q", 'QMM': "Q",
    'GLU':"E", '11W': "E", '3GL': "E", '3O3': "E", '5HP': "E", '9NE': "E", 'AR4': "E", 'B3E': "E", 'CGA': "E",
               'CGU': "E", 'CRU': "E", 'EME': "E", 'G01': "E", 'G8M': "E", 'GAU': "E", 'GHC': "E", 'GHW': "E",
               'GLJ': "E", 'GLK': "E", 'GLQ': "E", 'GMA': "E", 'GME': "E", 'GSU': "E", 'ILG': "E", 'LME': "E",
               'MEG': "E", 'MFN': "E", 'PCA': "E", 'X2W': "E",
    'GLX':"Z",
    'GLY':"G", '0AC': "G", '0YG': "G", '175': "G", '3XH': "G", '4F3': "G", '5PG': "G", '5ZA': "G", '9DS': "G",
               'AYG': "G", 'C12': "G", 'C99': "G", 'CCY': "G", 'CFY': "G", 'CH6': "G", 'CH7': "G", 'CHP': "G",
               'CJO': "G", 'CLV': "G", 'CQ1': "G", 'CQ2': "G", 'CQR': "G", 'CR0': "G", 'CR2': "G", 'CR5': "G",
               'CR7': "G", 'CR8': "G", 'CRF': "G", 'CRG': "G", 'CRK': "G", 'CRO': "G", 'CRQ': "G", 'CRU': "G",
               'CRW': "G", 'CRX': "G", 'CSH': "G", 'CSY': "G", 'CZO': "G", 'DYG': "G", 'EYG': "G", 'FGL': "G",
               'GEE': "G", 'GL3': "G", 'GLZ': "G", 'GSC': "G", 'GYC': "G", 'GYS': "G", 'IEY': "G", 'IGL': "G",
               'IIC': "G", 'IPG': "G", 'KWS': "G", 'LPG': "G", 'LVG': "G", 'M30': "G", 'MD6': "G", 'MDO': "G",
               'MEU': "G", 'MFC': "G", 'MGY': "G", 'MPQ': "G", 'MSA': "G", 'NLY': "G", 'NMC': "G", 'NRP': "G",
               'NRQ': "G", 'NYC': "G", 'NYG': "G", 'PGY': "G", 'PIA': "G", 'PRV': "G", 'QFG': "G", 'QLG': "G",
               'RC7': "G", 'SAR': "G", 'SHP': "G", 'SUI': "G", 'SWG': "G", 'UGY': "G", 'WCR': "G", 'X9Q': "G",
               'XXY': "G", 'XYG': "G",
    'HIS':"H", '2HF': "H", '2SO': "H", '3AH': "H", '56A': "H", 'B3U': "H", 'CR8': "H", 'CRG': "H", 'CSH': "H",
               'DDE': "H", 'HBN': "H", 'HHI': "H", 'HIA': "H", 'HIC': "H", 'HIP': "H", 'HIQ': "H", 'HS8': "H",
               'HS9': "H", 'HSO': "H", 'IEY': "H", 'IIC': "H", 'MH1': "H", 'MHS': "H", 'NEM': "H", 'NEP': "H",
               'NZH': "H", 'OHI': "H", 'OLD': "H", 'PSH': "H", 'PVH': "H", 'RC7': "H", 'XXY': "H", 'Z70': "H", 
    'ILE':"I", '7JA': "I", 'B2I': "I", 'BIU': "I", 'I2M': "I", 'IIL': "I", 'ILX': "I", 'IML': "I", 'QIL': "I",
               'TS9': "I",
    'LEU':"L", '0AG': "L", '2LU': "L", '2ML': "L", 'AN6': "L", 'BL2': "L", 'BLE': "L", 'BTA': "L", 'CLE': "L",
               'CR0': "L", 'DON': "L", 'EXY': "L", 'FLE': "L", 'HL2': "L", 'HLU': "L", 'L3O': "L", 'LED': "L",
               'LEF': "L", 'LEH': "L", 'LEM': "L", 'LEN': "L", 'LEX': "L", 'LNE': "L", 'LNM': "L", 'MHL': "L",
               'MK8': "L", 'MLE': "L", 'MLL': "L", 'MNL': "L", 'NLE': "L", 'NLN': "L", 'NLO': "L", 'NLP': "L",
               'NOT': "L", 'NRP': "L", 'QLG': "L", 'WLU': "L",
    'LYS':"K", '0A2': "K", '3QN': "K", '6CL': "K", '6HN': "K", 'ALY': "K", 'API': "K", 'APK': "K", 'AZK': "K",
               'B3K': "K", 'BLY': "K", 'BTK': "K", 'C1X': "K", 'CH7': "K", 'CLG': "K", 'CLH': "K", 'CR7': "K",
               'CYJ': "K", 'DLS': "K", 'DM0': "K", 'DNL': "K", 'DNS': "K", 'ELY': "K", 'FDL': "K", 'FH7': "K",
               'FHL': "K", 'FHO': "K", 'FZN': "K", 'GPL': "K", 'HHK': "K", 'I58': "K", 'IEL': "K", 'ILY': "K",
               'IT1': "K", 'KBE': "K", 'KCX': "K", 'KFP': "K", 'KGC': "K", 'KPI': "K", 'KPY': "K", 'KST': "K",
               'KYQ': "K", 'LA2': "K", 'LBY': "K", 'LCK': "K", 'LCX': "K", 'LDH': "K", 'LET': "K", 'LGY': "K",
               'LLO': "K", 'LLP': "K", 'LLY': "K", 'LLZ': "K", 'LMF': "K", 'LP6': "K", 'LSO': "K", 'LYF': "K",
               'LYK': "K", 'LYM': "K", 'LYN': "K", 'LYP': "K", 'LYR': "K", 'LYU': "K", 'LYX': "K", 'LYZ': "K",
               'M2L': "K", 'M3L': "K", 'M3R': "K", 'MCL': "K", 'ML3': "K", 'MLY': "K", 'MLZ': "K", 'MYK': "K",
               'OBS': "K", 'PE1': "K", 'PRK': "K", 'PYH': "K", 'SHR': "K", 'SLL': "K", 'SLZ': "K", 'TLY': "K",
               'TRG': "K", 'VB1': "K", 'XX1': "K",
    'MET':"M", '2FM': "M", '4CY': "M", 'CH6': "M", 'CRK': "M", 'CXM': "M", 'EPM': "M", 'ESC': "M", 'EYG': "M",
               'FME': "M", 'H1D': "M", 'HYI': "M", 'IZO': "M", 'KOR': "M", 'M2S': "M", 'ME0': "M", 'MHO': "M",
               'MME': "M", 'MSE': "M", 'MSL': "M", 'MSO': "M", 'MT2': "M", 'NRQ': "M", 'OMT': "M", 'SME': "M", 
    'PHE':"F", '0A9': "F", '0BN': "F", '1PA': "F", '200': "F", '23F': "F", '3CF': "F", '4AF': "F", '4CF': "F",
               '4PH': "F", '9NF': "F", 'B1F': "F", 'B2F': "F", 'BB8': "F", 'BIF': "F", 'BNN': "F", 'C99': "F",
               'CFY': "F", 'CLV': "F", 'DAH': "F", 'EHP': "F", 'F2F': "F", 'FC0': "F", 'FCL': "F", 'H14': "F",
               'HOX': "F", 'HPC': "F", 'HPE': "F", 'HPH': "F", 'HPQ': "F", 'MEA': "F", 'MHU': "F", 'NFA': "F",
               'PBF': "F", 'PCS': "F", 'PF5': "F", 'PFF': "F", 'PHA': "F", 'PHI': "F", 'PHL': "F", 'PHM': "F",
               'PM3': "F", 'PPN': "F", 'PSA': "F", 'QFG': "F", 'QPH': "F", 'SMF': "F", 'T11': "F", 'TEF': "F",
               'TFQ': "F", 'U3X': "F", 'WFP': "F", 'WPA': "F", 'X9Q': "F", 'ZCL': "F",
    'PRO':"P", '037': "P", '04U': "P", '04V': "P", '05N': "P", '0LF': "P", '0Y8': "P", '11Q': "P", '12L': "P",
               '12X': "P", '12Y': "P", '2MT': "P", '2P0': "P", '3PX': "P", '4FB': "P", 'DPL': "P", 'FP9': "P",
               'FPK': "P", 'H5M': "P", 'HY3': "P", 'HYP': "P", 'HZP': "P", 'LPD': "P", 'LWY': "P", 'MP8': "P",
               'N7P': "P", 'P2Y': "P", 'PCC': "P", 'PKR': "P", 'PLJ': "P", 'POM': "P", 'PR4': "P", 'PR7': "P",
               'PR9': "P", 'PRJ': "P", 'PRS': "P", 'PXU': "P", 'RT0': "P", 'TPJ': "P", 'TPK': "P", 'VH0': "P",
               'XPR': "P", 'ZYJ': "P", 'ZYK': "P",
    'PYL':"O", 
    'SEC':"U", 'PSW': "U", 'UOX': "U",
    'SER':"S", '0AH': "S", '175': "S", '1X6': "S", '8SP': "S", 'A9D': "S", 'AZS': "S", 'BG1': "S", 'BSE': "S",
               'CRW': "S", 'CRX': "S", 'CSH': "S", 'CSY': "S", 'CWR': "S", 'DBS': "S", 'DHA': "S", 'FGP': "S",
               'GFT': "S", 'GVL': "S", 'GYS': "S", 'HSE': "S", 'HSL': "S", 'IIC': "S", 'KWS': "S", 'LPS': "S",
               'MC1': "S", 'MDO': "S", 'MH6': "S", 'MIR': "S", 'MIS': "S", 'N10': "S", 'NC1': "S", 'OAS': "S",
               'OLZ': "S", 'OMH': "S", 'OSE': "S", 'PG1': "S", 'RVX': "S", 'RZ4': "S", 'S12': "S", 'S1H': "S",
               'SAC': "S", 'SBG': "S", 'SBL': "S", 'SDB': "S", 'SDP': "S", 'SEB': "S", 'SEE': "S", 'SEL': "S",
               'SEM': "S", 'SEN': "S", 'SEP': "S", 'SET': "S", 'SGB': "S", 'SOY': "S", 'SRZ': "S", 'SUN': "S",
               'SVA': "S", 'SVV': "S", 'SVW': "S", 'SVX': "S", 'SVY': "S", 'SVZ': "S", 'SWG': "S", 'SXE': "S",
               'TIS': "S", 'TNR': "S", 'UF0': "S",
    'THR':"T", '0E5': "T", '26B': "T", '28X': "T", '5ZA': "T", 'ALO': "T", 'B27': "T", 'BMT': "T", 'C12': "T",
               'C99': "T", 'CQ1': "T", 'CR0': "T", 'CRF': "T", 'CRG': "T", 'CTH': "T", 'DBU': "T", 'KWS': "T",
               'NYC': "T", 'OLT': "T", 'OTH': "T", 'TBM': "T", 'TH5': "T", 'TH6': "T", 'THC': "T", 'TMB': "T",
               'TMD': "T", 'TNY': "T", 'TPO': "T", 'XXY': "T", 'Z3E': "T", 'ZU0': "T",

    'TRP':"W", '0AF': "W", '0UO': "W", '1TQ': "W", '4AW': "W", '4DP': "W", '4FW': "W", '4HT': "W", '4IN': "W",
               '5ZA': "W", '6CW': "W", 'BTR': "W", 'CRF': "W", 'CTE': "W", 'FT6': "W", 'FTR': "W", 'HRP': "W",
               'HT7': "W", 'HTR': "W", 'KYN': "W", 'LTR': "W", 'NYC': "W", 'PAT': "W", 'R4K': "W", 'RE0': "W",
               'RE3': "W", 'SWG': "W", 'TCR': "W", 'TOQ': "W", 'TOX': "W", 'TPL': "W", 'TQI': "W", 'TQQ': "W",
               'TRF': "W", 'TRN': "W", 'TRO': "W", 'TRQ': "W", 'TRW': "W", 'TRX': "W", 'TRY': "W", 'TTQ': "W",
               'WRP': "W",
    'TYR':"Y", '0A1': "Y", '0EA': "Y", '0WZ': "Y", '0YG': "Y", '1OP': "Y", '1TY': "Y", '2TY': "Y", '3MY': "Y",
               '3NF': "Y", '3YM': "Y", '4BF': "Y", '4F3': "Y", '4HL': "Y", 'AGQ': "Y", 'AYG': "Y", 'AZY': "Y",
               'B3Y': "Y", 'C12': "Y", 'CCY': "Y", 'CFY': "Y", 'CH6': "Y", 'CH7': "Y", 'CJO': "Y", 'CQ1': "Y",
               'CQ2': "Y", 'CQR': "Y", 'CR2': "Y", 'CR7': "Y", 'CR8': "Y", 'CRK': "Y", 'CRO': "Y", 'CRQ': "Y",
               'CRU': "Y", 'CSY': "Y", 'CZO': "Y", 'DBY': "Y", 'DPQ': "Y", 'DYG': "Y", 'ESB': "Y", 'EYG': "Y",
               'FLT': "Y", 'FTY': "Y", 'GYC': "Y", 'GYS': "Y", 'IEY': "Y", 'IYR': "Y", 'MBQ': "Y", 'MDF': "Y",
               'MFC': "Y", 'MTY': "Y", 'NBQ': "Y", 'NIY': "Y", 'NRP': "Y", 'NRQ': "Y", 'NTR': "Y", 'NTY': "Y",
               'NYG': "Y", 'OMX': "Y", 'OMY': "Y", 'P2Q': "Y", 'P3Q': "Y", 'PAQ': "Y", 'PIA': "Y", 'PTH': "Y",
               'PTM': "Y", 'PTR': "Y", 'RC7': "Y", 'STY': "Y", 'T0I': "Y", 'TCQ': "Y", 'TPQ': "Y", 'TTS': "Y",
               'TY1': "Y", 'TY2': "Y", 'TY3': "Y", 'TY5': "Y", 'TY8': "Y", 'TY9': "Y", 'TYB': "Y", 'TYI': "Y",
               'TYJ': "Y", 'TYN': "Y", 'TYO': "Y", 'TYQ': "Y", 'TYS': "Y", 'TYT': "Y", 'TYW': "Y", 'TYY': "Y",
               'U2X': "Y", 'WCR': "Y", 'XYG': "Y", 'YOF': "Y", 'YPZ': "Y",
    'VAL':"V", '033': "V", '0AA': "V", '0AB': "V", '2VA': "V", '9NV': "V", 'A8E': "V", 'B2V': "V", 'BUG': "V",
               'DHN': "V", 'FVA': "V", 'HVA': "V", 'LE1': "V", 'LVN': "V", 'MNV': "V", 'MVA': "V", 'NVA': "V",
               'TBG': "V", 'VAD': "V", 'VAF': "V", 'VAH': "V", 'VAI': "V", 'WVL': "V",

    'XLE':"J",
    'UNK': "X", 'ACE': "X", '3FG': "X"
}

aminoacids1to3 = dict([[v,k] for k,v in list(aminoacids3to1.items())])
aminoacids1to3['A'] = "ALA"
aminoacids1to3['N'] = "ASN"
aminoacids1to3['R'] = "ARG"
aminoacids1to3['D'] = "ASP"
aminoacids1to3['C'] = "CYS"
aminoacids1to3['Q'] = "GLN"
aminoacids1to3['E'] = "GLU"
aminoacids1to3['G'] = "GLY"
aminoacids1to3['H'] = "HIS"
aminoacids1to3['I'] = "ILE"
aminoacids1to3['J'] = "XLE"
aminoacids1to3['L'] = "LEU"
aminoacids1to3['K'] = "LYS"
aminoacids1to3['M'] = "MET"
aminoacids1to3['F'] = "PHE"
aminoacids1to3['O'] = "PYL"
aminoacids1to3['P'] = "PRO"
aminoacids1to3['S'] = "SER"
aminoacids1to3['T'] = "THR"
aminoacids1to3['U'] = "SEC"
aminoacids1to3['W'] = "TRP"
aminoacids1to3['Y'] = "TYR"
aminoacids1to3['V'] = "VAL"

"""
REGULAR AMINOACIDS IDENTIFICATION
"""
aminoacids_main3 = set(['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLX','GLY','HIS','ILE', 'LEU',
                        'LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','SEC','PYL','XLE'])

aminoacids_main1 = set(['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'])

"""
PROPERTIES
"""
aminoacids_surface = {
    'A': 115, 'C': 149, 'D': 170, 'E': 207, 'F': 230, 'G': 86,  'H': 206,
    'I': 187, 'K': 222, 'L': 192, 'M': 210, 'N': 184, 'P': 140, 'Q': 208,
    'R': 263, 'S': 140, 'T': 164, 'V': 161, 'W': 269, 'Y': 257,
}

aminoacids_polarity_boolean = {
    'A': False, 'C': False, 'D': True,  'E': True,  'F': False, 'G': False, 'H': True,
    'I': False, 'K': True,  'L': False, 'M': False, 'N': True,  'P': False, 'Q': True,
    'R': True,  'S': True,  'T': True,  'V': False, 'W': False, 'Y': True
}

aminoacids_acceptors = {
    'N' : ['OD1', 'OD1'],
    'D' : ['OD1', 'OD1', 'OD2', 'OD2'],
    'Q' : ['OE1', 'OE1'],
    'E' : ['OE1', 'OE1', 'OE2', 'OE2'],
    'H' : ['ND1', 'NE1'],
    'S' : ['OG' , 'OG'],
    'T' : ['OG1', 'OG1'],
    'Y' : ['OH']
}

aminoacids_donors = {
    'R' : ['NE' , 'NH1', 'NH1', 'NH2', 'NH2'],
    'N' : ['ND2', 'ND2'],
    'Q' : ['NE2', 'NE2'],
    'H' : ['ND1', 'NE1'],
    'K' : ['NZ' , 'NZ',  'NZ'],
    'S' : ['OG'],
    'T' : ['OG1'],
    'W' : ['NE1'],
    'Y' : ['OH']
}

"""
#
# DNA/RNA
#
"""

"""
CODING TRANSFORMATION
http://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/
"""
nucleic3to1  = {
    'A': "A", 'DA': "A", 'ADE': "A", '+A': "A",
              '00A': "A", '12A': "A", '1MA': "A", '26A': "A", '2MA': "A", '5FA': "A", '6IA': "A", '6MA': "A", '6MC': "A",
              '6MP': "A", '6MT': "A", '6MZ': "A", '8AN': "A", 'A23': "A", 'A2L': "A", 'A2M': "A", 'A39': "A", 'A3P': "A",
              'A44': "A", 'A5O': "A", 'A6A': "A", 'A9Z': "A", 'ADI': "A", 'ADP': "A", 'AET': "A", 'AMD': "A", 'AMO': "A",
              'AP7': "A", 'AVC': "A", 'G3A': "A", 'LCA': "A", 'MA6': "A", 'MAD': "A", 'MGQ': "A", 'MIA': "A", 'MTU': "A",
              'N79': "A", 'P5P': "A", 'PPU': "A", 'PR5': "A", 'PU': "A", 'RIA': "A", 'SRA': "A", 'T6A': "A", 'TBN': "A",
              'TXD': "A", 'TXP': "A", 'V3L': "A", 'ZAD': "A", '0AM': "A", '0AV': "A", '0SP': "A", '1AP': "A", '2AR': "A",
              '2BU': "A", '2A': "A", '3A': "A", '5AA': "A", '6HA': "A", '7A': "A", '8BA': "A", 'A34': "A", 'A35': "A",
              'A38': "A", 'A3A': "A", 'A40': "A", 'A43': "A", 'A47': "A", 'A5L': "A", 'ABR': "A", 'ABS': "A", 'AD2': "A",
              'AF2': "A", 'AS': "A", 'DZM': "A", 'E': "A", 'E1X': "A", 'EA': "A", 'FA2': "A", 'MA7': "A", 'PRN': "A",
              'R': "A", 'RMP': "A", 'S4A': "A", 'SMP': "A", 'TCY': "A", 'TFO': "A", 'XAD': "A", 'XAL': "A", 'XUA': "A",
              'Y': "A",
    'C': "C", 'DC': "C", 'CYT': "C", '+C': "C",
              '10C': "C", '1SC': "C", '4OC': "C", '5IC': "C", '5MC': "C", 'A5M': "C", 'A6C': "C", 'C25': "C", 'C2L': "C",
              'C31': "C", 'C43': "C", 'C5L': "C", 'CBV': "C", 'CCC': "C", 'CH': "C", 'CSF': "C", 'IC': "C", 'LC': "C",
              'M4C': "C", 'M5M': "C", 'N5M': "C", 'OMC': "C", 'PMT': "C", 'RPC': "C", 'S4C': "C", 'ZBC': "C", 'ZCY': "C",
              '0AP': "C", '0R8': "C", '1CC': "C", '1FC': "C", '47C': "C", '4PC': "C", '4PD': "C", '4PE': "C", '4SC': "C",
              '5CM': "C", '5FC': "C", '5HC': "C", '5NC': "C", '5PC': "C", '6HC': "C", 'B7C': "C", 'C2S': "C", 'C32': "C",
              'C34': "C", 'C36': "C", 'C37': "C", 'C38': "C", 'C42': "C", 'C45': "C", 'C46': "C", 'C49': "C", 'C4S': "C",
              'CAR': "C", 'CB2': "C", 'CBR': "C", 'CDW': "C", 'CFL': "C", 'CFZ': "C", 'CMR': "C", 'CP1': "C", 'CSL': "C",
              'CX2': "C", 'CT': "C", 'DFC': "C", 'DNR': "C", 'DOC': "C", 'EXC': "C", 'GCK': "C", 'I5C': "C", 'IMC': "C",
              'MCY': "C", 'ME6': "C", 'PVX': "C", 'SC': "C", 'TC1': "C", 'TPC': "C", 'XCL': "C", 'XCR': "C", 'XCT': "C",
              'XCY': "C", 'YCO': "C", 'Z': "C", 
    'G': "G", 'DG': "G", 'GUA': "G", '+G': "G",
              '0AD': "G", '0UH': "G", '2PR': "G", '5CG': "G", '63G': "G", '63H': "G", '6HG': "G", '6OG': "G", '6PO': "G",
              '7GU': "G", '8AG': "G", '8FG': "G", '8MG': "G", '8OG': "G", 'AFG': "G", 'BGM': "G", 'C6G': "G", 'DCG': "G",
              'DG': "G", 'DFG': "G", 'G8': "G", 'GI': "G", 'GP': "G", 'EFG': "G", 'EHG': "G", 'FG': "G", 'FMG': "G",
              'FOX': "G", 'G2S': "G", 'G31': "G", 'G32': "G", 'G33': "G", 'G36': "G", 'G38': "G", 'G42': "G", 'G47': "G",
              'G49': "G", 'GDR': "G", 'GF2': "G", 'GFL': "G", 'GMS': "G", 'GN7': "G", 'GS': "G", 'GSR': "G", 'GSS': "G",
              'GX1': "G", 'HN0': "G", 'HN1': "G", 'IGU': "G", 'LCG': "G", 'LGP': "G", 'M1G': "G", 'MG1': "G", 'MRG': "G",
              'OGX': "G", 'P': "G", 'PG7': "G", 'PGN': "G", 'PPW': "G", 'S4G': "G", 'S6G': "G", 'SG': "G", 'TGP': "G",
              'X': "G", 'XGL': "G", 'XGR': "G", 'XGU': "G", 'XUG': "G", '102': "G", '18M': "G", '1MG': "G", '23G': "G",
              '2EG': "G", '2MG': "G", '7MG': "G", 'A6G': "G", 'CG1': "G", 'G1G': "G", 'G25': "G", 'G2L': "G", 'G3A': "G",
              'G46': "G", 'G48': "G", 'G7M': "G", 'GAO': "G", 'GDO': "G", 'GDP': "G", 'GH3': "G", 'GNG': "G", 'GOM': "G",
              'GRB': "G", 'GTP': "G", 'IG': "G", 'IMP': "G", 'KAG': "G", 'LG': "G", 'M2G': "G", 'MGT': "G", 'MGV': "G",
              'N6G': "G", 'O2G': "G", 'OMG': "G", 'PGP': "G", 'QUO': "G", 'TPG': "G", 'XTS': "G", 'YG': "G", 'YYG': "G",
              'ZGU': "G",
    'I': "I", 'DI': "I", 'INO': "I", '+I': "I",
              '2BD': "I", 'OIP': "I", 
    'T': "T", 'DT': "T", 'THY': "T", '+T': "T",
              '2AT': "T", '2BT': "T", '2T': "T", '2GT': "T", '2NT': "T", '2OT': "T", '2ST': "T", '5AT': "T", '5HT': "T",
              '5IT': "T", '5PY': "T", '64T': "T", '6CT': "T", '6HT': "T", 'ATD': "T", 'ATL': "T", 'ATM': "T", 'BOE': "T",
              'CTG': "T", 'D3T': "T", 'D4M': "T", 'DPB': "T", 'DRT': "T", 'EIT': "T", 'F3H': "T", 'F4H': "T", 'JT': "T",
              'MMT': "T", 'MTR': "T", 'NMS': "T", 'NMT': "T", 'P2T': "T", 'PST': "T", 'S2M': "T", 'SPT': "T", 'T32': "T",
              'T36': "T", 'T37': "T", 'T39': "T", 'T3P': "T", 'T48': "T", 'T49': "T", 'T4S': "T", 'T5S': "T", 'TA3': "T",
              'TAF': "T", 'TCP': "T", 'TDY': "T", 'TED': "T", 'TFE': "T", 'TFF': "T", 'TFT': "T", 'TLC': "T", 'TP1': "T",
              'TTD': "T", 'TTM': "T", 'US3': "T", 'XTF': "T", 'XTH': "T", 'XTL': "T", 'XTR': "T", 
    'U': "U", 'DU': "U", 'URA': "U", '+U': "U",
              '0AU': "U", '18Q': "U", '5HU': "U", '5IU': "U", '5SE': "U", 'BRU': "U", 'BVP': "U", 'DDN': "U", 'DRM': "U",
              'UZ': "U", 'GMU': "U", 'HDP': "U", 'HEU': "U", 'NDN': "U", 'NU': "U", 'OHU': "U", 'P2U': "U", 'PU': "U",
              'T5O': "U", 'TLN': "U", 'TTI': "U", 'U2N': "U", 'U33': "U", 'UBI': "U", 'UBR': "U", 'UCL': "U", 'UF2': "U",
              'UFR': "U", 'UFT': "U", 'UMS': "U", 'UMX': "U", 'UPE': "U", 'UPS': "U", 'URX': "U", 'US1': "U", 'US2': "U",
              'USM': "U", 'UVX': "U", 'ZU': "U", '125': "U", '126': "U", '127': "U", '1RN': "U", '2AU': "U", '2MU': "U",
              '2OM': "U", '3AU': "U", '3ME': "U", '3MU': "U", '3TD': "U", '4SU': "U", '5BU': "U", '5FU': "U", '5MU': "U",
              '70U': "U", 'A6U': "U", 'CNU': "U", 'DHU': "U", 'FHU': "U", 'FNU': "U", 'H2U': "U", 'IU': "U", 'LHU': "U",
              'MEP': "U", 'MNU': "U", 'OMU': "U", 'ONE': "U", 'PSU': "U", 'PYO': "U", 'RSQ': "U", 'RUS': "U", 'S4U': "U",
              'SSU': "U", 'SUR': "U", 'T31': "U", 'U25': "U", 'U2L': "U", 'U2P': "U", 'U31': "U", 'U34': "U", 'U36': "U",
              'U37': "U", 'U8U': "U", 'UAR': "U", 'UBB': "U", 'UBD': "U", 'UD5': "U", 'UPV': "U", 'UR3': "U", 'URD': "U",
              'US5': "U", 'UZR': "U", 'ZBU': "U"
}

nucleic1to3 = dict([[v,k] for k,v in list(nucleic3to1.items())])
nucleic1to3['A'] = "DA"
nucleic1to3['C'] = "DC"
nucleic1to3['G'] = "DG"
nucleic1to3['I'] = "DI"
nucleic1to3['T'] = "DT"
nucleic1to3['U'] = "DU"

"""
REGULAR NUCLEOTIDES IDENTIFICATION
"""
nucleic_main3 = set(['ADE','CYT','GUA','INO','THY','URA'])

nucleic_main2 = set(['DA','DC','DG','DI','DT','DU'])

nucleic_main1 = set(['A','C','G','I','T','U'])

"""
PROPERTIES
"""
nitrogenous_bases = {
    'A': "U", 'C': "Y", 'G': "U", 'I': "U", 'T': "Y", 'U': "Y"
}

dna_complementary = {
    'A': "T", 'C': "G", 'G': "C", 'T': "A", 'U': "A", 'I': "C", 'N': "N",
    'a': "t", 'c': "g", 'g': "c", 't': "a", 'u': "a", 'i': "c", 'n': "n",
}

rna_complementary = {
    'A': "U", 'C': "G", 'G': "C", 'U': "A", 'T': "A", 'I': "C", 'N': "N",
    'a': "u", 'c': "g", 'g': "c", 'u': "a", 't': "a", 'i': "c", 'n': "n",
}

nucleotides_acceptors = {
    'A' : ['N3', 'N7'],
    'C' : ['O2'],
    'G' : ['N3', 'N7', 'O6'],
    'T' : ['O2', 'O4'],
}

nucleotides_donors = {
    'A' : ['N6'],
    'C' : ['C5', 'N4'],
    'G' : ['N2'],
    'T' : ['C5M'],
}


"""
CRYSTALOGRAPHIC METHODS
"""
crystal_method_has_resolution = set(['X-RAY DIFFRACTION','ELECTRON MICROSCOPY','NEUTRON DIFFRACTION',
                                     'FIBER DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY'])

crystal_method_not_resolution = set(['SOLUTION NMR','POWDER DIFFRACTION','SOLUTION SCATTERING','SOLID-STATE NMR',
                                     'INFRARED SPECTROSCOPY','FLUORESCENCE TRANSFER'])

crystal_method                = crystal_method_has_resolution.union(crystal_method_not_resolution)

"""
ATOMS
radius from: https://en.wikipedia.org/wiki/Covalent_radius
"""
class Element(object):
    def __init__(self, number, symbol, name, radius=None):
        self.number = number
        self.symbol = symbol
        self.name   = name
        self.radius = radius


element_dic = {'H':  Element(  1, 'H',  'Hydrogen', 31),    'He': Element(  2, 'He', 'Helium'),
               'Li': Element(  3, 'Li', 'Lithium'),         'Be': Element(  4, 'Be', 'Beryllium'),
               'B':  Element(  5, 'B',  'Boron'),           'C':  Element(  6, 'C',  'Carbon', 69),
               'N':  Element(  7, 'N',  'Nitrogen', 71),    'O':  Element(  8, 'O',  'Oxygen', 66),
               'F':  Element(  9, 'F',  'Fluorine'),        'Ne': Element( 10, 'Ne', 'Neon'),
               'Na': Element( 11, 'Na', 'Sodium'),          'Mg': Element( 12, 'Mg', 'Magnesium'),
               'Al': Element( 13, 'Al', 'Aluminium'),       'Si': Element( 14, 'Si', 'Silicon'),
               'P':  Element( 15, 'P',  'Phosphorus', 107), 'S':  Element( 16, 'S',  'Sulfur', 105),
               'Cl': Element( 17, 'Cl', 'Chlorine'),        'Ar': Element( 18, 'Ar', 'Argon'),
               'K':  Element( 19, 'K',  'Potassium'),       'Ca': Element( 20, 'Ca', 'Calcium'),
               'Sc': Element( 21, 'Sc', 'Scandium'),        'Ti': Element( 22, 'Ti', 'Titanium'),
               'V':  Element( 23, 'V',  'Vanadium'),        'Cr': Element( 24, 'Cr', 'Chromium'),
               'Mn': Element( 25, 'Mn', 'Manganese'),       'Fe': Element( 26, 'Fe', 'Iron'),
               'Co': Element( 27, 'Co', 'Cobalt'),          'Ni': Element( 28, 'Ni', 'Nickel'),
               'Cu': Element( 29, 'Cu', 'Copper'),          'Zn': Element( 30, 'Zn', 'Zinc'),
               'Ga': Element( 31, 'Ga', 'Gallium'),         'Ge': Element( 32, 'Ge', 'Germanium'),
               'As': Element( 33, 'As', 'Arsenic'),         'Se': Element( 34, 'Se', 'Selenium', 120),
               'Br': Element( 35, 'Br', 'Bromine'),         'Kr': Element( 36, 'Kr', 'Krypton'),
               'Rb': Element( 37, 'Rb', 'Rubidium'),        'Sr': Element( 38, 'Sr', 'Strontium'),
               'Y':  Element( 39, 'Y',  'Yttrium'),         'Zr': Element( 40, 'Zr', 'Zirconium'),
               'Nb': Element( 41, 'Nb', 'Niobium'),         'Mo': Element( 42, 'Mo', 'Molybdenum'),
               'Tc': Element( 43, 'Tc', 'Technetium'),      'Ru': Element( 44, 'Ru', 'Ruthenium'),
               'Rh': Element( 45, 'Rh', 'Rhodium'),         'Pd': Element( 46, 'Pd', 'Palladium'),
               'Ag': Element( 47, 'Ag', 'Silver'),          'Cd': Element( 48, 'Cd', 'Cadmium'),
               'In': Element( 49, 'In', 'Indium'),          'Sn': Element( 50, 'Sn', 'Tin'),
               'Sb': Element( 51, 'Sb', 'Antimony'),        'Te': Element( 52, 'Te', 'Tellurium'),
               'I':  Element( 53, 'I',  'Iodine'),          'Xe': Element( 54, 'Xe', 'Xenon'),
               'Cs': Element( 55, 'Cs', 'Caesium'),         'Ba': Element( 56, 'Ba', 'Barium'),
               'La': Element( 57, 'La', 'Lanthanum'),       'Ce': Element( 58, 'Ce', 'Cerium'),
               'Pr': Element( 59, 'Pr', 'Praseodymium'),    'Nd': Element( 60, 'Nd', 'Neodymium'),
               'Pm': Element( 61, 'Pm', 'Promethium'),      'Sm': Element( 62, 'Sm', 'Samarium'),
               'Eu': Element( 63, 'Eu', 'Europium'),        'Gd': Element( 64, 'Gd', 'Gadolinium'),
               'Tb': Element( 65, 'Tb', 'Terbium'),         'Dy': Element( 66, 'Dy', 'Dysprosium'),
               'Ho': Element( 67, 'Ho', 'Holmium'),         'Er': Element( 68, 'Er', 'Erbium'),
               'Tm': Element( 69, 'Tm', 'Thulium'),         'Yb': Element( 70, 'Yb', 'Ytterbium'),
               'Lu': Element( 71, 'Lu', 'Lutetium'),        'Hf': Element( 72, 'Hf', 'Hafnium'),
               'Ta': Element( 73, 'Ta', 'Tantalum'),        'W':  Element( 74, 'W', 'Tungsten'),
               'Re': Element( 75, 'Re', 'Rhenium'),         'Os': Element( 76, 'Os', 'Osmium'),
               'Ir': Element( 77, 'Ir', 'Iridium'),         'Pt': Element( 78, 'Pt', 'Platinum'),
               'Au': Element( 79, 'Au', 'Gold'),            'Hg': Element( 80, 'Hg', 'Mercury'),
               'Tl': Element( 81, 'Tl', 'Thallium'),        'Pb': Element( 82, 'Pb', 'Lead'),
               'Bi': Element( 83, 'Bi', 'Bismuth'),         'Po': Element( 84, 'Po', 'Polonium'),
               'At': Element( 85, 'At', 'Astatine'),        'Rn': Element( 86, 'Rn', 'Radon'),
               'Fr': Element( 87, 'Fr', 'Francium'),        'Ra': Element( 88, 'Ra', 'Radium'),
               'Ac': Element( 89, 'Ac', 'Actinium'),        'Th': Element( 90, 'Th', 'Thorium'),
               'Pa': Element( 91, 'Pa', 'Protactinium'),    'U':  Element( 92, 'U', 'Uranium'),
               'Np': Element( 93, 'Np', 'Neptunium'),       'Pu': Element( 94, 'Pu', 'Plutonium'),
               'Am': Element( 95, 'Am', 'Americium'),       'Cm': Element( 96, 'Cm', 'Curium'),
               'Bk': Element( 97, 'Bk', 'Berkelium'),       'Cf': Element( 98, 'Cf', 'Californium'),
               'Es': Element( 99, 'Es', 'Einsteinium'),     'Fm': Element(100, 'Fm', 'Fermium'),
               'Md': Element(101, 'Md', 'Mendelevium'),     'No': Element(102, 'No', 'Nobelium'),
               'Lr': Element(103, 'Lr', 'Lawrencium'),      'Rf': Element(104, 'Rf', 'Rutherfordium'),
               'Db': Element(105, 'Db', 'Dubnium'),         'Sg': Element(106, 'Sg', 'Seaborgium'),
               'Bh': Element(107, 'Bh', 'Bohrium'),         'Hs': Element(108, 'Hs', 'Hassium'),
               'Mt': Element(109, 'Mt', 'Meitnerium')
}
