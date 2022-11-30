# -*-
#
# @author: jaumebonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.github.io
#
# @date:   2016-02-15 11:08:26
#
# @last modified by:   jaumebonet
# @last modified time: 2016-02-15 14:55:08
#
# -*-
from pynion import Multiton


class SetDict(dict):
    __metaclass__ = Multiton

    def __init__(self, name):
        self._representative = {}

    def __setitem__(self, key, value):
        for k in range(len(key)):
            if k == 0:
                self._representative[value] = key[k]
            super(SetDict, self).__setitem__(key[k], value)

    def representative(self, key):
        v = super(SetDict, self).__getitem__(key)
        return self._representative[v]


Aa_sdict = SetDict('AminoAcids')
a = ['ALA', 'AZT', 'CHA', 'HPH', 'NAL', 'AIB', 'BAL', 'DHA', 'BB9', 'ALM', 'AYA', 'BNN',
     'CHG', 'CSD', 'DAL', 'DNP', 'FLA', 'HAC', 'MAA', 'PRR', 'TIH', 'TPQ', 'BB9']
c = ['CYS', 'CYD', 'CYO', 'HCY', 'CSX', 'SMC', 'BCS', 'BUC', 'C5C', 'C6C', 'CCS', 'CEA',
     'CME', 'CSO', 'CSP', 'CSS', 'CSW', 'CY1', 'CY3', 'CYG', 'CYM', 'CYQ', 'DCY', 'OCS',
     'SOC', 'EFC', 'PR3', 'SCH', 'SCS', 'SCY', 'SHC', 'PEC']
d = ['ASP', 'ASZ', '2AS', 'ASA', 'ASB', 'ASK', 'ASL', 'ASQ', 'BHD', 'DAS', 'DSP']
r = ['ARG', 'ORN', 'ACL', 'ARM', 'AGM', 'HAR', 'HMR', 'DAR']
e = ['GLU', 'GLA', 'GLZ', 'PCA', '5HP', 'CGU', 'DGL', 'GGL', 'GMA']
f = ['PHE', 'DAH', 'HPQ', 'DPN', 'PHI', 'PHL']
g = ['GLY', 'GL3', 'GLZ', 'GSC', 'SAR', 'MPQ', 'NMC', 'MSA', 'DBU']
h = ['HIS', 'HSE', 'HSD', 'HI0', 'HIP', 'HID', 'HIE', '3AH', 'MHS', 'DHI', 'HIC', 'NEP',
     'NEM']
i = ['ILE', 'IIL', 'DIL']
k = ['LYS', 'LYZ', 'ALY', 'TRG', 'SHR', 'LYM', 'LLY', 'KCX', 'LLP', 'DLY', 'DM0']
l = ['LEU', 'NLE', 'LOV', 'NLN', 'NLP', 'MLE', 'BUG', 'CLE', 'DLE', 'MLU']
m = ['MET', 'MSE', 'CXM', 'FME', 'OMT']
n = ['ASN', 'MEN']
p = ['PRO', 'HYP', 'DPR', 'ECQ', 'POM', 'H5M']
q = ['GLN', 'DGN']
s = ['SER', 'HSE', 'STA', 'SVA', 'SAC', 'SEL', 'SEP', 'SET', 'OAS', 'DSN', 'MIS']
t = ['THR', 'PTH', 'ALO', 'TPO', 'BMT', 'DTH', 'CTH']
v = ['VAL', 'NVA', 'DVA', 'DIV', 'MVA']
w = ['TRP', 'TPL', 'TRO', 'DTR', 'HTR', 'LTR']
x = ['ACE', '3FG', 'UNK']
y = ['TYR', 'TYQ', 'TYS', 'TYY', 'TYB', 'STY', 'PTR', 'PAQ', 'DTY', 'IYR', 'GHP', 'D3P',
     'D4P', 'OMZ', 'OMY']
Aa_sdict[a]       = 'A'
Aa_sdict[['ASX']] = 'B'
Aa_sdict[c]       = 'C'
Aa_sdict[d]       = 'D'
Aa_sdict[e]       = 'E'
Aa_sdict[f]       = 'F'
Aa_sdict[g]       = 'G'
Aa_sdict[h]       = 'H'
Aa_sdict[i]       = 'I'
Aa_sdict[['XLE']] = 'J'
Aa_sdict[k]       = 'K'
Aa_sdict[l]       = 'L'
Aa_sdict[m]       = 'M'
Aa_sdict[n]       = 'N'
Aa_sdict[['PYL']] = 'O'
Aa_sdict[p]       = 'P'
Aa_sdict[q]       = 'Q'
Aa_sdict[r]       = 'R'
Aa_sdict[s]       = 'S'
Aa_sdict[t]       = 'T'
Aa_sdict[['SEC']] = 'U'
Aa_sdict[v]       = 'V'
Aa_sdict[w]       = 'W'
Aa_sdict[x]       = 'X'
Aa_sdict[y]       = 'Y'
Aa_sdict[['GLX']] = 'Z'
