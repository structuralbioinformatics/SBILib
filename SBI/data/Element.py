# -*-
#
# @author: jaumebonet
# @email:  jaume.bonet@gmail.com
# @url:    jaumebonet.github.io
#
# @date:   2016-02-12 13:20:21
#
# @last modified by:   jaumebonet
# @last modified time: 2016-02-15 10:59:52
#
# -*-
from pynion import Singleton

element_list = [
    ('H',  'Hydrogen'),    ('He', 'Helium'),       ('Li', 'Lithium'),      ('Be', 'Beryllium'),
    ('B',  'Boron'),       ('C',  'Carbon'),       ('N',  'Nitrogen'),     ('O',  'Oxygen'),
    ('F',  'Fluorine'),    ('Ne', 'Neon'),         ('Na', 'Sodium'),       ('Mg', 'Magnesium'),
    ('Al', 'Aluminium'),   ('Si', 'Silicon'),      ('P',  'Phosphorus'),   ('S',  'Sulfur'),
    ('Cl', 'Chlorine'),    ('Ar', 'Argon'),        ('K',  'Potassium'),    ('Ca', 'Calcium'),
    ('Sc', 'Scandium'),    ('Ti', 'Titanium'),     ('V',  'Vanadium'),     ('Cr', 'Chromium'),
    ('Mn', 'Manganese'),   ('Fe', 'Iron'),         ('Co', 'Cobalt'),       ('Ni', 'Nickel'),
    ('Cu', 'Copper'),      ('Zn', 'Zinc'),         ('Ga', 'Gallium'),      ('Ge', 'Germanium'),
    ('As', 'Arsenic'),     ('Se', 'Selenium'),     ('Br', 'Bromine'),      ('Kr', 'Krypton'),
    ('Rb', 'Rubidium'),    ('Sr', 'Strontium'),    ('Y',  'Yttrium'),      ('Zr', 'Zirconium'),
    ('Nb', 'Niobium'),     ('Mo', 'Molybdenum'),   ('Tc', 'Technetium'),   ('Ru', 'Ruthenium'),
    ('Rh', 'Rhodium'),     ('Pd', 'Palladium'),    ('Ag', 'Silver'),       ('Cd', 'Cadmium'),
    ('In', 'Indium'),      ('Sn', 'Tin'),          ('Sb', 'Antimony'),     ('Te', 'Tellurium'),
    ('I',  'Iodine'),      ('Xe', 'Xenon'),        ('Cs', 'Caesium'),      ('Ba', 'Barium'),
    ('La', 'Lanthanum'),   ('Ce', 'Cerium'),       ('Pr', 'Praseodymium'), ('Nd', 'Neodymium'),
    ('Pm', 'Promethium'),  ('Sm', 'Samarium'),     ('Eu', 'Europium'),     ('Gd', 'Gadolinium'),
    ('Tb', 'Terbium'),     ('Dy', 'Dysprosium'),   ('Ho', 'Holmium'),      ('Er', 'Erbium'),
    ('Tm', 'Thulium'),     ('Yb', 'Ytterbium'),    ('Lu', 'Lutetium'),     ('Hf', 'Hafnium'),
    ('Ta', 'Tantalum'),    ('W',  'Tungsten'),     ('Re', 'Rhenium'),      ('Os', 'Osmium'),
    ('Ir', 'Iridium'),     ('Pt', 'Platinum'),     ('Au', 'Gold'),         ('Hg', 'Mercury'),
    ('Tl', 'Thallium'),    ('Pb', 'Lead'),         ('Bi', 'Bismuth'),      ('Po', 'Polonium'),
    ('At', 'Astatine'),    ('Rn', 'Radon'),        ('Fr', 'Francium'),     ('Ra', 'Radium'),
    ('Ac', 'Actinium'),    ('Th', 'Thorium'),      ('Pa', 'Protactinium'), ('U',  'Uranium'),
    ('Np', 'Neptunium'),   ('Pu', 'Plutonium'),    ('Am', 'Americium'),    ('Cm', 'Curium'),
    ('Bk', 'Berkelium'),   ('Cf', 'Californium'),  ('Es', 'Einsteinium'),  ('Fm', 'Fermium'),
    ('Md', 'Mendelevium'), ('No', 'Nobelium'),     ('Lr', 'Lawrencium'),   ('Rf', 'Rutherfordium'),
    ('Db', 'Dubnium'),     ('Sg', 'Seaborgium'),   ('Bh', 'Bohrium'),      ('Hs', 'Hassium'),
    ('Mt', 'Meitnerium'),  ('Ds', 'Darmstadtium'), ('Rg', 'Roentgenium'),  ('Cn', 'Copernicium')
]


class Element(object):
    def __init__(self, number, symbol, name):
        self.number = int(number)
        self.symbol = symbol
        self.name   = name

    def __str__(self):
        return "{0.number:>3}\t{0.symbol:>2}\t{0.name}".format(self)


class PeriodicTable(object):
    __metaclass__ = Singleton

    def __init__(self):
        self._table = {}
        for i in range(len(element_list)):
            e = element_list[i]
            self._table.setdefault(e[0], Element(i + 1, e[0], e[1]))

    def __contains__(self, key):
        if isinstance(key, str) and key.isdigit():
            key = int(key)
        if isinstance( key, ( int, long ) ):
            key = key - 1
            if key > len(element_list) - 1 or key < 0:
                return False
            else:
                return True
        elif isinstance(key, str):
            return key.capitalize() in self._table
        else:
            return False

    def __getitem__(self, key):
        if isinstance(key, str) and key.isdigit():
            key = int(key)
        if isinstance( key, ( int, long ) ):
            key = key - 1
            if key > len(element_list) - 1 or key < 0:
                return KeyError('No element with this size is stored.')
            else:
                return self._table[element_list[key][0]]
        elif isinstance(key, str):
            return self._table[key.capitalize()]
        else:
            raise KeyError('KeyType should be integer or string')
