class CDhitHomolog(object):
    def __init__(self, name, length, homology):
        self._name     = name.lstrip('>').rstrip('.')
        self._length   = int(length.replace('aa','').rstrip(','))
        self._homology = homology.replace('%','')

    @property
    def name(self):        return self._name
    @property
    def length(self):      return self._length
    @property
    def homology(self): 
        if self.is_master: return self._homology
        else:              return int(self._homology)
    @property
    def is_master(self):   return self._homology == '*'

    def __repr__(self):
        if not self.is_master:
            return '{0.name}: {0.length:0004d} Aa with {0.homology:003d}%'.format(self)
        else:
            return '{0.name}: {0.length:0004d} Aa.'.format(self)
