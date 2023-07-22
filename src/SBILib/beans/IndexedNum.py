import re


class IndexedNum(object):
    '''
    {IndexedNum} has been created for PDB.
    Some PDB files contain residues in positions that are not only numeric, but also have an index.
    IndexedNum tries to paliate this by giving the properties of a number (as much as possible) to this kind of number-index combination.
    Thus, a 7 becomes a '7 ' but can still be used as a number.
    Comparissons, addition and int(), float() casts are implemented
    '''
    number_regex = re.compile('(\-*\d+)(\w*)')

    def __init__(self, value):
        '''
        {IndexedNum} has two simple attributes:
            - number (int)
            - index  (string)
        If it is initilized with X, it means that we have a "blank" position.
        To avoid mixups, the index is transformed to '_X_'

        @type  value: String
        @type  value: IndexedNum
        @type  value: int
        @param value: data to transform to a indexed number
        '''
        if isinstance(value, IndexedNum):
            self._num = value._num
            self._idx = value._idx
        else:
            value = str(value)
            if value != 'X':
                parts = re.search(self.number_regex, value)
                self._num = int(parts.group(1))
                self._idx = parts.group(2)
                if self._idx == '':
                    self._idx = ' '
            else:
                self._num = 0
                self._idx = '_X_'

    #
    # ATTRIBUTES
    #
    @property
    def number(self):
        return self._num

    @property
    def index(self):
        return self._idx

    #
    # BOOLEANS
    #
    @property
    def is_integer(self):
        return self._idx == ' '

    @property
    def is_blank(self):
        return self._idx == '_X_'

    #
    # DEFAULT COMPARATIVE FUNCTIONS
    #

    def __lt__(self, other):
        if isinstance(other, int):
            return self.number < other
        elif isinstance(other, IndexedNum):
            if self.number < other.number:
                return True
            elif self.number == other.number:
                return self.index < other.index
            else:
                return False
        return NotImplemented

    def __le__(self, other):
        if isinstance(other, int):
            return self.number <= other
        elif isinstance(other, IndexedNum):
            if self.number < other.number:
                return True
            elif self.number == other.number:
                return self.index <= other.index
            else:
                return False
        return NotImplemented

    def __eq__(self, other):
        if isinstance(other, int):
            return self.number == other
        elif isinstance(other, IndexedNum):
            return self.number == other.number and self.index == other.index
        elif isinstance(other, str):
            return str(self).strip() == other.strip()
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def __gt__(self, other):
        if isinstance(other, int):
            return self.number > other
        elif isinstance(other, IndexedNum):
            if self.number > other.number:
                return True
            elif self.number == other.number:
                return self.index > other.index
            else:
                return False
        return NotImplemented

    def __ge__(self, other):
        if isinstance(other, int):
            return self.number >= other
        elif isinstance(other, IndexedNum):
            if self.number > other.number:
                return True
            elif self.number == other.number:
                return self.index >= other.index
            else:
                return False
        return NotImplemented

    def __add__(self, other):
        if isinstance(other, int):
            return IndexedNum(str(self.number + other) + self.index)
        if isinstance(other, IndexedNum) and other.is_integer:
            return IndexedNum(str(self.number + other) + self.index)
        return NotImplemented

    def __iadd__(self, other):
        if isinstance(other, int):
            self.number += other
        if isinstance(other, IndexedNum) and other.is_integer:
            self.number += other
        return NotImplemented

    def __int__(self):
        return int(self.number)

    def __float__(self):
        return float(self.number)

    def __hash__(self):
        return repr(self)

    def __repr__(self):
        if self.is_blank:
            return 'X'
        else:
            return (str(self.number) + self.index).strip()

    def __str__(self):
        return repr(self)
