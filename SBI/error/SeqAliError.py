class SeqAliError(Exception):
    """
    Complementary Error Class: SeqAliError
    """
    def __init__(self, code = None, value = None):

        self.value = value
        self.code  = code

        if self.get_code() is not None:
            self._determine_message()

    def get_value(self):        return self.value
    def get_code(self):         return self.code

    def _determine_message(self):

        if self.get_code() ==  1:
            self.value = "Something has gone wrong with the sequence fragment detection\n"

        elif self.get_code() ==  2:
            self.value = "A multiple sequence alignment of {0} sequences cannot be evaluated with Rost's Twilight Zone\n".format(self.value)

    def __str__(self):
        error_str = "\n\n ###################\n[SeqAliError Type: %d]\n ###################\n%s\n" %(self.get_code(), self.get_value())
        return error_str