class BlastError(Exception):
    """
    Complementary Error Class: BlastError   -> raises when there are no hits for a given query
                                            -> raises when there is no path for blastpgp
                                            -> raises when there is an error in blast execution
    """
    def __init__(self, code = None, value = None):

        self.value = value
        self.code  = code

        if self.get_code() is not None:
            self._determine_message()

    def get_value(self):        return self.value
    def get_code(self):         return self.code

    def _determine_message(self):

        if self.get_code()   ==  0:
            self.value = "There is something wrong with the XML blast output file %s\n" %self.get_value()

        elif self.get_code() == -1:
            self.value = "blast execution failed. Returns error:\n\n%s" %self.get_value()

        elif self.get_code() == -2:
            self.value = "blast is not in $PATH enviroment variable.\n"
            self.value = self.get_value() + "The path can be added through the 'executablePath' option.\n"

        elif self.get_code() == -3:
            self.value = "Input fasta file %s does not exist.\n" %self.get_value()

        elif self.get_code() == -4:
            self.value = "Input fasta file %s is a multi-fasta file\n" %self.get_value() 

        elif self.get_code() == -5:
            self.value = "Query database %s does not exist\n" %self.get_value()

        elif self.get_code() == -6:
            if len(self.get_value()) == 3:
                text = "In order to perform a blastpgp search, the database MUST be formated (see blast's user guide).\n"
                self.value = text
            else:
                text = "Some of the format files for the database generated with formatdb are missing.\n"
                text = text + "Please check %s" %(",".join(self.get_value()))
                self.value = text

        elif self.get_code() == -7:
            self.value = "The specifyed weight matrix is not supported.\nSupported matrices are: %s\n" %self.get_value()

        elif self.get_code() == -8:
            self.value = "For the given matrix only certain gap opening and extension values are allowed:\n %s\n" %self.get_value()

        elif self.get_code() == -9:
            self.value = "The given path %s does not lead to the blastpgp executable.\n" %self.get_value()

        elif self.get_code() == -10:
            self.value = "The kind of blast can only be 'prot' or 'nucl'."

        elif self.get_code() == -11:
            pass
        
        elif self.get_code() ==  1:
            pass

        elif self.get_code() ==  2:
            self.value = "The BlastHandler has been loaded with a blast XML output file, no need to execute blastpgp\n"

        elif self.get_code() == 3:
            pass

    def __str__(self):
        error_str = "\n\n ###################\n[BlastError Type: %d]\n ###################\n%s\n" %(self.get_code(), self.get_value())
        if self.get_code() < 0:
            error_str = error_str + "\nBlastErrors with negative type indicate that the blast search has not been executed.\n"
        elif self.get_code() > 0:
            error_str = error_str + "\nBlastErrors with positive type are related to problems during the parsing of the blast XML output files.\n"
        return error_str