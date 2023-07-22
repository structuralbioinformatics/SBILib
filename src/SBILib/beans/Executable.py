"""
Executable

author: jbonet
date:   03/2013

@oliva's lab
"""

import os
import subprocess
from SBILib import SBIglobals


class Executable(object):
    """
    Checks the integrity of an executable program
    """
    def __init__(self, executable, path=None, variable_path=None):
        """
        @type  executable: String
        @param executable: Name of the executable file

        @type  path: String
        @param path: Putative path to the executable (if case is not in $PATH)

        @type  variable_path: String
        @param variable_path: Name of the enviroment variable with the path
        """

        self._exec = executable
        self._path = path

        if self._exec is None:
            raise AttributeError("The executable name MUST be specified")

        if variable_path is None:
            if path is not None:
                self._path = os.path.abspath(path)
            else:
                if not self._load_executable_path():
                    raise EnvironmentError("")
        else:
            self._load_variable_path(variable_path=variable_path)

        self._check_executable()

        self._command = []
        self._command.append(self.full_executable)

        self._outfile = None
        self._stdout  = None

    #
    # ATTRIBUTES
    #
    @property
    def executable(self):
        return self._exec

    @property
    def path(self):
        return self._path

    @property
    def command(self):
        return self._command

    @property
    def outfile(self):
        return self._outfile

    @outfile.setter
    def outfile(self, value):
        self._outfile = value

    @property
    def full_executable(self):
        return os.path.join(self._path, self._exec)

    #
    # METHODS
    #
    def add_attribute(self, attribute_value, attribute_id=None):
        if attribute_id is not None:
            self._command.append(attribute_id)
        self._command.append(str(attribute_value))

    def add_parameter(self, parameter):
        self._command.append(str(parameter))

    def clean_command(self):
        self._command = []
        self._command.append(self.full_executable)

    def execute(self, stdout=False, silent=False):
        """
        Executes the commands
        @type  stdout: Boolean
        @param stdout: determines if the output is through stdout
        """
        if self.command is None:
            return False

        if not stdout:
            stdoutPIPE = open('/dev/null', 'w')
        else:
            stdoutPIPE = subprocess.PIPE

        if silent:
            stderrPIPE = open('/dev/null', 'w')
        else:
            stderrPIPE = subprocess.PIPE

        SBIglobals.alert('debug', self, '\tExecuting command:\n\t{0}\n'.format(" ".join(self.command)))

        p = subprocess.Popen(self.command, stdout=stdoutPIPE,
                             stderr=stderrPIPE)

        out, err  = p.communicate()

        if stdout:
            self._stdout = out

        if not silent and err.strip() != b'':
            raise SystemError(err)

    #
    # PRIVATE METHODS
    #
    def _load_variable_path(self, variable_path):
        """
        Retrieves the path from a variable_path
        """
        try:
            self._path = os.environ[variable_path]
        except KeyError:
            raise EnvironmentError("The given Environment Variable %s is not defined" %variable_path)

    def _load_executable_path(self):
        """
        Retrieves the executable path in case self._path is None
        """
        if self._path is not None:
            return

        search = ["which", self.executable]
        p = subprocess.Popen(search, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out, err = p.communicate()

        if out != '':
            self._path = os.path.split(out.strip())[0]
            return True
        else:
            return False

    def _check_executable(self):
        """
        Checks that the final executable can be executed
        """
        if not os.path.isfile(self.full_executable):
            raise SystemError("The given executable file %s does not exist" %self.full_executable)
        if not os.access(self.full_executable, os.X_OK):
            raise SystemError("The given executable %s can not be executed" %self.full_executable)

    '''OVERLOAD DEFAULT FUNCTIONS'''
    def __repr__(self):
        return " ".join(self._command)

    def __str__(self):
        return repr(self)
