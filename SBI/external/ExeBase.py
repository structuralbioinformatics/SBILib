# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>

.. class:: ExeBase
"""
# Standard Libraries
import abc
import copy
import getpass
import os
import subprocess
import sys


# External Libraries
import numpy as np
import six

# This Library
import SBI.core as core

__all__ = ['ExeBase']


class ExeCommand( object ):
    """Unified management of the execution of external applications.

    :param str executable: Full path to the application's executable.

    :raises:
        :IOError: If the provided application cannot be found.
        :IOError: If the provided application cannot be executed by the user.
    """
    def __init__( self, executable ):
        # Attributes
        self._exec           = os.path.abspath(executable)
        self._command        = []
        self._backup_command = []

        # Check it exists
        if not os.path.isfile(self._exec):
            raise IOError('Application {} not found'.format(self._exec))
        if not os.access(self._exec, os.X_OK):
            raise IOError('User {} does not have access to execute {}'.format(getpass.getuser(), self._exec))

        self._command.append(self._exec)

    # METHODS
    def add_parameter( self, parameter ):
        """Add a new parameter to the command.

        :param parameter: New parameter to add. If the parameter requires a flag,
            it can be provided as a tuple flag-value.
        :type parameter: Union[:class:`str`, :class:`tuple`]
        """
        if isinstance(parameter, (six.string_types, six.integer_types, float, np.float)):
            self._command.append(str(parameter))
        elif isinstance(parameter, (list, tuple)):
            self._command.append(parameter[0])
            self._command.append(str(parameter[1]))
        else:
            raise ValueError('Unexpected object type for new parameter')

    def clean_command( self ):
        """Removes all parameters added to the command.
        """
        self._command = []
        self._command.append(self._exec)

    def backup_command( self ):
        """Store a copy of the command up to that point to retrieve import
        afterwards.
        """
        self._backup_command = copy.deepcopy(self._command)

    def restore_command( self ):
        """Retrieve the backup command into the working command.
        """
        self._command        = self._backup_command
        self._backup_command = []

    def execute( self ):
        """Execute the command.

        :raises:
            :SystemError: If an error occurs during execution.
        """
        verbose = core.get_option('system', 'verbose')
        PIPE    = subprocess.PIPE if verbose > 1 else open('/dev/null', 'w')

        command = " ".join(self._command)
        if verbose >= 1:
            sys.stdout.write('Executing command:\t{0}\n'.format(command))

        try:
            p = subprocess.Popen(self._command, stdout=PIPE, stderr=PIPE)
            if verbose > 1:
                for out in iter(p.stdout.readline, b''):
                    sys.stdout.write(out.strip() + '\n')
            p.communicate()
        except Exception as err:
            raise SystemError(err)

    # MAGIC METHODS
    def __repr__( self ):
        return " ".join(self._command)

    def __str__( self ):
        return repr(self)


@six.add_metaclass(abc.ABCMeta)
class ExeBase( object ):
    """Base class for all external execution classes.
    """
    def __init__( self ):
        self._EXE = None

    # PRIVATE METHODS
    def _set_default_executable( self, executable_id ):
        """Create the executable directly from the configuration system.
        """
        self._EXE = ExeCommand(core.get_option('bin', executable_id))

    @staticmethod
    def _set_dynamic_executable( self, executable ):
        """Create the class:`.Executable` from a provided path.

        :param str executable: Full path to the application's executable.
        """
        self._EXE = ExeCommand(executable)
