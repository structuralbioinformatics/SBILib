# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>
.. codeaffiliation:: Structural BioInformatics Lab <baldo.oliva@upf.edu>
"""

# Standard Libraries

# External Libraries
import six
import request
import pandas as pd

# This Library
import SBI.core as core

class PDBeChemConnect( object ):
    @staticmethod
    def retrieve_database():
        url = core.get_option('data', 'pdbechem')