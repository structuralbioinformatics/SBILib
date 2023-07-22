# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Structural BioInformatics Lab <sbi.upf.edu>
    Baldo Oliva <baldo.oliva@upf.edu>

.. class:: CATHLink
"""

# Standard Libraries
from tempfile import NamedTemporaryFile
import os

# External Libraries
import requests
import pandas as pd

# This Library
import SBI.core as core


__all__ = ['CATHLink']


class CATHLink( object ):
    """Controls data retrieval from the CATH_ database.

    .. note::
        This functionality depends on the :ref:`global configuration options <configuration>` ``db.cath``,
        that points to the top level url of the downloadable files of the database. The link will target
        the latest release of  the database.

    .. _CATH: http://www.cathdb.info/
    """
    def __init__( self ):
        self.__name__ = 'databases.CATHLink'

    def get_data( self ):
        """Retrieve data from the declared SCOP version.

        :return: :class:`~pandas.DataFrame`
        """
        url = '/'.join([core.get_option('db', 'cath'), 'cath-classification-data', 'cath-domain-boundaries.txt'])
        local_file = _retrieve_data(url)
        df = pd.read_fwf(local_file.name, comment='#', header=None, delimiter='~',
                         names=['chain_name', 'domain_count', 'frag_count', 'data'],
                         colspecs=[(0, 5), (6, 9), (9, 13), (13, None)])
        df[['domain_count']] = df['domain_count'].str.extract(r'(\d+)')
        df['domain_count'] = pd.to_numeric(df['domain_count'])
        df['data'] = df.apply(lambda row: process_CATH_data(row), axis=1)
        df = _split_dataframe_rows(df.drop(columns=['frag_count']), ['data']).rename(columns={'data': 'selectors'})
        df['chain_name'] = chain_rename(df['chain_name'])
        df['pdb'] = df['chain_name'].str[:4]
        df = df.rename(columns={'chain_name': 'domain_id'}).drop(columns=['domain_count'])
        df['selectors'] = df['selectors'].map(tuple)
        return df


def chain_rename(column):
    count = 1
    last = ''
    new_column = []
    for c in column:
        if c != last:
            last = c
            count = 1
        else:
            count += 1
        new_column.append(c.strip() + '{:02d}'.format(count))
    return new_column


def process_CATH_data( row ):
    total = row['domain_count']
    data = row['data'].split()
    domains = []
    segments = []
    status = 'count'
    new_status_count = 0
    positions = []
    for i, d in enumerate(data):
        if i == new_status_count:
            status = 'count'
            if len(domains) == total:
                break
        if status == 'count':
            domains.append([[], ] * int(d))
            segments.append((int(d), i))
            status = 0
            new_status_count = i + (6 * int(d)) + 1
            for x in range(int(d)):
                positions.append((len(domains) - 1, x))
            continue
    for x in reversed(segments):
        data.pop(x[1])
    segments = [x[0] for x in segments]
    for i, x in enumerate(range(0, len(data), 6)):
        if i >= sum(segments):
            break
        cdata = data[x:x + 6]
        if cdata[0] != cdata[3]:
            raise Exception('whaaat?')
        strdata = '{0}:{1}{2};{3}:{4}{5}'.format(*cdata).replace('-;', ';').rstrip('-')
        strdata = '{0[0]}:{0[1]}-{1[1]}'.format(*[x.split(':') for x in strdata.split(';')])
        domains[positions[i][0]][positions[i][1]] = strdata
    return domains


def _split_dataframe_rows( df, columns=None ):
    # modified from: https://gist.github.com/jlln/338b4b0b55bd6984f883#gistcomment-2321628
    def _split_list_to_rows_inner(row, row_accumulator, column_selectors):
        split_rows = {}
        max_split = 0
        for column_selector in column_selectors:
            split_row = row[column_selector]
            split_rows[column_selector] = split_row
            if len(split_row) > max_split:
                max_split = len(split_row)

        for i in range(max_split):
            new_row = row.to_dict()
            for column_selector in column_selectors:
                try:
                    new_row[column_selector] = split_rows[column_selector].pop(0)
                except IndexError:
                    new_row[column_selector] = ''
            row_accumulator.append(new_row)

    columns = list(df.columns) if columns is None else columns
    new_rows = []
    df.apply(_split_list_to_rows_inner, axis=1, args=(new_rows, columns))
    new_df = pd.DataFrame(new_rows, columns=df.columns)
    return new_df


def _retrieve_data( url ):
    """Download the url into a tempfile.
    """
    if core.get_option('system', 'verbose') == 2:
        print('Downloading from {} '.format(url))
    local_file = NamedTemporaryFile(delete=core.get_option('system', 'verbose') < 2)
    r = requests.get(url, stream=True)
    if r.status_code != requests.codes.ok:
        raise IOError('Remote path {} cannot be found'.format(url))
    for chunk in r.iter_content(chunk_size=1024):
        if chunk:
            local_file.write(chunk)
    if core.get_option('system', 'verbose') == 2:
        print('into temporary file {}'.format(local_file.name))
    local_file.flush()
    os.fsync(local_file.fileno())
    return local_file
