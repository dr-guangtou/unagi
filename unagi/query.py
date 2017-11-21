#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

"""
Original script provided by HSC:

    https://hsc-gitlab.mtk.nao.ac.jp/snippets/13

Before running, you must set in your terminal session(in, e.g., BASH):

    export HSC_SSP_CAS_USER=<your STARS username>
    export HSC_SSP_CAS_PASSWORD=<your STARS password>

You probably want to call with:

    python hsc_query.py -s query-test.sql

"""

# Standard library
from os import path
import os
import sys
import ssl
import time
import json
import getpass

# Third-party
from six.moves.urllib.request import Request, urlopen
from six.moves.urllib.error import HTTPError
import tqdm

client_version = 20160120.1
API_URL = 'https://hscdata.mtk.nao.ac.jp/datasearch/api/catalog_jobs/'


class HSCQueryError(Exception):
    pass


def ssl_http_post(url, data, headers):
    req = Request(url, data, headers)
    skipVerifying = None

    try:
        skipVerifying = ssl.SSLContext(ssl.PROTOCOL_TLSv1)
    except AttributeError:
        pass

    if skipVerifying:
        res = urlopen(req, context=skipVerifying)
    else:
        res = urlopen(req)

    return res


def ssl_http_post_json(url, data, headers=None):
    if headers is None:
        headers = dict()

    headers['Content-type'] = 'application/json'
    data['clientVersion'] = client_version
    json_data = json.dumps(data).encode('utf-8')
    return ssl_http_post(url, json_data, headers)


class HSCQuery(object):

    def __init__(self, release_version,
                 preview=False,
                 username=None,
                 password=None,
                 user_env='HSC_SSP_CAS_USER',
                 password_env='HSC_SSP_CAS_PASSWORD',
                 nomail=True,
                 skip_syntax_check=True,
                 delete_job_on_complete=True):

        # Login credentials
        if username is None:
            username = os.environ.get(user_env, None)
            if username is None:
                getpass.getpass('username: ')
        self.username = username

        if password is None:
            password = os.environ.get(password_env, None)
            if password is None:
                getpass.getpass('password: ')
        self.password = password
        self._credentials = {
            'account_name': self.username,
            'password': self.password
        }

        self.release_version = release_version
        self.nomail = nomail
        self.skip_syntax_check = skip_syntax_check
        self.delete_job_on_complete = delete_job_on_complete
        self._current_job = None

    def submit_job(self, sql_text, output_format='fits'):
        url = path.join(API_URL, 'submit')
        catalog_job = {
            'sql': sql_text,
            'out_format': output_format,
            'include_metainfo_to_body': True,
            'release_version': self.release_version,
        }

        post_data = {
            'credential': self._credentials,
            'catalog_job': catalog_job,
            'nomail': self.nomail,
            'skip_syntax_check': self.skip_syntax_check
        }
        res = ssl_http_post_json(url, post_data)
        job = json.loads(res.read().decode('utf-8'))

        return job

    def job_status(self, job_id):
        url = path.join(API_URL, 'status')
        post_data = {'credential': self._credentials, 'id': job_id}
        res = ssl_http_post_json(url, post_data)
        job = json.loads(res.read().decode('utf-8'))
        return job

    def cancel_job(self, job_id):
        url = path.join(API_URL, 'cancel')
        postData = {'credential': self._credentials, 'id': job_id}
        ssl_http_post_json(url, postData)

    def delete_job(self, job_id):
        url = path.join(API_URL, 'delete')
        postData = {'credential': self._credentials, 'id': job_id}
        ssl_http_post_json(url, postData)

    def download_job_results(self, job_id, output_file=None):
        url = path.join(API_URL, 'download')
        post_data = {'credential': self._credentials, 'id': job_id}
        res = ssl_http_post_json(url, post_data)

        if output_file is None:
            out = sys.stdout

        else:
            out = open(output_file, 'w')

        buffer_size = 64 * 1<<10  # 64k
        while True:
            buf = res.read(buffer_size)
            out.write(buf.decode(out.encoding))
            if len(buf) < buffer_size:
                break

        out.close()

    def run(self, sql_text, output_file=None, output_format='fits',
            overwrite=False):
        if (output_file is not None and
                path.exists(output_file) and not overwrite):
            raise IOError('File already exists ({}) -- use "overwrite".'
                          .format(output_file))

        self._current_job = self.submit_job(sql_text,
                                            output_format=output_format)

        # Block until the job finishes
        max_time = 600  # 10 minutes
        interval = 10   # check every 10 seconds

        run_time = 0
        for i in range(max_time // interval+1):
            for i in tqdm.tqdm(range(1, interval+1)):
                time.sleep(1)

            self._current_job = self.job_status(self._current_job['id'])

            if self._current_job['status'] == 'error':
                raise HSCQueryError('Query error: ' +
                                    self._current_job['error'])

            if self._current_job['status'] == 'done':
                break

            print('query still running...')
            run_time += interval

        self.download_job_results(self._current_job['id'],
                                  output_file=output_file)

        if self.delete_job_on_complete:
            self.delete_job(self._current_job['id'])

        self._current_job = None


def main(**kwargs):
    query = HSCQuery(kwargs['release_version'],
                     user_env=kwargs['user_env'],
                     password_env=kwargs['password_env'],
                     nomail=kwargs['nomail'],
                     skip_syntax_check=kwargs['skip_syntax_check'],
                     delete_job_on_complete=kwargs['delete_job'])

    with open(kwargs['sql_file'], 'r') as f:
        sql_text = f.read()

    try:
        query.run(sql_text, output_file=kwargs['output_file'],
                  output_format=kwargs['output_format'])
    except HTTPError as e:
        if e.code == 401:
            raise ValueError('invalid id or password.')
        if e.code == 406:
            raise ValueError(e.read())
        else:
            raise RuntimeError(e.read())

    except HSCQueryError as e:
        raise RuntimeError(e.read())

    except KeyboardInterrupt:
        if query._current_job is not None:
            query.cancel_job(query._current_job['id'])
        raise

    else:
        sys.exit(0)

    sys.exit(1)


if __name__ == "__main__":
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

    # Define parser object
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)

    parser.add_argument('--user-env', default='HSC_SSP_CAS_USER',
                        help='The name of the environment variable that '
                             'contains your STARS account username.')
    parser.add_argument('--password-env', default='HSC_SSP_CAS_PASSWORD',
                        help='The name of the environment variable that '
                             'contains your STARS account password.')
    parser.add_argument('--release-version', '-r',
                        choices='dr1 dr_early'.split(),
                        default='dr1', help='Data release version.')
    parser.add_argument('--nomail', '-M', action='store_true', default=True,
                        help='Suppress email notice.')
    parser.add_argument('--format', '-f', dest='output_format', default='csv',
                        choices=['csv', 'csv.gz', 'sqlite3', 'fits'],
                        help='The output file format.')
    parser.add_argument('--sql-file', '-s', dest='sql_file', required=True,
                        help='An SQL file containing the query to run.')
    parser.add_argument('--output', '-o', dest='output_file',
                        default=None, type=str,
                        help='The output name of the output file '
                             '(with or without an extension).')
    parser.add_argument('--skip-syntax-check', '-S', action='store_true',
                        default=True, help='Skip syntax checking the query.')
    parser.add_argument('--delete-job', '-d', action='store_true',
                        default=True, help='Delete the job after runnng.')

    args = parser.parse_args()

    main(**vars(args))
