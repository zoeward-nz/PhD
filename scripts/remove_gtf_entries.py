#!/usr/bin/env python
"""
NAME: test
=========

DESCRIPTION
===========

INSTALLATION
============

USAGE
=====

VERSION HISTORY
===============

0.0.3    201808  Removed the ref_gene_id extraction again.
0.0.2    201808  Get ref_gene_id if available instead of gene_id
0.0.1    2018    Initial version.

LICENCE
=======
2017, copyright Sebastian Schmeier, (s.schmeier@gmail.com), https://sschmeier.com

template version: 1.9 (2017/12/08)
"""
from signal import signal, SIGPIPE, SIG_DFL
import sys
import os
import os.path
import argparse
import csv
import collections
import gzip
import bz2
import zipfile
import time
import re

# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
signal(SIGPIPE, SIG_DFL)

__version__ = '0.0.3'
__date__ = '20180806'
__email__ = 's.schmeier@gmail.com'
__author__ = 'Sebastian Schmeier'

# For color handling on the shell
try:
    from colorama import init, Fore, Style
    # INIT color
    # Initialise colours for multi-platform support.
    init()
    reset=Fore.RESET
    colors = {'success': Fore.GREEN, 'error': Fore.RED, 'warning': Fore.YELLOW, 'info':''}
except ImportError:
    sys.stderr.write('colorama lib desirable. Install with "conda install colorama".\n\n')
    reset=''
    colors = {'success': '', 'error': '', 'warning': '', 'info':''}

def alert(atype, text, log, repeat=False):
    if repeat:
        textout = '{} [{}] {}\r'.format(time.strftime('%Y%m%d-%H:%M:%S'),
                                        atype.rjust(7),
                                        text)
    else:
        textout = '{} [{}] {}\n'.format(time.strftime('%Y%m%d-%H:%M:%S'),
                                        atype.rjust(7),
                                        text)

    
    log.write('{}{}{}'.format(colors[atype], textout, reset))
    if atype=='error': sys.exit()

def success(text, log=sys.stderr):
    alert('success', text, log)

def error(text, log=sys.stderr):
    alert('error', text, log)

def warning(text, log=sys.stderr):
    alert('warning', text, log)
    
def info(text, log=sys.stderr, repeat=False):
    alert('info', text, log)  

    
def parse_cmdline():
    """ Parse command-line args. """
    ## parse cmd-line -----------------------------------------------------------
    description = 'Read gtf-file and tx entries and remove all tx and exons from gtf.'
    version = 'version {}, date {}'.format(__version__, __date__)
    epilog = 'Copyright {} ({})'.format(__author__, __email__)

    parser = argparse.ArgumentParser(description=description, epilog=epilog)

    parser.add_argument('--version',
                        action='version',
                        version='{}'.format(version))

    parser.add_argument(
        'str_file',
        metavar='FILE',
        help=
        'GTF-file.')
    parser.add_argument(
        'str_file_tx',
        metavar='FILE',
        help=
        'TX entries to remove from GTF-file (also remove all exons of TX).  [if set to "-" or "stdin" reads from standard in]')
    parser.add_argument('-o',
                        '--out',
                        metavar='STRING',
                        dest='outfile_name',
                        default=None,
                        help='Out-file. [default: "stdout"]')

    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    return args, parser


def load_file(filename):
    """ LOADING FILES """
    if filename in ['-', 'stdin']:
        filehandle = sys.stdin
    elif filename.split('.')[-1] == 'gz':
        filehandle = gzip.open(filename, 'rt')
    elif filename.split('.')[-1] == 'bz2':
        filehandle = bz2.open(filename, 'rt')
    elif filename.split('.')[-1] == 'zip':
        filehandle = zipfile.ZipFile(filename)
    else:
        filehandle = open(filename)
    return filehandle


def main():
    """ The main funtion. """
    args, parser = parse_cmdline()

    if args.str_file in ["stdin", "-"]:
        error("Proper file expected for GTF-file. EXIT")
    try:
        fileobj_gtf = load_file(args.str_file)
    except:
        error('Could not load GTF-file {}. EXIT.'.format(args.str_file))

    try:
        fileobj_entries = load_file(args.str_file_tx)
    except:
        error('Could not load file {}. EXIT.'.format(args.str_file_tx))

    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wt')
    else:
        outfileobj = open(args.outfile_name, 'w')

    remove = {}
    reader = csv.reader(fileobj_entries, delimiter = '\t')
    for a in reader:
        remove[a[0]] = None
        
    regexp = re.compile('transcript_id\s"(.+?)";')

    i = 0
    for line in fileobj_gtf:
        if line[0] == "#":
            continue
        line = line.rstrip()
        res = regexp.search(line)
        assert res
        tx = res.groups()[0]
        if tx not in remove:
            outfileobj.write('{}\n'.format(line))
        else:
            i += 1

    info('{} entries removed from gtf-file.'.format(i))
    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == '__main__':
    sys.exit(main())

