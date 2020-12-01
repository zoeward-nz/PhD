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

0.0.1    2017    Initial version.

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
import pysam
import pybedtools


# When piping stdout into head python raises an exception
# Ignore SIG_PIPE and don't throw exceptions on it...
# (http://docs.python.org/library/signal.html)
signal(SIGPIPE, SIG_DFL)

__version__ = '0.0.1'
__date__ = '2018'
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

def alert(atype, text, log):
    textout = '%s [%s] %s\n' % (time.strftime('%Y%m%d-%H:%M:%S'),
                                atype.rjust(7),
                                text)
    log.write('%s%s%s' % (colors[atype], textout, reset))
    if atype=='error': sys.exit()
        
def success(text, log=sys.stderr):
    alert('success', text, log)
    
def error(text, log=sys.stderr):
    alert('error', text, log)
    
def warning(text, log=sys.stderr):
    alert('warning', text, log)
    
def info(text, log=sys.stderr):
    alert('info', text, log)
    
##----------------------------------------------------------


def parse_cmdline():
    # parse cmd-line -----------------------------------------------------------
    sDescription = 'Statistics for each tx. Output: num_elements, num_zero, num_non-zero, num_greaterequal1, num_unique, max, min, sum, mean, median, std, var' 
    sVersion='version %s, date %s' %(__version__,__date__)
    sEpilog = 'Copyright %s (%s)' %(__author__, __email__)

    parser = argparse.ArgumentParser(description=sDescription,
                                      epilog=sEpilog)
    parser.add_argument('--version',
                        action='version',
                        version='%s' % (sVersion))
    parser.add_argument('str_file1',
                         metavar='FILE',
                         help='Expression matrix.')
    parser.add_argument('-o',
                        '--out',
                        metavar='STRING',
                        dest='outfile_name',
                        default=None,
                        help='Out-file.')
    
    # if no arguments supplied print help
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    if not args.outfile_name:
        parser.error('No outputfile given, specify with -o')
        
    return args, parser


def main():
    args, parser = parse_cmdline()

    try:
        samfile = pysam.AlignmentFile(args.str_file1, "rb")
    except:
        error('Could not load file {}. EXIT.'.format(args.str_file1))

    scalingfactor = 1000000./samfile.mapped
    sys.stderr.write("# mapped: %i, scaling factor (1M/mapped)): %f\n" %(samfile.mapped, scalingfactor))
    samfile.close()
    try:
        bed = pybedtools.BedTool(args.str_file1)
    except:
        error('Could not load file {}. EXIT.'.format(args.str_file1))
    
    b = bed.genome_coverage(bga=True, split=True, scale=scalingfactor)
    b.saveas(args.outfile_name)
        
    return

if __name__ == '__main__':
    sys.exit(main())
