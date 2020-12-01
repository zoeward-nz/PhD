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

0.0.1 2018 Initial version. 

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

# When piping stdout into head python raises an exception Ignore SIG_PIPE and don't throw exceptions on it... (http://docs.python.org/library/signal.html)
signal(SIGPIPE, SIG_DFL) 
__version__ = '0.0.1' 
__date__ = '2018' 
__email__ = 's.schmeier@gmail.com' 
__author__ = 'Sebastian Schmeier'

# For color handling on the shell
try:
    from colorama import init, Fore, Style
    # INIT color Initialise colours for multi-platform support.
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
    description = 'Read delimited files and sub-select rows of file 1 if string in field 1 are found in field 2 of any row in file 2, i.e. file 2 is the reference.'
    version = 'version {}, date {}'.format(__version__, __date__)
    epilog = 'Copyright {} ({})'.format(__author__, __email__)
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('--version',
                        action='version',
                        version='{}'.format(version))
    parser.add_argument(
        'str_file1',
        metavar='FILE1',
        help=
        'Delimited file. [if set to "-" or "stdin" reads from standard in]')
    parser.add_argument(
        'str_file2',
        metavar='FILE2',
        help=
        'Delimited file.')
    parser.add_argument('-d1',
                        metavar='STRING',
                        dest='delimiter_str1',
                        default='\t',
                        help='Delimiter used in file 1.  [default: "tab"]')
    parser.add_argument('-d2',
                        metavar='STRING',
                        dest='delimiter_str2',
                        default='\t',
                        help='Delimiter used in file 2.  [default: "tab"]')
    parser.add_argument('-f1',
                        metavar='INT',
                        type=int,
                        dest='field1',
                        default=1,
                        help='Field number in file 1.  [default: 1]')
    parser.add_argument('-f2',
                        metavar='INT',
                        type=int,
                        dest='field2',
                        default=1,
                        help='Field number in file 2.  [default: 1]')
    parser.add_argument('--header',
                        action="store_true",
                        dest='header',
                        default=False,
                        help='Print header file 1.')
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
    try:
        fileobj1 = load_file(args.str_file1)
    except:
        error('Could not load file 1 {}. EXIT.'.format(args.str_file1))
    try:
        fileobj2 = load_file(args.str_file2)
    except:
        error('Could not load file 2 {}. EXIT.'.format(args.str_file2))
    assert args.field1 - 1 >= 0
    assert args.field2 - 1 >= 0
        
    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wt')
    else:
        outfileobj = open(args.outfile_name, 'w')

    # delimited file handler
    csv_reader_obj = csv.reader(fileobj2, delimiter=args.delimiter_str2)
    d = {}
    for a in csv_reader_obj:
        d[a[args.field2-1]] = None
    
    csv_reader_obj = csv.reader(fileobj1, delimiter=args.delimiter_str1)
    if args.header:
        header = next(csv_reader_obj)
        outfileobj.write("{}\n".format(args.delimiter_str1.join(header)))
    for a in csv_reader_obj:
        if a[args.field1-1] in d:
            outfileobj.write("{}\n".format(args.delimiter_str1.join(a)))
    outfileobj.close()
    return 

if __name__ == '__main__':
    sys.exit(main())
