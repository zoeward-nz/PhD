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

0.0.1 2017 Initial version. 

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
import numpy as np 
import pandas as pd

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
    # INIT color Initialise colours for multi-platform support.
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
    parser.add_argument('str_file2',
                         metavar='FILE',
                         help='GTF-file reference')
    parser.add_argument('-d', '--delimiter',
                         metavar='STRING',
                         dest='sDelim',
                         default='\t',
                         help='Delimiter for the columns of the expression matrix file. [default: tab]')
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
        filehandle = bz2.BZFile(filename)
    elif filename.split('.')[-1] == 'zip':
        filehandle = zipfile.Zipfile(filename)
    else:
        filehandle = open(filename)
    return filehandle 

def main():
    args, parser = parse_cmdline()
    try:
        fileobj1 = load_file(args.str_file1)
    except:
        error('Could not load file 1 {}. EXIT.'.format(args.str_file1))
    try:
        fileobj2 = load_file(args.str_file2)
    except:
        error('Could not load file 2 {}. EXIT.'.format(args.str_file2))
    
    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wt')
    else:
        outfileobj = open(args.outfile_name, 'w')
    
   
    # collect info from gffcompare file
    csv_reader = csv.reader(fileobj2, delimiter = '\t')
    header = next(csv_reader)
    d = {}
    for a in csv_reader:
        tx = a[4]
        num_exons = a[5]
        classTx = a[2]
        length = a[9]
        d[tx] = [length, classTx, num_exons]
   
    df = pd.read_table(fileobj1, header=0, index_col=0)
    df2 = pd.DataFrame()
    df2['no.SamplesNull'] = df[df==0].count(axis=1)
    df2['no.SamplesNotnull'] = df[df>0].count(axis=1)
    df2['no.SamplesGE1'] = df[df>=1].count(axis=1)
    df2['max'] = df.max(axis=1)
    df2['min'] = df.min(axis=1)
    df2['mean'] = df.mean(axis=1)
    df2['std'] = df.std(axis=1)
    df2['quant25%'] = df.quantile(q=0.25, axis=1)
    df2['quant50%'] = df.quantile(q=0.5, axis=1)
    df2['quant75%'] = df.quantile(q=0.75, axis=1)
    del(df)
    new_header = ["Tx", "txLength", "txClass", "no.Exons"] + list(df2.columns)
    outfileobj.write("{}\n".format("\t".join(new_header)))
    for index, row in df2.iterrows():
        if index not in d:
            warning("{} not in reference file.".format(index))
            res = [index, "NA","NA","NA"] + list(row)
        else:
            a = d[index]
            res = [index] + a + list(row)
        res = [str(x) for x in res]
        outfileobj.write("{}\n".format("\t".join(res)))
        
        
    return 

if __name__ == '__main__':
    sys.exit(main())
