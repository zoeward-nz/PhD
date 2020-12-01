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
    
## COMMON FUNCTIONS:
def sub_slice(list, list_with_indeces):
    """
    Using list comprehensions to build slice of list
    """
    return [list[i] for i in list_with_indeces] def sub_slice_comp(list, list_with_indeces):
    """
    Build slice of list
    Complement slice. All elements of list that are not in the list of indeces.
    """
    return [list[i] for i in range(len(list)) if i not in list_with_indeces] def uniq(aList, sort=False):
    """
    Make a list unique.
    @param aList: Python list
    """
    if sort:
        return sorted(list(set(aList)))
    else:
        return list(set(aList))
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
    parser.add_argument('sFile',
                         metavar='FILE',
                         help='Expression matrix.')
    parser.add_argument('sFile2',
                         metavar='FILE',
                         help='BED-file reference')
    parser.add_argument('sFileOut',
                         metavar='FILE',
                         help='Out-file')
    parser.add_argument('sFileOutStats',
                         metavar='FILE',
                         help='Out-file stats')
    parser.add_argument('-d', '--delimiter',
                         metavar='STRING',
                         dest='sDelim',
                         default='\t',
                         help='Delimiter for the columns of the expression matrix file. [default: tab]')
    
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

def open_file_write(filename):
    if filename.split('.')[-1] == 'gz':
        outfileobj = gzip.open(filename, 'wt')
    else:
        outfileobj = open(filename, 'w')
    return outfileobj 

def main():
    args, parser = parse_cmdline()
    try:
        fileobj2 = load_file(args.sFile2)
    except:
        error('Could not load file {}. EXIT.'.format(args.sFile2))
    out1 = open_file_write(args.sFileOut)
    out2 = open_file_write(args.sFileOutStats)
        
    # collect info from gtf
    csv_reader = csv.reader(fileobj2, delimiter = '\t')
    d = collections.OrderedDict()
    reg1 = re.compile('transcript_id "(.+?)";.+class_code "(.)";')
    reg2 = re.compile('transcript_id "(.+?)";')
    for a in csv_reader:
        if a[2] == "transcript":
            res = reg1.search(a[8])
            tx = res.group(1)
            classTx = res.group(2)
            length = int(a[4]) - int(a[3]) + 1
            if tx in d:
                error("{} already seen.".format(tx))
            else:
                d[tx] = [length, classTx, []]
        elif a[2] == "exon":
            res = reg2.search(a[8])
            tx = res.group(1)
            if tx not in d:
                error("{} parent transcript of exon not seen yet.".format(tx))
            else:
                length = int(a[4]) - int(a[3]) + 1
                d[tx][2].append(length)
    a = []
    for tx in d:
        exons = d[tx][2]
        maxE = max(exons)
        minE = min(exons)
        mean = np.mean(exons)
        std = np.std(exons, ddof=0) # unbiased N-1
        num = len(exons)
        a.append([tx, d[tx][0], d[tx][1], num, maxE, minE, mean, std])
    df_ref = pd.DataFrame(a, columns=["Tx", "txLength", "txClass", "#Exons", "maxExon", "minExon", "meanExons", "stdExons"])
    df_ref.index = df_ref["Tx"]
    del df_ref['Tx']
    # parse expression table
    try:
        fileobj = load_file(args.sFile)
    except:
        error('Could not load file {}. EXIT.'.format(args.sFile))
    df = pd.read_table(fileobj,header=0, index_col=0)
    df_subset = df.loc[df.index.intersection(list(df_ref.index))] ## subselect the ones in the gtf
    df2 = pd.DataFrame()
    df2['#null'] = df_subset[df_subset==0].count(axis=1)
    df2['#notnull'] = df_subset[df_subset>0].count(axis=1)
    df2['#ge1'] = df_subset[df_subset>=1].count(axis=1)
    df2['max'] = df_subset.max(axis=1)
    df2['min'] = df_subset.min(axis=1)
    df2['mean'] = df_subset.mean(axis=1)
    df2['std'] = df_subset.std(axis=1)
    df2['25%'] = df_subset.quantile(q=0.25, axis=1)
    df2['50%'] = df_subset.quantile(q=0.5, axis=1)
    df2['75%'] = df_subset.quantile(q=0.75, axis=1)
    df2 = pd.concat([df2, df_ref], axis=1, sort=False)
    # out
    df_subset.reset_index(level=0, inplace=True)
    df_subset.columns.values[0] = "Tx"
    df_subset.to_csv(out1, sep='\t', index=False, header=True)
    df2.reset_index(level=0, inplace=True)
    df2.columns.values[0] = "Tx"
    df2.to_csv(out2, sep='\t', header=True, index=False)
    out1.close()
    out2.close()
    return if __name__ == '__main__':
    main()
