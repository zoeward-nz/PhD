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
    description = 'Read gtf-file and extract tx to gene associations based on class codes for novels and genes. If a code is not in novels or genes, it is ignored. Extract tx to gene info from reference.'
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
        'GFFCOMPARE gtf file. [if set to "-" or "stdin" reads from standard in]')
    parser.add_argument(
        'str_file2',
        metavar='FILE',
        help=
        'GTF reference file, e.g. from gencode.')
    parser.add_argument('-n',
                        '--novel',
                        metavar='STRING',
                        dest='novels',
                        default="i,x,s,u,y",
                        help='Type codes to consider as novel. [default: "i,x,s,u,y"]')
    parser.add_argument('-g',
                        '--genes',
                        metavar='STRING',
                        dest='genes',
                        default="=,j,k,c,o",
                        help='Type codes to consider as novel. [default: "=,j,k,c,o"]')
    parser.add_argument('-o',
                        '--out',
                        metavar='STRING',
                        dest='outfile_name',
                        default=None,
                        help='Out-file. [default: "stdout"]')
    parser.add_argument('--gtf',
                        metavar='STRING',
                        dest='outfile_name_gtf',
                        default=None,
                        help='Out-file cleaned GTF.')

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
        fileobj = load_file(args.str_file)
    except:
        error('Could not load file. EXIT.')

    try:
        fileobj_ref = load_file(args.str_file2)
    except:
        error('Could not load ref file. EXIT.')

    # create outfile object
    if not args.outfile_name:
        outfileobj = sys.stdout
    elif args.outfile_name in ['-', 'stdout']:
        outfileobj = sys.stdout
    elif args.outfile_name.split('.')[-1] == 'gz':
        outfileobj = gzip.open(args.outfile_name, 'wt')
    else:
        outfileobj = open(args.outfile_name, 'w')

    outfileobjgtf = None
    if args.outfile_name_gtf:
        if args.outfile_name_gtf.split('.')[-1] == 'gz':
            outfileobjgtf = gzip.open(args.outfile_name_gtf, 'wt')
        else:
            outfileobjgtf = open(args.outfile_name_gtf, 'w')


    novels = [ s.strip() for s in args.novels.split(',') ]
    genes = [ s.strip() for s in args.genes.split(',') ]

    # parsing anno
    regRef = re.compile('gene_id\s"(.+?)";.+transcript_id\s"(.+?)";')
    dRef = {}
    for line in fileobj_ref:
        res = regRef.search(line)
        if not res:
            continue
        else:
            gene, tx = res.groups()
            dRef[tx] = gene
    fileobj_ref.close()      
    
    ## KNONW
    # GL000194.1      StringTie       transcript      55986   115057  .       -       .       transcript_id "MSTRG.52.3"; gene_id "MSTRG.52"; gene_name "MAFIP"; xloc "XLOC_000052"; cmp_ref "ENST00000400754.4"; class_code "j"; tss_id "TSS52";
    # GL000194.1      StringTie       exon    55986   62949   .       -       .       transcript_id "MSTRG.52.3"; gene_id "MSTRG.52"; exon_number "1";

    ## NOVEL
    # GL000008.2      StringTie       transcript      1       7898    .       +       .       transcript_id "MSTRG.4.1"; gene_id "MSTRG.4"; xloc "XLOC_000001"; class_code "u"; tss_id "TSS1";
    # GL000008.2      StringTie       exon    1       7898    .       +       .       transcript_id "MSTRG.4.1"; gene_id "MSTRG.4"; exon_number "1";

    
    regTxNovel = re.compile('transcript_id\s"(.+?)";.+gene_id\s"(.+?)";.+class_code "(.)"')
    regTxKnown = re.compile('transcript_id\s"(.+?)";.+gene_id\s"(.+?)";.+cmp_ref\s"(.+?)";.+class_code "(.)"')
    regExon = re.compile('transcript_id\s"(.+?)";.+gene_id\s"(.+?)";')

    dTX2NAME = {}
    d = {}
    dremove = {}
    
    reader = csv.reader(fileobj, delimiter='\t')
    for a in reader:
        tx, gene, cmp_ref, code = None, None, None, None

        if a[2] == "exon":
            res = regExon.search(a[8])
            assert res
            tx, gene = res.groups()
            if tx in dremove:
                continue
            else:
                if outfileobjgtf:
                    if tx in dTX2NAME:
                        # we need to adjust the gene_id of the exon if its tx is of a known gene
                        a[8] = re.sub('gene_id ".+?";', 'gene_id "{}";'.format(dTX2NAME[tx]), a[8])
                    outfileobjgtf.write('{}\n'.format('\t'.join(a)))
                    
        elif a[2] == "transcript":
            res = regTxKnown.search(a[8])
            if res:
                tx, gene, cmp_ref, code = res.groups()
            else:
                res = regTxNovel.search(a[8])
                assert res
                tx, gene, code = res.groups()
                cmp_ref = None

            # now based on known or novel codes
            if code in genes:
                # replace the gene_id with gene_id from refernce before printing gtf
                assert cmp_ref
                assert cmp_ref in dRef
                geneid = dRef[cmp_ref]
                d[(geneid,tx,code)] = None
                dTX2NAME[tx] = geneid  # if tx in here it is a known gene and we need to replace gene_id in exon too
                a[8] = re.sub('gene_id ".+?";', 'gene_id "{}";'.format(geneid), a[8])
            elif code in novels:
                # we are looking at a novel so we keep the stringtie asigned id
                d[(gene,tx,code)] = None
            else:
                dremove[tx] = None
                continue
                
            if outfileobjgtf:
                outfileobjgtf.write('{}\n'.format('\t'.join(a)))

        else:
            continue

    a = list(d.keys())
    a.sort()
    for t in a:
        gene,tx,code = t
        outfileobj.write('{}\t{}\t{}\n'.format(tx, gene, code))


    # ------------------------------------------------------
    outfileobj.close()
    return


if __name__ == '__main__':
    sys.exit(main())

