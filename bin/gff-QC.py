#! /usr/local/bin/python2.7
# Copyright (C) 2015  Mei-Ju Chen <arbula [at] gmail [dot] com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

import sys
import re
import logging
# try to import from project first
from os.path import dirname
if dirname(__file__) == '':
    lib_path = '../lib'
else:
    lib_path = dirname(__file__) + '/../lib'
sys.path.insert(1, lib_path)
from gff3_modified import Gff3
import single_feature
import inter_model
import intra_model

__version__ = '0.0.1'

if __name__ == '__main__':
    logger_stderr = logging.getLogger(__name__+'stderr')
    logger_stderr.setLevel(logging.INFO)
    stderr_handler = logging.StreamHandler()
    stderr_handler.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger_stderr.addHandler(stderr_handler)
    logger_null = logging.getLogger(__name__+'null')
    null_handler = logging.NullHandler()
    logger_null.addHandler(null_handler)
    import argparse
    from textwrap import dedent
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=dedent("""\
    
    Testing enviroment:
    1. Python 2.7

    Inputs:
    1. GFF3: reads from STDIN by default, may specify the file name with the -g argument
    2. fasta file: reads from STDIN by default, may specify the file name with the -f argument

    Outputs:
    1. Extract sequences from specific regions of genome based on gff file.

    """))
    parser.add_argument('-g', '--gff', type=str, help='Genome annotation file, gff3 format (default: STDIN)') 
    parser.add_argument('-f', '--fasta', type=str, help='Genome sequences, fasta format (default: STDIN)')
    parser.add_argument('-sc', '--species_code', type=str, help='If the species is hosted by I5K Workspace@NAL, you can give the I5K species code to have a url link to Web Apollo in I5K Workspace@NAL. eg. "lepdec", which was abbreviated from Leptinotarsa decemlineata. (default: STDIN)')
    parser.add_argument('-o', '--output_prefix', type=str, help='Prefix of output file name (default: STDIN)')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
    
    args = parser.parse_args()

    if args.gff:
        logger_stderr.info('Checking gff file (%s)...', args.gff)
    elif not sys.stdin.isatty(): # if STDIN connected to pipe or file
        args.gff = sys.stdin
        logger_stderr.info('Reading from STDIN...')
    else: # no input
        parser.print_help()
        sys.exit(1)

    if args.fasta:
        logger_stderr.info('Checking genome fasta (%s)...', args.fasta)
    elif not sys.stdin.isatty(): # if STDIN connected to pipe or file
        args.fasta = sys.stdin
        logger_stderr.info('Reading from STDIN...')
    else: # no input
        parser.print_help()
        sys.exit(1)

    if args.species_code:
        logger_stderr.info('Specifying species code: (%s)...', args.species_code)
    elif not sys.stdin.isatty(): # if STDIN connected to pipe or file
        args.species_code = sys.stdin
        logger_stderr.info('Reading from STDIN...')

    if args.output_prefix:
        logger_stderr.info('Specifying prefix of output file name: (%s)...\n', args.output_prefix)
        fname = '{0:s}_{1:s}.fa'.format(args.output_prefix, args.sequence_type)
        report_fh = open(fname, 'wb')


    ERROR_CODE = ['Esf0001', 'Ema0005', 'Emn0001'] 
    ERROR_TAG = ['pseudogene or not?', 'unusual child features in the type of pseudogene found', 'Duplicate transcripts found']
    ERROR_INFO = dict(zip(ERROR_CODE, ERROR_TAG))

    logger_stderr.info('Reading gff files: (%s)...\n', args.gff)
    gff3 = Gff3(gff_file=args.gff, logger=logger_null)
    logger_stderr.info('Checking missing attributes: (%s)...\n', 'single_feature.FIX_MISSING_ATTR()')
    single_feature.FIX_MISSING_ATTR(gff3, logger_stderr=logger_stderr)

    roots = [line for line in gff3.lines if line['line_type']=='feature' and not line['attributes'].has_key('Parent')]
    error_set=dict()

    logger_stderr.info('Checking whether pseudogene has weird child types: (%s)...', 'intra_model.pseudo_child_type()')
    trans_list = list()
    for root in roots:
        r = intra_model.pseudo_child_type(gff3, root)
        if not r == None:
            error_set = dict(error_set.items() + r.items())

        children = root['children']
        for child in children:
            trans_list.append(child)

    logger_stderr.info('Checking duplicate transcripts: (%s)...', 'inter_model.check_duplicate()')
    r = inter_model.check_duplicate(gff3, trans_list)
    if not r == None: 
        error_set = dict(error_set.items() + r.items())

    # Find out models with Note mentioned about pseudogene
    logger_stderr.info('Checking models with Note mentioned about pseudogene: (%s)...', 'single_feature.detect_pseudogene()')
    single_feature.FIX_PSEUDOGENE(gff3)
    for root in roots:
        r = single_feature.detect_pseudogene(gff3, root)
        if not r == None:
            error_set = dict(error_set.items() + r.items())
        children = root['children']
        for child in children:
            r = single_feature.detect_pseudogene(gff3, child)
            if not r == None:
                error_set = dict(error_set.items() + r.items())

    print '\n'
    print 'Transcript_ID\tError_code\tError_tag'
    for k,v in error_set.items():
        for e in v:
            print '{0:s}\t{1:s}\t[{2:s}]'.format(k, e, ERROR_INFO[e])
 
