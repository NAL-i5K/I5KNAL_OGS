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
import string
# try to import from project first
from os.path import dirname
if dirname(__file__) == '':
    lib_path = '../../lib'
else:
    lib_path = dirname(__file__) + '/../../lib'
sys.path.insert(1, lib_path)
from gff3_modified import Gff3
import function4gff
import intra_model
import single_feature

__version__ = '0.0.1'

COMPLEMENT_TRANS = string.maketrans('TAGCtagc', 'ATCGATCG')
def complement(seq):
    return seq.translate(COMPLEMENT_TRANS)

def get_subseq(gff, line):
    string = gff.fasta_external[line['seqid']]['seq'][(line['start']-1):(line['end'])-1]
    if line['strand'] == '-':
        string = complement(string[::-1])
    return string

def extract_start_end(gff, stype, dline):
    '''Extract seqeuces for a feature only use the Start and End information. The relationship between parent and children would be ignored.'''
    seq=dict()
    roots = [line for line in gff.lines if line['line_type'] == 'feature' and not line['attributes'].has_key('Parent')]
    if stype == 'pm':
        for root in roots:
            rid = 'NA'
            if root['attributes'].has_key('ID'):
                rid = root['attributes']['ID']
            children = root['children']
            for child in children:
                cid = 'NA'
                if child['attributes'].has_key('ID'):
                    cid = child['attributes']['ID']
                cname = cid
                if child['attributes'].has_key('Name'):
                    cname = child['attributes']['Name']
                defline='>{0:s}'.format(cid)
                if dline == 'complete':
                    defline = '>{0:s}:{1:d}..{2:d}:{3:s}|premature_transcript({4:s})|Parent={5:s}|ID={6:s}|Name={7:s}'.format(child['seqid'], child['start'], child['end'], child['strand'], child['type'], rid, cid, cname)
                seq[defline] = get_subseq(gff, child)
    elif stype == 'gene':
        for root in roots:
            rid = 'NA'
            if root['attributes'].has_key('ID'):
                rid = root['attributes']['ID']
            rname = rid
            if root['attributes'].has_key('Name'):
                rname = root['attributes']['ID']
            defline='>{0:s}'.format(rid)
            if dline == 'complete':
                defline = '>{0:s}:{1:d}..{2:d}:{3:s}|gene|ID={4:s}|Name={5:s}'.format(root['seqid'], root['start'], root['end'], root['strand'], rid, rname)
            seq[defline] = get_subseq(gff, root)
    elif stype == 'exon':
        exons = [line for line in gff.lines if line['type'] == 'exon']
        for exon in exons:
            eid = 'NA'
            if exon['attributes'].has_key('ID'):
                eid = exon['attributes']['ID']
            ename = eid
            if exon['attributes'].has_key('Name'):
                ename = exon['attributes']['Name']
            parents = exon['parents']
            plist = dict()
            for parent in parents:
                for p in parent:
                    plist[p['attributes']['ID']] = 1

            keys = plist.keys()
            pid = ','.join(keys)
            
            defline='>{0:s}'.format(eid)
            if dline == 'complete':
                defline = '>{0:s}:{1:d}..{2:d}:{3:s}|exon|Parent={5:s}|ID={6:s}|Name={7:s}'.format(exon['seqid'], exon['start'], exon['end'], exon['strand'], exon['type'], pid, eid, ename)

            seq[defline] = get_subseq(gff, exon)

    return seq
    
def main(gff_file=None, fasta_file=None, stype=None, dline=None):
    if not gff_file or not fasta_file or not stype:
        print('All of Gff file, fasta file, and type of extracted seuqences need to be specified')
        return
    logger_stderr.info('Reading files: {0:s}, {1:s}...'.format(gff_file, fasta_file))
    gff = Gff3(gff_file=gff_file, fasta_external=fasta_file, logger=logger_null)

    logger_stderr.info('Checking errors...')
    logger_stderr.warning('The extracted sequences might be wrong for the following features which have formatting errors...')
    gff.check_parent_boundary()
    gff.check_phase()
    gff.check_reference()
    print('ID\tError_Code\tError_Tag')
    single_feature.main(gff)
    intra_model.main(gff)
    error_set = function4gff.extract_internal_detected_errors(gff)
    

    seq=dict()
    if stype == 'pm' or stype == 'gene' or stype == 'exon':
        seq = extract_start_end(gff, stype, dline)        
    if len(seq):
        logger_stderr.info('Print out extracted sequences: {0:s}_{1:s}.fa...'.format(args.output_prefix, args.sequence_type))
        for k,v in seq.items():
            report_fh.write('{0:s}\n{1:s}\n'.format(k,v))

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
    Extract sequences from specific regions of genome based on gff file.
    Testing enviroment:
    1. Python 2.7

    Inputs:
    1. GFF3: reads from STDIN by default, may specify the file name with the -g argument
    2. fasta file: reads from STDIN by default, may specify the file name with the -f argument

    Outputs:
    1. Extract sequences from specific regions of genome based on gff file.

    """))
    parser.add_argument('-g', '--gff', type=str, help='Summary Report from Monica (default: STDIN)') 
    parser.add_argument('-f', '--fasta', type=str, help='File of typical errors (default: STDIN)')
    parser.add_argument('-st', '--sequence_type', type=str, help='Type of seuqences: please select from "g" - gene sequence for each record; "e" - exon sequence for each record; "pm" - premature transcripts; "m" - mature transcripts (only exons included); "cds"- coding sequences; "pep" - peptide seuqences.(default: STDIN)')
    parser.add_argument('-d', '--defline', type=str, help='"simple": only ID would be shown in the defline; "complete": complete information of the feature would be shown in the defline.')
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

    if args.sequence_type:
        logger_stderr.info('Specifying sequence type: (%s)...', args.sequence_type)
    elif not sys.stdin.isatty(): # if STDIN connected to pipe or file
        args.sequence_type = sys.stdin
        logger_stderr.info('Reading from STDIN...')
    else: # no input
        parser.print_help()
        sys.exit(1)

    if args.defline:
        logger_stderr.info('Defline format: (%s)...', args.defline)
    elif not sys.stdin.isatty(): # if STDIN connected to pipe or file
        args.defline = sys.stdin
        logger_stderr.info('Reading from STDIN...')
    else: # no input
        parser.print_help()
        sys.exit(1)

    if args.output_prefix:
        logger_stderr.info('Specifying prefix of output file name: (%s)...', args.output_prefix)
        fname = '{0:s}_{1:s}.fa'.format(args.output_prefix, args.sequence_type)
        report_fh = open(fname, 'wb')
    else:
        parser.print_help()
        sys.exit(1)

    main(args.gff, args.fasta, args.sequence_type, args.defline)
