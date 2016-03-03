#! /usr/local/bin/python2.7
# -*- coding: utf-8 -*-
# Copyright (C) 2016  Mei-Ju May Chen <arbula [at] gmail [dot] com>

"""
QC functions for processing multiple features between models (inter-model) in GFF3 file.
"""
from __future__ import print_function

#from collections import OrderedDict # not available in 2.6
from collections import defaultdict
from itertools import groupby
try:
    from urllib import quote, unquote
except ImportError:
    from urllib.parse import quote, unquote
from textwrap import wrap
import sys
import re
import logging
logger = logging.getLogger(__name__)
#log.basicConfig(level=logging.DEBUG, format='%(levelname)-8s %(message)s')
logger.setLevel(logging.INFO)
if not logger.handlers:
    lh = logging.StreamHandler()
    lh.setFormatter(logging.Formatter('%(levelname)-8s %(message)s'))
    logger.addHandler(lh)
from os.path import dirname
if dirname(__file__) == '':
    lib_path = '../../lib'
else:
    lib_path = dirname(__file__) + '/../../lib'
sys.path.insert(1, lib_path)
from gff3_modified import Gff3
import intra_model 

__version__ = '0.0.1'

def FIX_MISSING_ATTR(gff): 
    features = [line for line in gff.lines if line['line_type']=='feature']
    flag = 0
    for f in features:
        if not f['attributes'].has_key('owner'):
            f['attributes']['owner'] = 'Unassigned'
        if not f['attributes'].has_key('ID'):
            IDrequired = ['gene', 'pseudogene', 'mRNA', 'pseudogenic_transcript', 'exon', 'pseudogenic_exon']
            if f['type'] in IDrequired:
                logger_stderr.error('[Missing ID] A model needs to have a unique ID, but not. Please fix it first.\n{0:s}'.format(f['line_raw']))
                flag += 1
            else:
                if len(f['parents'])== 1 and len(f['parents'][0]) == 1:
                    tid = f['parents'][0][0]['attributes']['ID'] + '-' + f['type']
                    f['attributes']['ID'] = tid
                else:
                    logger_stderr.error('[Missing ID] The program try to automatically generate ID for this model, but failed becuase this model has multiple parent features.\n{0:s}'.format(f['line_raw']))
    if flag != 0:
        sys.exit()

def check_duplicate(gff, linelist):
    '''
    This function assumes that,
    1. Each gnee is unique
    2. Children features such as Exons/CDSs do not contain multiple Parent IDs

    Note: If there are additional transcript type in the input gfff, then you should go to intra_model.featureSort, and add the new transcript type to the dict of  FEATURECODE.
    '''

    eCode = 'Emn0001'
    result = dict()

    pairs = list()
    for i in range(len(linelist)-1):
        for j in range(i+1, len(linelist)):
            source, target = linelist[i], linelist[j]
            if source['seqid'] == target['seqid']:
                s7 = '{0:s}\t{1:s}\t{2:s}\t{3:d}\t{4:d}\t{5:s}\t{6:s}'.format(source['seqid'], source['source'], source['type'], source['start'], source['end'], source['score'], source['strand'], source['phase'])
                t7 = '{0:s}\t{1:s}\t{2:s}\t{3:d}\t{4:d}\t{5:s}\t{6:s}'.format(target['seqid'], target['source'], target['type'], target['start'], target['end'], target['score'], target['strand'], target['phase'])
                if s7 == t7:
                    pairs.append({'source':source, 'target':target})

    for pair in pairs:
        same_target = False
        if pair['source'].has_key('children') and pair['target'].has_key('children'):
            schildren = pair['source']['children']
            tchildren = pair['target']['children']
            if len(schildren) == len(tchildren):
                sort_schildren = intra_model.featureSort(schildren, reverse=True if pair['source']['strand'] == '-' else False)
                sort_tchildren = intra_model.featureSort(tchildren, reverse=True if pair['source']['strand'] == '-' else False)
                for i in range(len(sort_schildren)):
                    s7 = '{0:s}\t{1:s}\t{2:s}\t{3:d}\t{4:d}\t{5:s}\t{6:s}'.format(sort_schildren[i]['seqid'], sort_schildren[i]['source'], sort_schildren[i]['type'], sort_schildren[i]['start'], sort_schildren[i]['end'], sort_schildren[i]['score'], sort_schildren[i]['strand'], sort_schildren[i]['phase'])
                    t7 = '{0:s}\t{1:s}\t{2:s}\t{3:d}\t{4:d}\t{5:s}\t{6:s}'.format(sort_tchildren[i]['seqid'], sort_tchildren[i]['source'], sort_tchildren[i]['type'], sort_tchildren[i]['start'], sort_tchildren[i]['end'], sort_tchildren[i]['score'], sort_tchildren[i]['strand'], sort_tchildren[i]['phase'])
                    if s7 == t7:
                        same_target=True
                    else:
                        same_target=False
                        break
        if same_target:
            tmp = [pair['source']['attributes']['ID'], pair['target']['attributes']['ID']]
            sort_tmp = sorted(tmp)
            key = '{0:s},{1:s}'.format(sort_tmp[0], sort_tmp[1])
            if not result.has_key(key):
                result[key] = []
            result[key].append(eCode)
            pair['source']['line_errors'].append(eCode)
            pair['target']['line_errors'].append(eCode)
                   
    if len(result):
        return result


def main(gff):
    FIX_MISSING_ATTR(gff3)

    ERROR_CODE = ['Emn0001']
    ERROR_TAG = ['Duplicate transcripts found']
    ERROR_INFO = dict(zip(ERROR_CODE, ERROR_TAG))

    roots = [line for line in gff.lines if line['line_type']=='feature' and not line['attributes'].has_key('Parent')]
    error_set=dict()
    trans_list = list()
    for root in roots:
        children = root['children']
        for child in children:
            trans_list.append(child)

    r = check_duplicate(gff, trans_list)
    if not r == None:
        error_set = dict(error_set.items() + r.items())

    for k,v in error_set.items():
        for e in v:
            print (k, e, ERROR_INFO[e])



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
    QC functions for processing multiple features between models (inter-model) in GFF3 file.
    
    Testing enviroment:
    1. Python 2.7

    Inputs:
    1. GFF3: reads from STDIN by default, may specify the file name with the -g argument

    Outputs:
    1. GFF3: fixed GFF file

    """))
    parser.add_argument('-g', '--gff', type=str, help='Summary Report from Monica (default: STDIN)') 
    parser.add_argument('-o', '--output', type=str, help='Output file name (default: STDIN)')
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

    if args.output:
        logger_stderr.info('Specifying output file name: (%s)...\n', args.output)
        report_fh = open(args.output, 'wb')
    
    gff3 = Gff3(gff_file=args.gff, logger=logger_null)
    main(gff3)
