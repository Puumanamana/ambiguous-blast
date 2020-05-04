# This scripts is a wrapper around Biopython's BLAST commandline and parsing
# It implements a simple approch to blast primers with ambiguous nucleotides
# by making all possible sequences from an ambiguous one, blasting everything
# and retrieving any sequence matching any if those possibilities.

import sys
from tempfile import mkdtemp
from pathlib import Path
import argparse
from itertools import product

import pandas as pd
import numpy as np

from Bio import SeqIO, SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import (NcbiblastnCommandline,
                                    NcbiblastpCommandline,
                                    NcbiblastxCommandline,
                                    NcbitblastnCommandline)

NUCLEOTIDES = ['A', 'C', 'G', 'T']

AMBIGUOUS_CODES = {
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['C', 'G'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'T'],
    'B': ['C', 'G', 'T'],
    'V': ['A', 'C', 'G'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'N': ['A', 'C', 'G', 'T'],        
}

TABULAR_BLAST_FIELDS = {'qseqid': str, 'sseqid':str, 'staxid': int,
                        'qstart': int, 'qend': int, 'sstart': int, 'send': int,
                        'evalue': float, 'pident': float, 'qcovhsp': float}

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--cmd', type=str, choices=['blastn', 'blastp', 'blastx', 'tblastn'])
    parser.add_argument('--query', type=str)    
    parser.add_argument('--db', type=str, default='refseq')
    parser.add_argument('--primers', action='store_true', default=False)
    parser.add_argument('--min-id', type=float, default=.8)
    parser.add_argument('--min-cov', type=float, default=.5)
    parser.add_argument('--num-threads', type=int, default=30)
    args = parser.parse_args()

    args.query = check_input(args.query)
    args.db = check_input(args.db)

    args.task = '{}-short'.format(args.cmd) if args.primers else None

    return args

def check_input(filename):
    if Path(filename).parent.is_dir():
        return Path(filename)
    raise FileNotFoundError

def blast(cmd, query, db, **kwargs):

    outfmt = "'6 {}'".format(' '.join(TABULAR_BLAST_FIELDS.keys()))
    ext = '.tsv'

    out_file = Path(query.parent, 'result{}'.format(ext))

    if cmd == 'blastn':
        blast_on_db = NcbiblastnCommandline(query=str(query), db=str(db),
                                            out=str(out_file), outfmt=outfmt, **kwargs)
    elif cmd == 'blastp':
        blast_on_db = NcbiblastpCommandline(query=str(query), db=str(db),
                                            out=str(out_file), outfmt=outfmt, **kwargs)
    elif cmd == 'blastx':
        blast_on_db = NcbiblastxCommandline(query=str(query), db=str(db),
                                            out=str(out_file), outfmt=outfmt, **kwargs)
    elif cmd == 'tblastn':
        blast_on_db = NcbitblastnCommandline(query=str(query), db=str(db),
                                             out=str(out_file), outfmt=outfmt, **kwargs)
        
    else:
        sys.exit(f'Unknown command: {cmd}')
    print(blast_on_db)
    blast_on_db()

    return out_file

def convert_input(path):
    '''
    Parse input fasta and convert all sequences to multiple unambiguous sequences
    by computing all potential combinations
    '''

    tmpdir = mkdtemp()

    records = []
    for record in SeqIO.parse(path, 'fasta'):
        record_no_ambig = convert_ambiguous_sequence(record)
        records += record_no_ambig

    output = Path(tmpdir, '{}_unambiguous.fasta'.format(path.stem))
    SeqIO.write(records, output, 'fasta')

    return output

def convert_ambiguous_sequence(record):
    ambig_pos = [(i, AMBIGUOUS_CODES[nucl]) for (i, nucl) in enumerate(record.seq)
                 if nucl not in NUCLEOTIDES]
    (positions, ambig_nucl) = list(zip(*ambig_pos))
    positions = list(positions)

    sequence = np.array(list(record.seq))
    
    all_combinations = []

    for combination in product(*ambig_nucl):        
        sequence[positions] = combination
        record_no_ambig = SeqRecord(
            id="{}_combination={}".format(record.id, "+".join(combination)),
            seq=Seq(''.join(list(sequence))),
            description=''
        )
        all_combinations.append(record_no_ambig)

    return all_combinations

def parse_blast_tabular_output(blast_output, min_cov=-1, min_id=-1):

    results = []

    parser = SearchIO.parse(blast_output, 'blast-tab', fields=list(TABULAR_BLAST_FIELDS.keys()))

    def aln_filter(hsp):
        return (hsp.query_coverage>=100*min_cov) and (hsp.ident_pct>=100*min_id)

    for query_hits in map(lambda x: x.hsp_filter(aln_filter), parser):
        all_hits = [
            [hsp.query_id, hsp.hit_id, hsp.ident_pct, hsp.query_coverage]
            for hsp in query_hits.hsps
        ]

        results += all_hits

    results = pd.DataFrame(results, columns=['query', 'subject', 'identity', 'coverage'])

    return results

def get_matches(results, taxids, pe=False):

    if pe:
        results['primer'] = results['query'].str.split('_').str[0]
        both_hits = results.groupby('subject').primer.agg(len)
        both_hits = both_hits[both_hits > 1]
        results = results.loc[both_hits.index]

    results['taxids'] = taxids.loc[results.subject]

    return results

def main():
    '''
    '''

    args = parse_args()

    if args.cmd in {'blastn', 'blastx'}:
        unambiguous_query = convert_input(args.query)
    else:
        unambiguous_query = Path(args.query)

    blast_kwargs = {'num_threads': args.num_threads}
    if args.task is not None:
        blast_kwargs['task'] = args.task

    print('Unambiguous query stored in {}'.format(unambiguous_query))
    blast_file = blast(args.cmd, unambiguous_query, args.db, **blast_kwargs)
    print('Blast result stored in {}'.format(blast_file))
    parsed_results = parse_blast_tabular_output(blast_file, min_cov=args.min_cov, min_id=args.min_id)

    taxids = pd.read_csv(blast_file, sep='\t', index_col='subject', names=['subject', 'taxids'],
                                usecols=[1, 2]).taxids.drop_duplicates()

    matches = get_matches(parsed_results, taxids)
    
    import ipdb;ipdb.set_trace()

    return matches

if __name__ == '__main__':
    main()
