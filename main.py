import sys
from tempfile import mkdtemp
from pathlib import Path
from itertools import product
import argparse

import pandas as pd
from Bio import SeqIO
from Bio.Blast.Applications import (NcbiblastnCommandline,
                                    NcbiblastpCommandline,
                                    NcbiblastxCommandline,
                                    NcbitblastnCommandline)
from Bio.Blast import NCBIXML

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

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--cmd', type=str, choices=['blastn', 'blastp', 'blastx', 'tblastn'])
    parser.add_argument('--query', type=str)    
    parser.add_argument('--db', type=str, default='refseq')
    parser.add_argument('--min-id', type=float, default=0)
    parser.add_argument('--min-cov', type=float, default=0.5)    
    args = parser.parse_args()

    args.query = check_input(args.query)
    args.db = check_input(args.db)

    return args

def check_input(filename):
    if Path(filename).is_file():
        return Path(filename)
    if Path('databases', filename).is_file():
        return Path('databases', filename)
    raise FileNotFoundError

def blast(cmd, query, db):

    tmpdir = mkdtemp()
    out_file = Path(tmpdir, 'result.xml')

    if cmd == 'blastn':
        blast_on_db = NcbiblastnCommandline(query=query, db=db, out=str(out_file), outfmt=5)
    elif cmd == 'blastp':
        blast_on_db = NcbiblastpCommandline(query=query, db=db, out=str(out_file), outfmt=5)
    elif cmd == 'blastx':
        blast_on_db = NcbiblastxCommandline(query=query, db=db, out=str(out_file), outfmt=5)
    elif cmd == 'tblastn':
        blast_on_db = NcbitblastnCommandline(query=query, db=db, out=str(out_file), outfmt=5)
        
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
    with SeqIO.open(path, 'r') as handle:
        for record in handle:
            record_no_ambig = convert_ambiguous_sequence(record)
            records += record_no_ambig

    output = Path(tmpdir, '{}_unambiguous.fasta'.format(path.stem))
    SeqIO.write(record_no_ambig, output, 'fasta')

    return output

def convert_ambiguous_sequence(record):
    ambig_pos = [(i, AMBIGUOUS_CODES[nucl]) for (i, nucl) in enumerate(record.seq)
                 if nucl not in NUCLEOTIDES]
    (positions, ambig_nucl) = list(zip(*ambig_pos))

    all_combinations = []
    for combination in product(ambig_nucl):
        record_no_ambig = record.copy()

        for (position, nucl) in zip(positions, combination):
            record_no_ambig.seq[position] = nucl

        record_no_ambig.id = "{}_combination={}".format(record_no_ambig.id, "+".join(combination))
        all_combinations.append(record_no_ambig)

    return all_combinations

def parse_blast_output(blast_output, min_cov=-1, min_id=-1):

    results = {}
    
    for record in NCBIXML.parse(open(blast_output)):

        if record.alignments:
            record_id = record.query.split(' ')[0]
            results[record_id] = []

            for hit in record.alignments:
                # 1 alignment is one hit (with potentially multiple HSPs)

                matches = sum(hsp.identities for hsp in hit.hsps)
                covered = sum(hsp.align_length for hsp in hit.hsps)

                pct_identity = matches / covered
                pct_coverage = covered / record.query_length

                if pct_identity > min_id and pct_coverage > min_cov:
                    hit_id = hit.hit_def.split(' ')[0]

                    entry = [[hit_id, hsp.query_start, hsp.query_end, hsp.identities]
                             for hsp in hit.hsps]
                    results[record_id] += entry

    results = pd.DataFrame(results, columns=['source', 'start', 'end', 'identity'])

    return results

def main():
    '''
    '''

    args = parse_args()

    unambiguous_query = convert_input(args.query)
    print('Unambiguous query stored in {}'.format(unambiguous_query))
    blast_file = blast(args.cmd, unambiguous_query, args.db)
    print('Blast result stored in {}'.format(blast_file))
    parsed_results = parse_blast_output(blast_file, min_cov=args.min_cov, min_id=args.min_id)
    print(parsed_results.head(50))

if __name__ == '__main__':
    main()
