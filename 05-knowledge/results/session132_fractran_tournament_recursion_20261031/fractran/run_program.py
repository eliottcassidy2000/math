#!/usr/bin/env python3
"""Run the committed ordered fractions directly, without invoking the compiler.

Examples:
  python3 run_program.py a000568 3
  python3 run_program.py literal16 6
This intentionally uses literal integer divisibility and full ordered scans.
"""
import argparse
import json
from pathlib import Path
import literal16

ROOT=Path(__file__).resolve().parent


def run(rows,start,output_prime,step_limit):
    value=start; steps=0; peak=value.bit_length()
    while steps<step_limit:
        for numerator,denominator in rows:
            if value%denominator==0:
                value=value//denominator*numerator
                break
        else:
            output=0
            while value%output_prime==0:
                output+=1; value//=output_prime
            return {'output':output,'steps':steps,'max_integer_bits':peak}
        steps+=1;peak=max(peak,value.bit_length())
    raise RuntimeError(f'nonterminal after {step_limit} steps (not a nontermination claim)')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('program',choices=['a000568','literal16'])
    parser.add_argument('n',type=int)
    parser.add_argument('--max-steps',type=int,default=1_000_000)
    args=parser.parse_args()
    if args.program=='literal16':
        literal16.expected(args.n)
        rows=literal16.FRACTIONS
        start=306*13**(args.n-3); p=2
    else:
        if args.n<0: raise ValueError('n must be nonnegative')
        rows=[]
        for path in sorted(ROOT.glob('a000568.fractions.*.txt')):
            rows.extend(tuple(map(int,line.strip().split('/'))) for line in path.read_text().splitlines())
        if not rows: raise RuntimeError('ordered fraction file missing')
        start=157*5**args.n; p=3
    print(json.dumps({'program':args.program,'n':args.n,**run(rows,start,p,args.max_steps)},indent=2))

if __name__=='__main__': main()
