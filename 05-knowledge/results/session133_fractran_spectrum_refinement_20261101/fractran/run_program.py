#!/usr/bin/env python3
"""Read fixed fractions from disk, then perform literal first-eligible FRACTRAN.

No source evaluation, compiler call, formula, or sequence table is used here.
"""
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def require(test, message):
    if not test:
        raise RuntimeError(message)


def load(name):
    meta = json.loads((ROOT / f'{name}.meta.json').read_text())
    rows = []
    for part in meta['fraction_parts']:
        for line in (ROOT / part).read_text().splitlines():
            a,b = map(int,line.split('/'))
            require(a>0 and b>0,'positive fraction required')
            rows.append((a,b))
    states = []
    for part in meta['state_prime_parts']:
        states += json.loads((ROOT / part).read_text())
    return meta,rows,states


def run(name,n,limit=1_000_000,bit_limit=1_000_000):
    require(type(n) is int and n>=0,'nonnegative input required')
    meta,rows,states = load(name)
    primes = meta['data_primes']
    value = states[0] * primes[meta['input']]**n
    steps = 0
    maxbits = value.bit_length()
    while steps < limit:
        for a,b in rows:
            if value%b==0:
                value=value//b*a
                break
        else:
            result = {'steps':steps,'max_integer_bits':maxbits}
            for field in ('output','error'):
                p=primes[meta[field]]
                count=0
                while value%p==0:
                    count+=1
                    value//=p
                result[field]=count
            require(value%states[-1]==0,'terminal is not HALT')
            return result
        steps+=1
        maxbits=max(maxbits,value.bit_length())
        require(maxbits<=bit_limit,'integer bit limit reached')
    raise RuntimeError('literal FRACTRAN step limit reached')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('program')
    parser.add_argument('n',type=int)
    parser.add_argument('--steps',type=int,default=1_000_000)
    args=parser.parse_args()
    print(json.dumps(run(args.program,args.n,args.steps),sort_keys=True))
