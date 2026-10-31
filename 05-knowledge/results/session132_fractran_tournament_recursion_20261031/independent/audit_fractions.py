#!/usr/bin/env python3
"""Independent audit of emitted FRACTRAN files; no imports from their compiler.

Usage: python audit_fractions.py ../fractran/a000568 --through 5
Checks literal ordered fractions against every counter instruction, primality,
disjointness, reduction, guardedness, ordered zero tests, and no self-transition.
Then runs actual fraction exponents, not instruction dispatch. Source-control
indexing removes only provably ineligible rows; retained rows retain file order.
The n=0 run is additionally checked by scanning every actual rational fraction
against a literal integer. Finite replay supplements, not replaces, the proof.
"""
import argparse
from collections import Counter
from hashlib import sha256
from math import gcd, isqrt
from pathlib import Path
import json


def require(test, message):
    if not test:
        raise RuntimeError(message)


def prime(p):
    return p >= 2 and all(p % d for d in range(2,isqrt(p)+1))


def factor(v):
    result = Counter()
    d = 2
    while d*d <= v:
        while v % d == 0:
            result[d] += 1
            v //= d
        d += 1
    if v > 1:
        result[v] += 1
    return result


def load(stem):
    stem = Path(stem)
    paths = [stem.with_suffix('.ir.json')] + sorted(stem.parent.glob(stem.name+'.fractions.*.txt'))
    require(len(paths)>1,'No emitted fraction files')
    ir = json.loads(paths[0].read_text())
    ds, qs = ir['data_primes'],ir['state_primes']
    allp = ds+qs
    require(len(set(allp))==len(allp),'Data/control prime overlap or duplicates')
    require(all(prime(p) for p in allp),'Nonprime register or state label')
    require(len(ds)==ir['register_count'],'Register count mismatch')
    require(len(qs)==len(ir['instructions']),'State count mismatch')
    rows = []
    for path in paths[1:]:
        for line in path.read_text().splitlines():
            a,b = map(int,line.split('/'))
            require(a>0 and b>0 and gcd(a,b)==1,'Nonpositive or unreduced fraction')
            rows.append((a,b,factor(a),factor(b)))
    expected = []
    for pc,ins in enumerate(ir['instructions']):
        op = ins[0]
        if op=='HALT':
            require(len(ins)==1,'Invalid halt')
            continue
        if op=='INC':
            r,to = ins[1:]
            require(to!=pc,'Self-transition cancels guard')
            expected.append((ds[r]*qs[to],qs[pc]))
        elif op=='JUMP':
            to, = ins[1:]
            require(to!=pc,'Self-transition cancels guard')
            expected.append((qs[to],qs[pc]))
        elif op=='DECJZ':
            r,yes,no = ins[1:]
            require(yes!=pc and no!=pc,'Self-transition cancels guard')
            expected.extend(((qs[yes],ds[r]*qs[pc]),(qs[no],qs[pc])))
        else:
            raise RuntimeError('Unknown instruction '+op)
    require([(a,b) for a,b,_,_ in rows]==expected,'Emitted fractions differ from exact ordered translation')
    qset = set(qs)
    byq = {q:[] for q in qs}
    for i,(a,b,num,den) in enumerate(rows):
        qin = [(p,e) for p,e in den.items() if p in qset]
        qout = [(p,e) for p,e in num.items() if p in qset]
        require(len(qin)==len(qout)==1 and qin[0][1]==qout[0][1]==1,'Not one-hot guarded')
        require(qin[0][0]!=qout[0][0],'Control cancellation')
        require(set(num)|set(den)<=set(allp),'Unknown prime')
        # This index is built from actual factorization, not PC metadata.
        byq[qin[0][0]].append(i)
    return ir,rows,byq,{str(p.name):sha256(p.read_bytes()).hexdigest() for p in paths}


def exponent_run(ir,rows,byq,n,maxsteps):
    exps = Counter({ir['data_primes'][ir['input']]:n,ir['state_primes'][0]:1})
    controls = set(ir['state_primes'])
    q = ir['state_primes'][0]
    # Factorizations cached from actual integers; no IR instruction use here.
    encoded=[]
    for _,_,pos,neg in rows:
        nextq = next(p for p in pos if p in controls)
        encoded.append((tuple(neg.items()),tuple(pos.items()),nextq))
    step = 0
    while step < maxsteps:
        for idx in byq[q]:
            neg,pos,nxt = encoded[idx]
            if all(exps[p]>=e for p,e in neg):
                for p,e in neg: exps[p]-=e
                for p,e in pos: exps[p]+=e
                q=nxt
                break
        else:
            require(sum(exps[p] for p in controls)==1 and exps[q]==1,'Lost one-hot invariant')
            return dict(n=n,output=exps[ir['data_primes'][ir['output']]],
                        error=exps[ir['data_primes'][ir['error']]],steps=step,halt_prime=q)
        step+=1
    raise RuntimeError('Exceeded step bound '+str(maxsteps))


def literal_run(ir,rows,n,maxsteps=100000):
    value = ir['state_primes'][0]*ir['data_primes'][ir['input']]**n
    peak = value.bit_length()
    for steps in range(maxsteps+1):
        for a,b,_,_ in rows:
            if value % b == 0:
                value = value//b*a
                peak=max(peak,value.bit_length())
                break
        else:
            p = ir['data_primes'][ir['output']]
            out=0
            while value%p==0:
                out+=1
                value//=p
            return dict(n=n,output=out,steps=steps,max_integer_bits=peak)
    raise RuntimeError('Literal replay exceeded step bound')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('stem',type=Path)
    parser.add_argument('--through',type=int,default=3)
    parser.add_argument('--maxsteps',type=int,default=100000000)
    args=parser.parse_args()
    ir,rows,byq,hashes=load(args.stem)
    runs=[]
    for n in range(args.through+1):
        row=exponent_run(ir,rows,byq,n,args.maxsteps)
        require(row['error']==0,'Reachable compiler ERROR')
        runs.append(row)
        print(json.dumps(row),flush=True)
    lit=literal_run(ir,rows,0)
    require(lit['output']==runs[0]['output'] and lit['steps']==runs[0]['steps'],'Literal/exponent mismatch')
    report=dict(status='PASS',fraction_count=len(rows),registers=len(ir['data_primes']),states=len(ir['state_primes']),
                source_hashes=hashes,exponent_runs=runs,literal_run=lit)
    print(json.dumps(report,indent=2,sort_keys=True),flush=True)


if __name__=='__main__':
    main()
