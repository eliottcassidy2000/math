"""Audit actual integer FRACTRAN runs against independent matrix multiplication."""
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import random
import sys
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('linear_stream',ROOT/'distribution/stream.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def raw_check(a,v,rounds,clock_offset=0):
    p,_=module.compile_linear(a)
    d=len(a);Q=p['checkpoint'];clock=p['clock']
    N=Q*clock**clock_offset*math.prod(prime**e for prime,e in zip(p['x'],v))
    expected=v[:];elapsed=0;counts=[];bitmax=N.bit_length()
    for turn in range(rounds):
        old=expected[:]
        expected=[sum(coef*e for coef,e in zip(row,old)) for row in a]
        previous=elapsed
        while True:
            for numerator,denominator in p['fractions']:
                need(math.gcd(numerator,denominator)==1,'cancelled guard')
                if N%denominator==0:
                    N=N//denominator*numerator;elapsed+=1
                    break
            else:raise RuntimeError('unexpected halt')
            bitmax=max(bitmax,N.bit_length())
            need(elapsed<2_000_000,'literal step limit')
            if N%Q==0:break
        expected_integer=Q*clock**(turn+1+clock_offset)*math.prod(prime**e for prime,e in zip(p['x'],expected))
        need(N==expected_integer,('checkpoint mismatch',a,v,turn))
        count=elapsed-previous
        need(count==1+2*d+2*sum(old)+2*sum(expected),'round cost identity')
        counts.append(count)
    return {'matrix':a,'initial':v,'rounds':rounds,'last':expected,
            'steps_by_round':counts,'max_integer_bits':bitmax}

def main():
    cases=[]
    examples=[([[0]],[0],3),([[0]],[3],3),([[1]],[2],4),([[2]],[1],7),
              ([[1,1],[1,0]],[1,0],8),([[0,1],[1,0]],[2,3],6),
              ([[0,0,2],[1,0,0],[0,3,0]],[1,2,3],3),
              ([[0,0],[0,0]],[2,3],4)]
    rng=random.Random(133)
    for _ in range(24):
        d=rng.randint(1,4)
        examples.append(([[rng.randint(0,2) for _ in range(d)] for _ in range(d)],
                         [rng.randint(0,3) for _ in range(d)],3))
    for a,v,r in examples:cases.append(raw_check(a,v,r))
    rejected=0
    for a in [[],[[1,2]],[[1],[2]],[[-1]],[[1.0]],[[True]]]:
        try:module.compile_linear(a)
        except ValueError:rejected+=1
        else:raise RuntimeError(('invalid matrix accepted',a))
    atlas=json.loads((ROOT/'spectrum/atlas.json').read_text())
    atoms=[]
    from collections import Counter
    c=Counter((r['n'],r['H']) for r in atlas['classes'] if r['strong'] and r['H']<=21)
    atoms=sorted((n,h,k) for (n,h),k in c.items())
    need(atoms==sorted(module.ATOMS),'atlas/stream atom mismatch')
    # Inspect every matrix entry against the independently specified recurrence.
    H=module.H;L=module.L;a=module.matrix();d=len(H)*L
    target=[[0]*d for _ in range(d)]
    for hi,h in enumerate(H):
        target[hi*L][hi*L]+=1
        for n,q,k in atoms:
            if h%q==0:target[hi*L][H.index(h//q)*L+n-1]+=k
        for lag in range(1,L):target[hi*L+lag][hi*L+lag-1]=1
    need(a==target,'recurrence/matrix mismatch')
    def closed(n):
        return [1,n-2,n-3,0,(n*n-5*n+4)//2,n-4,n-4,
                (n-3)*(n-4),2*n-10,2*n-10,0]
    tail=module.tail_matrix()
    vector=closed(6)
    need(vector==[module.ordinary_table(6)[6][h] for h in H],'tail start')
    for n in range(7,1001):
        vector=[sum(c*v for c,v in zip(row,vector)) for row in tail]
        need(vector==closed(n),('tail closed form',n))
    tail_literal=raw_check(tail,closed(6),6,clock_offset=6)
    tail_program=json.loads((ROOT/'distribution/tail_program.json').read_text())
    p,_=module.compile_linear(tail)
    seed=p['checkpoint']*p['clock']**6*math.prod(q**e for q,e in zip(p['x'],closed(6)))
    need(int(tail_program['initial_integer'])==seed,'serialized tail seed')
    actual=[tuple(map(int,line.split('/'))) for line in (ROOT/'distribution/tail21.fractran').read_text().splitlines()]
    need(actual==p['fractions'],'serialized tail fractions')
    compact=module.compact_tail_matrix();compact_h=module.COMPACT_H
    compact_v=[closed(6)[H.index(h)] for h in compact_h]
    compact_literal=raw_check(compact,compact_v,6,clock_offset=6)
    cp,_=module.compile_linear(compact)
    compact_artifact=json.loads((ROOT/'distribution/compact_tail_program.json').read_text())
    zero=compact_artifact['unused_zero_prime']
    compact_seed=cp['checkpoint']*cp['clock']**6*math.prod(q**e for q,e in zip(cp['x'],compact_v))
    need(int(compact_artifact['initial_integer'])==compact_seed,'compact seed')
    compact_rows=[tuple(map(int,line.split('/'))) for line in (ROOT/'distribution/compact_tail21.fractran').read_text().splitlines()]
    need(compact_rows==cp['fractions'] and len(compact_rows)==43,'compact emitted list')
    need(compact_seed%zero!=0 and all(num%zero and den%zero for num,den in compact_rows),'unused zero prime not invariant')
    cv=compact_v
    for n in range(6,1001):
        values=dict(zip(cp['x'],cv));values[zero]=0
        actual=[values[compact_artifact['output_primes'][str(h)]] for h in H]
        need(actual==closed(n),('compact aliases',n))
        cv=[sum(c*v for c,v in zip(row,cv)) for row in compact]
    print(json.dumps({'status':'PASS','literal_matrix_cases':cases,
          'invalid_matrices_rejected':rejected,'strong_atoms_verified':atoms,
          'recurrence_matrix_dimension':d,
          'tail_closed_form_checks':'orders6..1000','tail_literal':tail_literal,
          'compact_alias_checks':'orders6..1000','compact_literal':compact_literal,
          'compact_unused_zero_prime':zero,
          'scope':'All-order recurrence correctness is conditional on completeness of the bounded strong inventory; literal generic compiler tests are finite.',
          'source_sha256':{str(q.relative_to(ROOT)):hashlib.sha256(q.read_bytes()).hexdigest()
                           for q in [ROOT/'distribution/stream.py',ROOT/'spectrum/atlas.json']}},indent=2))

if __name__=='__main__':main()
