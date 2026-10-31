#!/usr/bin/env python3
"""Independent exact audit; deliberately imports no producer code.

Construct the 884 labelled balanced inputs on thirteen primes, quotient only the
three spectator labels, and evaluate ordered rational multiplication directly.
All validation is always active, including under -O.
"""
import csv
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parent
P=(2,3,5,7,11,13,17,19,23,29)
B=P+(31,37,41)
Q=tuple(Fraction(x) for x in ('17/91','78/85','19/51','23/38','29/33',
    '77/29','95/23','77/19','1/17','11/13','13/11','15/2','1/7','55/1'))
COUNT=Counter()


def check(condition, message, group):
    COUNT[group]+=1
    if not condition:
        raise RuntimeError(message)


def exponents(n):
    result=[]
    for p in B:
        e=0
        while n%p==0:
            n//=p
            e+=1
        result.append(e)
    check(n==1,'unexpected factor outside basis','factorization')
    return tuple(result)


def orbit_key(e):
    return e[:10],tuple(sorted((k for k in e[10:] if k),reverse=True))


def value(e, primes=B):
    return prod(p**k for p,k in zip(primes,e))


def multiply_step(n):
    for i,q in enumerate(Q,1):
        candidate=n*q
        if candidate.denominator==1:
            return candidate.numerator,i
    raise RuntimeError('unconditional last fraction missing')


def counts(e):
    if not any(e):
        return None
    # Independent divisor enumeration by bounded products, not tau/omega formula.
    divs=[1]
    for p,k in zip(B,e):
        divs=[d*p**j for d in divs for j in range(k+1)]
    n=value(e)
    divs=[d for d in divs if d not in (1,n)]
    return (len(divs),sum(d in B for d in divs),
            sum(all(d%(p*p)!=0 for p in B) for d in divs))


def is_balanced(e):
    f=counts(e)
    return f is not None and f[0]==sum(f[1:])


def encoded_counts(e):
    result=counts(e)
    return 'undefined' if result is None else ','.join(map(str,result))


def main():
    labelled=set(B)|{p**3 for p in B}
    for p,q,r in combinations(B,3):
        labelled.update((p*p*q*r,p*q*q*r,p*q*r*r))
    wanted={orbit_key(exponents(n)) for n in labelled}
    raw=(ROOT/'transitions.tsv').read_bytes()
    rows=list(csv.DictReader(raw.decode().splitlines(),delimiter='\t'))
    check(len(labelled)==884 and len(wanted)==538,'universe counts','census')
    check(len(rows)==538,'table row count','census')
    observed=set()
    edges={}
    histograms={name:Counter() for name in ('active','all','active_valid','all_valid')}
    orbit_histograms={name:Counter() for name in histograms}
    for row in rows:
        n=int(row['n']);m=int(row['m']);i=int(row['rule'])
        e=exponents(n);f=exponents(m)
        key=orbit_key(e)
        check(key not in observed,'duplicate orbit','census')
        observed.add(key)
        check(is_balanced(e),'input not balanced','row')
        check(multiply_step(n)==(m,i),'integer/rational replay disagreement','row')
        check(row['input_fus']==encoded_counts(e) and row['output_fus']==encoded_counts(f),
              'direct divisor count mismatch','row')
        outcome='unit' if m==1 else 'balanced' if is_balanced(f) else 'unbalanced'
        check(outcome==row['outcome'],'outcome mismatch','row')
        check(row['active_only']==str(int(not any(e[10:]))),'active flag mismatch','row')
        valid=sum(e[4:10])<=1
        check(row['valid_control']==str(int(valid)),'valid flag mismatch','row')
        check(not valid or sum(f[4:10])<=1,'valid carrier not invariant','row')
        # An independent spectator reassignment keeps the fixed program unchanged.
        new_basis=P+(101,103,107)
        check(multiply_step(value(e,new_basis))==(value(f,new_basis),i),
              'spectator reassignment mismatch','row')
        edges[key]=orbit_key(f) if outcome=='balanced' else None
        universes=['all']
        if not any(e[10:]): universes.append('active')
        if valid:
            universes.append('all_valid')
            if not any(e[10:]): universes.append('active_valid')
        for name in universes:
            histograms[name][outcome]+=1
        row['_key']=key
        row['_universes']=universes
    check(observed==wanted,'missing or extra spectator orbit','census')
    cycle_sets=set()
    for row in rows:
        key=row['_key'];path=[];cursor=key
        while cursor is not None and cursor not in path:
            path.append(cursor)
            check(cursor in edges,'balanced output missing from complete table','graph')
            cursor=edges[cursor]
        if cursor is None:
            result=f'exit_at_{len(path)}'
            check(len(path)<=3,'exit bound violated','graph')
            if row['valid_control']=='1':
                check(len(path)<=2,'valid exit bound violated','graph')
        else:
            cycle=path[path.index(cursor):]
            check(len(cycle)==2 and len(path)==2,'nontrivial tail or non-2-cycle','graph')
            cycle_sets.add(tuple(sorted(cycle)))
            result='cycle_2'
        check(row['restricted_orbit']==result,'restricted orbit receipt mismatch','graph')
        for name in row['_universes']:
            orbit_histograms[name][result]+=1
    check(len(cycle_sets)==8,'spectator orbit count of cycles','graph')
    # Separate Fraction replay of the entire finite seed-2 balance window.
    # Tail exclusion is a symbolic consequence of the macro token bound.
    window=json.loads((ROOT/'seed_balance_window.json').read_text())
    n=2
    hits=[]
    for t in range(146):
        e=exponents(n)
        if is_balanced(e):
            hits.append((t,n))
        if t==145:
            check(n==5**5*7**4*13,'seed boundary C(5,4)','seed_window')
        else:
            n=multiply_step(n)[0]
    check(hits==[(r['step'],r['n']) for r in window['balanced_visits']],
          'complete seed balance window mismatch','seed_window')
    check(len(hits)==16 and hits[-1]==(70,60),'seed visit summary','seed_window')
    result={'status':'PASS','labelled_inputs':len(labelled),'spectator_orbits':len(wanted),
            'outcomes':{k:dict(v) for k,v in histograms.items()},
            'restricted_orbits':{k:dict(v) for k,v in orbit_histograms.items()},
            'balanced_two_cycle_types':len(cycle_sets),
            'seed_balance_visits':len(hits),'last_seed_balance_visit':hits[-1],
            'checks':dict(sorted(COUNT.items())),'checks_total':sum(COUNT.values()),
            'table_sha256':sha256(raw).hexdigest()}
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__=='__main__':
    main()
