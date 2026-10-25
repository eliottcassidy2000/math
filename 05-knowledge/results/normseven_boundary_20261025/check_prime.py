#!/usr/bin/env python3
"""Prime and composite carry conductors, independent from the binary checker."""
from collections import Counter
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(why)


def line(v, p):
    return frozenset(tuple(a*x % p for x in v) for a in range(p))


def check(parts, p, n):
    vectors = list(product(range(p), repeat=n))
    e = {}
    for selected in product(*parts):
        sums = tuple(sum(v[i] for v in selected) for i in range(n))
        z = tuple(s % p for s in sums)
        need(z not in e, 'direct sum, no repeated character')
        e[z] = tuple(s//p for s in sums)
    need(len(e) == p**n, 'full character coverage')
    active = [sum(any(v[i] for v in part) for part in parts) for i in range(n)]
    a = [sum(max(v[i] for v in part) for part in parts)-p+1 for i in range(n)]
    need(min(a) >= 0, 'nonnegative conductor exponent')
    if p in (2,3,5,7):
        need(a == [(t-1)*(p-1) for t in active], 'prime simplification')
    for z in vectors:
        calculated = [0]*n
        for y in vectors:
            w = tuple((x+t) % p for x,t in zip(z,y))
            for i in range(n):
                calculated[i] = max(calculated[i], e[w][i]-(z[i]+y[i])//p)
        predicted = [(a[i]+p-1-z[i])//p for i in range(n)]
        need(calculated == predicted, 'carry conductor formula')
        # Independent principal-generator multiplication coefficient.
        principal = [(a[i]+((z[i]-a[i]) % p))//p for i in range(n)]
        need(principal == calculated, 'principal conductor consequence')
    return tuple(a)


def main():
    census = []
    for p in (2,3,5,7):
        n=2
        lines = {line(v,p) for v in product(range(p), repeat=n) if any(v)}
        count=0
        exponents=Counter()
        for parts in combinations(lines,2):
            if len(parts[0] & parts[1]) != 1:
                continue
            exponents[check(parts,p,n)] += 1
            count += 1
        census.append({'prime':p,'dimension':n,'packets':2,'decompositions':count,
                       'exponent_vectors': {str(k):v for k,v in sorted(exponents.items())}})
    # All bases of projective lines in F_3^3, modulo line order and line scaling.
    p,n=3,3
    lines={line(v,p) for v in product(range(p),repeat=n) if any(v)}
    count=0
    exponents=Counter()
    for parts in combinations(lines,3):
        sums={tuple(sum(v[i] for v in chosen)%p for i in range(n))
              for chosen in product(*parts)}
        if len(sums) != p**n:
            continue
        exponents[check(parts,p,n)] += 1
        count += 1
    census.append({'prime':p,'dimension':n,'packets':3,'decompositions':count,
                   'exponent_vectors':{str(k):v for k,v in sorted(exponents.items())}})
    # A block of dimension two, not just one-dimensional factors.
    c = {tuple((a*x+b*y)%5 for x,y in zip((1,0,1),(0,1,1)))
         for a,b in product(range(5),repeat=2)}
    d=line((1,1,1),5)
    block = check((c,d),5,3)
    need(block == (4,4,4), 'two-plus-one dimensional packet')
    # Composite-modulus hostile: nonzero coordinate maps need not be onto.
    # Z/4 lines (1,2),(0,1) are complementary, but the second coordinate
    # in the first line has image {0,2}, not all four residues.
    composite_parts=(line((1,2),4),line((0,1),4))
    e={tuple(sum(v[i] for v in selected)%4 for i in range(2)):
       tuple(sum(v[i] for v in selected)//4 for i in range(2))
       for selected in product(*composite_parts)}
    composite={}
    for z in product(range(4),repeat=2):
        calculated=[0,0]
        for y in product(range(4),repeat=2):
            w=tuple((z[i]+y[i])%4 for i in range(2))
            for i in range(2):
                calculated[i]=max(calculated[i],e[w][i]-(z[i]+y[i])//4)
        composite[str(z)]=calculated
    need(composite['(0, 2)'][1] == 0, 'non-surjective composite hostile')
    need(check(composite_parts,4,2) == (0,2), 'repaired composite conductor')
    composite_counts = {}
    for m in (4,6,8):
        lines={line(v,m) for v in product(range(m),repeat=2)
               if len(line(v,m)) == m}
        count=0
        for parts in combinations(lines,2):
            if len(parts[0] & parts[1]) != 1:
                continue
            check(parts,m,2)
            count += 1
        composite_counts[m]=count
    # Nonfree cyclic-module summands: Z/6 = C2 direct_sum C3.
    need(check(({(0,),(3,)}, {(0,),(2,),(4,)}),6,1) == (2,),
         'cusp semigroup <2,3> has conductor 2')
    from math import gcd
    gap_checks=[]
    for a,b in combinations(range(2,13),2):
        if gcd(a,b) != 1:
            continue
        conductor=(a-1)*(b-1)
        elements={a*x+b*y for x,y in product(range(a*b+1),repeat=2)
                  if a*x+b*y <= a*b}
        gaps=[v for v in range(a*b) if v not in elements]
        need(max(gaps)+1 == conductor, 'coprime numerical semigroup conductor')
        need(len(gaps)*2 == conductor, 'exact gap count')
        gap_checks.append([a,b,conductor,len(gaps)])
    print(json.dumps({'status':'FINITE-EXACT','gates':GATES,'census':census,
                     'dimension_two_plus_one_p5_conductor_exponents':block,
                     'composite_dimension_two_decompositions':composite_counts,
                     'coprime_semigroup_checks':gap_checks,
                     'composite_hostile':{'modulus':4,'generators':[[1,2],[0,1]],
                                          'tested_character':[0,2],
                                          'actual_b2_exponent':0,'naive_prime_formula':1},
                     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
                    indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
