#!/usr/bin/env python3
"""Exact mixed-cluster Smith audit; all universal claims have a paper proof.

Universe: p=3,5,7,11; e=1,2,3; duplicate depth d=1,2,3,4;
three residue-lift systems (including nonzero duplicated residues, negative
lifts, unit multipliers and translations). Rational DVR and independent
modular DVR paths compute complete Smith lists. Four rank-eight cases also
enumerate EVERY minor by fraction-free determinants. Confluent residual
rank checks retain the divided duplicate row and all degree columns.
No project imports, floating point, or optimization-disabled assertions.
"""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb
import json
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def require(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


def vp(value, p):
    value = Fraction(value)
    require(value != 0, "valuation of nonzero input")
    a, b, v = abs(value.numerator), value.denominator, 0
    while a % p == 0:
        a //= p
        v += 1
    while b % p == 0:
        b //= p
        v -= 1
    return v


def matrix(nodes):
    return [[comb(q, r)*x**(q-r) if q >= r else 0
             for q in range(2*len(nodes))] for x in nodes for r in (0, 1)]


def rational_smith(entries, p):
    a = [[Fraction(v) for v in row] for row in entries]
    n = len(a)
    answer = []
    for k in range(n):
        i, j = min(((i, j) for i in range(k, n) for j in range(k, n) if a[i][j]),
                   key=lambda pair: vp(a[pair[0]][pair[1]], p))
        a[k], a[i] = a[i], a[k]
        for row in a:
            row[k], row[j] = row[j], row[k]
        pivot = a[k][k]
        answer.append(vp(pivot, p))
        for i in range(k+1, n):
            factor = a[i][k]/pivot
            require(not factor or vp(factor,p) >= 0, "DVR row operation is integral")
            a[i] = [left-factor*right for left,right in zip(a[i],a[k])]
            require(a[i][k] == 0, "DVR pivot column cleared")
    require(answer == sorted(answer), "DVR output ordered")
    return tuple(answer)


def modular_smith(nodes, p):
    total = 4*sum(vp(y-x,p) for i,x in enumerate(nodes) for y in nodes[i+1:])
    precision = total+1
    modulus = p**precision
    n = 2*len(nodes)
    a = [[pow(x,q,modulus) if r == 0 else
          (q*pow(x,q-1,modulus) % modulus if q else 0)
          for q in range(n)] for x in nodes for r in (0,1)]

    def val(x):
        if not x:
            return precision
        v = 0
        while x % p == 0:
            x //= p
            v += 1
        return v

    result = []
    for k in range(n):
        i,j = min(((i,j) for i in range(k,n) for j in range(k,n)),
                  key=lambda ij:val(a[ij[0]][ij[1]]))
        a[k],a[i] = a[i],a[k]
        for row in a:
            row[k],row[j] = row[j],row[k]
        v = val(a[k][k])
        require(v < precision, "independent precision exceeds every invariant")
        result.append(v)
        power = p**v
        reduced = modulus//power
        inverse = pow(a[k][k]//power,-1,reduced)
        for i in range(k+1,n):
            require(a[i][k] % power == 0,"modular divisibility")
            factor = (a[i][k]//power)*inverse % reduced
            a[i] = [(left-factor*right) % modulus for left,right in zip(a[i],a[k])]
            require(a[i][k] == 0,"modular pivot cleared")
    require(sum(result) == total,"independent determinant mass")
    return tuple(result)


def bareiss(entries):
    a = [list(row) for row in entries]
    n = len(a)
    if not n:
        return 1
    sign,old = 1,1
    for k in range(n-1):
        pivot = next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:
            return 0
        if pivot != k:
            a[k],a[pivot] = a[pivot],a[k]
            sign *= -1
        new = a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator = new*a[i][j]-a[i][k]*a[k][j]
                require(numerator % old == 0,"fraction-free division")
                a[i][j] = numerator//old
            a[i][k] = 0
        old = new
    return sign*a[-1][-1]


def all_minor_valuations(nodes,p):
    a = matrix(nodes)
    n = len(a)
    answer = [0]
    for h in range(1,n+1):
        least = None
        for rows in combinations(range(n),h):
            for columns in combinations(range(n),h):
                determinant = bareiss([[a[i][j] for j in columns] for i in rows])
                if determinant:
                    value = vp(determinant,p)
                    least = value if least is None else min(least,value)
        require(least is not None,"nonzero minor at every rank")
        answer.append(least)
    return tuple(answer)


def rank_mod(entries,p):
    a = [[x % p for x in row] for row in entries]
    rows,columns,rank = len(a),len(a[0]),0
    for j in range(columns):
        pivot = next((i for i in range(rank,rows) if a[i][j]),None)
        if pivot is None:
            continue
        a[rank],a[pivot] = a[pivot],a[rank]
        inv = pow(a[rank][j],-1,p)
        a[rank] = [inv*x % p for x in a[rank]]
        for i in range(rows):
            if i != rank:
                factor = a[i][j]
                a[i] = [(x-factor*y) % p for x,y in zip(a[i],a[rank])]
        rank += 1
        if rank == rows:
            break
    return rank


def full_residue_profile(p,e):
    return ((0,0)+tuple(e*j for j in range(1,p-1))+((p-1)*e+1,)
            +tuple(e*j for j in range(p+1,2*p-1))+((2*p-1)*e-1,))


def profile(p,e,d):
    m = min(e,d)
    return ((0,0)+tuple(e*j for j in range(1,p-1))+((p-1)*e+1,p*e+m)
            +tuple(e*j for j in range(p+2,2*p-1))
            +((2*p-1)*e-1,2*p*e+d-m,(2*p+1)*e+3*d))


def divisor_formula(p,e,d):
    m = min(e,d)
    return ((0,)+tuple(e*(h-1)*(h-2)//2 for h in range(1,p+1))
            +(e*p*(p-1)//2+1,)
            +tuple(e*(h*(h-1)//2-(p+1))+1+m for h in range(p+2,2*p))
            +(e*(2*p*p-2*p-1)+m,e*(2*p*p-1)+d,2*e*p*(p+1)+4*d))


def units_for(p,d,variant):
    if variant == 0:
        base = list(range(p))
        duplicate = p**d
    elif variant == 1:
        base = [i+p*(i*i-3*i-5) for i in range(p)]
        duplicate = base[0]+p**d*(1+2*p)
    else:
        base = [(2*i+1) % p+p*(2*i*i-i-7) for i in range(p)]
        j = p//2
        base[0],base[j] = base[j],base[0]
        duplicate = base[0]+p**d*(p*p+2)
    require(len(set(x % p for x in base)) == p,"all distinct primary residues")
    require(vp(duplicate-base[0],p) == d,"exact secondary depth")
    return tuple(base+[duplicate])


def consecutive_profile(n,p):
    result = []
    for residue in range(p):
        size = len(range(residue,n,p))
        require(size <= p+1,"consecutive universe")
        if size == p+1:
            result.extend(profile(p,1,1))
        elif size == p:
            result.extend(full_residue_profile(p,1))
        elif size:
            result.extend((0,0)+tuple(range(1,size))+tuple(range(size+1,2*size)))
    return tuple(sorted(result))


def main():
    require(sys.argv[1:] == [],"no arguments")
    cases,rank_cases,minor_cases,consecutive = [],[],[],[]
    for p in (3,5,7,11):
        for e in (1,2,3):
            for d in (1,2,3,4):
                expected = profile(p,e,d)
                require(len(expected) == 2*p+2,"profile length")
                require(expected == tuple(sorted(expected)),"profile ordering")
                running = [0]
                for value in expected:
                    running.append(running[-1]+value)
                require(tuple(running) == divisor_formula(p,e,d),"divisor differences")
                for variant in range(3):
                    units = units_for(p,d,variant)
                    nodes = tuple(7-p*d+p**e*u for u in units)
                    primary = rational_smith(matrix(nodes),p)
                    independent = modular_smith(nodes,p)
                    require(primary == expected,("mixed formula",p,e,d,variant,primary))
                    require(independent == primary,"independent full modular replay")
                    old = rational_smith(matrix(nodes[:-1]),p)
                    require(old == full_residue_profile(p,e),"inherited residue control")
                    require((primary[:2*p] == old) == (d >= e),"exact adjoining boundary")
                    cases.append((p,e,d,variant,primary))
        for d in (1,3):
            for variant in range(3):
                units = units_for(p,d,variant)
                delta = units[-1]-units[0]
                for h in range(p+2,2*p+1):
                    bank = [[q*u**(q-1) if q else 0 for q in range(h)] for u in units]
                    differences = [x-y for x,y in zip(bank[-1],bank[0])]
                    require(all(x % delta == 0 for x in differences),"integral divided row")
                    bank[-1] = [x//delta for x in differences]
                    required_rank = p if h < 2*p else p+1
                    require(rank_mod(bank,p) == required_rank,"Frobenius plus duplicate rank")
                    rank_cases.append((p,d,variant,h,required_rank))
                h = 2*p+1
                full = [[comb(q,r)*u**(q-r) if q >= r else 0
                         for q in range(h)] for u in units[:-1] for r in (0,1)]
                last = [(q*units[-1]**(q-1)-q*units[0]**(q-1))//delta if q else 0
                        for q in range(h)]
                require(rank_mod(full+[last],p) == h,"final mixed Hasse witness")
    for e,d,variant in ((1,1,0),(2,1,1),(1,2,2),(2,3,1)):
        units = units_for(3,d,variant)
        nodes = tuple(4+3**e*u for u in units)
        actual = all_minor_valuations(nodes,3)
        require(actual == divisor_formula(3,e,d),"every exact minor independent path")
        minor_cases.append((3,e,d,variant,actual))
    for p,nvalues in ((3,(9,10,11,12)),(5,(25,26,28,30))):
        for n in nvalues:
            nodes = tuple(range(n))
            primary = rational_smith(matrix(nodes),p)
            independent = modular_smith(nodes,p)
            require(primary == consecutive_profile(n,p),"new complete consecutive band")
            require(primary == independent,"independent consecutive path")
            consecutive.append((p,n,primary))
    hostile = rational_smith(matrix((0,2,4)),2)
    require(hostile == (0,0,2,2,5,7),"inherited dyadic hostile")
    require(sum(hostile[:5]) != 7+1,"odd-prime final rank fails at p=2")
    require(profile(3,1,1) == (0,0,1,3,4,4,6,10),"prime3 empty range")
    require(profile(3,2,1) == (0,0,2,5,7,9,12,17),"unchanged-prefix hostile")
    report = dict(mixed=cases,residual=rank_cases,all_minors=minor_cases,consecutive=consecutive)
    digest = sha256(json.dumps(report,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    print("PASS: odd-prime p+1-node mixed two-jet Smith law")
    print("arbitrary-lift mixed matrices:",len(cases),"; rational/modular agreement")
    print("divided derivative rank cases:",len(rank_cases))
    print("exhaustive all-minor rank-eight cases:",len(minor_cases))
    print("direct consecutive matrices:",len(consecutive))
    for p in (3,5,7):
        for e,d in ((1,1),(1,3),(2,1),(2,2),(2,3)):
            print("profile",p,e,d,profile(p,e,d))
    for row in minor_cases:
        print("all-minor divisors",row)
    print("dyadic hostile:",hostile)
    print("semantic_sha256:",digest)
    print("exact_gates:",GATES)
    print("SCOPE: odd prime; p+1 nodes; one duplicated residue; e,d>=1.")
    print("OPEN: multiple duplicates, deeper cluster trees, higher jets, grid extremality.")


if __name__ == "__main__":
    main()
