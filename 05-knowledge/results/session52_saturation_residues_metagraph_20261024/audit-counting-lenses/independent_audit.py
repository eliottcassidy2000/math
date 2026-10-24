#!/usr/bin/env python3
"""Independent GF(2), literal-class, and quotient-matrix audit; no producer imports."""
import argparse
import hashlib
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd
from pathlib import Path

HERE = Path(__file__).resolve().parent
NTEST = 0

def check(condition, label):
    global NTEST
    NTEST += 1
    if not condition:
        raise RuntimeError(label)

def partitions(total, maximum=None):
    if total == 0:
        yield ()
        return
    maximum = min(total, maximum if maximum is not None else total)
    for first in range(maximum, 0, -1):
        for rest in partitions(total-first, first):
            yield (first,) + rest

def permutation_from_cycles(cycles):
    permutation = list(range(sum(cycles)))
    offset = 0
    for size in cycles:
        for t in range(size):
            permutation[offset+t] = offset+(t+1)%size
        offset += size
    return permutation

def linear_fixed_count(p, anti):
    """Rank of literal affine GF(2) equations, not pair-cycle inventory."""
    pairs = list(combinations(range(len(p)), 2))
    index = {pair:i for i,pair in enumerate(pairs)}
    dimension = len(pairs)
    pivots = {}
    for source,(a,b) in enumerate(pairs):
        u,v = p[a],p[b]
        target = index[tuple(sorted((u,v)))]
        lhs = (1<<source) ^ (1<<target)
        rhs = int(u>v) ^ anti
        while lhs:
            leading = lhs.bit_length()-1
            if leading not in pivots:
                pivots[leading] = (lhs,rhs)
                break
            row,bit = pivots[leading]
            lhs ^= row
            rhs ^= bit
        if lhs == 0 and rhs:
            return 0
    return 1 << (dimension-len(pivots))

def expected_count(cycles,anti):
    if anti:
        if cycles.count(1)>1 or any(c!=1 and c%4!=2 for c in cycles):
            return 0
        even = [c for c in cycles if c!=1]
        exponent = sum(even)//2 + sum(gcd(a,b) for a,b in combinations(even,2)) + cycles.count(1)*len(even)
    else:
        if any(c%2==0 for c in cycles):
            return 0
        exponent = sum((c-1)//2 for c in cycles) + sum(gcd(a,b) for a,b in combinations(cycles,2))
    return 2**exponent

def denominator(cycles):
    result = 1
    for length,multiplicity in Counter(cycles).items():
        result *= length**multiplicity*factorial(multiplicity)
    return result

def literal_relabel(mask,p,pairs):
    arrows = {(a,b) if mask&(1<<i) else (b,a) for i,(a,b) in enumerate(pairs)}
    image = {(p[a],p[b]) for a,b in arrows}
    return sum(1<<i for i,(a,b) in enumerate(pairs) if (a,b) in image)

def matrix_product(a,b):
    n = len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]

def characteristic_polynomial(matrix):
    """Faddeev-LeVerrier over exact integers; descending coefficients."""
    size = len(matrix)
    b = [[int(i==j) for j in range(size)] for i in range(size)]
    coefficients = [1]
    for k in range(1,size+1):
        product = matrix_product(matrix,b)
        trace = sum(product[i][i] for i in range(size))
        check(trace%k==0, ('Newton divisibility',size,k))
        c = -trace//k
        coefficients.append(c)
        b = [[product[i][j]+c*int(i==j) for j in range(size)] for i in range(size)]
    return coefficients

def expected_characteristic(coefficients,d):
    poly = [1]
    for level,multiplicity in enumerate(coefficients):
        eigenvalue = d-2*level
        for _ in range(multiplicity):
            new = [0]*(len(poly)+1)
            for i,value in enumerate(poly):
                new[i] += value
                new[i+1] -= eigenvalue*value
            poly = new
    return poly

def literal_classes(n):
    pairs = list(combinations(range(n),2))
    pms = list(permutations(range(n)))
    d = len(pairs)
    universe = set(range(1<<d))
    orbits = []
    class_of = {}
    while universe:
        representative = min(universe)
        orbit = {literal_relabel(representative,p,pairs) for p in pms}
        for mask in orbit:
            class_of[mask] = len(orbits)
        universe.difference_update(orbit)
        orbits.append((representative,len(orbit)))
    matrix = [[0]*len(orbits) for _ in orbits]
    self_converse = 0
    full = (1<<d)-1
    for i,(representative,size) in enumerate(orbits):
        self_converse += class_of[representative^full]==i
        for edge in range(d):
            matrix[i][class_of[representative^(1<<edge)]] += 1
    return orbits,self_converse,matrix

def fixed_path_mask(n,free_bits):
    pairs = list(combinations(range(n),2))
    free = iter((free_bits>>j)&1 for j in range((n-1)*(n-2)//2))
    return sum((1 if b==a+1 else next(free))<<i for i,(a,b) in enumerate(pairs))

def class_key(n,mask,merged=False):
    pairs = list(combinations(range(n),2))
    values = [literal_relabel(mask,p,pairs) for p in permutations(range(n))]
    if merged:
        full = (1<<len(pairs))-1
        values += [value^full for value in values]
    return min(values)

def native_controls():
    keys = [class_key(3,fixed_path_mask(3,t)) for t in (1,0)]
    matrix = [[0,0],[0,0]]
    for t in (0,1):
        source = keys.index(class_key(3,fixed_path_mask(3,t)))
        target = keys.index(class_key(3,fixed_path_mask(3,t^1)))
        matrix[source][target] += 1
    check(matrix==[[0,1],[1,0]],'literal native n3')
    check(characteristic_polynomial(matrix)==[1,0,-1],'literal native n3 characteristic')
    first_failure = None
    for n in range(1,5):
        dimension = (n-1)*(n-2)//2
        targets = {}
        for t in range(1<<dimension):
            source = class_key(n,fixed_path_mask(n,t),True)
            target = class_key(n,fixed_path_mask(n,t^((1<<dimension)-1)),True)
            targets.setdefault(source,set()).add(target)
        if any(len(values)>1 for values in targets.values()):
            first_failure = n
            break
    check(first_failure==4,'literal minimal merged non-descent')
    check(class_key(4,fixed_path_mask(4,0))==class_key(4,fixed_path_mask(4,1)),'n4 sources isomorphic')
    check(class_key(4,fixed_path_mask(4,7),True)!=class_key(4,fixed_path_mask(4,6),True),'n4 merged targets unequal')
    return {'n3_matrix':matrix,'first_non_descent_order':first_failure,'n4_source_tilings':[0,1],'n4_target_tilings':[7,6]}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--producer',type=Path,default=HERE/'inputs'/'producer')
    ap.add_argument('--receipt',type=Path,default=HERE/'audit-normal.json')
    args = ap.parse_args()
    external = json.loads((HERE/'inputs'/'external-oeis-prefixes.json').read_text())
    producer = json.loads((args.producer/'receipt.json').read_text())
    producer_input = json.loads((args.producer/'inputs'/'counts-0-20.json').read_text())
    result = []
    for n in range(21):
        sums = [Fraction(),Fraction()]
        count = 0
        for cycles in partitions(n):
            permutation = permutation_from_cycles(cycles)
            for anti in range(2):
                value = linear_fixed_count(permutation,anti)
                check(value==expected_count(cycles,anti),('GF2 criterion/exponent',cycles,anti,value))
                sums[anti] += Fraction(value,denominator(cycles))
            count += 1
        check(all(v.denominator==1 for v in sums),('Burnside integrality',n))
        a,sc = [int(v) for v in sums]
        check(a==external['A000568'][n],('external A',n))
        check(sc==external['A002785_with_empty_n0'][n],('external SC',n))
        check(a==producer_input['ordinary'][n] and sc==producer_input['self_converse'][n],('producer frozen baseline',n))
        row = producer['table'][n]
        check((a,sc)==(row['ordinary'],row['self_converse']),('producer output',n))
        check((a+sc)%2==0 and row['merged']==(a+sc)//2,('merged count',n))
        item = {'n':n,'partitions':count,'ordinary':a,'self_converse':sc,'merged':(a+sc)//2}
        if n<=5:
            orbits,literal_sc,matrix = literal_classes(n)
            check((len(orbits),literal_sc)==(a,sc),('literal class census',n))
            charpoly = characteristic_polynomial(matrix)
            check(charpoly==expected_characteristic(row['signed_polynomial_coefficients'],n*(n-1)//2),('literal quotient spectrum',n))
            item.update(literal_orbit_sizes=[size for _,size in orbits],arc_quotient=matrix,characteristic_polynomial=charpoly)
        result.append(item)
    hostiles = []
    for cycles in ((1,1),(4,),(1,1,2),(3,),(2,),(2,2),(6,),(1,6)):
        fixed = linear_fixed_count(permutation_from_cycles(cycles),1)
        hostiles.append({'cycle_type':cycles,'anti_fixed':fixed})
    expected_hostiles = [0,0,0,0,2,16,8,16]
    check([h['anti_fixed'] for h in hostiles]==expected_hostiles,'hostile exact values')
    # At n=3 our literal ordering is transitive then cyclic.
    check(result[3]['arc_quotient']==[[2,1],[3,0]],'arc n3 explicit')
    native = native_controls()
    # The true arc eigenfunction is centered cyclicity; uncentered is not an eigenfunction.
    arc = result[3]['arc_quotient']
    centered = [-1,3]  # four times c(T)-C(3,3)/4.
    check([sum(arc[i][j]*centered[j] for j in range(2)) for i in range(2)]==[1,-3],'centered cyclicity')
    check([arc[0][1],arc[1][1]]!=[0,-1],'uncentered cyclicity hostile')
    out = {'status':'PASS','checks':NTEST,'method':'Independent literal GF(2) equation ranks through n=20; direct labelled-orbit class census and exact quotient characteristic polynomials through n=5; separately retrieved OEIS prefixes.',
           'external_source_urls':external['source_urls'],'table':result,'hostiles':hostiles,'native_controls':native,
           'scope_warning':'No native fixed-HP wiring spectrum, novel counting algorithm, or independent native U_11 prediction is certified.',
           'centered_cyclicity_precision':'The eigenfunction is c(T)-binom(n,3)/4; uncentered c(T) is not a nonconstant adjacency eigenfunction.',
           'producer_receipt_sha256':hashlib.sha256((args.producer/'receipt.json').read_bytes()).hexdigest()}
    args.receipt.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','checks':NTEST,'receipt':str(args.receipt)}))

if __name__=='__main__':
    main()
