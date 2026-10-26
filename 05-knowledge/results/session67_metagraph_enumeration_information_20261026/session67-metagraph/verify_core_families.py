#!/usr/bin/env python3
"""Independent exact family/core-size verifier, stdlib only; no repo imports.

Algorithmic counts: pair-free fixed permutations via signed ordered-list
inversion, then positive-weight cycle-index assembly. Independent controls:
arc-flip class BFS, actual homogeneous-pair contraction, actual automorphisms,
explicit weight orbits, and full-permutation isomorphism minima through six.
All validation remains active under python -O.
"""
from collections import Counter, defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import comb, factorial, gcd
from pathlib import Path
import argparse
import hashlib
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def relabel(rows, p):
    return tuple(sum(((rows[p[i]] >> p[j]) & 1) << j for j in range(len(p)))
                 for i in range(len(p)))


def invariant_cells(rows):
    n = len(rows)
    colors = tuple(r.bit_count() for r in rows)
    while True:
        labels = sorted(set(colors))
        sigs = [(colors[i], tuple(sum((rows[i] >> j) & 1
                for j in range(n) if colors[j] == c) for c in labels))
                for i in range(n)]
        rank = {s: k for k, s in enumerate(sorted(set(sigs)))}
        new = tuple(rank[s] for s in sigs)
        if new == colors:
            return tuple(tuple(i for i in range(n) if colors[i] == c)
                         for c in sorted(set(colors)))
        colors = new


def canonical(rows):
    return min(relabel(rows, tuple(v for block in bs for v in block))
               for bs in product(*(permutations(c) for c in invariant_cells(rows))))


def autos(rows):
    cells = invariant_cells(rows)
    ans = []
    for choices in product(*(permutations(c) for c in cells)):
        p = list(range(len(rows)))
        for cell, image in zip(cells, choices):
            for i, j in zip(cell, image):
                p[i] = j
        if relabel(rows, p) == rows:
            ans.append(tuple(p))
    return ans


def cycle_type(p):
    seen, lengths = set(), []
    for v in range(len(p)):
        if v in seen:
            continue
        j, k = v, 0
        while j not in seen:
            seen.add(j)
            k += 1
            j = p[j]
        lengths.append(k)
    return tuple(sorted(lengths))


def homogeneous_pairs(rows):
    n = len(rows)
    return [(i,j) for i,j in pairs(n)
            if ((rows[i] ^ rows[j]) & ~((1 << i) | (1 << j))) == 0]


def core_decomposition(rows):
    n = len(rows)
    parent = list(range(n))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for i,j in homogeneous_pairs(rows):
        parent[root(i)] = root(j)
    comps = defaultdict(list)
    for i in range(n):
        comps[root(i)].append(i)
    chains = []
    for block in comps.values():
        chain = sorted(block, key=lambda i: -sum((rows[i] >> j) & 1 for j in block))
        require(all((rows[chain[i]] >> chain[j]) & 1
                    for i,j in pairs(len(chain))), 'component is not transitive')
        require(all(len({(rows[v] >> z) & 1 for v in chain}) == 1
                    for z in range(n) if z not in chain), 'component is not module')
        chains.append(tuple(chain))
    qrows = tuple(sum(((rows[b[0]] >> c[0]) & 1) << j
                      for j,c in enumerate(chains) if j != i)
                  for i,b in enumerate(chains))
    require(not homogeneous_pairs(qrows), 'reduced quotient has homogeneous pair')
    return qrows, tuple(map(len, chains))


def inflate(rows, weights):
    owner = [i for i,w in enumerate(weights) for _ in range(w)]
    out = [0] * len(owner)
    for a,b in pairs(len(owner)):
        if owner[a] == owner[b] or (rows[owner[a]] >> owner[b]) & 1:
            out[a] |= 1 << b
        else:
            out[b] |= 1 << a
    return tuple(out)


def census(n):
    first = canonical(tuple(sum(1 << j for j in range(i+1,n)) for i in range(n)))
    found, queue = {first}, deque([first])
    while queue:
        rows = queue.popleft()
        for i,j in pairs(n):
            r = list(rows)
            r[i] ^= 1 << j
            r[j] ^= 1 << i
            nxt = canonical(tuple(r))
            if nxt not in found:
                found.add(nxt)
                queue.append(nxt)
    return sorted(found)


def odd_partitions(n, lower=1):
    if n == 0:
        yield ()
    for d in range(lower, n+1, 2):
        for tail in odd_partitions(n-d, d):
            yield (d,) + tail


def exponent(lam):
    return sum((d-1)//2 for d in lam) + sum(gcd(a,b)
            for i,a in enumerate(lam) for b in lam[i+1:])


def centralizer(lam):
    z = 1
    for d,m in Counter(lam).items():
        z *= d**m * factorial(m)
    return z


@lru_cache(None)
def fixed_core(lam):
    """Number of pair-free labelled tournaments fixed by cycle type lam."""
    mult = sorted(Counter(lam).items())
    if not mult:
        return 1
    value = 0
    for ks in product(*(range(1,m+1) for d,m in mult)):
        small, factor = [], 1
        for (d,m), k in zip(mult,ks):
            small.extend([d]*k)
            factor *= (-d)**(m-k) * factorial(m)//factorial(k) * comb(m-1,k-1)
        value += factor * (1 << exponent(small))
    require(value >= 0, ('negative fixed-core count',lam,value))
    return value


def denominator_coeffs(lam, bound):
    out = [1] + [0]*bound
    for d in lam:
        for n in range(d,bound+1):
            out[n] += out[n-d]
    return out


def arithmetic_triangle(bound):
    table = [[Fraction(0) for _ in range(bound+1)] for _ in range(bound+1)]
    table[0][0] = Fraction(1)
    for q in range(1,bound+1):
        for lam in odd_partitions(q):
            mass = Fraction(fixed_core(lam),centralizer(lam))
            den = denominator_coeffs(lam,bound-q)
            for n in range(q,bound+1):
                table[n][q] += mass*den[n-q]
    for row in table:
        require(all(x.denominator == 1 and x>=0 for x in row),'noninteger triangle')
    return [[int(x) for x in row] for row in table]


def direct_bivariate_triangle(bound):
    """Independent coefficient expansion of Z_T(u^d x^d/(1-x^d+u^d x^d))."""
    ans = [[Fraction(0) for q in range(bound+1)] for n in range(bound+1)]
    ans[0][0] = Fraction(1)
    for mass in range(1,bound+1):
        for lam in odd_partitions(mass):
            poly = {(0,0):1}
            for d,m in sorted(Counter(lam).items()):
                factor = defaultdict(int)
                for j in range(bound//d-m+1):
                    for h in range(j+1):
                        factor[(d*(m+j),d*(m+h))] += comb(m+j-1,j)*comb(j,h)*(-1)**h
                nxt = defaultdict(int)
                for (n,q),a in poly.items():
                    for (i,k),b in factor.items():
                        if n+i<=bound:
                            nxt[(n+i,q+k)] += a*b
                poly = {k:v for k,v in nxt.items() if v}
            coefficient = Fraction(1 << exponent(lam),centralizer(lam))
            for (n,q),v in poly.items():
                ans[n][q] += coefficient*v
    require(all(v.denominator==1 and v>=0 for row in ans for v in row),
            'direct bivariate coefficient integrality')
    return [[int(v) for v in row] for row in ans]


def davis(n):
    x = sum((Fraction(1 << exponent(lam),centralizer(lam))
             for lam in odd_partitions(n)), Fraction(0))
    require(x.denominator == 1, 'Davis noninteger')
    return int(x)


def compositions(n,q):
    if q == 0:
        if n == 0:
            yield ()
        return
    if q == 1:
        if n >= 1:
            yield (n,)
        return
    for k in range(1,n-q+2):
        for tail in compositions(n-k,q-1):
            yield (k,)+tail


def orbit_weight(w,group):
    return min(tuple(w[p[i]] for i in range(len(w))) for p in group)


def verify(bound, class_bound):
    triangle = arithmetic_triangle(bound)
    counts = [davis(n) for n in range(bound+1)]
    require([sum(row) for row in triangle] == counts, 'triangle sums disagree with Davis')
    known = [1,1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168]
    require(counts[:min(len(known),bound+1)] == known[:bound+1],'known count mismatch')
    direct_bound=min(bound,14)
    direct=direct_bivariate_triangle(direct_bound)
    require(direct == [row[:direct_bound+1] for row in triangle[:direct_bound+1]],
            'independent bivariate expansion mismatch')
    tail_checks=0
    uniform_checks=0
    for n in range(1,bound+1):
        for d in range(1,n):
            require(sum(triangle[n][:n-d+1]) <= comb(n-1,d)*counts[n-d],
                    ('core-tail inflation bound',n,d))
            tail_checks += 1
        for q in range(1,n+1):
            labelled_mass=Fraction(1 << comb(q,2),factorial(q))
            rigid_core_lower=(Fraction(counts[q])-(q-1)*counts[q-1]
                              -Fraction(3,2)*(counts[q]-labelled_mass))
            lower=max(Fraction(0),rigid_core_lower)*comb(n-1,q-1)
            upper=counts[q]*comb(n-1,q-1)
            require(lower <= triangle[n][q] <= upper,('uniform row bounds',n,q))
            uniform_checks += 1
    classes, cores, groups = {}, {}, {}
    structural = []
    cycle_mass_all = {}
    for n in range(1,class_bound+1):
        bank = census(n)
        classes[n] = set(bank)
        require(len(bank) == counts[n], ('class BFS count',n))
        actual = Counter()
        for r in bank:
            q,w = core_decomposition(r)
            actual[len(q)] += 1
            require(canonical(inflate(q,w)) == r, 'core round trip')
        require(dict(actual) == {q:v for q,v in enumerate(triangle[n]) if v},
                ('core triangle mismatch',n,actual))
        cores[n] = [r for r in bank if not homogeneous_pairs(r)]
        cycle_mass = defaultdict(Fraction)
        aut_mass = Fraction(0)
        for r in bank:
            g = autos(r)
            groups[r] = g
            aut_mass += Fraction(factorial(n),len(g))
            require(len(g)%2 == 1, 'even automorphism order')
            if not homogeneous_pairs(r):
                for p in g:
                    cycle_mass[cycle_type(p)] += Fraction(1,len(g))
        require(aut_mass == 1 << comb(n,2), 'labelled orbit mass')
        expected = {lam:Fraction(fixed_core(lam),centralizer(lam))
                    for lam in odd_partitions(n) if fixed_core(lam)}
        require(dict(cycle_mass) == expected, ('cycle-index inversion mismatch',n))
        if n <= 6:
            brute = {min(relabel(r,p) for p in permutations(range(n))) for r in bank}
            require(len(brute) == len(bank),'full-permutation audit')
        generated = set()
        allocations = 0
        for q in range(1,n+1):
            for r in cores[q]:
                orbits = {orbit_weight(w,groups[r]) for w in compositions(n,q)}
                allocations += len(orbits)
                for w in orbits:
                    child = canonical(inflate(r,w))
                    require(child not in generated, 'duplicate canonical family')
                    generated.add(child)
                    require(len(autos(child)) == sum(tuple(w[p[i]] for i in range(q)) == w
                            for p in groups[r]), 'weight-stabilizer automorphism law')
        require(generated == classes[n],'weighted family coverage')
        cycle_mass_all[str(n)] = {','.join(map(str,lam)):str(mass)
                                  for lam,mass in sorted(cycle_mass.items())}
        structural.append({'n':n,'classes':len(bank),'pair_free_cores':len(cores[n]),
                           'unique_weight_allocations':allocations,
                           'core_distribution':dict(sorted(actual.items())),
                           'labelled_orbit_mass':int(aut_mass)})
    # Canonical hostiles: identical group orders lose permutation cycle types.
    qa = (0,33,11,19,7,29)
    qb = (10,20,33,38,13,19)
    ga,gb = autos(qa),autos(qb)
    require(not homogeneous_pairs(qa) and not homogeneous_pairs(qb), 'hostile cores')
    require(len(ga)==len(gb)==3,'hostile group orders')
    wa = len({orbit_weight(w,ga) for w in compositions(7,6)})
    wb = len({orbit_weight(w,gb) for w in compositions(7,6)})
    require((wa,wb)==(4,2),'cycle-type hostile')
    # Plain pair homogeneity is not an equivalence; connected paths are required.
    tt3=(6,4,0)
    require(homogeneous_pairs(tt3)==[(0,1),(1,2)],'homogeneity nontransitivity')
    # General quotient weight lists overcount: TT2[1,2]=TT2[2,1]=TT3.
    require(canonical(inflate((2,0),(1,2)))==canonical(inflate((2,0),(2,1))),
            'nonreduced quotient merger hostile')
    # Family degree bound: largest core size can be generated by complement of TTn.
    full_partner=[]
    for n in range(3,31):
        rows=[0]*n
        for i,j in pairs(n):
            if j==i+1:
                rows[i]|=1<<j
            else:
                rows[j]|=1<<i
        b=len(core_decomposition(tuple(rows))[0])
        full_partner.append((n,b))
        if n>=5:
            require(b==n,'transitive partner full core')
    return {'status':'PASS','arithmetic_order':bound,'structural_order':class_bound,
            'direct_bivariate_order':direct_bound,'tail_bound_checks':tail_checks,
            'uniform_row_bound_checks':uniform_checks,'ordinary_counts':counts,
            'pair_free_counts':[triangle[n][n] for n in range(bound+1)],
            'core_size_triangle':triangle,'structural':structural,
            'actual_core_cycle_index':cycle_mass_all,
            'hostiles':{'same_group_order_different_weight7':[wa,wb],
                        'TT3_homogeneous_pairs':homogeneous_pairs(tt3),
                        'full_partner_core_size':full_partner}}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-order',type=int,default=20)
    ap.add_argument('--class-order',type=int,default=7)
    ap.add_argument('--out',type=Path)
    args=ap.parse_args()
    require(1 <= args.class_order <= args.max_order,'bad bounds')
    report=verify(args.max_order,args.class_order)
    text=json.dumps(report,sort_keys=True,indent=2)+'\n'
    if args.out:
        args.out.write_text(text)
        print(json.dumps({'status':'PASS','report':str(args.out),'sha256':hashlib.sha256(text.encode()).hexdigest(),
                          'arithmetic_order':args.max_order,'structural_order':args.class_order},sort_keys=True))
    else:
        print(text,end='')

if __name__=='__main__':
    main()
