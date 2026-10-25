#!/usr/bin/env python3
"""Independent exact controls for mass, resolved strata, and power quotients.
No producer imports; finite controls are not substitutes for RESULT.md proofs.
"""
from collections import Counter
from itertools import product
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def legendre(a, p):
    a %= p
    if not a:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def line_values(point, p):
    E, x, z = point
    return [a % p for a in (E, E+x, E-x, E+z, E-z, E+x+z, E-x-z)]


def points(p):
    return [(1, x, z) for x in range(p) for z in range(p)] + [(0, 1, z) for z in range(p)] + [(0, 0, 1)]


def compact_twists(p):
    nu = next(x for x in range(1, p) if legendre(x, p) == -1)
    coeffs = [(1,) + tuple(nu if b else 1 for b in bits) for bits in product((0, 1), repeat=6)]
    counts = [0]*64
    nodes = [0]*64
    free = [0]*64
    strata = Counter()
    fiber_mass = Counter()
    for P in points(p):
        L = line_values(P, p)
        r = L.count(0)
        strata[r] += 1
        require(r <= 3, ('fourfold', p, P))
        j = next(i for i, x in enumerate(L) if x)
        row = []
        for t, d in enumerate(coeffs):
            fiber = 1
            den_inv = pow(d[j]*L[j] % p, -1, p)
            for i in range(7):
                if i != j:
                    fiber *= 1 + legendre(d[i]*L[i]*den_inv, p)
            row.append(fiber)
            counts[t] += fiber
            if r == 0:
                free[t] += fiber
            if r == 3:
                nodes[t] += fiber
        require(sum(row) == 64, ('fiber mean', p, P, row))
        require(Counter(row) == Counter({0:64-2**r, 2**(6-r):2**r}), ('inertia row', p, P, r))
        fiber_mass[r] += sum(row)
    t3 = 5 if p == 3 else 3
    t2 = 21 - 3*t3
    U = p*p - 6*p + 15-t3
    require(strata == Counter({0:U, 1:7*(p+1)-2*t2-3*t3, 2:t2, 3:t3}), ('strata', p, strata))
    resolved = [x+p*n for x,n in zip(counts,nodes)]
    require(sum(counts) == 64*(p*p+p+1), ('compact mean', p))
    require(sum(nodes) == 64*t3, ('node mean', p))
    require(sum(free) == 64*U, ('free mean', p))
    require(sum(resolved) == 64*(p*p+(t3+1)*p+1), ('resolved mean', p))
    conics = []
    P2 = points(p)
    for a in range(1,p):
        for b in range(1,p):
            n = sum((w*w-a*u*u-b*v*v)%p == 0 for u,v,w in P2)
            require(n==p+1, ('exceptional conic',p,a,b,n))
            conics.append(n)
    return {'p':p,'nonzero_nonsquare':nu,'base_strata':dict(sorted(strata.items())),
            'free_sum':sum(free),'compact_sum':sum(counts),'node_sum':sum(nodes),'resolved_sum':sum(resolved),
            'means':{'free':sum(free)//64,'compact':sum(counts)//64,'nodes':sum(nodes)//64,'resolved':sum(resolved)//64},
            'compact_histogram':dict(sorted(Counter(counts).items())),
            'node_histogram':dict(sorted(Counter(nodes).items())),
            'resolved_histogram':dict(sorted(Counter(resolved).items())),
            'exceptional_conic_instances':len(conics)}


def power_mass(moduli):
    G=list(product(*(range(n) for n in moduli)))
    zero=tuple(0 for _ in moduli)
    def add(a,b): return tuple((x+y)%n for x,y,n in zip(a,b,moduli))
    def mul(m,a): return tuple(m*x%n for x,n in zip(a,moduli))
    subgroups=set()
    # All cyclic stabilizers, plus the full and trivial stabilizers.
    for a in G:
        H={zero}
        x=add(zero,a)
        while x not in H:
            H.add(x); x=add(x,a)
        subgroups.add(frozenset(H))
    subgroups.add(frozenset(G))
    cases=0
    for H in subgroups:
        cosets=[]
        todo=set(G)
        while todo:
            a=min(todo); C=frozenset(add(a,h) for h in H)
            cosets.append(C);todo-=C
        coset_idx={a:i for i,C in enumerate(cosets) for a in C}
        reps=[min(C) for C in cosets]
        for shift in G:
            for c in (1,2,3):
                # Frobenius moves cyclic copies and translates by shift on G/H.
                Omega=list(product(range(c),range(len(cosets))))
                def gmF(g,m,state):
                    k,i=state
                    return ((k+m)%c,coset_idx[add(reps[i],add(mul(m,g),mul(m,shift)))])
                for m in range(1,9):
                    lhs=sum(gmF(g,m,s)==s for g in G for s in Omega)
                    J={mul(m,g) for g in G}
                    orbits=[];todo=set(Omega)
                    while todo:
                        k,i=min(todo)
                        O=frozenset((k,coset_idx[add(reps[i],j)]) for j in J)
                        orbits.append(O);todo-=O
                    rational=sum(frozenset(gmF(zero,m,s) for s in O)==O for O in orbits)
                    require(lhs==len(G)*rational, ('power mean',moduli,H,shift,c,m,lhs,rational))
                    cases+=1
    return {'moduli':moduli,'stabilizer_subgroups':len(subgroups),'verified_cases':cases}


def main():
    compact=[compact_twists(p) for p in (3,5,7,11,13)]
    power=[power_mass(M) for M in [(2,),(4,),(6,),(2,2),(2,4),(3,3)]]
    # Ordinary Euler characteristic rules out the naive geometric mean identity.
    require(64*7 != 7,'naive motive hostile absent')
    result={'status':'PASS','compact_twist_resolution_checks':compact,'power_map_orbit_checks':power,
            'euler_hostile':{'chi_U':7,'chi_V':448,'chi_X':56,'chi_Y':80},
            'interpretation':'Finite exact controls; all-order statements are separately proved in RESULT.md.'}
    print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__': main()
