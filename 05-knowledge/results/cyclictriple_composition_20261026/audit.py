"""Independent coordinate arcs, walk inclusion-exclusion, and Walsh spectrum.

Imports no producer code or receipt. Standard library only.
"""
import itertools
import json
import math
from collections import Counter


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def arc(x, y, word):
    i, a = divmod(x, 4)
    j, b = divmod(y, 4)
    if i != j:
        return 0 < (j-i) % len(word) <= len(word)//2
    if a == b:
        return False
    if a == 3:
        return bool(word[i])
    if b == 3:
        return not word[i]
    return (b-a) % 3 == 1


def paths_by_walk_sieve(word):
    n = 4*len(word)
    neighbors = [[j for j in range(n) if arc(i, j, word)] for i in range(n)]
    total = 0
    # Sieve n-vertex walks of length n-1 by their visited set.
    for mask in range(1, 1 << n):
        vertices = [i for i in range(n) if mask >> i & 1]
        v = [int(i in vertices) for i in range(n)]
        for _ in range(n-1):
            v = [sum(v[j] for j in neighbors[i]) if mask >> i & 1 else 0 for i in range(n)]
        total += (-1)**(n-len(vertices))*sum(v)
    return total


def rotate_bits(x, m, k):
    return ((x << k) | (x >> (m-k))) & ((1 << m)-1)


def orbit(x, m):
    return sorted({rotate_bits(x, m, k) for k in range(m)})


def rank_mod(a, p=1000003):
    a = [[x % p for x in row] for row in a]
    r = 0
    for c in range(len(a[0])):
        pivot = next((j for j in range(r, len(a)) if a[j][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [x*inv % p for x in a[r]]
        for j in range(r+1, len(a)):
            t = a[j][c]
            a[j] = [(x-t*y) % p for x, y in zip(a[j], a[r])]
        r += 1
    return r


def spectral_audit():
    rows = []
    for m in [3, 5, 7, 9]:
        orbits = {min(orbit(x, m)): orbit(x, m) for x in range(1 << m)}
        reps = sorted(orbits)
        edges = [Counter(min(orbit(x ^ (1 << i), m)) for i in range(m)) for x in reps]
        # Walsh sums over subset-rotation orbits, evaluated on vertex orbits.
        f = [[sum((-1)**((x & y).bit_count()) for y in orbits[s]) for s in reps] for x in reps]
        index = {x:i for i, x in enumerate(reps)}
        for i in range(len(reps)):
            for j, s in enumerate(reps):
                lhs = sum(c*f[index[x]][j] for x, c in edges[i].items())
                require(lhs == (m-2*s.bit_count())*f[i][j], 'exact Walsh eigenvector')
        require(rank_mod(f) == len(reps), 'complete independent eigenbasis')
        weights = Counter(x.bit_count() for x in reps)
        for k in range(m+1):
            numerator = sum(sum(math.gcd(t,d)==1 for t in range(1,d+1))*math.comb(m//d,k//d)
                            for d in range(1,m+1) if m%d==0 and k%d==0)
            require(numerator == m*weights[k], 'weight-refined Burnside formula')
        primitive = {}
        for d in range(1,m+1):
            if m % d == 0:
                primitive[d] = 2**d-sum(v for e,v in primitive.items() if d%e==0)
        require(sum(primitive[d]//d for d in primitive) == len(reps), 'divisor recurrence independent of Mobius')
        require(len(edges[0]) != m, 'discarding move multiplicity destroys cube operator')
        rows.append({'m':m, 'matrix_size':len(reps), 'walsh_rank':len(reps),
                     'spectrum':{str(m-2*k):weights[k] for k in range(m+1)}})
    return rows


def emblem_audit():
    records = []
    triple_sets = []
    for word in [(0,0,0), (1,0,0)]:
        triples = {s for s in itertools.combinations(range(12),3)
                   if all(sum(arc(x,y,word) for y in s)==1 for x in s)}
        triple_sets.append(triples)
        aut = 0
        # The fixed strong-module tree forces these candidate permutations.
        # Unlike the producer, all 6 root permutations and all 6^3 triangle
        # permutations are tried, then the literal 132 directed arcs checked.
        for root in itertools.permutations(range(3)):
            for inner in itertools.product(list(itertools.permutations(range(3))), repeat=3):
                p = [4*root[i]+(inner[i][j] if j<3 else 3) for i in range(3) for j in range(4)]
                aut += all(arc(x,y,word)==arc(p[x],p[y],word) for x in range(12) for y in range(12))
        h = paths_by_walk_sieve(word)
        records.append({'word':list(word), 'automorphisms':aut, 'hamiltonian_paths':h})
    require(triple_sets[0] == triple_sets[1] and len(triple_sets[0])==67, 'literal C3 equality')
    require([x['automorphisms'] for x in records]==[81,27], 'independent stabilizers')
    require(all(x['hamiltonian_paths']==462915 for x in records), 'walk sieve versus exact target')
    return records


def paley_audit():
    # A noncyclic prime quotient, independent of all carousel routines.
    n = 7
    def edge(i,j): return (j-i)%7 in (1,2,4)
    group = [p for p in itertools.permutations(range(n))
             if all(edge(i,j)==edge(p[i],p[j]) for i in range(n) for j in range(n))]
    require(len(group)==21, 'Paley seven exact permutation group')
    for mask in range(1,(1<<n)-1):
        inside = [i for i in range(n) if mask>>i&1]
        outside = [i for i in range(n) if not mask>>i&1]
        module = all(len({edge(i,j) for j in inside})==1 for i in outside)
        require(not module or len(inside)==1, 'Paley seven prime')
    require(all(edge(i,j)==edge((-j)%7,(-i)%7) for i in range(n) for j in range(n)), 'self converse quotient')
    def move(x,p): return sum(1<<p[i] for i in range(n) if x>>i&1)
    orbits = {min(move(x,p) for p in group):sorted({move(x,p) for p in group}) for x in range(1<<n)}
    reps = sorted(orbits); index = {x:i for i,x in enumerate(reps)}
    canonical = {x:r for r,orb in orbits.items() for x in orb}
    c = [Counter(canonical[x^(1<<i)] for i in range(n)) for x in reps]
    f = [[sum((-1)**((x&y).bit_count()) for y in orbits[s]) for s in reps] for x in reps]
    for i,x in enumerate(reps):
        for y,b in c[i].items():
            require(len(orbits[x])*b==len(orbits[y])*c[index[y]][x], 'Paley flow balance')
        for j,s in enumerate(reps):
            require(sum(b*f[index[y]][j] for y,b in c[i].items())==(n-2*s.bit_count())*f[i][j], 'Paley Walsh eigenbasis')
    require(rank_mod(f)==len(reps), 'Paley full spectral rank')
    # Cycle-index polynomial without constructing any subset orbits.
    coeff = [0]*(n+1)
    for p in group:
        seen=set(); lengths=[]
        for i in range(n):
            if i in seen: continue
            j=i; length=0
            while j not in seen:
                seen.add(j);length+=1;j=p[j]
            lengths.append(length)
        poly=[1]+[0]*n
        for length in lengths:
            poly=[poly[k]+(poly[k-length] if k>=length else 0) for k in range(n+1)]
        coeff=[x+y for x,y in zip(coeff,poly)]
    weights=Counter(x.bit_count() for x in reps)
    require(all(coeff[k]==21*weights[k] for k in range(n+1)), 'Paley cycle index')
    aut=Counter(3**7*21//len(orb) for orb in orbits.values())
    require(aut=={45927:2,6561:6,2187:4}, 'Paley decorated automorphism strata')
    # Periodicity of words does not guarantee primitive prime divisors.
    primitive6={min(orbit(x,6)) for x in range(64) if len(orbit(x,6))==6}
    require(len(primitive6)==9 and 2**6-1==3**2*7, 'primitive-word / primitive-prime hostile')
    return {'quotient_vertices':7,'quotient_automorphisms':21,'classes':len(reps),
            'weight_classes':[weights[k] for k in range(n+1)],
            'decorated_automorphism_orders':dict(aut),'walsh_rank':len(reps),
            'primitive_necklaces_length6':len(primitive6)}


if __name__ == '__main__':
    print(json.dumps({'emblem':emblem_audit(), 'spectra':spectral_audit(), 'paley7':paley_audit()}, indent=2, sort_keys=True))
