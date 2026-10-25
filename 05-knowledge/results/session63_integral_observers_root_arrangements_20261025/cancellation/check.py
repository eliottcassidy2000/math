"""Portable exact arithmetic, subset-orbit, and tournament controls. Stdlib only."""
from collections import Counter
from itertools import product
from math import comb
from pathlib import Path
import hashlib
import json


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def primes_to(n):
    return [p for p in range(2, n+1) if all(p % q for q in range(2, int(p**0.5)+1))]


def prime_factors(n):
    return [p for p in primes_to(n) if n % p == 0]


def mul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                out[i+j] += x*y
    return out


def shift(poly):
    # Horner substitution t=x-1, deliberately different from the moment sum.
    out = [poly[-1]]
    for c in poly[-2::-1]:
        new = [0]*(len(out)+1)
        for i, a in enumerate(out):
            new[i] -= a
            new[i+1] += a
        new[0] += c
        out = new
    return out


def depth(poly, ell=0):
    # Repeated exact division by monic t+1, not shifted-coefficient scanning.
    q = [a % ell for a in poly] if ell else poly[:]
    d = 0
    while len(q) > 1:
        r = [q[0]]
        for a in q[1:-1]:
            r.append((a-r[-1]) % ell if ell else a-r[-1])
        if ((q[-1]-r[-1]) % ell if ell else q[-1]-r[-1]):
            break
        q = r
        d += 1
    return d


def elementary_poly(p, k, s=1):
    n = p**k
    degree = s*n
    numerator = [comb(degree, j) for j in range(degree+1)]
    for j in range(s*n//p+1):
        numerator[p*j] += (n-1)*comb(s*n//p, j)
    need(all(a % n == 0 for a in numerator), 'coefficientwise integral Burnside division')
    return [a//n for a in numerator]


def elementary_expected(p, k, s, ell=0):
    n = p**k
    if ell == 0:
        return s*n//p
    if ell == p:
        return s*n if s % p == 0 else s*n-p+1
    return s*n if (n-1) % ell == 0 else s*n//p


def regular(moduli, copies=1):
    points = list(product(*(range(p) for p in moduli)))
    index = {a:i for i,a in enumerate(points)}
    size = len(points)
    return [tuple(block*size+index[tuple((a[j]+v[j]) % moduli[j]
                                          for j in range(len(moduli)))]
                  for block in range(copies) for a in points) for v in points]


def cycles(perm):
    seen = set()
    result = []
    for i in range(len(perm)):
        if i in seen:
            continue
        cycle = []
        j = i
        while j not in seen:
            seen.add(j)
            cycle.append(j)
            j = perm[j]
        need(j == i, 'permutation cycle closes')
        result.append(cycle)
    return result


def act(mask, perm):
    return sum(1 << perm[i] for i in range(len(perm)) if mask >> i & 1)


def orbit_poly(group):
    n = len(group[0])
    seen = set()
    counts = [0]*(n+1)
    sizes = Counter()
    for mask in range(1 << n):
        if mask in seen:
            continue
        orbit = {act(mask, g) for g in group}
        need(not (orbit & seen), 'subset orbits partition')
        seen.update(orbit)
        counts[mask.bit_count()] += 1
        sizes[len(orbit)] += 1
    need(len(seen) == 1 << n, 'full declared subset universe visited')
    return counts, sizes


def burnside(group):
    n = len(group[0])
    total = [0]*(n+1)
    for g in group:
        poly = [1]
        for cycle in cycles(g):
            factor = [0]*(len(cycle)+1)
            factor[0] = factor[-1] = 1
            poly = mul(poly, factor)
        total = [a+b for a,b in zip(total, poly)]
    need(all(a % len(group) == 0 for a in total), 'generic Burnside divisibility')
    return [a//len(group) for a in total], total


def valuation(a, p):
    if a == 0:
        return None
    ans = 0
    while a % p == 0:
        a //= p
        ans += 1
    return ans


def rooted_refinement(adj):
    n = len(adj)
    colors = [int(i != 0) for i in range(n)]
    trace = [colors]
    while True:
        palette = sorted(set(colors))
        keys = [(colors[i], tuple(sum(adj[i][j] for j in range(n) if colors[j] == z)
                                  for z in palette)) for i in range(n)]
        rank = {key:i for i,key in enumerate(sorted(set(keys)))}
        new = [rank[key] for key in keys]
        trace.append(new)
        if len(set(new)) == len(set(colors)):
            return trace
        colors = new


def rooted_automorphisms(adj):
    # Exhaustive bijection search: only arc consistency and unused images.
    n = len(adj)
    answer = []
    nodes = 0
    def visit(mapping, used):
        nonlocal nodes
        nodes += 1
        if len(mapping) == n:
            answer.append(tuple(mapping[i] for i in range(n)))
            return
        choice = None
        for u in range(n):
            if u in mapping:
                continue
            candidates = [v for v in range(n) if v not in used and
                          all(adj[u][a] == adj[v][b] for a,b in mapping.items())]
            if not candidates:
                return
            if choice is None or len(candidates) < len(choice[1]):
                choice = (u, candidates)
        u, candidates = choice
        for v in candidates:
            mapping[u] = v
            visit(mapping, used | {v})
            del mapping[u]
    visit({0:0}, {0})
    return answer, nodes


def homogeneous_pairs(adj):
    n = len(adj)
    return [(u,v) for u in range(n) for v in range(u+1,n)
            if all(adj[u][w] == adj[v][w] for w in range(n) if w not in (u,v))]


def tournament50():
    points = list(product(range(5), repeat=2))
    S = {(0,1),(0,2),(1,1),(1,2),(1,3),(1,4),(2,0),(2,1),(2,2),(2,3),(2,4),(4,0)}
    neg = {tuple(-a % 5 for a in v) for v in S}
    need(not S & neg and len(S | neg) == 24, 'connection set chooses each antipodal pair once')
    adj = [[int(tuple((v[j]-u[j]) % 5 for j in range(2)) in S)
            for v in points] for u in points]
    need(all(sum(row) == 12 for row in adj), 'core regular of degree12')
    need(not homogeneous_pairs(adj), 'core has no homogeneous pair')
    trace = rooted_refinement(adj)
    need(len(set(trace[-1])) == 25, 'rooted color refinement discrete')
    stabilizer, nodes = rooted_automorphisms(adj)
    need(stabilizer == [tuple(range(25))], 'independent exact root stabilizer enumeration')
    translations = regular((5,5))
    need(len(set(translations)) == 25, '25 different translations')
    need(all(all(adj[u][v] == adj[g[u]][g[v]] for u in range(25) for v in range(25))
             for g in translations), 'all translations preserve every arc')
    expanded = [[int((u//2 == v//2 and u % 2 < v % 2) or
                     (u//2 != v//2 and adj[u//2][v//2])) for v in range(50)]
                for u in range(50)]
    pairs = homogeneous_pairs(expanded)
    need(pairs == [(2*i,2*i+1) for i in range(25)], 'exactly intended canonical clone pairs')
    need(all(expanded[i][j]+expanded[j][i] == 1 for i in range(50) for j in range(i+1,50)),
         'expanded object literally a tournament')
    polys, numerator = burnside(translations)
    need(polys == elementary_poly(5,2), 'literal gap-action cycle index equals elementary formula')
    return {'core_vertices':25, 'vertices':50, 'connection_set':sorted(S),
            'core_degree':12, 'core_homogeneous_pairs':0, 'clone_homogeneous_pairs':25,
            'rooted_refinement_colors':trace, 'refinement_class_counts':[len(set(c)) for c in trace],
            'independent_root_stabilizer_size':len(stabilizer), 'backtracking_nodes':nodes,
            'full_automorphisms':25, 'gap_orbits':1,
            'cycle_types':dict(Counter(','.join(map(str,sorted(map(len,cycles(g))))) for g in translations)),
            'orbit_polynomial':polys, 'shifted_polynomial':shift(polys),
            'depths':{'Q':depth(polys), **{str(p):depth(polys,p) for p in primes_to(31)}},
            'full_2_to_25_subset_census':False,
            'adjacency_sha256':hashlib.sha256(bytes(sum(adj,[]))).hexdigest()}


def main():
    arithmetic = []
    for p in (3,5,7,11,13,17):
        for k in range(1,6):
            if p**k > 343:
                continue
            n = p**k
            for s in range(1,7):
                if s*n > 500:
                    continue
                poly = elementary_poly(p,k,s)
                local = shift(poly)
                need(poly == poly[::-1], 'complement reciprocity')
                need(local == [sum(a*comb(i,j)*(-1)**(i-j) for i,a in enumerate(poly) if i>=j)
                               for j in range(len(poly))], 'Horner shift equals independent binomial moments')
                chars = sorted(set(primes_to(41)+prime_factors(n-1)+[p]))
                for ell in [0]+chars:
                    want = elementary_expected(p,k,s,ell)
                    got = depth(poly,ell)
                    need(got == want, f'depth spectrum p={p}, k={k}, s={s}, ell={ell}')
                    first = next(i for i,a in enumerate(local) if (a % ell if ell else a))
                    need(first == got, 'monic-division depth equals shifted coefficient depth')
                h = [0]+[comb(p,j)*(-1)**(p-j)//p for j in range(1,p)]
                reduced = [0]*(s*n+1)
                reduced[-1] = 1
                for j,a in enumerate(h):
                    reduced[s*n-p+j] = (reduced[s*n-p+j]-s*a) % p
                need([a % p for a in local] == reduced, 'full bad-prime polynomial congruence')
                need([a % 2 for a in poly] == [comb(s*n,j)%2 for j in range(s*n+1)], 'universal odd-action parity')
                arithmetic.append({'p':p,'k':k,'diagonal_copies':s,'points':s*n,
                                   'Q':depth(poly),'characteristic_p':depth(poly,p),'tested_primes':chars})
    small = []
    controls = [('empty',[()]), ('trivial3',[tuple(range(3))]),
                ('C3',regular((3,))), ('C5',regular((5,))),('C7',regular((7,))),
                ('C9',regular((9,))),('C3xC3',regular((3,3))),
                ('C3_diag2',regular((3,),2)),('C3_diag3',regular((3,),3)),
                ('C3_diag4',regular((3,),4)),('C5_diag2',regular((5,),2))]
    for label,group in controls:
        poly,sizes = orbit_poly(group)
        cyc,numerator = burnside(group)
        need(poly == cyc, 'literal subset orbits agree with generic Burnside')
        local = shift(poly)
        numer_local = shift(numerator)
        for ell in primes_to(31):
            e = valuation(len(group),ell)
            first = next(i for i,a in enumerate(numer_local) if a and valuation(a,ell) == e)
            need(first == depth(poly,ell), 'valuation-e integral averaging criterion')
        need([a%2 for a in poly] == [comb(len(poly)-1,j)%2 for j in range(len(poly))], 'all odd actions parity')
        small.append({'name':label,'points':len(group[0]),'group_order':len(group),
                      'subset_orbits':sum(poly),'orbit_size_histogram':dict(sorted(sizes.items()))})
    even, _ = orbit_poly([(0,1),(1,0)])
    need(even == [1,1,1] and depth(even,2) == 0, 'even transposition defeats universal dyadic law')
    need(depth(elementary_poly(3,1,3),3) == 9, 'diagonal three C3 orbits: depth9')
    independent = [1]
    for _ in range(3):
        independent = mul(independent,elementary_poly(3,1))
    need(depth(independent,3) == 3, 'independent three C3 actions: depth3')
    witness = tournament50()
    need(witness['depths']['Q'] == 5 and witness['depths']['5'] == 21 and
         witness['depths']['2'] == witness['depths']['3'] == 25, '50-vertex specialization')
    print(json.dumps({'status':'PASS','scope':'finite checks support, not replace, all-parameter proofs',
                      'arithmetic_cases':len(arithmetic),'arithmetic_universe':arithmetic,
                      'literal_subset_controls':small,'even_action_hostile':[1,1,1],
                      'diagonal_versus_product_C3_char3_depths':[9,3],
                      'tournament50':witness},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
