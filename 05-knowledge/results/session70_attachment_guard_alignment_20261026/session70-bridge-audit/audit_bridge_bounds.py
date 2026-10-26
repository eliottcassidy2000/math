#!/usr/bin/env python3
"""Independent, exact Session70 audit. No producer imports; explicit checks survive -O."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, factorial
import hashlib
import json

CHECKS = Counter()


def check(condition, label):
    CHECKS[label] += 1
    if not condition:
        raise RuntimeError(label)


def seed(mask):
    pairs = [(i, j) for i in range(7) for j in range(i+2, 7)]
    arcs = {(i, i+1) for i in range(6)}
    arcs.update((i, j) if mask & (1 << b) else (j, i)
                for b, (i, j) in enumerate(pairs))
    return arcs


def partner(arcs, n):
    hinges = {(i, i+1) for i in range(n-1)}
    return {(i, j) if (i, j) in hinges else (j, i) for i, j in arcs}


def transport(arcs, p):
    return {(p[i], p[j]) for i, j in arcs}


def path_guard(arcs, order):
    return all((u, v) in arcs for u, v in zip(order, order[1:]))


def strong(arcs, n):
    for e in (arcs, {(v, u) for u, v in arcs}):
        seen = {0}
        pending = [0]
        while pending:
            v = pending.pop()
            for a, b in e:
                if a == v and b not in seen:
                    seen.add(b)
                    pending.append(b)
        if len(seen) != n:
            return False
    return True


def native(bits, k, seeds):
    arcs = set()
    for i in range(k):
        arcs.update((7*i+a, 7*i+b) for a, b in seeds[(bits >> i) & 1])
        arcs.update((a, b) for a in range(7*i, 7*i+7)
                    for b in range(7*i+7, 7*k))
    arcs.update((a, 7*k) for a in range(7*k))
    arcs.add((7*k, 7*k+1))
    arcs.update((7*k+1, a) for a in range(7*k))
    return arcs


def rank_mod(matrix, p):
    rows = [[x % p for x in r] for r in matrix]
    rank = 0
    for c in range(len(rows[0])):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][c]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        inv = pow(rows[rank][c], -1, p)
        rows[rank] = [(x * inv) % p for x in rows[rank]]
        for i in range(rank+1, len(rows)):
            a = rows[i][c]
            if a:
                rows[i] = [(x-a*y) % p for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def audit_alignment():
    seeds = [seed(1056), seed(4160)]
    expected_rows = [(2,12,41,17,39,75,31), (2,20,9,51,37,71,31)]
    for e, rows in zip(seeds, expected_rows):
        check(tuple(sum(1 << v for u, v in e if u == i) for i in range(7)) == rows,
              'seed independently decoded from free mask')
        check(strong(e, 7), 'seed strong')
    probe = (0,1,2,3,5,6,4)
    check([path_guard(e, probe) for e in seeds] == [False, True], 'seed native guard')
    expected_maps = [(0,1,4,2,3,5,6), (0,2,4,1,5,3,6)]
    maps = []
    aut_counts = []
    for a, b in ((seeds[0], seeds[1]), (partner(seeds[0], 7), partner(seeds[1], 7))):
        isos, autos = [], []
        for p in permutations(range(7)):
            mapped = transport(a, p)
            if mapped == b:
                isos.append(p)
            if mapped == a:
                autos.append(p)
        check(len(isos) == len(autos) == 1, 'unique unrestricted seed isomorphism and rigidity')
        maps.append(isos[0])
        aut_counts.append(len(autos))
    check(maps == expected_maps, 'source target map orientation')
    inv_dst = {v:i for i,v in enumerate(maps[1])}
    defect = tuple(inv_dst[v] for v in maps[0])
    check(defect == (0,3,2,1,5,4,6), 'relative defect equals (1 3)(4 5)')
    check(transport(seeds[0], defect) != seeds[0], 'relative defect is not source automorphism')
    tables = []
    for k in range(1, 7):
        n, size = 7*k+2, 1 << k
        matrix = [[0] * size for _ in range(size)]
        base = native(0, k, seeds)
        for b in range(size):
            e = native(b, k, seeds)
            check(len(e) == comb(n, 2), 'native arc count')
            check(strong(e, n), 'native strong source')
            check(path_guard(e, tuple(range(n))), 'standard path legal')
            degree = Counter(u for u, v in e)
            check(degree[7*k] == 1 and degree[7*k+1] == 7*k and
                  all(2 <= degree[v] <= 7*k-1 for v in range(7*k)),
                  'unique wrapper degree strata')
            check(degree[0]+degree[n-1] != n-1, 'black endpoint degree obstruction')
            for local, initial, target in ((maps[0], base, e),
                    (maps[1], partner(base,n), partner(e,n))):
                p = tuple(7*i+(local[j] if (b >> i) & 1 else j)
                          for i in range(k) for j in range(7))+(7*k,7*k+1)
                check(p[0] == 0 and p[-1] == n-1 and transport(initial,p) == target,
                      'same jointly rooted endpoint types')
            for s in range(size):
                order = tuple(7*i+j for i in range(k)
                              for j in (probe if (s >> i) & 1 else range(7)))+(7*k,7*k+1)
                matrix[s][b] = int(path_guard(e, order))
                check(matrix[s][b] == int((s & b) == s), 'simultaneous native subset guard')
        # This inverse is constructed independently from native guards; integer
        # multiplication checks the all-field statement, including characteristic 2.
        inverse = [[((-1) ** ((s ^ b).bit_count())) if (b & s) == b else 0
                    for s in range(size)] for b in range(size)]
        for s in range(size):
            for t in range(size):
                check(sum(matrix[s][b]*inverse[b][t] for b in range(size)) == int(s == t),
                      'exact integral Mobius inverse')
        trunc = []
        for t in range(k+1):
            rows = [matrix[s] for s in range(size) if s.bit_count() <= t]
            d = sum(comb(k,j) for j in range(t+1))
            for p in (2,3,5):
                check(rank_mod(rows,p) == d, 'truncated rank in three characteristics')
            trunc.append([t,d,size-d])
        laws = [[Fraction(int(b.bit_count()%2 == parity), size//2)
                 for b in range(size)] for parity in (0,1)]
        moments = [[sum(matrix[s][b]*law[b] for b in range(size))
                    for s in range(size)] for law in laws]
        check(moments[0][:-1] == moments[1][:-1], 'normalized proper moments equal')
        check(moments[k%2][-1] == Fraction(1,size//2) and moments[1-k%2][-1] == 0,
              'normalized native full guard gap')
        for s in range(size-1):
            marginals = [Counter() for _ in (0,1)]
            for parity in (0,1):
                for b, weight in enumerate(laws[parity]):
                    marginals[parity][b & s] += weight
            check(marginals[0] == marginals[1], 'all proper full marginals equal')
        tables.append({'k':k,'n':n,'states':size,'guards':size,
                       'native_guard_tests':size**2,'rank_rows_t_d_missing':trunc,
                       'parity_gap':str(Fraction(1,size//2)),
                       'matrix_sha256':hashlib.sha256(bytes(x for r in matrix for x in r)).hexdigest()})
    return {'seed_automorphism_orders':aut_counts,'maps':maps,'relative_defect':defect,'tables':tables}


def edge(a,b):
    return tuple(sorted((a,b)))


def matchings(vertices):
    if not vertices:
        yield ()
    else:
        first, rest = vertices[0], vertices[1:]
        for b in rest:
            for pairing in matchings(tuple(x for x in rest if x != b)):
                yield (edge(first,b),)+pairing


def graph_facts(n, edges):
    adjacency = [set() for _ in range(n)]
    for a,b in edges:
        adjacency[a].add(b)
        adjacency[b].add(a)
    unseen, parts = set(range(n)), []
    while unseen:
        stack = [unseen.pop()]
        part = set(stack)
        while stack:
            for b in adjacency[stack.pop()]:
                if b in unseen:
                    unseen.remove(b)
                    part.add(b)
                    stack.append(b)
        parts.append(part)
    return sorted(map(len,parts)), [len(a) for a in adjacency], len(edges) == n-len(parts)


def connect(pairing, choices):
    n = 2*len(pairing)
    edges = set(pairing)
    for a,b in choices:
        edges.update((edge(a,n),edge(n,b)))
        n += 1
    return n, edges


def audit_pairings():
    rows = []
    for r in range(1,6):
        bank = list(matchings(tuple(range(2*r))))
        expected = factorial(2*r)//(2**r*factorial(r))
        check(len(bank) == expected, 'endpoint matching cardinality')
        profiles = {tuple(edge(a,b) not in p for a,b in combinations(range(2*r),2)) for p in bank}
        check(len(profiles) == expected, 'one connector profiles separate every matching')
        contexts = 0
        if r <= 4:
            for p in bank:
                mates = {v:u for a,b in p for v,u in ((a,b),(b,a))}
                for q in bank:
                    if p == q:
                        continue
                    a,b = min(set(q)-set(p))
                    choices = [(a,b)]
                    occupied = {a,b,mates[a],mates[b]}
                    last = mates[b]
                    for c,d in p:
                        if c not in occupied:
                            choices.append((last,c))
                            occupied.update((c,d))
                            last = d
                    check(len(choices) == r-1, 'context length r-1')
                    for j in range(1,len(choices)+1):
                        n, pe = connect(p,choices[:j])
                        parts,deg,forest = graph_facts(n,pe)
                        check(forest and max(deg) <= 2, 'every positive prefix remains a path forest')
                    n,pe = connect(p,choices)
                    _,qe = connect(q,choices)
                    parts,deg,forest = graph_facts(n,pe)
                    qparts,qdeg,qforest = graph_facts(n,qe)
                    check(parts == [n] and sorted(deg) == [1,1]+[2]*(n-2) and forest,
                          'final positive object is a spanning path')
                    check(not qforest and max(qdeg) <= 2, 'same final context has irrevocable negative cycle')
                    contexts += 1
        rows.append({'r':r,'states':expected,'ordered_final_contexts':contexts})
    return rows


def cycles(vertices, edges):
    first = min(vertices)
    result = []
    for tail in permutations(sorted(set(vertices)-{first})):
        if tail[0] > tail[-1]:
            continue
        word = (first,)+tail
        e = {edge(word[i],word[(i+1)%len(word)]) for i in range(len(word))}
        if e <= edges:
            result.append(e)
    return result


def audit_petersen():
    outer = {edge(i,(i+1)%5) for i in range(5)}
    inner = {edge(5+i,5+(i+2)%5) for i in range(5)}
    spokes = {(i,i+5) for i in range(5)}
    pet = outer|inner|spokes
    quotient = {edge(a%5,b%5) for a,b in outer|inner}
    check(quotient == set(combinations(range(5),2)), 'all-spoke contraction K5')
    qcycles = cycles(range(5),quotient)
    check(len(qcycles) == 12, 'K5 bare Hamilton cycles')
    two_factors = []
    for chosen in combinations(sorted(pet),10):
        if all(sum(v in e for e in chosen) == 2 for v in range(10)):
            two_factors.append(frozenset(chosen))
    check(len(two_factors) == 6, 'native Petersen two-factor count')
    check(all(graph_facts(10,e)[0] == [5,5] for e in two_factors), 'all native two-factors are two pentagons')
    lifted = []
    for mask in range(1<<10):
        chosen = {e for i,e in enumerate(sorted(quotient)) if mask & (1<<i)}
        lifted_nonspokes = {e for e in outer|inner if edge(e[0]%5,e[1]%5) in chosen}
        left = [sum(i in e for e in lifted_nonspokes) for i in range(5)]
        right = [sum(i+5 in e for e in lifted_nonspokes) for i in range(5)]
        lawful = all(left[i] == right[i] and left[i] in (1,2) for i in range(5))
        if lawful:
            e = frozenset(lifted_nonspokes|{(i,i+5) for i in range(5) if left[i] == 1})
            lifted.append(e)
            check(all(sum(v in a for a in e) == 2 for v in range(10)), 'local port rule lifts exactly')
    check(set(lifted) == set(two_factors), 'quotient port rule neither loses nor invents 2-factors')
    balanced = [e for e in qcycles if all(sum(i in a for a in e&outer) == 1 for i in range(5))]
    check(not balanced, 'no K5 Hamilton cycle uses balanced spoke ports')
    # Delete 0 and suppress exactly its three degree-two neighbours.
    missing = {e for e in pet if 0 not in e}
    chains = [(2,1,6),(3,4,9),(7,5,8)]
    branch = {2,3,6,7,8,9}
    reduced = {e for e in missing if set(e) <= branch}
    marks = {edge(p[0],p[-1]) for p in chains}
    for chain in chains:
        check(all(edge(a,b) in missing for a,b in zip(chain,chain[1:])), 'actual suppression chains')
    reduced |= marks
    parts = [(set(s),branch-set(s)) for s in combinations(sorted(branch),3)
             if reduced == {edge(a,b) for a in s for b in branch-set(s)}]
    check(len(parts) == 2, 'suppressed graph exact K33')
    small_cycles = cycles(branch,reduced)
    native_deleted_cycles = cycles(range(1,10),missing)
    expanded = []
    for c in small_cycles:
        e = set(c)
        for path in chains:
            key = edge(path[0],path[-1])
            if key in e:
                e.remove(key)
                e.update(edge(a,b) for a,b in zip(path,path[1:]))
        expanded.append(e)
        check(set(v for a in e for v in a) == set(range(1,10)) if marks <= c else
              set(v for a in e for v in a) != set(range(1,10)), 'K33 lift spanning iff all subdivision marks used')
    spanning_expanded = [e for c,e in zip(small_cycles,expanded) if marks <= c]
    check(set(map(frozenset,spanning_expanded)) == set(map(frozenset,native_deleted_cycles)),
          'K33 marked response reconstructs all deleted-vertex Hamilton cycles')
    return {'Petersen_two_factors':len(two_factors),'Petersen_component_sizes':[[5,5]]*len(two_factors),
            'K5_Hamilton_cycles':len(qcycles),'port_lifts':len(lifted),'balanced_K5_cycles':len(balanced),
            'K33_parts':[sorted(s) for s in parts[0]],'K33_subdivision_chains':chains,
            'K33_Hamilton_cycles':len(small_cycles),
            'K33_Hamilton_cycles_using_all_subdivision_marks':len(spanning_expanded),
            'Petersen_minus_0_Hamilton_cycles':len(native_deleted_cycles)}


def main():
    result = {'alignment':audit_alignment(),'endpoint_pairing':audit_pairings(),'Petersen':audit_petersen()}
    # Non-assert negative control confirms an intentionally forgotten port is rejected.
    try:
        check(result['Petersen']['balanced_K5_cycles'] == 12, 'hostile bare quotient lift rejected')
    except RuntimeError:
        result['negative_control'] = 'EXPECTED REJECTION'
    else:
        raise RuntimeError('hostile control unexpectedly accepted')
    result.update(status='PASS',producer_imports=False,checks=dict(sorted(CHECKS.items())),
                  total_checks=sum(CHECKS.values()))
    print(json.dumps(result,sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
