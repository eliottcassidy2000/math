"""Independent literal controls for the closed-middle endpoint-three proof.

No imports from proof/probe packages. Complete small standard-blue universes
and deterministic larger sources; no graph-isomorphism or mixedness filters.
"""
from collections import Counter
from itertools import combinations
from pathlib import Path
import hashlib
import json
import random


def need(p, message):
    if not p:
        raise RuntimeError(message)


def arc(t, a, b):
    return bool(t[a] >> b & 1)


def path(t, p):
    return len(set(p)) == len(t) and all(arc(t, a, b) for a, b in zip(p, p[1:]))


def blue(t, p):
    return all(arc(t, p[i], p[j]) == arc(t, p[-1-j], p[-1-i])
               for i in range(len(p)) for j in range(i+1, len(p)))


def partner(t, p):
    retained = {frozenset((a, b)) for a, b in zip(p, p[1:])}
    result = [0]*len(t)
    for a, b in combinations(range(len(t)), 2):
        value = arc(t, a, b) ^ (frozenset((a, b)) not in retained)
        result[a if value else b] |= 1 << (b if value else a)
    return tuple(result)


def middle_orbits(m):
    seen, free, fixed = set(), [], []
    for a, b in combinations(range(m), 2):
        if (a, b) in seen:
            continue
        orbit = {(a, b), (m-1-b, m-1-a)}
        seen |= orbit
        (fixed if b-a == 1 else free).append(orbit)
    return fixed, free


def source(n, incoming, bits):
    m = n-2
    fixed, free = middle_orbits(m)
    t = [0]*n
    for index, orbit in enumerate(fixed+free):
        value = True if index < len(fixed) else bool(bits >> (index-len(fixed)) & 1)
        for a, b in orbit:
            a, b = a+1, b+1
            t[a if value else b] |= 1 << (b if value else a)
    outgoing = {n-1-u for u in incoming}
    for v in range(1, n-1):
        t[v if v in incoming else 0] |= 1 << (0 if v in incoming else v)
        t[n-1 if v in outgoing else v] |= 1 << (v if v in outgoing else n-1)
    t[-1] |= 1
    t = tuple(t)
    p = tuple(range(n))
    need(path(t, p) and blue(t, p), 'constructed blue source')
    need(t[0].bit_count() == n-4 and t[-1].bit_count() == 3, 'endpoint scores')
    return t


def full_blue(n, bits):
    fixed, free = middle_orbits(n)
    t = [0]*n
    for index, orbit in enumerate(fixed+free):
        value = True if index < len(fixed) else bool(bits >> (index-len(fixed)) & 1)
        for a,b in orbit:
            t[a if value else b] |= 1 << (b if value else a)
    return tuple(t)


def words(n, incoming, epsilon):
    a, b, x, c, z = 0, 1, 2, n-2, n-1
    w = tuple(range(2, n-2))
    if c in incoming and n-3 not in incoming:
        h = next(v for v in incoming if v != c)
        k = n-1-h
        l,r = tuple(range(1,k)),tuple(range(k,n-2))
        return 'movable_last_cut', ((a,)+l+(z,)+r+(c,), l+(z,)+r+(a,c), l+r+(a,c,z))
    if set(incoming) == {x,c} and not epsilon:
        d,y = tuple(range(3,n-3)),n-3
        return 'paired_last_second', ((a,b,x)+d+(z,y,c), (b,x)+d+(z,y,a,c), (b,x)+d+(y,a,c,z))
    if not epsilon:
        return 'epsilon0', ((c,)+w+(b,z,a), (z,c)+w+(b,a), (a,z,b)+w+(c,))
    if x not in incoming:
        return 'epsilon1', ((c,b,z,a)+w, (z,c,b,a)+w, (a,z,b)+w+(c,))
    if 3 not in incoming:
        u, y = tuple(range(3, n-3)), n-3
        return 'second_port', ((a,b,x)+u+(z,y,c), (a,x)+u+(z,y,b,c), (a,x)+u+(y,b,c,z))
    u, e, y = tuple(range(3, n-4)), n-4, n-3
    return 'two_second_ports', ((a,b,x)+u+(z,e,y,c), (a,x)+u+(z,e,y,b,c), (a,x)+u+(e,y,b,c,z))


def replay(t, incoming, stats, strict):
    n = len(t)
    if set(incoming) == {n-3,n-2} or arc(t, 1, n-2):
        return
    name, qs = words(n, incoming, arc(t, 2, n-2))
    if name == 'two_second_ports' and n < 8:
        return
    now = t
    flags = []
    for q in qs:
        need(path(now, q), ('directed ledger', n, name, incoming, q, t))
        flags.append(blue(now, q))
        now = partner(now, q)
    need(now == partner(t, tuple(range(n))), 'exact labelled target')
    if strict or name in {'paired_last_second','movable_last_cut'}:
        need(not any(flags), ('claimed three-black theorem', n, name, incoming, flags, t))
    stats[(n, name, 'blue' if any(flags) else 'all_black')] += 1


def cycles(p, q):
    n = len(p)
    rp, rq = dict(zip(p, p[::-1])), dict(zip(q, q[::-1]))
    g, seen, answer = [rq[rp[i]] for i in range(n)], set(), []
    for i in range(n):
        if i in seen:
            continue
        orbit, j = [], i
        while j not in seen:
            seen.add(j)
            orbit.append(j)
            j = g[j]
        answer.append(set(orbit))
    return answer


def main():
    stats, rng = Counter(), random.Random(20261020)
    small = {}
    for n in range(6, 10):
        count = 0
        for bits in range(1 << len(middle_orbits(n-2)[1])):
            for incoming in combinations(range(2, n-1), 2):
                t = source(n, incoming, bits)
                count += 1
                replay(t, incoming, stats, strict=n >= 11)
        small[n] = count
    larger = Counter()
    for n in [11, 12, 13, 14, 15, 16, 19, 20, 25, 32, 41, 60]:
        choices = list(combinations(range(2, n-2), 2))
        for incoming in choices[:20] + [(2,3),(2,n-4),(4,5),(4,n-3)] + [(h,n-2) for h in range(2,n-3)]:
            for _ in range(12):
                bits = rng.getrandbits(len(middle_orbits(n-2)[1]))
                t = source(n, incoming, bits)
                replay(t, incoming, stats, strict=True)
                larger[n] += 1
    for n in range(8, 201):
        p, w = tuple(range(n)), tuple(range(2, n-2))
        q0 = (n-2,)+w+(1,n-1,0)
        cf = (0,n-1,1)+w+(n-2,)
        relocated = (0,)+w+(1,n-2,n-1)
        for q in [q0, cf]:
            c = cycles(p, q)
            need(sorted(map(len,c)) == ([n] if n%2 else [n//2]*2), 'full orbit sizes')
        c = cycles(p, relocated)
        need(sorted(map(len,c)) == ([1,1,n-2] if n%2 else [1,1,(n-2)//2,(n-2)//2]), 'middle orbit sizes')
    degree_free = []
    for n in range(6,10):
        count, degree_hist = 0, Counter()
        total = 1 << len(middle_orbits(n)[1])
        for bits in range(total):
            t = full_blue(n,bits)
            if not (arc(t,n-1,0) and arc(t,n-2,1) and arc(t,n-1,1)
                    and arc(t,0,n-3) and t[-1].bit_count()>=3):
                continue
            k = next(v for v in range(2,n-1) if arc(t,n-1,v))
            l,r = tuple(range(1,k)),tuple(range(k,n-2))
            qs = ((0,)+l+(n-1,)+r+(n-2,), l+(n-1,)+r+(0,n-2), l+r+(0,n-2,n-1))
            now = t
            for q in qs:
                need(path(now,q) and not blue(now,q), 'degree-free three-black corollary')
                now = partner(now,q)
            need(now == partner(t,tuple(range(n))), 'degree-free exact target')
            count += 1
            degree_hist[t[-1].bit_count()] += 1
        degree_free.append({'n':n,'all_standard_blue_presentations':total,'qualifying_sources':count,'last_degree_counts':dict(degree_hist)})
    receipt = {
        'status': 'FINITE-EXACT independent literal controls, separate from the symbolic proof',
        'universe': 'All standard-blue endpoint-(n-4,3) presentations n6..9; no strong-partner or mixedness filters; proof guards applied only before ledger replay',
        'complete_small_counts': small,
        'larger_deterministic_generated_sources': dict(larger),
        'ledger_counts': [{'n':n,'branch':branch,'outcome':outcome,'count':count} for (n,branch,outcome),count in sorted(stats.items())],
        'permutation_orders_checked': [8,200],
        'degree_free_corollary_complete_controls': degree_free,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'reproduce': 'python 05-knowledge/results/endpoint_three_repairs_20261020h/theory_verify.py; repeat with python -O',
    }
    Path(__file__).with_suffix('.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'complete_small_counts':small, 'larger_generated':sum(larger.values()),'replayed':sum(stats.values()),'degree_free_qualifying':sum(r['qualifying_sources'] for r in degree_free),'status':'PASS'}))


if __name__ == '__main__':
    main()
