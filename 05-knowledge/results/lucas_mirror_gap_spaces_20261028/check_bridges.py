"""Exact controls for Lucas mirrors, biplane frames, and threshold spaces.

Run normally and with -O. No floating arithmetic or external packages.
The finite checks support the accompanying proofs; they are not LRC tests.
"""

from collections import Counter, deque
from fractions import Fraction as Q
from itertools import combinations, product
from math import factorial
import json
from pathlib import Path

CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def words(n):
    return [b for b in range(1 << n)
            if all(not (b >> i & 1 and b >> ((i + 1) % n) & 1)
                   for i in range(n))]


def reader(b, n, p, a):
    return sum((b >> i & 1) * pow(a, i, p) for i in range(n)) % p


def mirror(b, n):
    return sum((b >> i & 1) << ((-i) % n) for i in range(n))


def compose(a, b):
    return tuple(a[b[x]] for x in range(len(a)))


def cycles(g):
    left, result = set(range(len(g))), []
    while left:
        x, c = min(left), []
        while x in left:
            left.remove(x)
            c.append(x)
            x = g[x]
        result.append(c)
    return result


def closure(generators):
    identity = tuple(range(len(generators[0])))
    seen, todo = {identity}, deque([identity])
    while todo:
        g = todo.popleft()
        for s in generators:
            h = compose(s, g)
            if h not in seen:
                seen.add(h)
                todo.append(h)
    return seen


# Z[zeta_5], in basis 1,z,z^2,z^3; z^4=-1-z-z^2-z^3.
def reduce5(a):
    a = list(a) + [0] * max(0, 4 - len(a))
    for k in range(len(a) - 1, 3, -1):
        c = a[k]
        for j in range(4):
            a[k - 4 + j] -= c
    return tuple(a[:4])


def mul5(a, b):
    c = [0] * 7
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return reduce5(c)


def conjugate5(a):
    c = [0] * 5
    for i, x in enumerate(a):
        c[-i % 5] += x
    return reduce5(c)


def square_distance(a, b):
    d = tuple(x - y for x, y in zip(a, b))
    return mul5(d, conjugate5(d))


def main():
    W = words(5)
    label = {b: reader(b, 5, 11, 4) for b in W}
    points = {b: reduce5([(b >> i) & 1 for i in range(5)]) for b in W}
    require(len(set(label.values())) == len(set(points.values())) == 11, 'faithful labels')
    edges = {frozenset((a, b)) for a, b in combinations(W, 2)
             if (a ^ b).bit_count() == 1}
    require(len(edges) == 15, 'Lucas edges')
    for a, b in combinations(W, 2):
        require((square_distance(points[a], points[b]) == (1, 0, 0, 0))
                == (frozenset((a, b)) in edges), 'exact unit relation')
        # Independent route: a degree<=4 polynomial vanishes at zeta_5
        # iff all five coefficients are equal (a multiple of Phi_5).
        eps = [(a >> i & 1)-(b >> i & 1) for i in range(5)]
        coeff = [sum(eps[i]*eps[(i-k) % 5] for i in range(5)) for k in range(5)]
        coeff[0] -= 1
        require((len(set(coeff)) == 1) == (frozenset((a, b)) in edges),
                'independent cyclic-autocorrelation unit check')
    degrees = Counter(sum(b in e for e in edges) for b in W)
    require(degrees == {5: 1, 3: 5, 2: 5}, 'degree strata')
    nonzero_edges = [e for e in edges if 0 not in e]
    require(len(nonzero_edges) == 10 and
            all(sum(b in e for e in nonzero_edges) == 2 for b in W if b), 'outer ten-cycle')
    H = {1, 3, 4, 5, 9}
    arcs = {(label[a], label[b]) for a in W for b in W
            if (a ^ b).bit_count() == 1 and a.bit_count() < b.bit_count()}
    paley = {(x, y) for x in range(11) for y in range(11) if (y-x) % 11 in H}
    translated = {((x+t) % 11, (y+t) % 11) for x, y in arcs for t in range(11)}
    require(arcs <= paley and translated == paley and len(paley) == 55, 'translation closure')
    R = [0] * 11
    for b in W:
        R[label[b]] = label[mirror(b, 5)]
    R = tuple(R)
    require(all(R[x] == (3*pow(x, 9, 11)-2*pow(x, 4, 11)) % 11
                for x in range(11)), 'reflection polynomial')
    require([x for x in range(11) if R[x] == x] == [0, 1, 7], 'three fixed labels')
    require({(R[x], R[y]) for x, y in arcs} == arcs, 'mirror preserves guarded arcs')
    require(sum((R[x], R[y]) in paley for x, y in paley) == 25, 'mixed Paley signs')

    I = tuple(range(11))
    T = tuple((x+1) % 11 for x in I)
    M = tuple(4*x % 11 for x in I)
    blocks = [frozenset((x+h) % 11 for h in H) for x in I]
    lookup = {b: i for i, b in enumerate(blocks)}
    require(all(len(a & b) == 2 for a, b in combinations(blocks, 2)), 'biplane intersections')

    def on_blocks(g):
        return tuple(lookup[frozenset(g[y] for y in b)] for b in blocks)

    def on_frame(g, f):
        bg, out = on_blocks(g), [None] * 11
        for x in I:
            out[g[x]] = bg[f[x]]
        return tuple(out)

    G = closure((T, R))
    require(len(G) == 660 and M in G, 'independent BFS group order')
    diagonal = [g for g in G if on_blocks(g) == g]
    require(len(diagonal) == 55 and set(diagonal) == closure((T, M)), 'matching stabilizer')
    f0 = on_frame(R, I)
    frames = {'inf': I}
    f = f0
    for b in I:
        frames[b] = f
        f = on_frame(T, f)
    require(len(set(frames.values())) == 12, 'twelve distinct frames')
    require({on_frame(g, I) for g in G} == set(frames.values()), 'whole matching orbit')
    for b, f in frames.items():
        require(on_frame(T, f) == frames['inf' if b == 'inf' else (b+1) % 11], 'projective T')
        require(on_frame(M, f) == frames['inf' if b == 'inf' else 4*b % 11], 'projective M')
        rb = 0 if b == 'inf' else 'inf' if b == 0 else 6*pow(b, -1, 11) % 11
        require(on_frame(R, f) == frames[rb], 'projective R')
        if b != 'inf':
            require([x for x in I if f[x] == x] == [b], 'unique fixed point gives faithful frame action')
        require(all((y in blocks[f[x]]) != (x in blocks[f[y]])
                    for x, y in combinations(I, 2)), 'transported matching decodes a tournament')
        require(all(x not in blocks[f[x]] for x in I), 'no loops in decoded tournament')
    scores = Counter(sum(y in blocks[f[x]] for x, y in arcs) for f in frames.values())
    require(scores == {15: 2, 6: 10}, 'Lucas observations leave exactly one bit')
    require(2 in blocks[I[1]] and 2 not in blocks[f0[1]], 'probe 1->2 resolves bit')
    inverse_label = {x: b for b, x in label.items()}
    probe_norm = square_distance(points[inverse_label[1]], points[inverse_label[2]])
    require(probe_norm != (1, 0, 0, 0), 'probe is not a unit edge')
    # phi = 1+z+z^-1 = -z^2-z^3, hence 4-phi = (4,0,1,1).
    require(probe_norm == (4, 0, 1, 1), 'probe norm 4-phi')
    require(factorial(11)//20 == 1995840, 'signed-copy pairing counts')
    require(sum(frozenset((mirror(a, 5), mirror(b, 5))) == frozenset((a, b))
                for a, b in map(tuple, edges)) == 1, 'Lucas edge reflection trace')
    k6_vertices = list(range(6))
    k6_mirror = [0, 4, 3, 2, 1, 5]
    require(sum({k6_mirror[a], k6_mirror[b]} == {a, b}
                for a, b in combinations(k6_vertices, 2)) == 3, 'K6 trace hostile')

    # Next odd Lucas prime: a short certificate, independent of group-order software.
    W7 = words(7)
    label7 = {b: reader(b, 7, 29, 24) for b in W7}
    require(len(set(label7.values())) == 29, 'period seven reader')
    R7 = [0] * 29
    for b in W7:
        R7[label7[b]] = label7[mirror(b, 7)]
    require(sum(x == R7[x] for x in range(29)) == 5, 'period seven mirror')
    exponents = [19, 17, 22, 23, 18, 1, 18, 21, 6, 1, 15, 17]
    g = tuple(range(29))
    for a in exponents:
        g = compose(R7, compose(tuple((x+a) % 29 for x in range(29)), g))
    require(sorted(map(len, cycles(g))) == [1, 3, 25], 'short Jordan witness')
    h = tuple(range(29))
    for _ in range(25):
        h = compose(g, h)
    require([c for c in cycles(h) if len(c) > 1] == [[0, 20, 23]], 'single three-cycle')

    # Exact rational probes of the general threshold-space homotopies.
    tested_gaps = 0
    for n in range(3, 9):
        a = [Q(1, n)] * n
        legal = set(words(n))
        for mask in legal - {0}:
            support = [i for i in range(n) if mask >> i & 1]
            k = len(support)
            gap = [Q(n+k, 2*n*k) if i in support else Q(1, 2*n) for i in range(n)]
            require(sum(gap) == 1 and sum((gap[i] >= a[i]) << i for i in range(n)) == mask,
                    'every safe mask has a rational gap source')
            for stretch in (Q(1, 5), Q(1, 2), Q(1)):
                g = [a[i]+stretch*(gap[i]-a[i]) for i in range(n)]
                e = [max(g[i]-a[i], 0) for i in range(n)]
                r = [x/sum(e) for x in e]
                section = [(a[i]+r[i])/2 for i in range(n)]
                for t in (Q(0), Q(1, 7), Q(1, 2), Q(1)):
                    homotopy = [(1-t)*g[i]+t*section[i] for i in range(n)]
                    safe = sum((homotopy[i] >= a[i]) << i for i in range(n))
                    require(sum(homotopy) == 1 and min(homotopy) > 0 and safe in legal-{0},
                            'gap homotopy stays in bad space')
                tested_gaps += 1
    # Boundary equality, nonuniform thresholds, and arbitrary forbidden-pair complexes.
    for n in range(3, 7):
        a = [Q(i+1, n*(n+1)//2) for i in range(n)]
        for numerators in product(range(5), repeat=n):
            if sum(numerators) != 4:
                continue
            g = [Q(x, 4) for x in numerators]
            safe = {i for i in range(n) if g[i] >= a[i]}
            if {0, 1} <= safe:
                continue
            e = [max(g[i]-a[i], 0) for i in range(n)]
            require(sum(e) > 0, 'nonzero positive part including boundary')
            r = [x/sum(e) for x in e]
            section = [(a[i]+r[i])/2 for i in range(n)]
            for t in (Q(0), Q(1, 2), Q(1)):
                h = [(1-t)*g[i]+t*section[i] for i in range(n)]
                require(not (h[0] >= a[0] and h[1] >= a[1]), 'marked-target homotopy')

    signed_counts = {}
    for n in range(3, 18):
        ws = words(n)
        signed = sum((-1)**b.bit_count() for b in ws)
        require(signed == [2, 1, -1, -2, -1, 1][n % 6], 'six-period signed count')
        signed_counts[n] = [len(ws), signed]
    receipt = {
        'status': 'PASS; exact finite controls, not a major-conjecture proof',
        'checks': CHECKS, 'labels': sorted(label.items()), 'reflection': R,
        'group_order_BFS': len(G), 'matching_stabilizer': len(diagonal),
        'frames': {str(k): v for k, v in frames.items()},
        'Lucas_arc_agreement': sorted(scores.items()), 'probe_norm_zeta_basis': probe_norm,
        'period7_Jordan_word_exponents': exponents, 'period7_Jordan_power': 25,
        'period7_three_cycle': [0, 20, 23], 'rational_gap_sources': tested_gaps,
        'word_and_signed_counts_n3_to17': signed_counts,
    }
    Path(__file__).with_name('receipt.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({k: receipt[k] for k in ['status', 'checks', 'group_order_BFS',
          'matching_stabilizer', 'Lucas_arc_agreement', 'rational_gap_sources']}, indent=2))


if __name__ == '__main__':
    main()
