#!/usr/bin/env python3
"""Independent finite checks of pair-state, cancellation-fibre and path-mass formulas.

No external packages. Explicit exceptions retain all checks under python -O.
The complete literal universe is every loopless labelled digraph through order4.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
from math import factorial, gcd, lcm
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def need(test, message):
    if not test:
        raise RuntimeError(message)


def parts(n, lo=1):
    if n == 0:
        yield ()
    for a in range(lo, n + 1):
        for tail in parts(n - a, a):
            yield (a,) + tail


def centralizer(lam):
    out = 1
    for length, count in Counter(lam).items():
        out *= length ** count * factorial(count)
    return out


def cycle_permutation(lam):
    p = []
    start = 0
    for length in lam:
        p.extend(range(start + 1, start + length))
        p.append(start)
        start += length
    return tuple(p)


def pair_orbits_literal(p):
    n = len(p)
    unseen = {(i, j) for i in range(n) for j in range(i + 1, n)}
    preserving, reversing = [], []
    while unseen:
        initial = min(unseen)
        a, b = initial
        length = 0
        while True:
            unseen.discard(tuple(sorted((a, b))))
            a, b = p[a], p[b]
            length += 1
            if tuple(sorted((a, b))) == initial:
                break
        (preserving if (a, b) == initial else reversing).append(length)
    return sorted(preserving), sorted(reversing)


def pair_orbits_formula(lam):
    preserving, reversing = [], []
    for i, a in enumerate(lam):
        preserving.extend([a] * ((a - 1) // 2))
        if a % 2 == 0:
            reversing.append(a // 2)
        for b in lam[:i]:
            preserving.extend([lcm(a, b)] * gcd(a, b))
    return sorted(preserving), sorted(reversing)


def multiply(poly, terms):
    out = defaultdict(int)
    for (a, b, c), coeff in poly.items():
        for (da, db, dc), value in terms:
            out[a + da, b + db, c + dc] += coeff * value
    return dict(out)


def weighted_fixed(lam):
    poly = {(0, 0, 0): 1}
    preserving, reversing = pair_orbits_formula(lam)
    for ell in preserving:
        poly = multiply(poly, [((ell, 0, 0), 1), ((0, ell, 0), 1), ((0, 0, ell), 2)])
    for ell in reversing:
        poly = multiply(poly, [((ell, 0, 0), 1), ((0, ell, 0), 1)])
    return poly


def orbit_polynomial(n):
    out = defaultdict(Fraction)
    for lam in parts(n):
        z = centralizer(lam)
        for degrees, coeff in weighted_fixed(lam).items():
            out[degrees] += Fraction(coeff, z)
    need(all(v.denominator == 1 for v in out.values()), 'Integral orbit coefficients')
    return {k: int(v) for k, v in out.items()}


def palette_count(n, q, t, loops=False):
    out = Fraction(0)
    for lam in parts(n):
        a, b = pair_orbits_formula(lam)
        out += Fraction(q ** len(a) * t ** len(b) * (2 ** len(lam) if loops else 1), centralizer(lam))
    need(out.denominator == 1, 'Integral palette count')
    return int(out)


def diagonal_polynomial(n):
    out = defaultdict(Fraction)
    for lam in parts(n):
        preserving, reversing = pair_orbits_formula(lam)
        current = {0: 1}
        for ell in lam:
            nxt = defaultdict(int)
            for k, coeff in current.items():
                nxt[k] += coeff
                nxt[k + ell] += coeff
            current = nxt
        for k, coeff in current.items():
            out[k] += Fraction(coeff * 4 ** len(preserving) * 2 ** len(reversing), centralizer(lam))
    need(all(v.denominator == 1 for v in out.values()), 'Integral diagonal polynomial')
    return [int(out[k]) for k in range(n + 1)]


def pair_state(rows):
    absent = both = single = 0
    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            total = ((rows[i] >> j) & 1) + ((rows[j] >> i) & 1)
            absent += total == 0
            single += total == 1
            both += total == 2
    return absent, both, single


def relabel(rows, p):
    out = [0] * len(rows)
    for i, row in enumerate(rows):
        for j in range(len(rows)):
            if (row >> j) & 1:
                out[p[i]] |= 1 << p[j]
    return tuple(out)


def canonical(rows, perms):
    return min(relabel(rows, p) for p in perms)


def cancel(rows):
    return tuple(sum(1 << j for j in range(len(rows))
                     if i != j and (rows[i] >> j) & 1 and not ((rows[j] >> i) & 1))
                 for i in range(len(rows)))


def literal_paths(rows):
    return sum(all((rows[p[i]] >> p[i + 1]) & 1 for i in range(len(rows) - 1))
               for p in permutations(range(len(rows))))


def subset_paths(rows):
    n = len(rows)
    if n == 0:
        return 1
    counts = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        counts[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for u in range(n):
            if counts[mask][u]:
                for v in range(n):
                    if not ((mask >> v) & 1) and (rows[u] >> v) & 1:
                        counts[mask | (1 << v)][v] += counts[mask][u]
    return sum(counts[-1])


def marked_formula(n, cancelled=False):
    need(n >= 1, 'Marked formula is stated for n>=1')
    poly = {(0, 0, 0): 1}
    for _ in range(n - 1):
        terms = [((0, 0, 1), 1)]
        if not cancelled:
            terms.append(((0, 1, 0), 1))
        poly = multiply(poly, terms)
    for _ in range(n * (n - 1) // 2 - n + 1):
        poly = multiply(poly, [((1, 0, 0), 1), ((0, 1, 0), 1), ((0, 0, 1), 2)])
    return poly


def cycle_lengths_on(objects, action):
    unseen = set(objects)
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while True:
            need(current in unseen, 'Orbit remains in prescribed finite set')
            unseen.remove(current)
            current = action(current)
            length += 1
            if current == start:
                break
        lengths.append(length)
    return lengths


def fibre_polynomial(rows, perms, loops=False):
    n = len(rows)
    missing = [(i, j) for i in range(n) for j in range(i + 1, n)
               if not ((rows[i] >> j) & 1 or (rows[j] >> i) & 1)]
    auts = [p for p in perms if relabel(rows, p) == rows]
    out = defaultdict(Fraction)
    for p in auts:
        current = {(0, 0): 1}
        for ell in cycle_lengths_on(missing, lambda ij: tuple(sorted((p[ij[0]], p[ij[1]])))):
            nxt = defaultdict(int)
            for (a, b), coeff in current.items():
                nxt[a, b] += coeff
                nxt[a + ell, b] += coeff
            current = nxt
        if loops:
            for ell in cycle_lengths_on(range(n), lambda i: p[i]):
                nxt = defaultdict(int)
                for (a, b), coeff in current.items():
                    nxt[a, b] += coeff
                    nxt[a, b + ell] += coeff
                current = nxt
        for degree, coeff in current.items():
            out[degree] += Fraction(coeff, len(auts))
    need(all(v.denominator == 1 for v in out.values()), 'Integral fibre coefficients')
    return {k: int(v) for k, v in out.items()}, len(auts), missing


def literal_fibre(rows, perms, loops=False):
    poly, _, missing = fibre_polynomial(rows, perms, loops)
    types = {}
    for bits in product((0, 1), repeat=len(missing) + (len(rows) if loops else 0)):
        raw = list(rows)
        for (i, j), bit in zip(missing, bits):
            if bit:
                raw[i] |= 1 << j
                raw[j] |= 1 << i
        if loops:
            for i, bit in enumerate(bits[len(missing):]):
                raw[i] |= bit << i
        key = canonical(raw, perms)
        degree = sum(bits[:len(missing)]), sum(bits[len(missing):])
        need(key not in types or types[key] == degree, 'Decoration degrees preserved')
        types[key] = degree
    hist = dict(Counter(types.values()))
    need(hist == poly, 'Literal versus Burnside cancellation fibre')
    return poly


def main():
    profile_checks = 0
    for n in range(11):
        for lam in parts(n):
            need(pair_orbits_formula(lam) == pair_orbits_literal(cycle_permutation(lam)), 'Pair-orbit lengths and reversal')
            profile_checks += 1
    censuses = []
    fibre_checks = 0
    for n in range(1, 5):
        perms = list(permutations(range(n)))
        arcs = [(i, j) for i in range(n) for j in range(n) if i != j]
        classes = Counter()
        for mask in range(1 << len(arcs)):
            rows = [0] * n
            for k, (i, j) in enumerate(arcs):
                if (mask >> k) & 1:
                    rows[i] |= 1 << j
            classes[canonical(rows, perms)] += 1
        hist = Counter()
        mass = Counter()
        cancelled_mass = Counter()
        fibres = defaultdict(Counter)
        oriented = []
        for rows, orbit_size in classes.items():
            aut = sum(relabel(rows, p) == rows for p in perms)
            need(orbit_size * aut == factorial(n), 'Orbit-stabilizer')
            h, hc = literal_paths(rows), literal_paths(cancel(rows))
            need(h == subset_paths(rows) and hc == subset_paths(cancel(rows)), 'Independent H counts')
            need(h % aut == 0 and hc % aut == 0, 'Free action on ordered Hamilton paths')
            degree = pair_state(rows)
            hist[degree] += 1
            mass[degree] += h // aut
            cancelled_mass[degree] += hc // aut
            fibres[canonical(cancel(rows), perms)][degree[1], 0] += 1
            if degree[1] == 0:
                oriented.append(rows)
        clean = lambda p: {k: v for k, v in p.items() if v}
        need(dict(hist) == orbit_polynomial(n), 'Full weighted orbit polynomial')
        need(clean(mass) == marked_formula(n), 'Ordinary marked path polynomial')
        need(clean(cancelled_mass) == marked_formula(n, True), 'Cancelled marked path polynomial')
        loop_classes = 0
        loop_row = Counter()
        for rows in oriented:
            predicted, _, _ = fibre_polynomial(rows, perms)
            need(dict(fibres[rows]) == predicted, 'Aggregated actual cancellation fibre')
            literal_fibre(rows, perms)
            loop_poly = literal_fibre(rows, perms, True)
            loop_classes += sum(loop_poly.values())
            for (_, loops), coeff in loop_poly.items():
                loop_row[loops] += coeff
            fibre_checks += 2
        need(loop_classes == palette_count(n, 4, 2, True), 'Looped matrix classes from fibre sum')
        need([loop_row[k] for k in range(n + 1)] == diagonal_polynomial(n), 'Diagonal polynomial from literal fibres')
        censuses.append(dict(n=n, labelled_matrices=sum(classes.values()), raw_types=len(classes),
                             semantic_oriented_types=len(oriented), looped_raw_types=loop_classes,
                             ordinary_path_marked_types=sum(mass.values()),
                             cancelled_path_marked_types=sum(cancelled_mass.values())))
    c4 = (2, 4, 8, 1)
    p4 = list(permutations(range(4)))
    c4poly, c4aut, c4missing = fibre_polynomial(c4, p4)
    c4loop = literal_fibre(c4, p4, True)
    c4raw = []
    seen = set()
    for bits in product((0, 1), repeat=len(c4missing)):
        rows = list(c4)
        for (i, j), bit in zip(c4missing, bits):
            if bit:
                rows[i] |= 1 << j
                rows[j] |= 1 << i
        key = canonical(rows, p4)
        if key not in seen:
            seen.add(key)
            c4raw.append(dict(rows=list(key), digons=sum(bits), ordinary_H=literal_paths(key), cancelled_H=literal_paths(cancel(key))))
    sources = [Path(__file__).resolve(), ROOT/'01-canon/theorems/THM-5041-matrix-tiling-pair-dictionary.md',
               ROOT/'01-canon/theorems/THM-4720-k-valued-tournaments-by-burnside-the-parity-dichotomy-of-the-tie-state-and-three-hypertournaments.md',
               ROOT/'05-knowledge/results/session239_sequence_engines_20261110/burnside/PALETTE_PROOF.md',
               ROOT/'05-knowledge/results/session234_digraphs_20261109/PROOFS.md']
    result = dict(status='PASS', scope='Exact pair-state, cancellation-fibre and path-marked identities; finite independent universe through order4; no Hamilton spectrum closure.',
                  pair_orbit_profile_checks=profile_checks, cancellation_fibre_polynomial_checks=fibre_checks,
                  complete_censuses=censuses,
                  sequence_prefixes={name:[palette_count(n,q,t,loops) for n in range(11)]
                                     for name,q,t,loops in [('tournaments_A000568',2,0,False),('oriented_A001174',3,1,False),
                                         ('digraphs_A000273',4,2,False),('graphs_A000088',2,2,False),('looped_digraphs_A000595',4,2,True)]},
                  diagonal_refinement_A384105=[diagonal_polynomial(n) for n in range(7)],
                  weighted_orbit_coefficients={str(n):[dict(missing=k[0],digons=k[1],single=k[2],types=v)
                                                      for k,v in sorted(orbit_polynomial(n).items())] for n in range(1,6)},
                  c4_fibre=dict(effective_rows=list(c4),effective_H=4,automorphism_order=c4aut,
                                missing_pairs=c4missing,labelled_preimages=4,unlabelled_preimages=sum(c4poly.values()),
                                by_digons={str(k[0]):v for k,v in sorted(c4poly.items())},
                                looped_preimage_types=sum(c4loop.values()),raw_types=sorted(c4raw,key=lambda a:a['digons'])),
                  source_sha256={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sources})
    print(json.dumps(result,indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
