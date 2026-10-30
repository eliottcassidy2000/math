"""Exact ramified affine cycles, native attachments and path complements."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import random
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent
GATES = 0


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


native = load('native_parent', RESULTS/'regular_parent_symmetry_20261029/verify.py')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def mul(x, y, q, a):
    u, v = x
    w, z = y
    return ((u*w+v*z) % q, (u*z+v*w+a*v*z) % q)


def clock(p, k, a):
    q = p**k
    g = ((a*a+1) % q, (a*(a*a+2)) % q)
    powers = [(1, 0)]
    for _ in range(q-1):
        powers.append(mul(powers[-1], g, q, a))
    need(mul(powers[-1], g, q, a) == (1, 0) and len(set(powers)) == q, 'exact clock order')
    return q, powers


def cycle_type(perm):
    seen = bytearray(len(perm))
    counts = Counter()
    for i in range(len(perm)):
        if seen[i]:
            continue
        v, size = i, 0
        while not seen[v]:
            seen[v] = 1
            size += 1
            v = perm[v]
        need(v == i, 'permutation closes at cycle source')
        counts[size] += 1
    return tuple(sorted(counts.items()))


def formula_inventory(p, k):
    n = p**(2*k)
    out = Counter({((1, n),): 1})
    for m in range(1, k+1):
        r = 2*k-2*m+1
        cyc = {1: p**r}
        for j in range(1, m):
            cyc[p**j] = (p**(r+2*j)-p**(r+2*j-2))//p**j
        cyc[p**m] = (n-p**(2*k-1))//p**m
        need(sum(cyc.values()) == 2*p**(2*k-m)-p**(2*k-2*m), 'fixed-point type total cycles')
        out[tuple(sorted(cyc.items()))] += (p-1)*p**(3*m-2)
        out[((p**m, n//p**m),)] += (p-1)*(p*p+1)*p**(3*m-3)
    need(sum(out.values()) == p**(3*k), 'group inventory mass')
    need(all(sum(length*count for length, count in typ) == n for typ in out), 'point mass per type')
    return out


def affine_checks():
    reports = []
    for p, k, a in ((5, 1, 1), (5, 2, 1), (13, 1, 3), (17, 1, 8)):
        q, powers = clock(p, k, a)
        points = list(product(range(q), repeat=2))
        hist = Counter()
        for h in powers:
            images = [mul(h, x, q, a) for x in points]
            for u, v in points:
                perm = [((x+u) % q)*q+(y+v) % q for x, y in images]
                hist[cycle_type(perm)] += 1
        expected = formula_inventory(p, k)
        need(hist == expected, 'complete literal affine cycle census')
        reports.append({'p': p, 'k': k, 'a': a, 'vertices': q*q,
                        'all_group_elements': q**3, 'inventory': [(list(t), c) for t, c in sorted(hist.items())]})
    # Larger depths: specified samples, not complete group censuses.
    rng = random.Random(20261030)
    for p, k, a in ((5, 3, 1), (13, 2, 3)):
        q, powers = clock(p, k, a)
        points = list(product(range(q), repeat=2))
        allowed = formula_inventory(p, k)
        samples = [(0, 0, 0)]
        for j in range(k):
            samples.extend([(p**j, 0, 0), (p**j, 1, 0), (0, p**j, 0)])
        samples.extend((rng.randrange(q), rng.randrange(q), rng.randrange(q)) for _ in range(24))
        for t, u, v in samples:
            h = powers[t]
            perm = []
            for x in points:
                w, z = mul(h, x, q, a)
                perm.append(((w+u) % q)*q+(z+v) % q)
            need(cycle_type(perm) in allowed, 'larger-depth literal cycle type')
        reports.append({'p': p, 'k': k, 'a': a, 'vertices': q*q, 'sample_elements': len(samples)})
    recurrence = []
    for p in (5, 13, 17):
        for k in range(1, 4):
            inv = formula_inventory(p, k)
            old = formula_inventory(p, k-1)
            for colors in (2, 3):
                phi = sum(mult*colors**sum(c for _, c in typ) for typ, mult in inv.items())
                rhs = sum(mult*colors**(p*p*sum(c for _, c in typ)) for typ, mult in old.items())
                rhs += (p-1)*p**(3*k-2)*colors**(2*p**k-1)
                rhs += (p-1)*(p*p+1)*p**(3*k-3)*colors**(p**k)
                need(phi == rhs and phi % p**(3*k) == 0, 'all-order Burnside recurrence')
                n = p**(2*k)
                d = n*(p-1)**2//(p*p)
                others = [(typ, mult) for typ, mult in inv.items() if typ != ((1, n),)]
                need(max(sum(c for _, c in typ) for typ, _ in others) == n-d, 'exact minimum index')
                lead = [(typ, mult) for typ, mult in others if sum(c for _, c in typ) == n-d]
                need(len(lead) == 1 and lead[0][1] == p*(p-1), 'leading multiplicity')
                need(all(sum(c for _, c in typ) <= n//p for typ, _ in others if typ != lead[0][0]), 'exponential remainder gap')
                recurrence.append({'p': p, 'k': k, 'Q': colors, 'orbit_count_bits': (phi//p**(3*k)).bit_length(),
                                   'orbit_count_sha256': hashlib.sha256((phi//p**(3*k)).to_bytes((phi.bit_length()+7)//8, 'big')).hexdigest()})
    return {'literal': reports, 'recurrences': recurrence}


def matrix(rows):
    n = len(rows)
    return sp.Matrix([[0 if i == j else (1 if rows[i] >> j & 1 else -1) for j in range(n)] for i in range(n)])


def norm2_small(rows):
    n = len(rows)
    if n <= 3:
        return sp.Integer(n*(n-1)//2)
    a = n*(n-1)//2
    b = sum(int(matrix([sum(((rows[i] >> j) & 1) << t for t, j in enumerate(sub)) for i in sub]).det())
            for sub in combinations(range(n), 4))
    return (a+sp.sqrt(a*a-4*b))/2


def small_checks():
    borders = paths = 0
    comparisons = set()
    for n in range(1, 5):
        for mask in range(1 << (n*(n-1)//2)):
            rows = native.from_mask(n, mask)
            old = norm2_small(rows)
            for cut in range(1 << n):
                val = norm2_small(native.add_vertex(rows, cut))
                comparisons.add((old, val, n))
                borders += 1
    for old, val, n in comparisons:
        need(sp.simplify(val-old).is_nonnegative is True, 'exact spectral interlacing')
        need(sp.simplify(old+n-val).is_nonnegative is True, 'exact squared border budget')
    for n in range(2, 6):
        for mask in range(1 << (n*(n-1)//2)):
            rows = native.from_mask(n, mask)
            def visit(path, used):
                nonlocal paths
                if len(path) == n:
                    target = actual_complement(rows, path)
                    need(actual_complement(target, path) == rows, 'native complement involution')
                    paths += 1
                    return
                for v in range(n):
                    if not (used >> v & 1) and (not path or rows[path[-1]] >> v & 1):
                        visit(path+[v], used | 1 << v)
            visit([], 0)
    need(borders == 1098 and paths == 7886, 'exhaustive small universes')
    prime_cases = []
    for name, rows in (('C3', native.carousel(3)), ('R5', native.carousel(5)), ('P7', native.paley7())):
        n = len(rows)
        group = native.iso_maps(rows, rows)
        bad = {0, (1 << n)-1} | set(rows) | {r | 1 << v for v, r in enumerate(rows)}
        need(len(bad) == 2+2*n, 'prime-extension forbidden masks')
        prime = rigid_prime = self_converse_forbidden = 0
        for cut in range(1 << n):
            child = native.add_vertex(rows, cut)
            yes = is_prime(child)
            need(yes == (cut not in bad), 'native all-module prime extension test')
            st = sum(native.permute_mask(cut, g) == cut for g in group)
            need(len(native.iso_maps(child, child)) == st, 'full unmarked child group')
            reverse = [((1 << (n+1))-1) ^ (1 << i) ^ r for i, r in enumerate(child)]
            self_converse = bool(native.iso_maps(child, reverse, stop=True))
            need(not (yes and self_converse), 'prime regular-parent child is chiral')
            self_converse_forbidden += self_converse
            prime += yes
            rigid_prime += yes and st == 1
        need(self_converse_forbidden == 2*n, 'self-converse clone hostile in these controls')
        prime_cases.append({'parent': name, 'prime_cuts': prime, 'rigid_prime_cuts': rigid_prime,
                            'self_converse_forbidden_cuts': self_converse_forbidden})
    return {'all_border_cases': borders, 'all_legal_complement_cases': paths, 'exact_radical_cases': len(comparisons),
            'prime_extension_controls': prime_cases}


def is_prime(rows):
    n = len(rows)
    for mask in range(1, (1 << n)-1):
        if mask.bit_count() < 2:
            continue
        if all((mask >> v & 1) or (rows[v] & mask) in (0, mask) for v in range(n)):
            return False
    return True


def actual_complement(rows, path):
    n = len(rows)
    need(sorted(path) == list(range(n)), 'path uses every named vertex')
    need(all(rows[u] >> v & 1 for u, v in zip(path, path[1:])), 'legal Hamilton path guard')
    target = [((1 << n)-1) ^ (1 << i) ^ row for i, row in enumerate(rows)]
    for u, v in zip(path, path[1:]):
        target[u] ^= 1 << v
        target[v] ^= 1 << u
    return target


def hamilton(rows, order):
    path = []
    for v in order:
        pos = next((i for i, u in enumerate(path) if rows[v] >> u & 1), len(path))
        path.insert(pos, v)
    return path


def refinement(rows, parts):
    rounds = [list(map(len, parts))]
    while True:
        masks = [sum(1 << v for v in cell) for cell in parts]
        nxt = []
        for cell in parts:
            groups = {}
            for v in cell:
                key = tuple((rows[v] & mask).bit_count() for mask in masks)
                groups.setdefault(key, []).append(v)
            nxt.extend(groups[key] for key in sorted(groups))
        if len(nxt) == len(parts):
            return nxt, [Counter(r) for r in rounds]
        parts = nxt
        rounds.append(list(map(len, parts)))


def large_native():
    p, q, a = 5, 25, 1
    points = list(product(range(q), repeat=2))
    left = set(points)-{(0, 0)}
    orbits, owner = [], {}
    while left:
        v = first = min(left)
        orb = []
        while v not in orb:
            orb.append(v)
            v = mul((2, 3), v, q, a)
        need(v == first, 'literal connection orbit')
        i = len(orbits)
        for x, y in orb:
            owner[x, y] = (i, 1)
            owner[-x % q, -y % q] = (i, -1)
        left -= set(orb) | {(-x % q, -y % q) for x, y in orb}
        orbits.append(orb)
    source = json.loads((RESULTS/'source_cutoffs_balanced_lifts_20261029/lifts.json').read_text())['level25']
    eps = source['signs']
    need(len(eps) == len(orbits) == 24, 'pinned 625 source and orientation gauge')
    connection = {v for v, (i, s) in owner.items() if eps[i]*s == 1}
    rows = [sum(1 << (((x+u) % q)*q+(y+v) % q) for u, v in connection) for x, y in points]
    n = len(rows)
    need(all(row.bit_count() == 312 for row in rows), 'regular actual parent')
    allbits = (1 << n)-1
    # Independent all-root-pair module closures; translations cover every pair.
    closure_rounds = Counter()
    for v in range(1, n):
        module = 1 | 1 << v
        rounds = 0
        while module != allbits:
            forced = sum(1 << w for w in range(n) if not (module >> w & 1) and (rows[w] & module) not in (0, module))
            need(forced != 0, 'no proper module contains the named pair')
            module |= forced
            rounds += 1
        closure_rounds[rounds] += 1
    rng = random.Random(31052026)
    cut = rng.getrandbits(n)
    need(cut not in ({0, allbits} | set(rows) | {r | 1 << i for i, r in enumerate(rows)}), 'prime extension guard')
    child = native.add_vertex(rows, cut)
    need(native.regular_cards(child) == [n], 'unique intrinsic regular-deletion root')
    cells, rounds = refinement(child, [[n], [v for v in range(n) if cut >> v & 1], [v for v in range(n) if not (cut >> v & 1)]])
    need(len(cells) == n+1, 'native full child rigidity by canonical refinement')
    # Fourier enclosure replays every character, not just orbit representatives.
    pins = load('pinned_sines', RESULTS/'pinned_discrepancy_reset_tails_20261029/check_pins.py')
    sines = pins.sine_centres(q)
    peaks = []
    for u, v in points:
        peaks.append(abs(2*sum(sines[(u*x+v*y) % q] for x, y in connection)))
    peak = max(peaks)
    lo, hi = F(peak-2*(n-1), pins.Q), F(peak+2*(n-1), pins.Q)
    need(lo > 0 and lo*lo > p*n, 'literal inherited lower scale')
    need([str(lo), str(hi)] == source['best_norm_interval'], 'independent full Fourier replay')
    # cos x >=1-x^2/2, sin x <=x, pi<22/7 give a rational transitive norm floor.
    transitive_lower = F(7*n, 11)-F(11, 14*n)
    distance_lower = (transitive_lower-hi)/4
    certified_moves = -((-distance_lower.numerator)//distance_lower.denominator)
    need(certified_moves == 81, '625 native complement-distance certificate')
    out = {'parent_vertices': n, 'child_vertices': n+1, 'source_signs': eps,
           'all_root_pair_closures': n-1, 'closure_rounds': dict(closure_rounds),
           'cut_hex': hex(cut), 'rigidity_refinement_histograms': rounds,
           'discrete_vertex_order': [cell[0] for cell in cells],
           'parent_norm_interval': [str(lo), str(hi)], 'child_squared_norm_upper': str(hi*hi+n),
           'parent_to_transitive_complements_at_least': certified_moves,
           'distance_rational_lower': str(distance_lower),
           'parent_relation_sha256': digest(list(map(hex, rows))), 'child_relation_sha256': digest(list(map(hex, child)))}
    current, history = child, []
    for step in range(3):
        ordering = list(range(n+1))
        rng.shuffle(ordering)
        path = hamilton(current, ordering)
        nxt = actual_complement(current, path)
        # K_next+K_current has exactly the aligned two-sided path support.
        edges = {frozenset((u, v)) for u, v in zip(path, path[1:])}
        for i, j in combinations(range(n+1), 2):
            need(((current[i] >> j) & 1) == ((nxt[i] >> j) & 1) if frozenset((i, j)) in edges
                 else ((current[i] >> j) & 1) != ((nxt[i] >> j) & 1), 'full actual complement relation')
        need(actual_complement(nxt, path) == current, 'large involution with common path frame')
        history.append({'path': path, 'target_relation_sha256': digest(list(map(hex, nxt)))})
        current = nxt
    out['actual_complements'] = history
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=HERE/'certificate.json')
    args = parser.parse_args()
    result = {'status': 'FINITE-EXACT controls; all-size proofs separate', 'affine': affine_checks(),
              'small_native': small_checks(), 'large_native': large_native()}
    before = native.GATES
    result['hostile_81_27'] = native.hostile_checks()
    result['inherited_hostile_gates'] = native.GATES-before
    result['gates'] = GATES
    text = json.dumps(result, indent=2, sort_keys=True)+'\n'
    args.output.write_bytes(text.encode())
    print(json.dumps({'gates': GATES, 'receipt_sha256': hashlib.sha256(text.encode()).hexdigest(), 'output': str(args.output)}))


if __name__ == '__main__':
    main()
