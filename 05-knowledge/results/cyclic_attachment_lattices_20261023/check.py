"""Declared exhaustive and sampled controls; all checks survive Python -O."""
import argparse
from collections import Counter
from itertools import combinations, permutations, product
import json
from pathlib import Path
import random
from math import gcd
from lattice import (block_record, charge_spectrum, contains, encode_relation,
                     in_solution, necklace_relation, rows_for_edge, shifts, solve)

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def matrix(mask, m, n):
    return [[(mask >> (i*n+j)) & 1 for j in range(n)] for i in range(m)]


def translate(A, a, b):
    m, n = len(A), len(A[0])
    return [[A[(i-a) % m][(j-b) % n] for j in range(n)] for i in range(m)]


def matching(p):
    return [[int(j == p[i]) for j in range(len(p))] for i in range(len(p))]


def exhaustive_blocks():
    hist = Counter()
    for mask in range(512):
        A = matrix(mask, 3, 3)
        for a, b in product(range(3), repeat=2):
            B = translate(A, a, b)
            R = shifts(A, B)
            record = encode_relation(3, 3, R)
            need(record is not None, 'lost translated block')
            need(len(R) == 9//(record['h']*record['v']), 'lattice index')
            need(all(contains(record, u, v) == ((u, v) in R)
                     for u, v in product(range(3), repeat=2)), 'block encoding')
        H = shifts(A, A)
        hist[len(H)] += 1
        B = [r[:] for r in A]
        B[0][0] ^= 1
        need(block_record(A, B) is None, 'changed edge count must reject')
    regular = Counter()
    choices = list(combinations(range(5), 2))
    for rows in product(choices, repeat=5):
        if any(sum(j in row for row in rows) != 2 for j in range(5)):
            continue
        A = [[int(j in rows[i]) for j in range(5)] for i in range(5)]
        H = shifts(A, A)
        rec = encode_relation(5, 5, H)
        regular[len(H)] += 1
        need(len(H) in (1, 5), 'regular prime block has nonfunctional kernel')
        need(len({a for a, b in H}) == len(H) == len({b for a, b in H}),
             'regular prime block is not functional')
        need(rec is not None, 'self comparison lost')
    necklaces = {}
    for n in (3, 5, 7):
        histn = Counter()
        for p in permutations(range(n)):
            d = tuple((p[(i+1) % n]-p[i]) % n for i in range(n))
            word = min(d[a:]+d[:a] for a in range(n))
            q = [0]
            for t in word[:-1]:
                q.append((q[-1]+t) % n)
            R = necklace_relation(p, q)
            need(R == shifts(matching(p), matching(q)), 'necklace phase decoder')
            need(len(R) in (1, n), 'prime matching relation size')
            histn[len(R)] += 1
        necklaces[n] = dict(histn)
    return dict(all_binary_3_by_3=512, translated_comparisons=4608,
                stabilizer_histogram=dict(hist), regular_5_by_5_degree2=dict(regular),
                owned_necklace_relation_histograms=necklaces)


def brute(moduli, rows):
    return [a for a in product(*(range(m) for m in moduli))
            if all((sum(c*x for c, x in zip(coeff, a))-rhs) % q == 0
                   for coeff, rhs, q in rows)]


def check_solution(moduli, rows, points=None):
    points = brute(moduli, rows) if points is None else points
    cert = solve(moduli, rows)
    need((cert is None) == (not points), 'feasibility mismatch')
    if cert is not None:
        need(cert['count'] == len(points), 'exact solution count')
        decoded = [x for x in product(*(range(m) for m in moduli)) if in_solution(cert, x)]
        need(decoded == points, 'whole affine solution lattice mismatch')
        need(tuple(cert['offset']) in points, 'actual witness missing')
    return cert, points


def triangle_bank():
    # Every coset of every subgroup of F3^2, plus the empty relation.
    full = tuple(product(range(3), repeat=2))
    sets = {full}
    for a, b in full:
        sets.add(((a, b),))
        for u, v in full[1:]:
            sets.add(tuple(sorted({((a+t*u) % 3, (b+t*v) % 3) for t in range(3)})))
    return [()] + sorted(sets)


def exhaustive_global():
    bank = triangle_bank()
    need(len(bank) == 23, 'affine relation bank coverage')
    universe = list(product(range(3), repeat=2))
    accepted = 0
    for mask in range(512):
        relation = [p for i, p in enumerate(universe) if (mask >> i) & 1]
        expected = tuple(relation) in bank
        try:
            encode_relation(3, 3, relation)
            actual = True
        except ValueError:
            actual = False
        need(actual == expected, 'arbitrary relation misclassified as one coset')
        accepted += actual
    need(accepted == 23, 'coset validation coverage')
    records = [encode_relation(3, 3, R) for R in bank]
    hist = Counter()
    for inds in product(range(len(bank)), repeat=3):
        rows = []
        for idx, (i, j) in zip(inds, ((0, 1), (1, 2), (0, 2))):
            rows.extend(rows_for_edge(records[idx], i, j, 3))
        direct = [x for x in product(range(3), repeat=3)
                  if (x[0], x[1]) in bank[inds[0]] and
                     (x[1], x[2]) in bank[inds[1]] and
                     (x[0], x[2]) in bank[inds[2]]]
        cert, pts = check_solution((3, 3, 3), rows, direct)
        hist[len(pts)] += 1
    # Composite nonfunctional true tournament blocks, with named cyclic cells.
    records9 = []
    A = [[int((j-i) % 3 == 0) for j in range(9)] for i in range(9)]
    for a, b in product(range(9), repeat=2):
        rec = block_record(A, translate(A, a, b))
        need(81//(rec['h']*rec['v']) == 27, 'three-lift relation size')
        records9.append(rec)
    rows0 = rows_for_edge(records9[0], 0, 1, 3)+rows_for_edge(records9[0], 1, 2, 3)
    flat = rows0+rows_for_edge(records9[0], 0, 2, 3)
    twisted = rows0+rows_for_edge(records9[1], 0, 2, 3)
    flat_cert, flat_pts = check_solution((9, 9, 9), flat)
    twisted_cert, _ = check_solution((9, 9, 9), twisted)
    need(flat_cert['count'] == 81 and twisted_cert is None, 'native composite triangle')
    # Deterministic mixed-modulus systems from arbitrary periodic linear rows.
    rng = random.Random(2301026)
    mixed = Counter()
    moduli = (3, 5, 9, 15)
    for _ in range(96):
        rows = []
        for _ in range(rng.randrange(1, 9)):
            q = rng.choice((3, 5, 9, 15, 45))
            coeff = [rng.randrange(q)*q//gcd(q, m) for m in moduli]
            rows.append((coeff, rng.randrange(q), q))
        cert, pts = check_solution(moduli, rows)
        mixed[len(pts)] += 1
    return dict(affine_F3_relations=len(bank), triangle_systems=len(bank)**3,
                arbitrary_binary_relations_validated=512,
                triangle_solution_histogram=dict(hist), composite_flat=flat_cert,
                composite_twist_feasible=twisted_cert is not None,
                translated_9_by_9_nonfunctional_relations=len(records9),
                mixed_modulus_systems=96, mixed_solution_histogram=dict(mixed))


def charges():
    rng = random.Random(2301027)
    tested = 0
    constant = 0
    for moduli in ((3, 3, 3), (3, 5, 9), (5, 5, 5), (9, 9, 3)):
        for _ in range(80):
            witness = tuple(rng.randrange(m) for m in moduli)
            rows = []
            for _ in range(rng.randrange(5)):
                q = rng.choice((3, 5, 9, 15))
                coeff = [rng.randrange(q)*q//gcd(q, m) for m in moduli]
                rhs = sum(c*t for c, t in zip(coeff, witness)) % q
                rows.append((coeff, rhs, q))
            cert, points = check_solution(moduli, rows)
            for weights in product((-1, 0, 1), repeat=3):
                spec = charge_spectrum(cert, weights)
                values = [sum(a*w for a, w in zip(t, weights)) for t in points]
                need(spec['constant'] == (min(values) == max(values)), 'exact charge constancy')
                from fractions import Fraction
                need(Fraction(spec['mean']) == Fraction(sum(values), len(values)), 'charge mean')
                tested += 1
                constant += spec['constant']
    # Testing the two generator lines alone misses mixed carries.
    cert = solve((3, 3, 3), [([1, -1, -1], 0, 3)])
    w = (1, -1, -1)
    u, v = (1, 1, 0), (1, 0, 1)
    for g in (u, v):
        need(all(sum(wi*((t*gi) % 3) for wi, gi in zip(w, g)) == 0
                 for t in range(3)), 'generator-line hostile setup')
    mixed = tuple((2*a+b) % 3 for a, b in zip(u, v))
    need(mixed == (0, 2, 1) and sum(wi*xi for wi, xi in zip(w, mixed)) == -3,
         'mixed carry hostile consequence')
    carry = charge_spectrum(cert, w)
    need(not carry['constant'], 'carry hidden by Fourier certificate')
    # A genuinely nonfunctional family with constant observed difference.
    paired = solve((9, 9, 9), [([1, -1, 0], 0, 9), ([0, 1, -1], 0, 3)])
    positive = charge_spectrum(paired, (1, -1, 0))
    need(paired['count'] == 27 and positive['constant'] and positive['mean'] == '0',
         'nonfunctional positive observation control')
    # Six actual recolour-position owners, retained in both leaf and debt ledgers.
    rows = []
    for i in range(6):
        row = [0]*12
        row[i], row[i+6] = 1, -1
        rows.append((row, 0, 2))
    same = solve((2,)*12, rows)
    separated = solve((2,)*12, [])
    balance = charge_spectrum(same, (-1,)*6+(1,)*6)
    false_balance = charge_spectrum(separated, (-1,)*6+(1,)*6)
    need(same['count'] == 64 and balance['constant'], 'shared zero-slot charges')
    need(separated['count'] == 4096 and not false_balance['constant'], 'erased charge owners')
    return dict(weighted_tests=tested, constant_tests=constant,
                zero_slot_shared=balance, zero_slot_erased=false_balance,
                mixed_carry_hostile=carry, mixed_carry_point=mixed,
                nonfunctional_positive_control=positive)


def large_family():
    # Alternating exact and modulo-nine attachments, on64 actual cyclic cells.
    n = 64
    moduli = (27,)*n
    equations = []
    for i in range(n):
        j = (i+1) % n
        row = [0]*n
        row[i], row[j] = -1, 1
        equations.append((row, 0, 27 if i % 2 == 0 else 9))
    cert = solve(moduli, equations)
    need(cert['count'] == 9*3**32, 'large-family exact count')
    need(all(0 <= b <= moduli[i] for i, row in enumerate(cert['basis']) for b in row), 'Hermite certificate size')
    for coeff, rhs, q in equations:
        need((sum(c*r for c, r in zip(coeff, cert['offset']))-rhs) % q == 0, 'large source witness')
        need(all(sum(coeff[i]*cert['basis'][i][j] for i in range(n)) % q == 0
                 for j in range(n)), 'large homogeneous directions')
    spectrum = charge_spectrum(cert, (1, -1)*32)
    need(spectrum['constant'] and spectrum['mean'] == '0', 'large shared-pair cancellation')
    pinned = list(equations)
    row = [0]*n
    row[0] = 1
    pinned.append((row, 0, 27))
    one_owner = solve(moduli, pinned)
    need(one_owner['count'] == 3**31, 'one actual owner')
    for i in range(2, n, 2):
        row = [0]*n
        row[i] = 1
        pinned.append((row, 0, 27))
    all_owners = solve(moduli, pinned)
    need(all_owners['count'] == 1, 'full pair-owner base')
    return dict(cells=n, cell_modulus=27, actual_source_vertices=1728,
                compatible_frames=cert['count'], independent_pair_lifts=32,
                one_actual_owner_frames=one_owner['count'], all_32_owners_frames=1,
                lattice_entries=n*n, max_lattice_entry=max(map(max, cert['basis'])),
                spectrum=spectrum)


if __name__ == '__main__':
    pa = argparse.ArgumentParser()
    pa.add_argument('--part', choices=('blocks', 'global', 'charges', 'large'), required=True,
                    help='select one independently frozen receipt')
    pa.add_argument('--write', action='store_true')
    args = pa.parse_args()
    output = {'status': 'FINITE-EXACT controls; proof scopes in README.md'}
    for name, fn in (('blocks', exhaustive_blocks), ('global', exhaustive_global), ('charges', charges), ('large', large_family)):
        if args.part == name:
            output[name] = fn()
    output['checks'] = CHECKS
    output = json.loads(json.dumps(output))
    path = HERE/('receipt-'+args.part+'.json')
    if args.write:
        path.write_text(json.dumps(output, indent=2)+'\n')
    else:
        need(json.loads(path.read_text()) == output, 'frozen receipt mismatch')
    print(json.dumps(output, indent=2))
