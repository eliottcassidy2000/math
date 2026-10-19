"""Exact weighted mod15/mod35 interface; no orbit or ROOT oracle.

Universe: five stated literal source-charge tables; every support graph on
3+3 vertices; full K3,7 rectangle bases; declared native all-one prefixes.
The original source is retained. Marginal orthogonality is not arrival.
"""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys

import root_clock_channels_20261018o as inherited

ROOT = Path(__file__).resolve().parents[2]
CHECKS = 0
SCALES = (4096, 5997, 11002, 16384, 65536)


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def integer(x, minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError("literal integer outside domain")
    return x


def rational(x):
    if type(x) not in (int, F):
        raise ValueError("literal rational required")
    return F(x)


def table(source):
    if type(source) not in (tuple, list) or not source:
        raise ValueError("nonempty rectangular table required")
    if any(type(row) not in (tuple, list) or not row for row in source):
        raise ValueError("nonempty rows required")
    n = len(source[0])
    if any(len(row) != n for row in source):
        raise ValueError("ragged table")
    out = tuple(tuple(rational(x) for x in row) for row in source)
    if any(x < 0 for row in out for x in row):
        raise ValueError("negative source mass")
    return out


def analyse(source):
    M = table(source)
    nr, nc = len(M), len(M[0])
    R = tuple(map(sum, M))
    C = tuple(sum(M[i][j] for i in range(nr)) for j in range(nc))
    T = sum(R)
    if not T:
        return dict(mass=T, rows=R, columns=C, delta=F(0), independent=True,
                    edges=0, vertices=0, components=0, interaction_dimension=0)
    delta = sum(M[i][j]**2/(R[i]*C[j]) for i in range(nr)
                for j in range(nc) if R[i] and C[j]) - 1
    defect = sum((M[i][j]-R[i]*C[j]/T)**2/(R[i]*C[j])
                 for i in range(nr) for j in range(nc) if R[i] and C[j])
    check(delta == defect and delta >= 0, "exact rational Pearson identity")
    independent = all(M[i][j]*T == R[i]*C[j]
                      for i in range(nr) for j in range(nc))
    check((delta == 0) == independent, "zero defect iff product table")
    adjacency = defaultdict(set)
    edges = 0
    for i, j in product(range(nr), range(nc)):
        if M[i][j]:
            adjacency[i].add(nr+j)
            adjacency[nr+j].add(i)
            edges += 1
    unseen = set(adjacency)
    components = 0
    while unseen:
        components += 1
        stack = [unseen.pop()]
        while stack:
            for v in adjacency[stack.pop()]:
                if v in unseen:
                    unseen.remove(v)
                    stack.append(v)
    dim = edges-len(adjacency)+components
    return dict(mass=T, rows=R, columns=C, delta=delta, independent=independent,
                edges=edges, vertices=len(adjacency), components=components,
                interaction_dimension=dim)


def mixed_receipt(source, left, right):
    M = table(source)
    if type(left) not in (list, tuple) or type(right) not in (list, tuple):
        raise ValueError("literal feature vectors required")
    f, g = tuple(map(rational, left)), tuple(map(rational, right))
    if len(f) != len(M) or len(g) != len(M[0]):
        raise ValueError("feature dimensions")
    a = analyse(M)
    R, C = a['rows'], a['columns']
    if sum(x*y for x, y in zip(R, f)) or sum(x*y for x, y in zip(C, g)):
        raise ValueError("features must be conditionally centered")
    cross = sum(M[i][j]*f[i]*g[j] for i in range(len(f)) for j in range(len(g)))
    ef = sum(w*x*x for w, x in zip(R, f))
    eg = sum(w*x*x for w, x in zip(C, g))
    check(cross*cross <= a['delta']*ef*eg, "rational mixed-term certificate")
    return dict(cross=cross, left_energy=ef, right_energy=eg, delta=a['delta'],
                bound_squared=a['delta']*ef*eg)


def qrank(matrix):
    a = [list(map(F, row)) for row in matrix]
    if not a:
        return 0
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        z = a[rank][col]
        a[rank] = [x/z for x in a[rank]]
        for i in range(len(a)):
            if i != rank and a[i][col]:
                z = a[i][col]
                a[i] = [x-z*y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def _projection(values, weights, keys):
    totals = defaultdict(F)
    masses = defaultdict(F)
    for v, w, key in zip(values, weights, keys):
        totals[key] += w*v
        masses[key] += w
    return tuple(totals[key]/masses[key] for key in keys)


def cycle_controls():
    # Independent incidence-matrix rank for every 3 by 3 support.
    for mask in range(512):
        edges = [(i, j) for i, j in product(range(3), repeat=2)
                 if mask >> (3*i+j) & 1]
        M = [[int((i, j) in edges) for j in range(3)] for i in range(3)]
        a = analyse(M)
        incidence = [[int(i == v)-int(3+j == v) for i, j in edges]
                     for v in range(6)]
        rank = qrank(incidence) if edges else 0
        check(len(edges)-rank == a['interaction_dimension'], "cycle-space nullity")
        # All positive support atoms, including disconnected zero-cell cases.
        if edges:
            weights = [F(1)]*len(edges)
            rowkeys, colkeys = zip(*edges)
            product_is_mean = True
            for k in range(len(edges)):
                v = [F(int(j == k)) for j in range(len(edges))]
                result = _projection(_projection(v, weights, colkeys), weights, rowkeys)
                product_is_mean &= all(x == F(1, len(edges)) for x in result)
            check(product_is_mean == a['independent'], "conditional projection iff product")
    basis = []
    for i, j in product(range(1, 3), range(1, 7)):
        row = [0]*21
        for a, b, sign in [(0, 0, 1), (i, j, 1), (0, j, -1), (i, 0, -1)]:
            row[7*a+b] = sign
        check(all(sum(row[7*a+b] for b in range(7)) == 0 for a in range(3)),
              "rectangle row conservation")
        check(all(sum(row[7*a+b] for a in range(3)) == 0 for b in range(7)),
              "rectangle column conservation")
        basis.append(row)
    a = analyse([[1]*7 for _ in range(3)])
    check(qrank(basis) == a['interaction_dimension'] == 12, "K3,7 basis complete")
    check(5*12 == 60 and 105-60 == 45, "five-fibre visible and hidden dimensions")
    return dict(basis=basis, per_fibre_dimension=12, full_105_hidden=60,
                full_105_marginal_visible=45, support_universe=512)


def native_rectangle(length):
    integer(length, 1)
    modulus = 1 << (length+1)
    # One literal positive source per active CRT corner, at fixed mod5=1.
    sources = []
    for i, j in product((0, 1), repeat=2):
        residue = next(x for x in range(105) if (x%3, x%7, x%5) == (i, j, 1))
        k = ((modulus-1-residue)*pow(105, -1, modulus)) % modulus
        n = residue+105*k
        check(n > 1 and n % modulus == modulus-1, "native source CRT")
        orbit = [n]
        for _ in range(length):
            numerator = 3*orbit[-1]+1
            check(numerator%4 == 2, "exact odd valuation one")
            orbit.append(numerator//2)
        check(orbit[-1] == 3**length*(n+1)//2**length-1 > n,
              "actual growing injective branch")
        sources.append(dict(source=n, mod15=n%15, mod35=n%35, sign=(-1)**(i+j),
                            endpoint=orbit[-1], odd_word=[1]*length,
                            shortcut_clock=length, ordinary_clock=2*length))
    for modulus in (15, 35):
        bins = defaultdict(int)
        for row in sources:
            bins[row['source']%modulus] += row['sign']
        check(not any(bins.values()), "both signed marginals vanish")
    check(len({r['endpoint'] for r in sources}) == 4, "no actual terminal cancellation")
    a = analyse([[1, 1], [1, 1]])
    check(a['delta'] == 0 and a['interaction_dimension'] == 1,
          "independent measure with unpaid four-cycle")
    return dict(length=length, rows=sources, source_l1=4, terminal_l1=4,
                source_energy=4, marginal_energy=0, interaction_energy=4, delta=F(0))


def independent_charge(s, n):
    """Separate shell formula: map each residue weight directly to (a,b)."""
    integer(s, 1)
    integer(n, 1)
    if s%2 == 0 or s%8 == 5 or s >= 4*n:
        raise ValueError("not a canonical source")
    # Inherited three-real-row map inverted once, not called in this path.
    inner = ((1, 2), (-2, -1), (1, -1))
    outer = ((0, 1), (-1, -1), (1, 0))
    a = b = 0
    x = s
    while x < 4*n:
        u, v = (inner if x < 2*n else outer)[x%3]
        a += u
        b += v
        x = 4*x+1
    return a, b


def source_table(n):
    integer(n, 1)
    if n not in SCALES:
        raise ValueError("source-table consumer is scoped to the five declared scales")
    M = [[[0]*7 for _ in range(3)] for _ in range(5)]
    first = [[[0]*7 for _ in range(3)] for _ in range(5)]
    records = []
    second = 0
    positive = 0
    for s in range(1, 4*n, 2):
        if s%8 == 5:
            continue
        a, b = inherited.charge(s, n)
        check((a, b) == independent_charge(s, n), "independent source charge")
        mu = a*a-a*b+b*b
        check(mu == inherited.norm((a, b)) and mu in (0, 1, 3), "source norm")
        records.append(f'{s}:{a}:{b}:{mu}')
        if mu:
            positive += 1
        c, i, j = s%5, s%3, s%7
        M[c][i][j] += mu
        first[c][i][j] += mu*s
        second += mu*s*s
    stats = [analyse(t) for t in M]
    # Exact weighted ANOVA for the diagnostic observable F(s)=s.
    # This is a source-coordinate test, not a ROOT amplitude estimate.
    m_energy = x_energy = y_energy = z_energy = cross = within = F(0)
    cross_receipts = []
    for c in range(5):
        t, v, stat = M[c], first[c], stats[c]
        R, C, T = stat['rows'], stat['columns'], stat['mass']
        if not T:
            continue
        mean = F(sum(map(sum, v)), T)
        amean = [F(sum(v[i]), R[i]) if R[i] else F(0) for i in range(3)]
        bmean = [F(sum(v[i][j] for i in range(3)), C[j]) if C[j] else F(0)
                 for j in range(7)]
        x = [a-mean if R[i] else F(0) for i, a in enumerate(amean)]
        y = [b-mean if C[j] else F(0) for j, b in enumerate(bmean)]
        rec = mixed_receipt(t, x, y)
        cross_receipts.append(rec)
        m_energy += T*mean*mean
        x_energy += rec['left_energy']
        y_energy += rec['right_energy']
        cross += rec['cross']
        for i, j in product(range(3), range(7)):
            h = amean[i]+bmean[j]-mean
            z_energy += -2*h*v[i][j]+h*h*t[i][j]
            if t[i][j]:
                within -= F(v[i][j]**2, t[i][j])
    z_energy += second
    within += second
    check(F(second) == m_energy+x_energy+y_energy+z_energy-2*cross,
          "exact source ANOVA with interaction retained")
    check(z_energy >= within >= 0, "physical within-cell residual retained")
    check(all(s['interaction_dimension'] == 12 for s in stats), "actual tables full support")
    signed_rows = '\n'.join(records)+'\n'
    maxdelta = max(s['delta'] for s in stats)
    reciprocal_bound = maxdelta.denominator//maxdelta.numerator
    check(maxdelta <= F(1, reciprocal_bound), "exact reciprocal defect bound")
    return dict(N=n, source_rule='1 <= s < 4N, s odd, s mod8 != 5',
                source_columns=records, source_sha256=sha256(signed_rows.encode()).hexdigest(),
                source_count=len(records), positive_source_count=positive,
                zero_source_count=len(records)-positive, source_norm=sum(s['mass'] for s in stats),
                mass_tables=M, tables=stats, maximum_delta=maxdelta,
                delta_upper_reciprocal_denominator=reciprocal_bound,
                hidden_cell_dimension=60, within_cell_dimension=positive-105,
                source_coordinate_diagnostic=dict(total_energy=second, mean_energy=m_energy,
                    mod15_detail_energy=x_energy, mod35_detail_energy=y_energy,
                    mixed_term=cross, interaction_energy=z_energy,
                    within_cell_energy=within, per_fibre_mixed_receipts=cross_receipts))


def controls():
    f, g = (1, -1, 0), (1, -1, 0, 0, 0, 0, 0)
    M = [[2+f[i]*g[j] for j in range(7)] for i in range(3)]
    sharp = mixed_receipt(M, f, g)
    check(sharp == dict(cross=F(4), left_energy=F(28), right_energy=F(12),
                        delta=F(1, 21), bound_squared=F(16)), "sharp positive table")
    check(analyse([[0, 0], [0, 0]])['delta'] == 0, "empty mass convention")
    check(analyse([[1, 0], [0, 1]])['delta'] == 1, "disconnected commutation hostile")
    rejected = 0
    for bad in (True, 2.0, -1, '2'):
        try:
            native_rectangle(bad)
        except ValueError:
            rejected += 1
    for bad in ([[True]], [[1.0]], [[-1]], [[1], [1, 2]], []):
        try:
            analyse(bad)
        except ValueError:
            rejected += 1
    for bad in (True, 4096.0, 1):
        try:
            source_table(bad)
        except ValueError:
            rejected += 1
    check(rejected == 12, "strict public types and declared source universe")
    rectangles = [native_rectangle(k) for k in (1, 2, 3, 8, 16, 32, 64)]
    check([r['source'] for r in rectangles[1]['rows']] == [231, 351, 511, 631],
          "literal rectangle source binding")
    return dict(sharp_positive_table=M, sharp=sharp, native_rectangles=rectangles,
                support_cycles=cycle_controls(), malformed_rejections=rejected)


def jsonable(x):
    if type(x) is F:
        return str(x)
    if type(x) is dict:
        return {k: jsonable(v) for k, v in x.items()}
    if type(x) in (tuple, list):
        return [jsonable(v) for v in x]
    return x


def main():
    global CHECKS
    CHECKS = 0
    c = controls()
    rows = [source_table(n) for n in SCALES]
    output = dict(status='PROVED weighted interface; FINITE-EXACT source tables; arrival OPEN',
                  scales=SCALES, controls=c, source_tables=rows, checks=CHECKS,
                  source_dependency='root_clock_channels_20261018o.py:charge/norm only',
                  oracle_calls=0)
    encoded = json.dumps(jsonable(output), indent=2, sort_keys=True)+'\n'
    path = ROOT/'05-knowledge/results/weighted_crt_20261018g.json'
    if '--verify' in sys.argv:
        check(path.read_text(encoding='utf-8') == encoded, 'saved receipt exact match')
    else:
        path.write_text(encoded, encoding='utf-8')
    print(json.dumps(jsonable(dict(checks=output['checks'], sha256=sha256(encoded.encode()).hexdigest(),
        source_tables=[dict(N=r['N'], sources=r['source_count'], norm=r['source_norm'],
                            delta=str(r['maximum_delta'])) for r in rows])), indent=2))


if __name__ == '__main__':
    main()
