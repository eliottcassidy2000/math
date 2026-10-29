"""Source-preserving C7 primary charts and minimum injective pairing completions.

No assertions: all proof gates remain active under python -O.
"""
from collections import Counter
from itertools import combinations_with_replacement, product
from pathlib import Path
import hashlib
import json

from three_run_glue import three_runs, mm_mod, hostile_singleton_slices
from paley_pencil import local_layers

HERE = Path(__file__).resolve().parent
GATES = 0


def need(ok, why):
    global GATES
    GATES += 1
    if not ok:
        raise ValueError(why)


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def inverse_mod(a, q, p):
    n = len(a)
    b = [list(row)+unit for row, unit in zip(a, identity(n))]
    for j in range(n):
        r = next((r for r in range(j, n) if b[r][j] % p), None)
        need(r is not None, 'unit matrix over the local ring')
        b[r], b[j] = b[j], b[r]
        u = pow(b[j][j], -1, q)
        b[j] = [x*u % q for x in b[j]]
        for r in range(n):
            if r != j:
                c = b[r][j]
                b[r] = [(x-c*y) % q for x, y in zip(b[r], b[j])]
    return [row[n:] for row in b]


class Ring:
    """Unramified O/(p^K), with a supplied monic factor (low coefficients)."""
    def __init__(self, p, k, f):
        self.p, self.k, self.q, self.f = p, k, p**k, tuple(f)
        self.d = len(f)-1
        self.zero = (0,)*self.d
        self.one = (1,)+(0,)*(self.d-1)

    def add(self, a, b):
        return tuple((x+y) % self.q for x, y in zip(a, b))

    def sub(self, a, b):
        return tuple((x-y) % self.q for x, y in zip(a, b))

    def mul(self, a, b):
        if a == self.zero or b == self.zero:
            return self.zero
        c = [0]*(2*self.d-1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                c[i+j] += x*y
        for r in range(len(c)-1, self.d-1, -1):
            for j in range(self.d):
                c[r-self.d+j] -= c[r]*self.f[j]
        return tuple(x % self.q for x in c[:self.d])

    def scalar(self, a):
        return (a % self.q,)+(0,)*(self.d-1)

    def times(self, a):
        cols = [self.mul(a, tuple(int(i == j) for i in range(self.d)))
                for j in range(self.d)]
        return [list(row) for row in zip(*cols)]

    def inv(self, a):
        return tuple(row[0] for row in inverse_mod(self.times(a), self.q, self.p))

    def val(self, a):
        for v in range(self.k):
            if any(x % (self.p**(v+1)) for x in a):
                return v
        return self.k

    def matmul(self, a, b):
        bt = list(zip(*b))
        out = []
        for row in a:
            line = []
            for col in bt:
                z = self.zero
                for x, y in zip(row, col):
                    z = self.add(z, self.mul(x, y))
                line.append(z)
            out.append(line)
        return out


def local_smith_with_sources(M, ring):
    """O-linear row/column operations, retaining the quotient source map U."""
    n = len(M)
    a = [row[:] for row in M]
    I = [[ring.scalar(int(i == j)) for j in range(n)] for i in range(n)]
    U, Uinv = [row[:] for row in I], [row[:] for row in I]
    exponents = []
    for r in range(n):
        v, i, j = min((ring.val(a[i][j]), i, j)
                      for i in range(r, n) for j in range(r, n))
        need(v < ring.k, 'precision exceeds every actual elementary divisor')
        exponents.append(v)
        a[r], a[i] = a[i], a[r]
        U[r], U[i] = U[i], U[r]
        for row in Uinv:
            row[r], row[i] = row[i], row[r]
        for row in a:
            row[r], row[j] = row[j], row[r]
        pv = ring.p**v
        unit = tuple(x//pv for x in a[r][r])
        inv = ring.inv(unit)
        a[r] = [ring.mul(inv, x) for x in a[r]]
        U[r] = [ring.mul(inv, x) for x in U[r]]
        for row in Uinv:
            row[r] = ring.mul(unit, row[r])
        need(a[r][r] == ring.scalar(pv), 'normalized DVR pivot')
        for i in range(r+1, n):
            need(all(x % pv == 0 for x in a[i][r]), 'row quotient integral')
            c = tuple(x//pv for x in a[i][r])
            a[i] = [ring.sub(x, ring.mul(c, y)) for x, y in zip(a[i], a[r])]
            U[i] = [ring.sub(x, ring.mul(c, y)) for x, y in zip(U[i], U[r])]
            for row in Uinv:
                row[r] = ring.add(row[r], ring.mul(c, row[i]))
        for j in range(r+1, n):
            need(all(x % pv == 0 for x in a[r][j]), 'column quotient integral')
            c = tuple(x//pv for x in a[r][j])
            for i in range(n):
                a[i][j] = ring.sub(a[i][j], ring.mul(c, a[i][r]))
    need(a == [[ring.scalar(ring.p**exponents[i]) if i == j else ring.zero
                for j in range(n)] for i in range(n)], 'literal local diagonal')
    need(ring.matmul(U, Uinv) == I, 'source chart has an explicit inverse')
    UM = ring.matmul(U, M)
    for e, row in zip(exponents, UM):
        need(all(x % (ring.p**e) == 0 for a in row for x in a),
             'every original relation maps to zero at its exact quotient level')
    return exponents, U, Uinv


def split_factors(p, k):
    q = p**k
    if p == 2:
        eta = 0
        for h in range(1, k):
            eta = next(x for x in (eta, eta+2**h) if (x*x+x+2) % (2**(h+1)) == 0)
        eta_bar = (-1-eta) % q
        factors = [[-1, eta_bar, -eta, 1], [-1, eta, -eta_bar, 1]]
        f, g = factors
        coeffs = [sum(f[i]*g[j] for i in range(4) for j in range(4) if i+j == r) % q
                  for r in range(7)]
        need(coeffs == [1]*7, 'Hensel factors multiply to Phi7 at full precision')
        need(all(sum(c*x**i for i, c in enumerate(f)) % 2 for f in factors for x in (0, 1)),
             'both cubic factors irreducible modulo2')
    else:
        need(p == 3 and [pow(3, i, 7) for i in range(1, 7)] == [3, 2, 6, 4, 5, 1],
             'Phi7 irreducible over F3 by Frobenius degree6')
        factors = [[1]*7]
    return [[-1, 1]]+factors


def paired_partition(partition):
    a = sorted(partition, reverse=True)
    return [x for e in a[::2] for x in (e, e)]


def common_partition(a, b):
    a, b = sorted(a, reverse=True), sorted(b, reverse=True)
    return [max(x, y) for x, y in zip(a+[0]*max(0, len(b)-len(a)),
                                      b+[0]*max(0, len(a)-len(b))) ]


def trace_dual_certificate(plus, minus, exponent):
    """Actual reciprocal trace pairing, including the dyadic alternating guard."""
    q, d = plus.p**exponent, plus.d
    t = (0, 1)+(0,)*(d-2)
    tinv = plus.inv(t)
    pos, neg = [plus.one], [plus.one]
    for _ in range(d-1):
        pos.append(plus.mul(pos[-1], t))
        neg.append(plus.mul(neg[-1], tinv))
    T = [[sum(plus.times(plus.mul(x, y))[i][i] for i in range(d)) % q
          for y in neg] for x in pos]
    inverse_mod(T, q, plus.p)
    A, B = plus.times(t), minus.times(t)
    need(mm_mod(mm_mod(list(map(list, zip(*A))), T, q), B, q) == T,
         'reciprocal trace pairing preserves the literal two clocks')
    J = [[0]*(2*d) for _ in range(2*d)]
    for i, j in product(range(d), repeat=2):
        J[i][d+j], J[d+j][i] = T[i][j], -T[i][j] % q
    need(all(J[i][i] == 0 for i in range(2*d)) and
         all((J[i][j]+J[j][i]) % q == 0 for i, j in product(range(2*d), repeat=2)),
         'genuine alternating pairing, not merely skew in characteristic2')
    inverse_mod(J, q, plus.p)
    return {'exponent':exponent, 'cross_trace_matrix':T, 'alternating_matrix':J}


def minimum_tests():
    # Complete bounded search over ambient partitions, not the predicted formula.
    parts = [list(reversed(a)) for n in range(6)
             for a in combinations_with_replacement(range(1, 5), n)]
    paired = [[x for e in a for x in (e, e)] for n in range(4)
              for a in combinations_with_replacement(range(4, 0, -1), n)]
    def embeds(a, b):
        return len(a) <= len(b) and all(x <= y for x, y in zip(a, b))
    for a in parts:
        candidates = [sum(b) for b in paired if embeds(a, b)]
        need(sum(paired_partition(a)) == min(candidates), 'exhaustive abstract embedding minimum')
    small = [a for a in parts if len(a) <= 3]
    for a, b in product(small, repeat=2):
        expected = min(sum(c) for c in small if embeds(a, c) and embeds(b, c))
        need(sum(common_partition(a, b)) == expected, 'exhaustive reciprocal embedding minimum')
    # Explicit hostile: Z/2 + Z/4 embeds in (Z/4)^2 with index2,
    # whereas a direct-sum repair has index8.
    images = {(2*a % 4, b) for a in range(2) for b in range(4)}
    need(len(images) == 8, 'root-adjoining map is injective')
    need(all(all((x[0]*y[1]-x[1]*y[0]) % 4 == 0 for y in product(range(4), repeat=2))
             == (x == (0, 0)) for x in product(range(4), repeat=2)), 'perfect ambient pairing')
    # No alternative choice of source generator helps with ordinary reduction:
    # every order-dividing2 element of (Z/4)^2 reduces to zero in (Z/2)^2.
    order_two = [x for x in product(range(4), repeat=2) if all(2*y % 4 == 0 for y in x)]
    need(len(order_two) == 4 and all(all(y % 2 == 0 for y in x) for x in order_two),
         'pointwise minimum completions cannot form an injective quotient tower')
    return {'abstract_partitions':len(parts), 'reciprocal_pairs':len(small)**2,
            'hostile_source_orders':[2, 4], 'embedding':'(a,b)->(2a,b) in (Z/4)^2',
            'embedding_index':2, 'direct_sum_index':8,
            'ordinary_reduction_tower_hostile':{'order_two_elements_mod4':order_two,
                                               'all_reduce_to_zero_mod2':True}}


def primary_case(p, k, M, ports, orbits):
    q = p**k
    rings = [Ring(p, k, f) for f in split_factors(p, k)]
    powers = []
    for ring in rings:
        t = ring.scalar(1) if ring.d == 1 else (0, 1)+(0,)*(ring.d-2)
        row = [ring.one]
        for _ in range(7):
            row.append(ring.mul(row[-1], t))
        need(row[7] == ring.one, 'actual clock has order dividing7 over full local ring')
        powers.append(row)
    F = [list(row) for ring, ps in zip(rings, powers) for row in zip(*ps[:7])]
    Finv = inverse_mod(F, q, p)
    need(mm_mod(F, Finv, q) == identity(7), 'complete source Fourier chart, no lost coefficient')
    nb = len(orbits)
    operators = [[[ring.zero for _ in orbits] for _ in orbits] for ring in rings]
    for i, j in product(range(nb), repeat=2):
        B = [[M[u][v] % q for v in orbits[j]] for u in orbits[i]]
        kernel = B[0]
        expected = [[0]*7 for _ in range(7)]
        offset = 0
        for h, (ring, ps) in enumerate(zip(rings, powers)):
            value = ring.zero
            for r, c in enumerate(kernel):
                value = ring.add(value, ring.mul(ring.scalar(c), ps[(-r) % 7]))
            operators[h][i][j] = value
            block = ring.times(value)
            for u, v in product(range(ring.d), repeat=2):
                expected[offset+u][offset+v] = block[u][v]
            offset += ring.d
        need(mm_mod(mm_mod(F, B, q), Finv, q) == expected,
             'every named 7x7 block has the claimed factor action and inverse phase')
    records = []
    all_layers = [0]*k
    partitions = []
    for ring, A in zip(rings, operators):
        exps, U, Uinv = local_smith_with_sources(A, ring)
        expanded = [[0]*(nb*ring.d) for _ in range(nb*ring.d)]
        for i, j in product(range(nb), repeat=2):
            block = ring.times(A[i][j])
            for u, v in product(range(ring.d), repeat=2):
                expanded[i*ring.d+u][j*ring.d+v] = block[u][v]
        actual = local_layers(expanded, p, k)
        counts = Counter(exps)
        need(actual == [ring.d*counts[e] for e in range(k)],
             'independent integer local elimination verifies every O-module exponent')
        all_layers = [a+b for a, b in zip(all_layers, actual)]
        positive = sorted((e for e in exps if e), reverse=True)
        partitions.append(positive)
        records.append({'factor_low_coefficients':[x % q for x in ring.f],
                        'degree':ring.d, 'O_module_partition':positive,
                        'all_30_pivot_valuations':exps, 'Zp_layers':actual,
                        'source_row_map_U':U, 'source_row_map_inverse':Uinv})
    known = {2:[165, 31, 2, 3, 9, 0], 3:[201, 2, 6, 0, 1, 0]}[p]
    need(all_layers == known, 'independent complete integer Smith profile recovered')
    trivial = paired_partition(partitions[0])
    if p == 2:
        dual = common_partition(partitions[1], partitions[2])
        targets = [trivial, dual, dual]
    else:
        targets = [trivial, partitions[1]]
    source_order = sum(r.d*sum(a) for r, a in zip(rings, partitions))
    ambient_order = sum(r.d*sum(a) for r, a in zip(rings, targets))
    maps = []
    for record, source, target in zip(records, partitions, targets):
        need(len(source) <= len(target) and all(a <= b for a, b in zip(source, target)),
             'all original cyclic module generators inject into the ambient group')
        source_rows = sorted((i for i, e in enumerate(record['all_30_pivot_valuations']) if e),
                             key=lambda i:(-record['all_30_pivot_valuations'][i], i))
        maps.append([{'source_row':source_rows[j], 'source_exponent':a,
                      'ambient_slot':j, 'ambient_exponent':b,
                      'multiply_by':p**(b-a)} for j, (a, b) in enumerate(zip(source, target))])
    direct = sum(r.d*sum(e*(m % 2) for e, m in Counter(partitions[0]).items())
                 for r in rings[:1])
    if p == 2:
        c1, c2 = Counter(partitions[1]), Counter(partitions[2])
        direct += 3*sum(e*abs(c1[e]-c2[e]) for e in c1.keys() | c2.keys())
    layer_costs = []
    for h in range(1, 1+max(e for a in partitions for e in a)):
        ranks = [sum(e >= h for e in a) for a in partitions]
        cost = ranks[0] % 2 + (3*abs(ranks[1]-ranks[2]) if p == 2 else 0)
        truncated = [[min(e, h) for e in a] for a in partitions]
        min_target = [paired_partition(truncated[0])]
        min_target += ([common_partition(*truncated[1:])]*2 if p == 2 else [truncated[1]])
        actual_cost = sum(r.d*(sum(b)-sum(a)) for r, a, b in zip(rings, truncated, min_target))
        need(actual_cost == sum(x['added_log_index'] for x in layer_costs)+cost,
             'all-level minimum order recurrence agrees with independent partition matching')
        layer_costs.append({'height':h, 'isotype_ranks':ranks, 'added_log_index':cost,
                            'cumulative_log_index':actual_cost})
    trace_certificates = ([trace_dual_certificate(rings[1], rings[2], e)
                           for e in sorted(set(targets[1]))] if p == 2 else [])
    return {'prime':p, 'precision':k, 'modulus':q, 'source_orbit_chart':F,
            'source_orbit_chart_inverse':Finv, 'factors':records,
            'combined_integer_layers':all_layers, 'ambient_partitions':targets,
            'source_to_ambient_maps':maps, 'source_log_order':source_order,
            'minimum_ambient_log_order':ambient_order,
            'minimum_embedding_log_index':ambient_order-source_order,
            'minimum_direct_sum_log_index':direct, 'layer_recurrence':layer_costs,
            'reciprocal_pairing_certificates':trace_certificates}


def main():
    rows = [sum(1 << ((u+s) % 7) for s in (1, 2, 4)) for u in range(7)]
    _, ports, M = three_runs(rows)
    index = {x:i for i, x in enumerate(ports)}
    move = [index[tuple((x+1) % 7 for x in t)] for t in ports]
    orbits, seen = [], set()
    for i in range(len(ports)):
        if i in seen:
            continue
        orbit = [i]
        for _ in range(6):
            orbit.append(move[orbit[-1]])
        need(len(set(orbit)) == 7 and move[orbit[-1]] == i, 'free C7 source orbit')
        orbits.append(orbit)
        seen.update(orbit)
    need(len(orbits) == 30 and len(seen) == 210, 'complete named source partition')
    need(all(M[move[i]][move[j]] == M[i][j] for i, j in product(range(210), repeat=2)),
         'native translation is an actual source and target clock')
    cases = []
    for p in (2, 3):
        case = primary_case(p, 6, M, ports, orbits)
        cases.append(case)
        print(json.dumps({'prime':p, 'partitions':[r['O_module_partition'] for r in case['factors']],
                          'embedding_log_index':case['minimum_embedding_log_index'],
                          'direct_sum_log_index':case['minimum_direct_sum_log_index']}), flush=True)
    tests = minimum_tests()
    hostile = hostile_singleton_slices()
    out = {'status':'FINITE-EXACT source modules; PROVED scoped embedding theorem in PRIMARY_CLOCK_COMPLETION.md',
           'ports':ports, 'translation_orbits':orbits,
           'native_matrix_sha256':hashlib.sha256(json.dumps(M, separators=(',', ':')).encode()).hexdigest(),
           'primary_cases':cases, 'minimum_tests':tests, 'hostile_81_27':hostile, 'gates':GATES,
           'scope':'C7-equivariant injective completion; additional affine/run clocks separate'}
    (HERE/'primary_clock_layers.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'gates':GATES, 'hostile_81_27':hostile}), flush=True)


if __name__ == '__main__':
    main()
