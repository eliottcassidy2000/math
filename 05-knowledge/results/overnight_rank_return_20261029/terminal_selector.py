"""Exact selector closure of the inherited periodic complement observer.

The operation is permutation of the three path names / terminal assignments,
not arbitrary vertex re-marking or another actual complement.
"""
from collections import Counter, deque
from itertools import combinations, permutations
import json

from flint import fmpz_mat
from sympy import factorint

from complement_kernel import HERE, need, rank_mod
from defect_transfer import interval_target, query, subset_query
from paley_pencil import local_layers, valuation
from periodic_defect import initial, bulk, step
from two_run import subset_paths

PERMS = list(permutations(range(3)))
SIGNS = [1, -1, -1, 1, 1, -1]


def terminal(state, perm):
    weights = {state: 1}
    last = {0: 3, 1: 2, 2: 5, 3: 4, 4: -1, 5: -1}
    for v, back in ((2, {1}), (3, {0}), (4, {3}), (5, {2})):
        weights = step(weights, v, back, last,
                       end=perm[v-3] if v >= 3 else None)
    return weights.get((7, ()), 0)


def action(state, perm):
    done, components = state
    return (sum(1 << perm[i] for i in range(3) if done >> i & 1),
            tuple(sorted((-(perm[-h-2]+2) if h < -1 else h,
                          -(perm[-t-2]+2) if t < -1 else t)
                         for h, t in components)))


def mv(B, v, modulus=None):
    out = [sum(a*b for a, b in zip(row, v)) for row in B]
    return [a % modulus for a in out] if modulus else out


def polynomial(B, v, coefficients, modulus=None):
    out = [0]*len(v)
    for c in coefficients:
        out = [a+c*b for a, b in zip(out, v)]
        v = mv(B, v, modulus)
    return [a % modulus for a in out] if modulus else out


def krylov(B, v, length):
    rows = []
    for _ in range(length):
        rows.append(v)
        v = mv(B, v)
    return rows


def lattice(rows):
    return [list(map(int, row)) for row in fmpz_mat(rows).hnf().table()
            if any(row)]


def brute_covers(n, defects, perm=(0, 1, 2)):
    inverse = [perm.index(i) for i in range(3)]
    ends = [n-3+inverse[i] for i in range(3)]
    for middle in permutations(range(3, n-3)):
        for a in range(n-5):
            for b in range(a, n-5):
                paths = ((0,)+middle[:a]+(ends[0],),
                         (1,)+middle[a:b]+(ends[1],),
                         (2,)+middle[b:]+(ends[2],))
                if all((u < v) != ((min(u, v), max(u, v)) in defects)
                       for path in paths for u, v in zip(path, path[1:])):
                    yield paths


def run():
    old = json.loads((HERE/'periodic_defect.json').read_text())
    states = [(d, tuple(map(tuple, pairs))) for d, pairs in old['states']]
    index = {s: i for i, s in enumerate(states)}
    B, alpha = old['transfer'], old['initial']
    init = initial()
    need(alpha == [init.get(s, 0) for s in states], 'inherited initial law')
    seen, pending = set(init), deque(init)
    while pending:
        s = pending.popleft()
        transitions = bulk(s)
        need(B[index[s]] == [transitions.get(t, 0) for t in states],
             'rebuild every bulk transition')
        for t in transitions:
            need(t in index, 'closed 35-state domain')
            if t not in seen:
                seen.add(t)
                pending.append(t)
    need(seen == set(states) and len(states) == 35, 'exact reachable universe')
    betas = [[terminal(s, perm) for s in states] for perm in PERMS]
    need(betas[0] == old['terminal'], 'identity terminal specializes')
    transports = []
    for perm, beta in zip(PERMS, betas):
        transport = [index[action(s, perm)] for s in states]
        transports.append(transport)
        need(all(B[transport[i]][transport[j]] == B[i][j]
                 for i in range(35) for j in range(35)), 'common S3 bulk transport')
        inverse = tuple(perm.index(i) for i in range(3))
        need(beta == [betas[0][index[action(s, inverse)]] for s in states],
             'terminal assignment is the inverse common-name transport')
    coefficients = old['recurrence_low_first']
    for beta in betas:
        need(not any(polynomial(B, beta, coefficients)), 'same full vector polynomial')
    Os = [krylov(B, beta, 18) for beta in betas]
    O = [row for rows in Os for row in rows]
    need(Os[0] == old['observable'], 'old observer retained literally')
    need(fmpz_mat(O).rank() == rank_mod(O) == 24, 'exact joint rank 24')
    total = [sum(v[i] for v in betas) for i in range(35)]
    alternating = [sum(s*v[i] for s, v in zip(SIGNS, betas)) for i in range(35)]
    standard = [[6*v[i]-total[i]-s*alternating[i] for i in range(35)]
                for s, v in zip(SIGNS, betas)]
    f11 = [-576, -1232, 459, 3197, 2072, -1048, -1170, -382, -268, -16, -13, 1]
    f6 = [-4, 5, -7, 4, -6, -1, 1]
    sectors = []
    for name, vectors, poly, expected in (
            ('trivial', [total], f11, 11),
            ('sign', [alternating], [1, 1], 1),
            ('standard', standard, f6, 12)):
        need(all(not any(polynomial(B, v, poly)) for v in vectors),
             'exact sector annihilation')
        rows = [row for v in vectors for row in krylov(B, v, len(poly)-1)]
        need(fmpz_mat(rows).rank() == rank_mod(rows) == expected,
             'sector dimension with independent modular lower bound')
        sectors.append(dict(name=name, rank=expected, polynomial=poly))

    H = lattice(O)
    D = fmpz_mat(O).snf()
    ds = [abs(int(D[i, i])) for i in range(35) if D[i, i]]
    expected_ds = [1]*9+[2]*3+[4]*6+[12]*3+[24, 1344, 147676566374644608]
    need(ds == expected_ds, 'full integral selector observation Smith factors')
    layers = {}
    primes = sorted({int(p) for value in ds for p in factorint(value)})
    for p in primes:
        counts = Counter(valuation(value, p) for value in ds)
        precision = max(counts)+1
        layers[str(p)] = local_layers(O, p, precision)
        need(layers[str(p)] == [counts[j] for j in range(precision)],
             'independent primary Smith elimination')
    banks = []
    for size in range(1, 7):
        for rest in combinations(range(1, 6), size-1):
            bank = (0,)+rest
            rows = [row for i in bank for row in Os[i]]
            banks.append(dict(selectors=bank, rational_rank=fmpz_mat(rows).rank(),
                              binary_rank=rank_mod(rows, 2), ternary_rank=rank_mod(rows, 3),
                              integral_complete=lattice(rows) == H))
    minima = dict(rational=min(len(b['selectors']) for b in banks if b['rational_rank'] == 24),
                  binary=min(len(b['selectors']) for b in banks if b['binary_rank'] == 9),
                  ternary=min(len(b['selectors']) for b in banks if b['ternary_rank'] == 18),
                  integral=min(len(b['selectors']) for b in banks if b['integral_complete']))
    need(minima == dict(rational=2, binary=4, ternary=2, integral=4), 'sharp fixed-bank minima')
    need(next(b for b in banks if b['selectors'] == (0, 1, 2, 3))['integral_complete'],
         'explicit universal four-selector bank')

    witness = [0]*35
    for i, value in ((2, -1), (4, -1), (11, 1), (18, 1)):
        witness[i] = value
    need(all(sum(x*y for x, y in zip(row, witness)) == 0 for row in Os[0]),
         'equal old observations for every future bulk insertion')
    responses = [sum(x*y for x, y in zip(beta, witness)) for beta in betas]
    need(responses == [0, -1, 1, 0, 0, 0] and sum(witness) == 0,
         'two nonnegative mass-two laws separated by one terminal swap')
    witness_counts = [[sum(beta[i] for i in (11, 18)), sum(beta[i] for i in (2, 4))]
                      for beta in betas]

    modular = []
    q3 = old['ternary_recurrence_low_first']
    for p, poly, primary in ((2, [1, 0, 0, 0, 1], [([1, 1], 4)]),
                            (3, q3, [([1, 1], 2), ([-1, 0, -1, -1, 1, 1], 2)])):
        need(all(not any(polynomial(B, beta, poly, p)) for beta in betas),
             'full modular joint-clock annihilator')
        nullities = []
        d = rank_mod(O, p)
        for factor, exponent in primary:
            rows = O
            dims = []
            for _ in range(exponent):
                rows = [polynomial(B, row, factor, p) for row in rows]
                dims.append(d-rank_mod(rows, p))
            nullities.append(dict(factor=factor, dimensions=dims))
        expected = [[4, 7, 8, 9]] if p == 2 else [[2, 3], [10, 15]]
        need([v['dimensions'] for v in nullities] == expected, 'primary Jordan layer ranks')
        modular.append(dict(prime=p, rank=d, polynomial=poly, nullities=nullities))

    values = [[sum(a*b for a, b in zip(alpha, row)) for row in krylov(B, beta, 54)]
              for beta in betas]
    parity = [v[0] % 2 for v in values]
    need(sum(coefficients) % 2 == 0 and parity == [0, 0, 0, 1, 1, 1]
         and all(all(x % 2 == c for x in v) for v, c in zip(values, parity)),
         'constant actual parity for each selector; monic recurrence closes all sizes')
    state = tuple(tuple(x % 3 for x in v[:18]) for v in values)
    seen = {}
    while state not in seen:
        seen[state] = len(seen)
        state = tuple(v[1:]+((-sum(c*x for c, x in zip(coefficients, v)) % 3),)
                      for v in state)
    need(seen[state] == 0 and len(seen) == 726, 'joint physical clock period remains 726')
    controls = []
    for r in range(8):
        n, path, defects = interval_target([2]*(r+4))
        if r < 3:
            rows = [sum(1 << v for v in range(n) if u != v
                        and ((u < v) != ((min(u, v), max(u, v)) in defects)))
                    for u in range(n)]
            table = subset_paths(rows)
        for j, perm in enumerate(PERMS):
            ends = tuple(n-3+perm.index(i) for i in range(3))
            count, stats = query(n, defects, (0, 1, 2), ends)
            need(count == values[j][r], 'independent general frontier query')
            if r < 3:
                need(subset_query(table, (0, 1, 2), ends) == count,
                     'independent subset query with selected matching')
            controls.append(dict(n=n, permutation=perm, value=count))
    reflection = []
    for n in (8, 10, 12):
        _, _, defects = interval_target([2]*(n//2))
        covers = list(brute_covers(n, defects))
        universe = set(covers)
        fixed = []
        for cover in covers:
            partner = tuple(tuple(n-1-v for v in reversed(path)) for path in reversed(cover))
            need(partner in universe, 'literal reflection preserves the selected covers')
            if partner == cover:
                fixed.append(cover)
        need(len(covers) == values[0][(n-8)//2], 'independent full cover enumeration')
        reflection.append(dict(n=n, covers=len(covers), fixed=len(fixed), witness=fixed[0]))
    need([r['fixed'] for r in reflection] == [2, 6, 26], 'natural reflection is not free')
    _, _, defects8 = interval_target([2]*4)
    need([len(list(brute_covers(8, defects8, perm))) for perm in PERMS]
         == [v[0] for v in values], 'six independent size-eight matching counts')
    return dict(permutations=PERMS, transports=transports, betas=betas,
                joint_rank=24, sectors=sectors, smith=ds, smith_layers=layers,
                row_lattice=H, banks=banks, minima=minima,
                witness=witness, witness_states=[states[i] for i in (2, 4, 11, 18)],
                witness_counts=witness_counts, modular=modular,
                first_values=[v[:8] for v in values], parity=parity,
                physical_ternary_period=726, independent_queries=controls,
                reflection=reflection)


if __name__ == '__main__':
    import sys
    receipt = json.loads(json.dumps(run()))
    path = HERE/'terminal_selector.json'
    if '--write' in sys.argv:
        path.write_text(json.dumps(receipt, separators=(',', ':'))+'\n')
    else:
        need(receipt == json.loads(path.read_text()), 'unchanged selector certificate')
    print('PASS terminal-selector closure: rank24; minimal banks Q2 / F2 4 / F3 2 / Z4')
