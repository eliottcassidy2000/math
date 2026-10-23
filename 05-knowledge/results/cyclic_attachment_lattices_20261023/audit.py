"""Separate direct tournament and small-group paths, without lattice imports."""
from collections import Counter
from itertools import permutations, product
from pathlib import Path
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
CHECKS = 0


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def cyclic(m):
    return [[int(1 <= (j-i) % m <= (m-1)//2) for j in range(m)] for i in range(m)]


def native(m, twist):
    inside = cyclic(m)
    n = 3*m
    T = [[0]*n for _ in range(n)]
    for c in range(3):
        for i, j in product(range(m), repeat=2):
            T[c*m+i][c*m+j] = inside[i][j]
    for c, d in ((0, 1), (1, 2), (0, 2)):
        phase = twist if (c, d) == (0, 2) else 0
        for i, j in product(range(m), repeat=2):
            val = int((j-i-phase) % 3 == 0)
            T[c*m+i][d*m+j] = val
            T[d*m+j][c*m+i] = 1-val
    return T


def actual_graphs():
    automorphisms = {}
    for m in (3, 5, 7, 9):
        T = cyclic(m)
        maps = [p for p in permutations(range(m))
                if all(T[i][j] == T[p[i]][p[j]] for i in range(m) for j in range(m))]
        need(set(maps) == {tuple((i+a) % m for i in range(m)) for a in range(m)},
             'full internal automorphism group differs')
        automorphisms[m] = len(maps)
    rows = []
    for m in (3, 9, 15):
        A, B = native(m, 0), native(m, 1)
        n = 3*m
        need(all(A[i][j]+A[j][i] == int(i != j)
                 and B[i][j]+B[j][i] == int(i != j)
                 for i in range(n) for j in range(n)), 'not tournaments')
        degrees = [(m-1)//2+2*m//3, (m-1)//2+m, (m-1)//2+4*m//3]
        need([sum(row) for row in A] == [d for d in degrees for _ in range(m)], 'cell degree separation')
        need([sum(row) for row in B] == [d for d in degrees for _ in range(m)], 'target degree separation')
        counts = Counter()
        for a in product(range(m), repeat=3):
            perm = [c*m+(i+a[c]) % m for c in range(3) for i in range(m)]
            for label, target in (('flat', A), ('twisted', B)):
                valid = all(A[i][j] == target[perm[i]][perm[j]] for i in range(n) for j in range(n))
                counts[label] += valid
        need(counts['flat'] == m**3//9 and counts['twisted'] == 0, 'full-arc native count')
        rows.append(dict(cell_size=m, tournament_order=n, degree_classes=degrees,
                         maps=dict(counts), source_sha256=hashlib.sha256(bytes(sum(A, []))).hexdigest(),
                         target_sha256=hashlib.sha256(bytes(sum(B, []))).hexdigest()))
    return dict(internal_full_permutation_tests=automorphisms, native_comparisons=rows)


def small_groups():
    def mul(x, y):
        a, b, c = x
        d, e, f = y
        return ((a+d) % 3, (b+e) % 3, (c+f+a*e) % 3)
    def inv(x):
        a, b, c = x
        return (-a % 3, -b % 3, (a*b-c) % 3)
    H = list(product(range(3), repeat=3))
    need(all(mul(mul(a, b), c) == mul(a, mul(b, c)) for a, b, c in product(H, repeat=3)), 'UT3 associativity')
    need(all(mul(a, inv(a)) == (0, 0, 0) for a in H), 'UT3 inverse')
    X, Y = (1, 0, 0), (0, 1, 0)
    comm = mul(mul(mul(X, Y), inv(X)), inv(Y))
    need(comm == (0, 0, 1), 'lost odd central charge')
    centre = [a for a in H if all(mul(a, b) == mul(b, a) for b in H)]
    need(centre == [(0, 0, c) for c in range(3)], 'UT3 centre')

    def qm(x, y):
        s, a = x
        t, b = y
        if a == 0:
            return (s*t, b)
        if b == 0:
            return (s*t, a)
        if a == b:
            return (-s*t, 0)
        sign = 1 if (a, b) in ((1, 2), (2, 3), (3, 1)) else -1
        return (s*t*sign, a ^ b)
    Q = list(product((1, -1), range(4)))
    need(all(qm(qm(a, b), c) == qm(a, qm(b, c)) for a, b, c in product(Q, repeat=3)), 'quaternion associativity')
    qcomm = qm(qm(qm((1, 1), (1, 2)), (-1, 1)), (-1, 2))
    need(qcomm == (-1, 0), 'quaternion lift hostile')

    # Toy common-C2 product; not an enumeration of the Monster or He.2.
    A = list(product(range(1, 7), range(7)))
    C = list(product(range(5), range(2)))
    def am(x, y):
        return (x[0]*y[0] % 7, (x[1]+x[0]*y[1]) % 7)
    def eps(x):
        return int(pow(x[0], 3, 7) == 6)
    need(all(eps(am(x, y)) == (eps(x)+eps(y)) % 2 for x, y in product(A, repeat=2)), 'sign homomorphism')
    fibre = [(a, b) for a in A for b in C if eps(a) == b[1]]
    need(len(fibre) == 210 and len(A)*len(C) == 420, 'common-sign fibre size')
    need({a for a, b in fibre} == set(A) and {b for a, b in fibre} == set(C), 'common-sign projections')
    need(all(eps(am(a, c)) == (b[1]+d[1]) % 2 for (a, b), (c, d) in product(fibre, repeat=2)), 'fibre closure')
    plus, minus = 1+4371+96255, 96256
    need(plus+minus == 196883 and plus-minus == 4371, 'cited character arithmetic')
    return dict(UT3_order=len(H), UT3_central_commutator=comm, quaternion_commutator=qcomm,
                common_sign_toy=dict(product_size=420, compatible_size=210),
                cited_dimension_arithmetic=dict(plus=plus, minus=minus, total=plus+minus, trace=plus-minus))


def phase_orbits():
    output = []
    for d, v, edges in ((3, 4, [(i, j) for i in range(4) for j in range(i+1, 4)]),
                        (9, 3, [(0, 1), (1, 2), (0, 2)]),
                        (5, 3, [(0, 1), (1, 2)])):
        unseen = set(product(range(d), repeat=len(edges)))
        sizes = Counter()
        while unseen:
            seed = min(unseen)
            seen = {seed}
            queue = [seed]
            while queue:
                state = queue.pop()
                for k in range(v):
                    target = tuple((c+int(j == k)-int(i == k)) % d for c, (i, j) in zip(state, edges))
                    if target not in seen:
                        seen.add(target)
                        queue.append(target)
            need(seen <= unseen, 'gauge orbits overlap')
            unseen -= seen
            sizes[len(seen)] += 1
        beta = len(edges)-v+1
        need(dict(sizes) == {d**(v-1): d**beta}, 'all-size phase formula control')
        output.append(dict(quotient_modulus=d, cells=v, edges=edges,
                           labelled_phase_systems=d**len(edges), orbit_size_counts=dict(sizes)))
    return output


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--write', action='store_true')
    args = p.parse_args()
    output = dict(status='FINITE-EXACT independent direct paths; same-agent audit',
                  actual_graphs=actual_graphs(), groups=small_groups(),
                  phase_orbits=phase_orbits(), checks=CHECKS)
    output = json.loads(json.dumps(output))
    path = HERE/'audit.json'
    if args.write:
        path.write_text(json.dumps(output, indent=2)+'\n')
    else:
        need(json.loads(path.read_text()) == output, 'audit receipt mismatch')
    print(json.dumps(output, indent=2))
