#!/usr/bin/env python3
"""Independent final audit: binary row operators and pencil/parity Fano generation.
No producer imports. The finite universe ends at n=5 and 30 labelled systems.
"""
import hashlib
import itertools as it
import json
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
CHECKS = 0
PINS = {
    'bridges/RESULT.md': '29fc9ca464171bd1f5ef30602786cac2bbe222bda2b5b7cab9deb34e7e9699b8',
    'bridges/CORRECTION.md': 'dfa1cee6f976feb6b3873957eceac8d46edcae14e215f3a1a91e234944d20c63',
    'bridges/check_bridges.py': '89fae3be3a595d5de0b7427485092c24a5099e9cba3894718a08746775f3da95',
    'bridges/normal.json': '1d98df0f4bbc99b7f4dc6969eb8f6291a22b6ac21a6998c40b4e5ca3b28a82f9',
    'bridges/optimized.json': '1d98df0f4bbc99b7f4dc6969eb8f6291a22b6ac21a6998c40b4e5ca3b28a82f9',
    'bridges/FANO_ORBITS.md': 'cd7bb7e6d9c8d79c995e835f8b374c4d819e554eba86f03b34431acb027f28cb',
    'bridges/check_fano_orbits.py': '8effd5070d1b6d79718a205d969fa22b3ae54ea408d8d1c3d8d5606da6a504cd',
    'bridges/fano.normal.json': '2d0525d93012115e5ce251e76d9e73f8702e29b219f9943eeab47f7ecd4518e7',
    'bridges/fano.optimized.json': '2d0525d93012115e5ce251e76d9e73f8702e29b219f9943eeab47f7ecd4518e7',
    'incoming/check_stencil.py': 'e5d932483fc4e52fd381398ad262310bbf56a6b5ff13bebf1e2c8f1f6653b63a',
}

def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)

def parity(x):
    return x.bit_count() & 1

def mat(n, code):
    a = [0] * n
    for k, (i, j) in enumerate(it.combinations(range(n), 2)):
        u, v = (i, j) if code >> k & 1 else (j, i)
        a[u] |= 1 << v
    return tuple(r | (parity(r) << i) for i, r in enumerate(a))

def apply(a, x):
    return sum(parity(r & x) << i for i, r in enumerate(a))

def phi(a):
    rows = []
    for i, r in enumerate(a):
        square = 0
        for j, s in enumerate(a):
            if r >> j & 1:
                square ^= s
        rows.append(square ^ r ^ (1 << i))
    return tuple(rows)

def matchings(xs):
    if not xs:
        yield ()
    else:
        a = xs[0]
        for b in xs[1:]:
            for tail in matchings(tuple(x for x in xs[1:] if x != b)):
                yield ((a, b),) + tail

def system(lines):
    return frozenset(tuple(sorted(t)) for t in lines)

def act(F, p):
    return system(tuple(p[v] for v in t) for t in F)

def main():
    for name, digest in PINS.items():
        need(hashlib.sha256((PACKET / name).read_bytes()).hexdigest() == digest, 'frozen bytes: ' + name)
    matrix_rows = []
    parity_hostile = None
    for n in range(1, 6):
        full = (1 << (n * (n - 1) // 2)) - 1
        for code in range(full + 1):
            a, b = mat(n, code), mat(n, full ^ code)
            p, q = phi(a), phi(b)
            need(p == q, 'same labelled Phi under converse')
            need(all((p[i] >> j & 1) == (p[j] >> i & 1) for i in range(n) for j in range(n)), 'symmetric Phi')
            need(apply(p, (1 << n) - 1) == (1 << n) - 1, 'Phi fixes ones')
            for x in range(1 << n):
                if apply(p, x):
                    continue
                need(parity(x) == 0, 'kernel in ones orthogonal')
                need(apply(b, x) ^ apply(a, x) == (x if n & 1 else 0), 'exact embedding change')
                y = apply(a, x) ^ x
                need(apply(a, x) ^ x ^ y == 0 and x ^ apply(a, y) == 0, 'literal Monsky block kernel')
                if n == 3 and x and parity_hostile is None:
                    parity_hostile = {'n': n, 'code': code, 'x': x, 'second_half': y, 'converse_second_half': apply(b, x) ^ x}
            need(p == q and apply(p, 0) == 0, 'kernel equality follows same operator')
        matrix_rows.append({'n': n, 'all_labelled': full + 1})
    need(parity_hostile is not None and parity_hostile['second_half'] != parity_hostile['converse_second_half'], 'odd embedding really can change')
    witness_kernels = {c: [x for x in range(32) if apply(phi(mat(5, c)), x) == 0] for c in (8, 10)}
    need(witness_kernels == {8: [0, 15, 17, 30], 10: [0]}, 'same-Q witness kernel certificates')

    # Distinct from the producer's exact-cover and vector relabelling generators:
    # choose one of 15 pencils through 0 and either transversal parity.
    pencils = list(matchings(tuple(range(1, 7))))
    need(len(pencils) == 15, '15 pencils')
    fano = set()
    for pencil in pencils:
        for eps in (0, 1):
            lines = [(0, *pair) for pair in pencil]
            lines += [tuple(pencil[j][bits[j]] for j in range(3))
                      for bits in it.product((0, 1), repeat=3) if sum(bits) % 2 == eps]
            F = system(lines)
            need(len(F) == 7 and Counter(pair for t in F for pair in it.combinations(t, 2)) == Counter(it.combinations(range(7), 2)), 'literal Fano predicate')
            fano.add(F)
    need(len(fano) == 30, 'complete pencil census')
    group = [tuple((a*x+b) % 7 for x in range(7)) for a in (1, 2, 4) for b in range(7)]
    plus = system(tuple((x+b) % 7 for x in (1, 2, 4)) for b in range(7))
    minus = system(tuple((x+b) % 7 for x in (3, 5, 6)) for b in range(7))
    need(plus != minus and plus in fano and minus in fano, 'distinct difference systems')
    todo, orbit_sizes = set(fano), []
    while todo:
        F = min(todo, key=lambda z: tuple(sorted(z)))
        orbit = {act(F, p) for p in group}
        need(orbit <= fano, 'literal action preserves census')
        orbit_sizes.append(len(orbit))
        todo -= orbit
    need(sorted(orbit_sizes) == [1, 1, 7, 7, 7, 7], 'Fano profile')
    fixed_counts = [sum(act(F, p) == F for F in fano) for p in group]
    need(fixed_counts == [30] + [2]*6 + [6]*14, 'element fixed counts')
    need(sum(fixed_counts) == 126, 'Burnside six orbits')
    addresses = {}
    for F in fano - {plus, minus}:
        stab = [p for p in group if act(F, p) == F]
        centres = [c for c in range(7) if all(p[c] == c for p in stab)]
        need(len(stab) == 3 and len(centres) == 1, 'unique stabilizer centre')
        typ = (len(F & plus), len(F & minus))
        key = (*typ, centres[0])
        need(key not in addresses, 'type centre injective')
        addresses[key] = F
    types = {(0, 1), (1, 0), (1, 3), (3, 1)}
    need(set(addresses) == {(*t, c) for t in types for c in range(7)}, 'four affine copies')
    for (u, v, c), F in addresses.items():
        for p in group:
            need(act(F, p) == addresses[u, v, p[c]], 'complete address equivariance')
    displayed = {
        (0, 1): '123 034 025 145 016 246 356',
        (1, 0): '124 034 025 135 016 236 456',
        (1, 3): '023 124 015 345 136 046 256',
        (3, 1): '013 234 125 045 026 146 356',
    }
    for typ, row in displayed.items():
        need(system(tuple(map(int, t)) for t in row.split()) == addresses[typ + (0,)], 'displayed centre-zero inverse')
    larger = (0, 1, 4, 3, 2, 6, 5)
    need(act(plus, larger) == plus and act(minus, larger) not in (plus, minus), 'explicit larger-group failure')
    aut_count = sum(act(plus, p) == plus for p in it.permutations(range(7)))
    need(aut_count == 168, 'full Fano automorphism group')
    print(json.dumps({'status': 'FINITE-EXACT independent audit, no status promotion', 'checks': CHECKS,
                      'matrix_universe': matrix_rows, 'odd_embedding_hostile': parity_hostile,
                      'same_Q_witness_kernels': witness_kernels, 'Fano_pencils': len(pencils),
                      'Fano_systems': len(fano), 'orbit_profile': sorted(orbit_sizes),
                      'fixed_count_multiplicities': dict(Counter(fixed_counts)),
                      'address_types': sorted(types), 'full_Fano_automorphisms': aut_count,
                      'source_pins': PINS, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, sort_keys=True, indent=2))

if __name__ == '__main__':
    main()
