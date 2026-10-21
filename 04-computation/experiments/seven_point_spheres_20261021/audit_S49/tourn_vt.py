"""HYP-9353 predictions 2 and 3: vertex-transitive tournaments on 21 vertices (repo d6 files, parsed as data) and
Cayley tournaments on Z3^3, Z9xZ3, Z27 (all connection sets).  Random search for an Aut+- -regular subset; anything
without one is tested for isomorphism with the family member of that order.  Also: all 16 Cayley(Z3^2) ~ C3[C3]?"""
import sys, time, itertools
import numpy as np
from tlib import *

t0 = time.time()
rng = np.random.default_rng(21)
C3 = circ(3, {1}); T5 = circ(5, {1, 2}); P7 = circ(7, {1, 2, 4}); C33 = lexprod(C3, C3)


def parse_d6(line):
    s = line.strip()
    assert s[0] == '&'
    n = ord(s[1]) - 63
    assert 0 < n < 63
    bits = []
    for ch in s[2:]:
        v = ord(ch) - 63
        for k in range(5, -1, -1):
            bits.append((v >> k) & 1)
    rows = []
    for i in range(n):
        r = 0
        for j in range(n):
            if bits[i * n + j]:
                r |= 1 << j
        rows.append(r)
    return from_rows(rows)


# Z3^2 Cayley all ~ C3[C3]
Z32 = [(a, b) for a in range(3) for b in range(3)]
idx = {g: i for i, g in enumerate(Z32)}
pairs = [((1, 0), (2, 0)), ((0, 1), (0, 2)), ((1, 1), (2, 2)), ((1, 2), (2, 1))]
isoc = 0
for bits in range(16):
    Sset = [pairs[i][(bits >> i) & 1] for i in range(4)]
    rows = [sum(1 << idx[((g[0] + s[0]) % 3, (g[1] + s[1]) % 3)] for s in Sset) for g in Z32]
    isoc += iso_exists(from_rows(rows), C33)
print(f'Cayley tournaments on Z3^2 isomorphic to C3[C3]: {isoc} of 16')

# VT 21
base = sys.argv[1]
fam21 = lexprod(P7, C3)
for fn in ('circtourn21.d6', 'othervttourn21.d6'):
    with open(base + '/' + fn) as fh:
        Ts = [parse_d6(l) for l in fh if l.strip()]
    nonfound = []
    for i, T in enumerate(Ts):
        assert is_tournament(T) and len(T) == 21
        S, t = random_regular_search(T, 4000, rng)
        if S is None:
            nonfound.append(i)
    print(f'{fn}: {len(Ts)} tournaments; without a found regular subset: {nonfound}; '
          f'iso to P7[C3]: {[iso_exists(Ts[i], fam21) for i in nonfound]}  [{time.time() - t0:.0f}s]', flush=True)


# Cayley tournaments of order 27
def cayley_all(elems, addf, negf):
    nz = [g for g in elems if g != elems[0]]
    reps = []; seen = set()
    for g in nz:
        if g in seen:
            continue
        seen.add(g); seen.add(negf(g)); reps.append(g)
    ix = {g: i for i, g in enumerate(elems)}
    out = []
    for bits in range(1 << len(reps)):
        S = [r if (bits >> j) & 1 else negf(r) for j, r in enumerate(reps)]
        rows = [sum(1 << ix[addf(g, s)] for s in S) for g in elems]
        out.append(from_rows(rows))
    return out


fam27 = lexprod(C3, C33)
groups = {
    'Z3^3': ([(a, b, c) for a in range(3) for b in range(3) for c in range(3)],
             lambda x, y: tuple((u + v) % 3 for u, v in zip(x, y)), lambda x: tuple((-u) % 3 for u in x)),
    'Z9xZ3': ([(a, b) for a in range(9) for b in range(3)],
              lambda x, y: ((x[0] + y[0]) % 9, (x[1] + y[1]) % 3), lambda x: ((-x[0]) % 9, (-x[1]) % 3)),
    'Z27': ([(a,) for a in range(27)], lambda x, y: ((x[0] + y[0]) % 27,), lambda x: ((-x[0]) % 27,)),
}
for nm, (el, addf, negf) in groups.items():
    Ts = cayley_all(el, addf, negf)
    nonfound = []
    tries_hist = {}
    for i, T in enumerate(Ts):
        S, t = random_regular_search(T, 400, rng)
        if S is None:
            nonfound.append(i)
        else:
            tries_hist[min(t, 10)] = tries_hist.get(min(t, 10), 0) + 1
    isos = [iso_exists(Ts[i], fam27) for i in nonfound]
    print(f'{nm}: {len(Ts)} Cayley tournaments; tries histogram {dict(sorted(tries_hist.items()))}; without regular after 400 tries: '
          f'{len(nonfound)}, of which iso to C3^[3]: {sum(isos)}  [{time.time() - t0:.0f}s]', flush=True)
    bad = [nonfound[j] for j in range(len(nonfound)) if not isos[j]]
    for i in bad[:5]:
        S, t = random_regular_search(Ts[i], 20000, rng)
        print(f'   non-family {nm}[{i}]: regular found with 20000 more tries: {S is not None}', flush=True)
