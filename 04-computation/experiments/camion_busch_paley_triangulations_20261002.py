"""Paley face-cyclic triangulations (companion to 05-knowledge/results/camion_busch_gaps_polyhedra_collatz_20261002.md,
section 6.6).

For a prime p = 1 (mod 6), p = 3 (mod 4), let P_p be the Paley tournament on Z_p (i -> j iff j - i is a nonzero square).
We look for Z_p-invariant Steiner triple systems S on Z_p (translates of (p-1)/6 base blocks {0, a, b} whose difference
classes {+-a, +-b, +-(b-a)} partition the (p-1)/2 classes) all of whose blocks are directed triangles of P_p.  For a pair
{S1, S2} with no common block, orient S1-faces along P_p and S2-faces against it.  Every edge lies in one block of each
system and is traversed oppositely by the two faces, so the faces glue to an oriented pseudo-surface on which P_p is the
tournament of the face orientation and every face is a directed triangle.  It is a genuine surface iff the rotation at
vertex 0 (and so, by translation, at every vertex) is a single (p-1)-cycle.

Usage: python camion_busch_paley_triangulations_20261002.py [p ...]   (default: 7 19 31 43)
"""
import itertools
import sys
from collections import Counter


def run(p):
    QR = {(x * x) % p for x in range(1, p)}
    beats = lambda x, y: (y - x) % p in QR

    def directed(t):
        a, b, c = t
        if beats(a, b) and beats(b, c) and beats(c, a):
            return (a, b, c)
        if beats(a, c) and beats(c, b) and beats(b, a):
            return (a, c, b)
        return None

    dclass = lambda d: min(d % p, (-d) % p)
    orbit = lambda t: frozenset(tuple(sorted((x + s) % p for x in t)) for s in range(p))
    base, seen = [], set()
    for a, b in itertools.combinations(range(1, p), 2):
        cls = frozenset({dclass(a), dclass(b), dclass(b - a)})
        if len(cls) == 3 and directed((0, a, b)):
            o = orbit((0, a, b))
            if o not in seen:
                seen.add(o); base.append((cls, o))
    allcls = frozenset(range(1, (p - 1) // 2 + 1))
    systems = set()

    def rec(used, chosen):
        if used == allcls:
            systems.add(frozenset(chosen)); return
        first = min(allcls - used)
        for j, (cls, o) in enumerate(base):
            if first in cls and not (cls & used):
                rec(used | cls, chosen + [j])
    rec(frozenset(), [])
    systems = sorted(sorted(s) for s in systems)
    faces = [frozenset().union(*[base[j][1] for j in s]) for s in systems]

    def rotation_cycles(F1, F2):
        rho = {}
        for f in F1 | F2:
            if 0 not in f:
                continue
            cyc = directed(f)
            if f in F2:
                cyc = cyc[::-1]
            i = cyc.index(0)
            rho[cyc[(i + 1) % 3]] = cyc[(i + 2) % 3]
        assert len(rho) == p - 1 and sorted(rho.values()) == sorted(rho)
        seen_, n = set(), 0
        for s in rho:
            if s in seen_:
                continue
            n += 1; x = s
            while x not in seen_:
                seen_.add(x); x = rho[x]
        return n

    stats, sols = Counter(), []
    for i, j in itertools.combinations(range(len(systems)), 2):
        F1, F2 = faces[i], faces[j]
        if F1 & F2:
            continue
        # every pair of points lies in one block of each system (Steiner property) -> each edge in exactly 2 faces
        nc = rotation_cycles(F1, F2)
        stats[nc] += 1
        if nc == 1:
            sols.append(F1 | F2)
    QRl = sorted(QR)
    canon = {min(tuple(sorted(tuple(sorted((q * v) % p for v in f)) for f in F)) for q in QRl) for F in sols}
    genus = (p - 3) * (p - 4) // 12
    print(f'p = {p}: {len(base)} orbits of usable P_p-directed triangles; {len(systems)} Z_p-invariant Steiner triple '
          f'systems with P_p-directed blocks; {sum(stats.values())} disjoint pairs, rotation-cycle counts '
          f'{dict(sorted(stats.items()))}; genuine surfaces {len(sols)} (genus {genus}), {len(canon)} up to multipliers')
    return len(systems), sum(stats.values()), len(sols), len(canon)


if __name__ == '__main__':
    ps = [int(a) for a in sys.argv[1:]] or [7, 19, 31, 43]
    expect = {7: (2, 1, 1, 1), 19: (8, 4, 0, 0), 31: (192, 9056, 1280, 112), 43: (1024, 93696, 0, 0)}
    for p in ps:
        r = run(p)
        if p in expect:
            assert r == expect[p], (p, r)
    print('ALL CHECKS PASSED')
