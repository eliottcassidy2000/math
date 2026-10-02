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


def zero_sum_partitions(p):
    """N_p: partitions of the quadratic residues mod p into triples with zero sum.  (For p = 3 mod 4 each class
    {+-d} has exactly one residue, and the three steps of a P_p-directed triangle are residues summing to 0, so an
    all-directed Z_p-invariant STS is such a partition plus a cyclic order on each triple: N_p 2^((p-1)/6) systems.)"""
    QR = sorted({(x * x) % p for x in range(1, p)})
    S = set(QR)
    triples = [t for t in itertools.combinations(QR, 3) if sum(t) % p == 0]
    by = {x: [t for t in triples if x in t] for x in QR}
    count = 0

    def rec(left):
        nonlocal count
        if not left:
            count += 1; return
        x = min(left)
        for t in by[x]:
            if all(y in left for y in t):
                rec(left - set(t))
    rec(frozenset(S))
    return count, len(triples)


def verify(p, S1, S2):
    """Literal check of a proposed Paley face-cyclic triangulation given by base blocks of S1 (along P_p) and S2."""
    QR = {(x * x) % p for x in range(1, p)}
    beats = lambda x, y: (y - x) % p in QR
    dclass = lambda d: min(d % p, (-d) % p)

    def directed(t):
        a, b, c = t
        if beats(a, b) and beats(b, c) and beats(c, a):
            return (a, b, c)
        if beats(a, c) and beats(c, b) and beats(b, a):
            return (a, c, b)
        return None
    for S in (S1, S2):
        assert all(directed(t) for t in S)
        cls = [dclass(d) for (a, b, c) in [tuple(t) for t in S] for d in (b - a, c - a, c - b)]
        assert sorted(cls) == list(range(1, (p - 1) // 2 + 1))
    faces1 = {tuple(sorted(((x + s) % p) for x in t)) for t in S1 for s in range(p)}
    faces2 = {tuple(sorted(((x + s) % p) for x in t)) for t in S2 for s in range(p)}
    assert not (faces1 & faces2)
    darts = Counter()
    for F, sgn in ((faces1, 1), (faces2, -1)):
        for f in F:
            c = directed(f)
            if sgn < 0:
                c = c[::-1]
            for k in range(3):
                darts[(c[k], c[(k + 1) % 3])] += 1
    assert len(darts) == p * (p - 1) and set(darts.values()) == {1}
    for f in faces1:                     # the induced tournament (edges along S1-faces) is P_p
        c = directed(f)
        assert all(beats(c[k], c[(k + 1) % 3]) for k in range(3))
    rho = {}
    for F, sgn in ((faces1, 1), (faces2, -1)):
        for f in F:
            if 0 in f:
                c = directed(f)
                if sgn < 0:
                    c = c[::-1]
                i = c.index(0)
                rho[c[(i + 1) % 3]] = c[(i + 2) % 3]
    x, L = 1, 0
    while True:
        x = rho[x]; L += 1
        if x == 1:
            break
    return L == p - 1


EXAMPLES = {   # base blocks found by the independent audit (audit24, paley_verify_example.out); re-checked here
    67: ([(0, 1, 20), (0, 3, 33), (0, 4, 28), (0, 5, 16), (0, 6, 46), (0, 7, 9), (0, 10, 32), (0, 12, 25), (0, 14, 31),
          (0, 15, 44), (0, 18, 26)],
         [(0, 1, 5), (0, 2, 15), (0, 6, 32), (0, 9, 42), (0, 10, 27), (0, 11, 29), (0, 12, 19), (0, 14, 51), (0, 20, 23),
          (0, 21, 45), (0, 28, 36)]),
    79: ([(0, 1, 56), (0, 2, 7), (0, 4, 35), (0, 8, 33), (0, 9, 29), (0, 10, 28), (0, 11, 63), (0, 12, 26), (0, 13, 34),
          (0, 17, 32), (0, 19, 41), (0, 30, 36), (0, 37, 40)],
         [(0, 1, 53), (0, 2, 12), (0, 4, 54), (0, 5, 47), (0, 6, 9), (0, 8, 59), (0, 11, 24), (0, 14, 31), (0, 15, 45),
          (0, 16, 60), (0, 18, 41), (0, 21, 57), (0, 33, 40)]),
    139: ([(0, 1, 130), (0, 2, 55), (0, 4, 110), (0, 5, 88), (0, 6, 102), (0, 7, 114), (0, 8, 30), (0, 11, 82), (0, 12, 31),
           (0, 13, 58), (0, 14, 35), (0, 15, 41), (0, 16, 50), (0, 17, 65), (0, 20, 87), (0, 24, 103), (0, 28, 70),
           (0, 38, 92), (0, 39, 66), (0, 40, 63), (0, 44, 90), (0, 59, 77), (0, 61, 64)],
          [(0, 1, 53), (0, 2, 16), (0, 4, 58), (0, 5, 84), (0, 6, 119), (0, 7, 90), (0, 9, 115), (0, 11, 74), (0, 12, 30),
           (0, 13, 104), (0, 17, 57), (0, 19, 34), (0, 21, 29), (0, 22, 45), (0, 25, 103), (0, 28, 95), (0, 31, 97),
           (0, 32, 71), (0, 37, 75), (0, 41, 88), (0, 43, 46), (0, 50, 77), (0, 59, 69)]),
    163: ([(0, 1, 75), (0, 2, 10), (0, 4, 18), (0, 5, 53), (0, 6, 68), (0, 9, 127), (0, 12, 85), (0, 15, 79), (0, 16, 42),
           (0, 21, 125), (0, 22, 82), (0, 23, 34), (0, 24, 70), (0, 25, 108), (0, 29, 49), (0, 30, 96), (0, 31, 58),
           (0, 32, 35), (0, 33, 120), (0, 37, 56), (0, 39, 116), (0, 40, 109), (0, 41, 92), (0, 44, 61), (0, 50, 57),
           (0, 52, 65), (0, 63, 91)],
          [(0, 1, 66), (0, 2, 14), (0, 4, 30), (0, 5, 16), (0, 6, 110), (0, 7, 15), (0, 9, 92), (0, 10, 50), (0, 13, 58),
           (0, 20, 90), (0, 21, 106), (0, 22, 86), (0, 24, 67), (0, 25, 125), (0, 27, 69), (0, 28, 47), (0, 29, 60),
           (0, 32, 84), (0, 33, 82), (0, 34, 89), (0, 35, 122), (0, 36, 75), (0, 37, 54), (0, 44, 62), (0, 46, 102),
           (0, 48, 51), (0, 68, 91)]),
}


if __name__ == '__main__':
    ps = [int(a) for a in sys.argv[1:]] or [7, 19, 31, 43]
    expect = {7: (2, 1, 1, 1), 19: (8, 4, 0, 0), 31: (192, 9056, 1280, 112), 43: (1024, 93696, 0, 0)}
    for p in ps:
        r = run(p)
        if p in expect:
            assert r == expect[p], (p, r)
        Np, ntr = zero_sum_partitions(p)
        assert r[0] == Np * 2 ** ((p - 1) // 6)
        print(f'  reduction: N_{p} = {Np} partitions of QR({p}) into zero-sum triples ({ntr} such triples); '
              f'systems = N_p 2^((p-1)/6) = {Np * 2 ** ((p - 1) // 6)}')
    if not sys.argv[1:]:
        print('N_67 =', zero_sum_partitions(67)[0])
        for p, (S1, S2) in EXAMPLES.items():
            ok = verify(p, S1, S2)
            print(f'p = {p} (= {p % 24} mod 24): the audit example is a genuine surface of genus {(p - 3) * (p - 4) // 12}, '
                  f'every face a directed triangle of P_{p}: {ok}')
            assert ok
    print('ALL CHECKS PASSED')
