#!/usr/bin/env python3
"""
dihedral_quaternion_tournament_census_boxeph_S150.py  (HYP-8205)

THE DIHEDRAL QUESTION FOR TOURNAMENTS (owner S150: tournaments <-> dihedral).

For a tournament T, converse = complement = T^op (reverse all arcs).  An
ANTI-AUTOMORPHISM is a permutation sigma with sigma(T) = T^op.  Facts:
  * Aut+-(T) := Aut(T) u {anti-autos} is a group; [Aut+- : Aut] <= 2 with
    equality iff T is self-converse (= self-complementary, SC).
  * Every anti-automorphism has EVEN order (odd order would put sigma in Aut).
  * SPLIT/DIHEDRAL question: does an SC tournament always admit an INVOLUTIVE
    anti-automorphism (sigma^2 = id)?  If yes, Aut+- = Aut x| Z2 (dihedral
    type).  If some SC tournament has all anti-autos of order divisible by 4,
    call it QUATERNIONIC (non-split extension witness at the coset level).
  * The staircase bridge: a grid-symmetric tiling (isGridSym, reflection
    (x,y) -> (n-y+1, n-x+1)) realizes sigma(v) = n+1-v as an involutive
    anti-automorphism (reflecting the base path reverses it, i.e. conjugates
    T to T^op).  The verified fact "transpose-self classes are NEVER pure
    black" (n <= 7) therefore PREDICTS: no quaternionic tournaments at n <= 7.
    This census verifies the prediction independently (full n <= 6) by direct
    group computation, and verifies the circulant instance (negation) n <= 13.

OUTPUT per n <= 6: #iso classes, #SC classes, per-SC-class |Aut|, # anti-autos,
minimal anti-auto order, DIHEDRAL (has involutive anti) vs QUATERNIONIC.
Plus: circulant tournaments n = 3,5,7,9,11,13: negation x -> -x is always an
involutive anti-automorphism (proved trivially; verified).

boxeph-2026-07-20-S150.
"""

from itertools import permutations

def all_tournaments(n, pairs):
    for bits in range(1 << len(pairs)):
        yield bits

def arc_matrix(bits, n, pairs):
    adj = [[0] * n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if (bits >> k) & 1: adj[i][j] = 1
        else: adj[j][i] = 1
    return adj

def apply_perm(adj, p, n):
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if adj[i][j]: out[p[i]][p[j]] = 1
    return out

def op(adj, n):
    return [[adj[j][i] for j in range(n)] for i in range(n)]

def key(adj, n):
    return tuple(tuple(r) for r in adj)

def canon(adj, n, perms):
    best = None
    for p in perms:
        k = key(apply_perm(adj, p, n), n)
        if best is None or k < best: best = k
    return best

def perm_order(p, n):
    seen = [False] * n
    from math import lcm
    o = 1
    for s in range(n):
        if seen[s]: continue
        l, v = 0, s
        while not seen[v]:
            seen[v] = True; v = p[v]; l += 1
        o = lcm(o, l)
    return o

print("=" * 96)
print("FULL CENSUS n <= 6: SC classes and the dihedral/quaternion split")
print("=" * 96)
summary = []
for n in range(3, 7):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    perms = list(permutations(range(n)))
    seen = {}
    for bits in all_tournaments(n, pairs):
        adj = arc_matrix(bits, n, pairs)
        c = canon(adj, n, perms)
        if c not in seen:
            seen[c] = adj
    classes = list(seen.values())
    sc_rows = []
    quat = 0
    for adj in classes:
        T, Top = key(adj, n), key(op(adj, n), n)
        autos, antis = [], []
        for p in perms:
            k = key(apply_perm(adj, p, n), n)
            if k == T: autos.append(p)
            if k == Top: antis.append(p)
        if antis:
            orders = sorted({perm_order(p, n) for p in antis})
            has_inv = 2 in orders
            if not has_inv: quat += 1
            sc_rows.append((len(autos), len(antis), orders, has_inv))
    print("n=%d: classes=%d  SC=%d  quaternionic=%d" % (n, len(classes), len(sc_rows), quat))
    for aut, na, orders, inv in sorted(sc_rows):
        print("    |Aut|=%-3d #anti=%-3d anti-orders=%s  %s"
              % (aut, na, orders, "DIHEDRAL(split)" if inv else "*** QUATERNIONIC ***"))
    summary.append((n, len(classes), len(sc_rows), quat))

print("\nsummary (n, classes, SC, quaternionic):", summary)
print("A000568 check: classes should be 2, 4, 12, 56 for n=3..6")
print("prediction from never-pure-black (n<=7): quaternionic = 0 everywhere")

print("\n" + "=" * 96)
print("CIRCULANT INSTANCE: negation x -> -x is an involutive anti-automorphism (n odd)")
print("=" * 96)
for n in (3, 5, 7, 9, 11, 13):
    half = [s for s in range(1, n) if s <= (n - 1) // 2]
    import itertools
    ok_all = True
    cnt = 0
    for choice in itertools.product([1, -1], repeat=len(half)):
        S = set()
        for s, e in zip(half, choice):
            S.add(s % n if e == 1 else (-s) % n)
        # T_S: arc i -> j iff (j - i) mod n in S ; negation p(x) = -x mod n
        adj = [[1 if (j - i) % n in S else 0 for j in range(n)] for i in range(n)]
        for i in range(n):
            adj[i][i] = 0
        p = [(-x) % n for x in range(n)]
        img = apply_perm(adj, p, n)
        if key(img, n) != key(op(adj, n), n): ok_all = False
        cnt += 1
    o = perm_order([(-x) % n for x in range(n)], n)
    print("  n=%2d: all %4d circulant tournaments: negation anti-auto %s ; order(negation)=%d"
          % (n, cnt, "OK" if ok_all else "FAIL", o))
print("\n=> the dihedral group D_n = <rotation, negation> acts on every circulant")
print("   tournament with rotations as automorphisms and reflections as")
print("   anti-automorphisms; the reflection is involutive: circulants are always")
print("   DIHEDRAL-split.  The staircase reflection sigma(v) = n+1-v plays the")
print("   same role for grid-symmetric tilings (base-path dihedral).")
print("DONE.")
