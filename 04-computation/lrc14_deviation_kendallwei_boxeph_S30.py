#!/usr/bin/env python3
"""
S30 twin deliverables (boxeph-2026-07-16-S30):

PART 1 -- N2: THE DEVIATION IDENTITY for the invariant mean (THM-892 addendum):
    <Q_s> = (pi^2/3) M (1 + 1/P)  +  (pi^2/(3 P^2)) sum_{q | P, q > 1} J_2(q) delta_{P/q},
    delta_g = A_g - PM/(P-1)  (class imbalance),  using sum_{q|P} J_2(q) = P^2.
  The universal term is EXACT (not asymptotic); the correction is the J_2-weighted
  class imbalance of the coincidence spectrum. Referee: exact on four clusters.

PART 2 -- THE KENDALL-WEI / PERRON ATLAS per iso class, n = 4..7 (owner-fed claims
  verified): R o R counts 2-paths; iterating (Kendall-Wei) converges to the Perron
  eigenvector -- every vertex's strength defined through all others simultaneously.
  Claims: (i) lambda = 0 iff transitive (classical: adjacency nilpotent iff acyclic
  iff transitive -- proof in write-up); (ii) corr(lambda, -x) ~ 0.93 with x = score
  variance; (iii) lambda has nonzero spread inside every x-level (the all-frames
  invariant strictly refines the frame-dependent first moment).
  Bonus: SC detection per class => A000568(n) == #SC(n) (mod 2) (complementation
  involution), n = 4..7.
"""

import sys
from math import gcd, lcm, pi, sin
from itertools import permutations, combinations

sys.path.insert(0, '04-computation')

# ---------------------------------------------------------------- PART 1

def part1():
    from lrc14_invariant_mean_boxeph_S29 import (cluster_data, S_spec, Ncoin,
                                                 divisors, mu, J2, phi)
    print("=" * 74)
    print("PART 1 -- N2: the deviation identity")
    for E, s, name in [([1, 2, 3, 4, 5, 6, 60], 0, "family{1..6,60}"),
                       ([1, 2, 3, 4, 5, 36, 60], 0, "two-owner"),
                       ([12, 15, 20, 21, 28, 30, 35], 0, "balanced"),
                       ([1, 2, 3, 4, 5, 6, 120], 0, "family{1..6,120}")]:
        P, pos, sgn = cluster_data(E, s)
        M = len(pos)
        if M == 0:
            continue
        # Jordan identity
        assert sum(J2(q) for q in divisors(P)) == P * P
        spec = S_spec(pos, sgn, P)
        # class averages A_g and the closed mean (from THM-892)
        total = 0.0
        dev_term = 0.0
        Abar = P * M / (P - 1)
        for g in divisors(P):
            if g == P:
                continue
            q = P // g
            # A_g via direct spectrum (referee-grade)
            cls = [m for m in range(1, P) if gcd(m, P) == g]
            Ag = sum(spec[m] for m in cls) / len(cls)
            total += (J2(q) / 3) * Ag
            dev_term += J2(q) * (Ag - Abar)
        closed = (pi * pi / (P * P)) * total
        identity_rhs = (pi * pi / 3) * M * (1 + 1 / P) + \
                       (pi * pi / (3 * P * P)) * dev_term
        err = abs(closed - identity_rhs) / closed
        print(f"  [{name}] P={P} M={M}: <Q_s> = {closed:.6f}; "
              f"universal (pi^2/3)M(1+1/P) = {(pi*pi/3)*M*(1+1/P):.6f}; "
              f"deviation = {closed - (pi*pi/3)*M*(1+1/P):+.6f}; identity err {err:.1e}")

# ---------------------------------------------------------------- PART 2

def all_tournament_classes(n):
    """orbit-marking enumeration of iso classes; returns list of adjacency
    bitmask representatives."""
    pairs = list(combinations(range(n), 2))
    npairs = len(pairs)
    seen = bytearray(1 << npairs)
    reps = []
    perms = list(permutations(range(n)))
    pair_index = {p: i for i, p in enumerate(pairs)}
    for code in range(1 << npairs):
        if seen[code]:
            continue
        reps.append(code)
        # mark the whole orbit
        for sigma in perms:
            img = 0
            for i, (a, b) in enumerate(pairs):
                bit = (code >> i) & 1          # 1 means a -> b
                sa, sb = sigma[a], sigma[b]
                if sa < sb:
                    j = pair_index[(sa, sb)]
                    if bit:
                        img |= 1 << j
                else:
                    j = pair_index[(sb, sa)]
                    if not bit:
                        img |= 1 << j
            seen[img] = 1
        # also mark complement orbit? no -- complement is a DIFFERENT class
    return reps, pairs

def adj_from_code(code, pairs, n):
    A = [[0] * n for _ in range(n)]
    for i, (a, b) in enumerate(pairs):
        if (code >> i) & 1:
            A[a][b] = 1
        else:
            A[b][a] = 1
    return A

def perron(A, n, iters=800):
    v = [1.0] * n
    lam = 0.0
    for _ in range(iters):
        w = [sum(A[i][j] * v[j] for j in range(n)) for i in range(n)]
        s = sum(w)
        if s < 1e-300:
            return 0.0
        lam = s / sum(v)
        v = [x / s * n for x in w]
    return lam

def complement_code(code, npairs):
    return ((1 << npairs) - 1) ^ code

def canonical(code, pairs, pair_index, perms, n):
    best = None
    for sigma in perms:
        img = 0
        for i, (a, b) in enumerate(pairs):
            bit = (code >> i) & 1
            sa, sb = sigma[a], sigma[b]
            if sa < sb:
                if bit:
                    img |= 1 << pair_index[(sa, sb)]
            else:
                if not bit:
                    img |= 1 << pair_index[(sb, sa)]
        if best is None or img < best:
            best = img
    return best

def part2(nmax=7):
    print()
    print("=" * 74)
    print("PART 2 -- the Kendall-Wei / Perron atlas")
    for n in range(4, nmax + 1):
        reps, pairs = all_tournament_classes(n)
        npairs = len(pairs)
        perms = list(permutations(range(n)))
        pair_index = {p: i for i, p in enumerate(pairs)}
        rows = []
        n_sc = 0
        canon_set = set(canonical(c, pairs, pair_index, perms, n) for c in reps)
        for code in reps:
            A = adj_from_code(code, pairs, n)
            scores = [sum(r) for r in A]
            xvar = sum((s - (n - 1) / 2) ** 2 for s in scores) / n
            lam = perron(A, n)
            trans = sorted(scores) == list(range(n))
            # SC: complement isomorphic to itself?
            sc = canonical(complement_code(code, npairs), pairs, pair_index,
                           perms, n) == canonical(code, pairs, pair_index, perms, n)
            n_sc += sc
            rows.append((lam, xvar, trans, tuple(sorted(scores))))
        # (i) lambda = 0 iff transitive
        ok_i = all((abs(l) < 1e-9) == t for l, x, t, sc_ in rows)
        # (ii) correlation
        lams = [r[0] for r in rows]
        xs = [r[1] for r in rows]
        mL = sum(lams) / len(lams); mX = sum(xs) / len(xs)
        cov = sum((l - mL) * (x - mX) for l, x in zip(lams, xs))
        vL = sum((l - mL) ** 2 for l in lams) ** 0.5
        vX = sum((x - mX) ** 2 for x in xs) ** 0.5
        corr = cov / (vL * vX) if vL * vX > 0 else 0.0
        # (iii) spread within score-multiset levels
        from collections import defaultdict
        lev = defaultdict(list)
        for l, x, t, sm in rows:
            lev[sm].append(l)
        multi = {k: v for k, v in lev.items() if len(v) > 1}
        spreads = [max(v) - min(v) for v in multi.values()]
        all_spread = all(s > 1e-9 for s in spreads) if spreads else True
        print(f"  n={n}: classes {len(rows)} (A000568 ok: "
              f"{len(rows) in (4,12,56,456)}); SC = {n_sc} "
              f"(A000568 == SC mod 2: {len(rows) % 2 == n_sc % 2}); "
              f"(i) lambda=0 iff transitive: {ok_i}; "
              f"(ii) corr(lambda, -xvar) = {-corr:+.4f}; "
              f"(iii) multi-class levels {len(multi)}, all with lambda-spread > 0: "
              f"{all_spread}")

if __name__ == "__main__":
    part1()
    nmax = 7 if len(sys.argv) < 2 else int(sys.argv[1])
    part2(nmax)
    print("done")
