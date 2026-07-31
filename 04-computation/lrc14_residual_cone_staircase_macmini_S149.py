#!/usr/bin/env python3
"""THM-2589 referee — the residual cone and the triangular staircase
(mac-mini-2026-07-27-S149).

KEY: THM-2588's Lemma 2 (snap-fold) is RATIO-FREE.  So for any family V of
m distinct speeds, sorting v_(1) > ... > v_(m) and writing the consecutive
ratios rho_i = v_(i)/v_(i+1), the cascade bound holds UNCONDITIONALLY:

    M(V) >= B(V) := max_{k=0..m-1} [ F(m-k) - sum_{i<=k} 1/(2(rho_i+1)) ]

with F(j) = 1/(j+1) the settled LRC floor for j <= 13 speeds (F(1) = 1/2
exact for a single speed).  Consequences for a 13-speed GAP family
(1/14 < M < 3/41):

  CONE: for every k = 1..12,  sum_{i<=k} 1/(2(rho_i+1)) > delta_k,
        delta_k := 1/(14-k) - 3/41.
  STAIRCASE: delta_k - delta_{k-1} = 1/((14-k)(15-k)) = 1/(2 C(15-k,2)):
        the marginal ratio profile is rho*_k = C(15-k,2) - 1 — TRIANGULAR
        NUMBERS (edges of complete graphs on the residual runner count).
        k=1 rung modulus 13*14 = 182 = the deep-well far speed; and
        rho*_1 + 1 = 91 = C(14,2) = the 91-line.
  PREFIX-MIN: min(rho_1..rho_k) < c_k := k/(2 delta_k)  (c_1 = 533/4, ...,
        c_12 = 492/35 < 15).
  MASS: mu(V) := sum_{i=1}^{12} 1/(2(rho_i+1)) > 35/82.
  RATIO-14: if every consecutive ratio >= 14 then M >= 3/41 (12/30 = 2/5
        <= 35/82 fails? no: 2/5 < 35/82 CLEARS); ratio 13 does NOT clear
        (12/28 = 3/7 > 35/82): sharp at 14.

Referee parts:
 (1) unconditional bound B(V) <= exact M on random families (sizes 4..8
     exhaustively random, size 13 samples), including CLUSTERED ones;
 (2) exact delta/telescoping/staircase identities; prefix-min table;
 (3) canon families: AP {1..13}, GW {1..11,13,24}, deep well {1..12,182}:
     exact M, cone membership consistency;
 (4) sharpness probes: staircase-boundary towers close at predicted depth,
     just-inside towers do not close by cascade alone;
 (5) ratio-14 corollary exact; mass corollary exact;
 (6) ladder transversality: n(n+1) is never a prime power (n >= 2);
     escape-window/staircase ordering 48 < 90 = rho*_1.
"""

from fractions import Fraction as F
from itertools import combinations
from math import floor, gcd
import random

GAP_HI = F(3, 41)


def dist(x: F) -> F:
    fr = x - floor(x)
    return min(fr, 1 - fr)


def exact_M(W):
    W = sorted(W)
    g = 0
    for w in W:
        g = gcd(g, w)
    W = [w // g for w in W]
    mods = set()
    for u, v in combinations(W, 2):
        mods.add(u + v)
        mods.add(v - u)
    for w in W:
        mods.add(2 * w)
    mods.discard(0)
    best = F(0)
    for q in mods:
        for k in range(1, q):
            m = min(min((w * k) % q, q - (w * k) % q) for w in W)
            if F(m, q) > best:
                best = F(m, q)
    return best


def cascade_bound(V):
    """Unconditional: B(V) = max over peel depth k of floor(m-k) - losses."""
    V = sorted(V, reverse=True)
    m = len(V)
    best = F(1, m + 1)                      # k = 0: the settled floor itself
    loss = F(0)
    for k in range(1, m):
        h, w = V[k - 1], V[k]
        loss += F(w, 2 * (h + w))           # = 1/(2(rho_k+1))
        best = max(best, F(1, m - k + 1) - loss)
    return best


def ratios(V):
    V = sorted(V, reverse=True)
    return [F(V[i], V[i + 1]) for i in range(len(V) - 1)]


def main():
    rng = random.Random(149)

    print("== (1) unconditional bound B(V) <= exact M ==")
    tested = 0
    worst_slack = None
    for m in (4, 5, 6, 7, 8):
        for _ in range(60):
            style = rng.random()
            if style < 0.4:                          # clustered
                base = rng.randint(3, 40)
                V = [base]
                while len(V) < m:
                    V.append(V[-1] + rng.randint(1, max(1, V[-1] // 3)))
            elif style < 0.7:                        # mixed tower
                V = [rng.randint(1, 6)]
                while len(V) < m:
                    V.append(V[-1] * rng.randint(2, 9) + rng.randint(0, 3))
            else:                                    # random
                V = rng.sample(range(1, 90), m)
            V = sorted(set(V))
            if len(V) < m:
                continue
            B = cascade_bound(V)
            M = exact_M(V)
            assert B <= M, (V, B, M)
            tested += 1
            s = M - B
            if worst_slack is None or s < worst_slack[0]:
                worst_slack = (s, V)
    for _ in range(8):                               # size-13 spot checks
        V = sorted(rng.sample(range(1, 46), 13))
        B = cascade_bound(V)
        M = exact_M(V)
        assert B <= M, (V, B, M)
        tested += 1
    print(f"  {tested} families sizes 4..8 + 13 (clustered/mixed/random): "
          f"B(V) <= M(V) always; tightest slack {worst_slack[0]} on {worst_slack[1]}")

    print("\n== (2) exact staircase identities ==")
    deltas = [F(1, 14 - k) - GAP_HI for k in range(0, 13)]
    for k in range(1, 13):
        inc = deltas[k] - deltas[k - 1]
        rect = F(1, (14 - k) * (15 - k))
        assert inc == rect, (k, inc, rect)
    print("  delta_k - delta_{k-1} = 1/((14-k)(15-k)) exact, k = 1..12  OK")
    C = lambda n, r: F(1) * __import__('math').comb(n, r)
    star = [(k, (14 - k) * (15 - k) // 2 - 1) for k in range(1, 13)]
    assert all(F(s) == C(15 - k, 2) - 1 for k, s in star)
    print("  staircase rho*_k = C(15-k,2) - 1 =",
          ", ".join(str(s) for _, s in star))
    assert star[0][1] == 90 and 2 * (90 + 1) == 182
    print("  rho*_1 + 1 = 91 = C(14,2) [the 91-line]; k=1 rung modulus "
          "13*14 = 182 = the deep-well far speed  OK")
    print("  prefix-min ceilings c_k = k/(2 delta_k):")
    print("   ", ", ".join(f"k={k}: {F(k, 2 * deltas[k])}"
                           f"~{float(F(k, 2 * deltas[k])):.2f}"
                           for k in range(1, 13)))

    print("\n== (3) canon families ==")
    for name, W in [("AP {1..13}", list(range(1, 14))),
                    ("GW {1..11,13,24}", list(range(1, 12)) + [13, 24]),
                    ("deep well {1..12,182}", list(range(1, 13)) + [182])]:
        M = exact_M(W)
        B = cascade_bound(W)
        mu = sum(F(1, 2 * (r + 1)) for r in ratios(W))
        in_gap = F(1, 14) < M < GAP_HI
        cone_ok = all(
            sum(F(1, 2 * (r + 1)) for r in ratios(W)[:k]) > deltas[k]
            for k in range(1, 13))
        print(f"  {name}: M = {M}, B = {B}, mu = {mu} ~ {float(mu):.4f}, "
              f"in gap: {in_gap}, satisfies full cone system: {cone_ok}")
        assert not in_gap
        assert B <= M

    print("\n== (4) sharpness probes at the staircase boundary ==")
    # tower with ratios ~ rho*_k + 1 (just above): must close by cascade
    V = [1]
    for k in range(12, 0, -1):                     # build bottom-up
        V.append(V[-1] * ((14 - k) * (15 - k) // 2 + 1))
    V = sorted(set(V), reverse=True)
    B = cascade_bound(V)
    print(f"  just-ABOVE-staircase tower (13 speeds, v_max ~ 10^{len(str(V[0]))-1}): "
          f"B = {float(B):.6f} >= 3/41: {B >= GAP_HI}")
    assert B >= GAP_HI
    # tower with all ratios = 13 (just inside cone at k=12): must NOT close
    V = [13 ** i for i in range(13)]
    B = cascade_bound(V)
    print(f"  all-ratios-13 tower: B = {float(B):.6f} >= 3/41: {B >= GAP_HI} "
          f"(cascade alone must NOT close it)")
    assert B < GAP_HI
    # all ratios = 14: must close
    V = [14 ** i for i in range(13)]
    B = cascade_bound(V)
    print(f"  all-ratios-14 tower: B = {float(B):.6f} >= 3/41: {B >= GAP_HI}")
    assert B >= GAP_HI

    print("\n== (5) ratio-14 + mass corollaries, exact ==")
    assert F(12, 30) < F(35, 82) and F(12, 28) > F(35, 82)
    print("  12/(2*15) = 2/5 < 35/82 < 3/7 = 12/(2*14): ratio >= 14 closes "
          "(k = 12 peel to a single speed, floor 1/2); ratio 13 does not: SHARP.")
    print("  MASS: every gap family has mu(V) > 35/82 (k = 12 cone line); "
          "every gap family has some consecutive ratio < 492/35 < 15.")

    print("\n== (6) ladder transversality ==")
    for k in range(1, 13):
        n = 14 - k
        v = n * (n + 1)
        f = {}
        x = v
        d = 2
        while d * d <= x:
            while x % d == 0:
                f[d] = f.get(d, 0) + 1
                x //= d
            d += 1
        if x > 1:
            f[x] = f.get(x, 0) + 1
        assert len(f) >= 2, (k, v)
    print("  rung moduli (14-k)(15-k), k=1..12: all have >= 2 prime factors "
          "(consecutive integers are coprime) — never prime powers: the")
    print("  cascade rungs NEVER sit on unabsorbed-prime-power exits (49 = 7^2);")
    print("  ordering: escape window ceiling 48 < 90 = rho*_1 < 533/4 = c_1.")

    print("\nALL PARTS PASS — the residual cone + triangular staircase referee.")


if __name__ == "__main__":
    main()
