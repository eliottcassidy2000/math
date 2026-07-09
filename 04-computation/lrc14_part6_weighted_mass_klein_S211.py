#!/usr/bin/env python3
"""
lrc14_part6_weighted_mass_klein_S211.py

HYP-5761 corrected: THM-671 part 6 via WEIGHTED resonance mass (the binary
resolved/bad dichotomy of the first attempt is vacuous — every modulus has
thousands of low-height representations q = j·v, each of tiny weight; the
right object is the MASS).

Definitions (support ≤ 5 throughout — the quintic truncation sees only
≤5-tuples; weights w(j) = Π_l min(1/7, 1/(2|j_l|)) = the K̂ envelope):

  E_H(S)      = Σ_{j exact: j·v = 0, 0<‖j‖∞≤H} w(j)      [hits EVERY modulus]
  R_H(S, q)   = Σ_{j: j·v = m q, m≠0, 0<‖j‖∞≤H} w(j)      [modulus-specific]

LEMMA A' (mass dispersal; provable): Σ_{q∈(V,2V]} R_H(S,q) ≤ W₁(H) :=
  Σ_j (#divisors of |j·v| in range ≤ ‖j‖₁) · w(j) — V-INDEPENDENT.
  ⟹ avg_q R_H ≤ W₁(H)/V ⟹ for any δ, ≥ (1 − W₁/(δV)) of moduli have R_H ≤ δ.

LEMMA B' (deviation control; tested here): deficit(q) := 0.1221 − B5(S,q)/(q−1)
  should be predicted by the visible mass: deficit ≈ c·(E_H(S) + R_H(S,q)) + tail.
  Test: per-instance correlation/fit of deficit(q) vs R_H(q); cross-instance
  fit of the intercept vs E_H(S).

LEMMA C' (branch data): min over low-mass moduli of B5/(q−1) vs E_H(S);
  the E₀ threshold; are the over-budget covering sets near-dilated-intervals?
"""

import random
from math import gcd, comb
from itertools import combinations, product

random.seed(20260711)
QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


def coverage_hist(S, q):
    def r_safe(r):
        return 14 * r >= q and 14 * (q - r) >= q
    cov = bytearray(q)
    seen = set()
    ncl = 0
    for v in S:
        r = v % q
        key = min(r, (q - r) % q)
        if key in seen:
            continue
        seen.add(key)
        ncl += 1
        if r == 0:
            for p in range(q):
                cov[p] += 1
            continue
        g = gcd(r, q)
        rr, qq = r // g, q // g
        inv = pow(rr, -1, qq)
        for m in range(qq):
            s = m * g
            if not r_safe(s):
                p0 = (m * inv) % qq
                for t in range(g):
                    cov[p0 + t * qq] += 1
    hist = [0] * (ncl + 1)
    for p in range(1, q):
        hist[cov[p]] += 1
    return hist


def B5_of(S, q):
    hist = coverage_hist(S, q)
    return sum(n * sum((-1) ** d * comb(c, d) for d in range(6))
               for c, n in enumerate(hist))


def w_of(j):
    w = 1.0
    for a in j:
        w *= min(1 / 7, 1 / (2 * abs(a)))
    return w


def masses(S, V, H, max_supp=5):
    """returns (E_H, R: dict q -> mass, W1_actual = sum_q R[q])."""
    EH = 0.0
    R = {}
    coeffs = [c for c in range(-H, H + 1) if c != 0]
    for s in range(2, max_supp + 1):
        for U in combinations(range(13), s):
            vs = [S[i] for i in U]
            for j in product(coeffs, repeat=s):
                jv = sum(a * b for a, b in zip(j, vs))
                if jv == 0:
                    EH += w_of(j)
                    continue
                a = abs(jv)
                if a <= V:
                    continue
                mlo = max(1, -(-a // (2 * V)))
                mhi = a // (V + 1)
                for m in range(mlo, mhi + 1):
                    if a % m == 0:
                        q = a // m
                        if V < q <= 2 * V:
                            R[q] = R.get(q, 0.0) + w_of(j)
    return EH / 2, R  # /2 for ±j symmetry on the exact part


def gen_instance(V, style='slowheavy'):
    P = random.choice([(8, 9, 10, 12), (7, 9, 10, 11, 12), (11, 12, 13),
                       (10, 11, 12, 13), (9, 11, 13)])
    k = 13 - len(P)
    if k < 8:
        return None
    L = {V}
    missed = [q for q in QS if not any(p % q == 0 for p in P)]
    for q in missed:
        if any(u % q == 0 for u in L):
            continue
        lo, hi = -(-14 // q), V // q
        if lo > hi:
            return None
        L.add(q * random.randint(lo, hi))
    if style == 'slowheavy':
        for _ in range(3):
            if len(L) < k:
                L.add(random.randint(max(14, V // 14 + 1), max(16, 9 * V // 14 - 1)))
    while len(L) < k:
        L.add(random.randint(14, V))
    S = sorted(set(P) | L)
    if len(S) == 13 and is_covering(S):
        return S
    return None


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return float('nan')
    mx, my = sum(xs) / n, sum(ys) / n
    sx = (sum((a - mx) ** 2 for a in xs)) ** 0.5
    sy = (sum((b - my) ** 2 for b in ys)) ** 0.5
    if sx == 0 or sy == 0:
        return float('nan')
    return sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / (sx * sy)


def main():
    print("=" * 78)
    print("HYP-5761 (corrected): weighted resonance mass and the B5 deficit (klein-S211)")
    print("=" * 78)
    H = 2

    corpus = []
    for V in (120, 200):
        got = 0
        while got < 3:
            S = gen_instance(V)
            if S:
                corpus.append((S, f"rand-V{V}-{got}"))
                got += 1
    corpus += [
        ([9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91], "@91"),
        ([12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120], "adv-worst"),
        ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26], "dilAP-2x"),
        ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 28], "near-dilAP"),
    ]
    corpus = [(S, n) for (S, n) in corpus if is_covering(S) and len(set(S)) == 13]

    print(f"\nper-instance summary (H={H}, support<=5):")
    print(f"{'name':>14} {'V':>4} {'E_H':>7} {'sumR=W1act':>10} {'W1act/V':>8} "
          f"{'corr(def,R)':>11} {'minB5/q lowR':>12} {'minB5/q all':>11} {'E_H+minR->pred':>14}")
    table = []
    for S, name in corpus:
        V = max(S)
        EH, R = masses(S, V, H)
        defs, Rs = [], []
        lowmass_b5 = []
        min_all = None
        for q in range(V + 1, 2 * V + 1):
            b5 = B5_of(S, q) / (q - 1)
            rq = R.get(q, 0.0)
            defs.append(0.1221 - b5)
            Rs.append(rq)
            if rq <= 0.02:
                lowmass_b5.append(b5)
            min_all = b5 if min_all is None or b5 < min_all else min_all
        corr = pearson(Rs, defs)
        W1act = sum(R.values())
        minlow = min(lowmass_b5) if lowmass_b5 else float('nan')
        table.append((name, V, EH, W1act, corr, minlow, min_all))
        print(f"{name:>14} {V:>4} {EH:>7.4f} {W1act:>10.3f} {W1act/V:>8.4f} "
              f"{corr:>11.3f} {minlow:>12.4f} {min_all:>11.4f}")

    print("\n[Lemma A' verdict] sumR/V = the avg per-modulus m!=0 mass -- small and")
    print("V-shrinking => low-mass moduli abound (Markov).")

    print("\n[Lemma B'/C' verdict] intercept check: min-B5-over-low-mass-moduli vs E_H:")
    for name, V, EH, W1, corr, minlow, minall in sorted(table, key=lambda t: t[2]):
        print(f"   E_H={EH:>7.4f}  min-lowR-B5/q={minlow:>8.4f}  {name}")
    print("\n   If min-lowR-B5/q ≈ 0.1221 - c1*E_H with a stable c1, Lemma B' holds and")
    print("   E0 = 0.1221/c1 is the branch threshold for Lemma C'.")
    xs = [t[2] for t in table if t[5] == t[5]]
    ys = [0.1221 - t[5] for t in table if t[5] == t[5]]
    if len(xs) >= 3:
        n = len(xs)
        mx, my = sum(xs) / n, sum(ys) / n
        var = sum((a - mx) ** 2 for a in xs)
        c1 = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / var if var else float('nan')
        print(f"   fit: deficit_at_lowR ≈ {my - c1 * mx:+.4f} + {c1:.3f}·E_H   "
          f"(r = {pearson(xs, ys):.3f})")

    print("\nDONE.")


if __name__ == '__main__':
    main()
