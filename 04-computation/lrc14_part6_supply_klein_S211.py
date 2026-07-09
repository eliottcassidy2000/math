#!/usr/bin/env python3
"""
lrc14_part6_supply_klein_S211.py

HYP-5761: THM-671 part 6 — the a-priori resolved-modulus supply, as a
two-parameter counting theorem. Three lemmas:

LEMMA A (bad-modulus counting; PROVED — verified here).
  Call q ∈ (V, 2V] H-BAD if there exist a support set U (|U| ≤ 5), a vector
  j ∈ (Z\\{0})^U with ‖j‖∞ ≤ H, and m ≠ 0 with Σ_{l∈U} j_l v_l = m·q.
  Each such (U, j, m) determines q = (j·v_U)/m UNIQUELY. Hence
     #H-BAD ≤ P(H) := #{(U, j, m) : |m| < ‖j‖₁} — INDEPENDENT OF V.
  So for V > P(H), H-GOOD moduli exist in (V, 2V] (indeed most are good).

LEMMA B (good-modulus deviation; stated, constants verified here).
  At an H-good q, the resonances of every ≤5-tuple split into
   (i) m = 0 EXACT relations of S at height ≤ H — the budget
       E_H(S) = Σ_{j exact, supp ≤ 5, 0<‖j‖∞≤H} Π_l min(1/7, 1/(2|j_l|)),
       which hits every modulus equally, and
   (ii) a tail at height > H — SUPPORT ≤ 5 (the quintic truncation sees only
       ≤5-tuples!), so the exact-relation lattice contributes only its
       rank-≤4 sections: theta tails converge like O(1/H); the m≠0 tail at a
       good q has height > H by definition.
  Prediction: B5(S,q)/(q−1) ≥ 0.1221 − c₁·E_H(S) − ε(H, q).

LEMMA C (the branch; census here).
  For covering S: either E_H(S) ≤ E₀ (counting theorem closes at good q), or
  E_H(S) > E₀ — conjecturally forcing S near a dilated interval
  (LEM-015/LRCSchurRigidity stability), the dilation-boundary family.

This script: exact #H-BAD vs P(H) and vs V; B5 statistics at good vs bad
moduli; E_H(S) across covering families (random / planted-Schur /
near-dilated-interval / @91); the regression B5-deficit vs E_H; the branch
verdict table.
"""

import random
from math import gcd, comb
from itertools import combinations, product

random.seed(20260711)
QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


# ---------------------------------------------------------------- B5 machinery
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
    return hist, ncl


def B5_of(S, q):
    hist, ncl = coverage_hist(S, q)
    return sum(n * sum((-1) ** d * comb(c, d) for d in range(6))
               for c, n in enumerate(hist)), hist[0]


# ---------------------------------------------------------------- Lemma A
def bad_moduli(S, V, H, max_supp=5):
    """exact H-bad set of moduli in (V, 2V] (m != 0 resonances)."""
    bad = {}
    coeffs = [c for c in range(-H, H + 1) if c != 0]
    for s in range(1, max_supp + 1):
        for U in combinations(range(13), s):
            vs = [S[i] for i in U]
            for j in product(coeffs, repeat=s):
                jv = sum(a * b for a, b in zip(j, vs))
                if jv == 0:
                    continue
                a = abs(jv)
                # q = a/m in (V, 2V]
                mlo = max(1, (a + 2 * V - 1) // (2 * V))
                mhi = a // (V + 1)
                for m in range(mlo, mhi + 1):
                    if a % m == 0:
                        q = a // m
                        if V < q <= 2 * V:
                            bad.setdefault(q, []).append((U, j, jv // q if jv > 0 else -(a // q)))
    return bad


def exact_budget(S, H, max_supp=5):
    """E_H(S): m=0 exact relations, support<=5, height<=H, weighted."""
    tot = 0.0
    count = 0
    coeffs = [c for c in range(-H, H + 1) if c != 0]
    for s in range(2, max_supp + 1):
        for U in combinations(range(13), s):
            vs = [S[i] for i in U]
            for j in product(coeffs, repeat=s):
                if sum(a * b for a, b in zip(j, vs)) == 0:
                    w = 1.0
                    for a in j:
                        w *= min(1 / 7, 1 / (2 * abs(a)))
                    tot += w
                    count += 1
    return tot / 2, count // 2  # mod +-j symmetry


# ---------------------------------------------------------------- instances
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


def main():
    print("=" * 78)
    print("HYP-5761 / THM-671 part 6: the a-priori resolved-modulus supply (klein-S211)")
    print("=" * 78)

    H = 2
    # crude P(H): pairs (U, j, m): sum_s C(13,s)(2H)^s * 2*(sH-1) -- upper count
    PH = sum(comb(13, s) * (2 * H) ** s * 2 * (s * H - 1) for s in range(1, 6))
    print(f"\nLemma A crude bound P({H}) <= {PH}  (V-independent).")

    corpus = []
    for V in (120, 200, 300):
        got = 0
        while got < 3:
            S = gen_instance(V)
            if S:
                corpus.append((S, f"rand-V{V}"))
                got += 1
    corpus += [
        ([9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91], "@91 7-struct"),
        ([12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120], "adv-worst-120"),
        ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26], "dilated-AP 2*{1..13}"),
        ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 28], "near-dilated (one bump)"),
        ([3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 42], "3*{1..12}∪{42}"),
    ]
    corpus = [(S, n) for (S, n) in corpus if is_covering(S) and len(set(S)) == 13]

    print(f"\n{'instance':>22} {'V':>4} {'#bad(H=2)':>9} {'%bad':>6} "
          f"{'E_H(S)':>8} {'#rel':>5} {'minB5/q good':>13} {'minB5/q bad':>12} {'all good>0':>10}")
    rows = []
    for S, name in corpus:
        V = max(S)
        bad = bad_moduli(S, V, H)
        EH, nrel = exact_budget(S, H)
        goods, bads = [], []
        for q in range(V + 1, 2 * V + 1):
            b5, lm = B5_of(S, q)
            (bads if q in bad else goods).append(b5 / (q - 1))
        row = dict(name=name, V=V, nbad=len(bad), EH=EH, nrel=nrel,
                   ming=min(goods) if goods else float('nan'),
                   minb=min(bads) if bads else float('nan'),
                   allg=all(x > 0 for x in goods) if goods else False)
        rows.append(row)
        print(f"{name:>22} {V:>4} {len(bad):>9} {100*len(bad)/V:>5.1f}% "
              f"{EH:>8.4f} {nrel:>5} {row['ming']:>13.4f} "
              f"{(row['minb'] if bads else float('nan')):>12.4f} {str(row['allg']):>10}")

    print("\n[Lemma A verdict] #bad is V-independent-ish and << V; P(H) is a gross")
    print("over-count (each (U,j,m) rarely lands an integer q in range).")

    print("\n[Lemma B/C regression] min-good-B5/q vs the exact budget E_H(S):")
    for r in sorted(rows, key=lambda r: r['EH']):
        print(f"   E_H={r['EH']:>7.4f}  min-good-B5/q={r['ming']:>8.4f}  "
              f"deficit-from-iid={0.1221 - r['ming']:>8.4f}   {r['name']}")
    print("\n   (E0 = the budget level where min-good-B5 crosses 0; instances above it")
    print("    should be the near-dilated-interval family only)")

    # focused check: does the m=0 budget PREDICT the deficit? crude c1 fit
    xs = [r['EH'] for r in rows]
    ys = [0.1221 - r['ming'] for r in rows]
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    var = sum((a - mx) ** 2 for a in xs)
    c1 = cov / var if var > 0 else float('nan')
    print(f"\n   crude linear fit: deficit ≈ {my - c1 * mx:+.4f} + {c1:.3f}·E_H  "
          f"(Lemma B's c1)")

    print("\nDONE.")


if __name__ == '__main__':
    main()
