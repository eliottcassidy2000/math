#!/usr/bin/env python3
"""
lrc14_signed_box_klein_S212.py

HYP-5766: (C1) THE SIGNED LOW-RELATION BOX — closed-form per-relation
contributions to the aggregated B5 supply, verification against measured
deficits, and the LARGE-V BRANCH TEST.

THE CLOSED FORM. Kill window K_q = [−(c−1), c−1] mod q, c = ⌈q/14⌉, is
SYMMETRIC, so the discrete Fourier coefficient is real and explicit:
    K̂_q(j) = sin(π j (2c−1)/q) / (q sin(π j/q))     [Dirichlet kernel]
            = k(j) + O(|j|/q),   k(j) := sin(π j/7)/(π j)   (q-free limit).
Signs of k: + for j ≡ 1..6 (mod 14), 0 at 7|j, − for j ≡ 8..13 (mod 14).

For an exact relation j (support U, Σ_{l∈U} j_l v_l = 0), its first-order
contribution to avg_q B5(S,q)/(q−1) is

    Δ(j) = [Π_{l∈U} k(j_l)] · (−1)^{|U|} · T_{|U|},
    T_u := Σ_{e=0}^{5−u} (−1)^e C(13−u, e) 7^{−e}
    (T_2 = +[1−11/7+55/49−165/343] ≈ 0.0700, T_3 ≈ +0.4898,
     T_4 = 1−9/7 ≈ −0.2857, T_5 = 1),
because the relation enters S_d for every T ⊇ U with the free 13−u
coordinates carrying K̂(0) = (2c−1)/q ≈ 1/7, and the quintic alternating sum
over d ≥ u telescopes to (−1)^u T_u. Both ±j appear (equal contributions).

TESTS:
 [1] VERIFY the closed form: predicted avg-B5 deficit = −Σ_{relations, both
     signs, height ≤ H} Δ(j) vs MEASURED (0.1221-ref and exact-κ_q-ref) on
     the corpus. Tight match ⟹ (C1) is an equality-grade decomposition.
 [2] THE LARGE-V BRANCH TEST: S = 20·{1..13} with one gcd-breaking bump
     (primitive, covering, V = 260, E3 ≈ dilated level). Prediction from the
     box: avg-B5 ≈ 0.1221 − (Schur+midpoint budget ≈ 0.05–0.09) > 0.
     If measured avg-B5 > 0: the (near-)dilated branch is a SMALL-V artifact
     and (C2) rigidity is UNNECESSARY — the aggregated theorem needs only
     [closed form + box + dispersal] at V > V₀ with modest V₀.
 [3] The per-relation ledger on the dilated family: which relations carry
     the deficit; the exact budget Σ|Δ| vs the 0.1221 slack.
"""

import random
from math import gcd, comb, sin, pi
from itertools import combinations, product

random.seed(20260712)
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
    return hist, ncl


def B5_of(S, q):
    hist, ncl = coverage_hist(S, q)
    return sum(n * sum((-1) ** d * comb(c, d) for d in range(6))
               for c, n in enumerate(hist))


def k_cont(j):
    if j % 7 == 0:
        # sin(pi j/7)=0 exactly at multiples of 7
        return 0.0
    return sin(pi * j / 7) / (pi * j)


T_u = {}
for u in range(2, 6):
    T_u[u] = sum((-1) ** e * comb(13 - u, e) / 7 ** e for e in range(0, 6 - u))


def relations(S, H, max_supp=5):
    """all exact relations (U, j) with 0<||j||inf<=H, support 2..5, BOTH signs."""
    rels = []
    coeffs = [c for c in range(-H, H + 1) if c != 0]
    for s in range(2, max_supp + 1):
        for U in combinations(range(13), s):
            vs = [S[i] for i in U]
            for j in product(coeffs, repeat=s):
                if sum(a * b for a, b in zip(j, vs)) == 0:
                    rels.append((U, j))
    return rels


def predicted_deficit(S, H):
    """- sum of Delta(j) over all relations (both signs), first order."""
    tot = 0.0
    ledger = {}
    for (U, j) in relations(S, H):
        u = len(U)
        prod = 1.0
        for a in j:
            prod *= k_cont(a)
        d = prod * ((-1) ** u) * T_u[u]
        tot += d
        key = (u, tuple(sorted(abs(a) for a in j)))
        ledger[key] = ledger.get(key, 0.0) + d
    return -tot, ledger  # deficit = -(sum of Delta)


def measured_avg(S):
    V = max(S)
    vals = []
    iid_exact = []
    for q in range(V + 1, 2 * V + 1):
        b5 = B5_of(S, q) / (q - 1)
        vals.append(b5)
        c = -(-q // 14)
        kq = (2 * c - 1) / q
        # exact-kq iid reference at 13 classes (merges rare, ignore here)
        iid_exact.append(sum((-1) ** e * comb(13, e) * kq ** e for e in range(6)))
    return sum(vals) / len(vals), sum(iid_exact) / len(iid_exact)


def main():
    print("=" * 78)
    print("HYP-5766: (C1) the signed low-relation box (klein-S212)")
    print("=" * 78)
    print(f"\nT_u factors: " + ", ".join(f"T_{u}={T_u[u]:+.4f}" for u in sorted(T_u)))
    print(f"k(j) values: " + ", ".join(f"k({j})={k_cont(j):+.4f}" for j in range(1, 8)))
    print(f"unit costs: Schur (1,1,-1): {k_cont(1)**2*k_cont(1)*(-1)**3*T_u[3]:+.6f}"
          f"   midpoint (1,-2,1): {k_cont(1)*k_cont(2)*k_cont(1)*(-1)**3*T_u[3]:+.6f}"
          f"   ratio (2,-1): {k_cont(2)*k_cont(1)*T_u[2]:+.6f}")

    corpus = [
        ([12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120], "adv-worst-120"),
        ([62, 66, 69, 102, 109, 118, 120, 126, 130, 136, 159, 185, 200], "adv-200"),
        ([9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91], "@91"),
        ([10, 46, 48, 59, 66, 148, 177, 181, 208, 213, 236, 261, 280], "adv-280"),
        ([2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26], "dilAP-2x V=26"),
        # large-V near-dilations (THE BRANCH TEST): 20*{1..13} with one bump
        ([20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260], "20*{1..13} 40->41 V=260"),
        ([20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 261], "20*{1..13} 260->261 V=261"),
        ([21, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260], "20*{1..13} 20->21 V=260"),
    ]
    corpus = [(S, n) for (S, n) in corpus
              if len(set(S)) == 13 and is_covering(S)]
    from functools import reduce
    corpus = [(S, n) for (S, n) in corpus if reduce(gcd, S) == 1]

    print(f"\n[1]/[2] closed form vs measured (H=2 and H=3 predictions):")
    print(f"{'instance':>26} {'V':>4} {'predH2':>8} {'predH3':>8} "
          f"{'meas avg':>9} {'iid(kq)avg':>10} {'meas deficit':>12} {'avg>0':>6}")
    for S, name in corpus:
        V = max(S)
        pred2, _ = predicted_deficit(S, 2)
        pred3, ledger = predicted_deficit(S, 3)
        meas, iidref = measured_avg(S)
        print(f"{name:>26} {V:>4} {pred2:>8.4f} {pred3:>8.4f} "
              f"{meas:>9.4f} {iidref:>10.4f} {iidref - meas:>12.4f} "
              f"{'YES' if meas > 0 else 'NO':>6}")

    print("\n[3] per-relation-type ledger on the V=260 near-dilation (20*{..} 40->41):")
    S = [20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260]
    pred, ledger = predicted_deficit(S, 3)
    for key, val in sorted(ledger.items(), key=lambda kv: kv[1])[:10]:
        print(f"   support={key[0]} |j|={key[1]}: total {val:+.5f}")
    print(f"   TOTAL predicted deficit: {pred:+.4f} vs slack 0.1221")

    print("\n[4] does the CERTIFICATE still fire on the V=260 branch? (max B5 check)")
    for S, name in [([20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260], "40->41"),
                    ([21, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260], "20->21")]:
        V = max(S)
        best = None
        npos = 0
        for q in range(V + 1, 2 * V + 1):
            b5 = B5_of(S, q)
            if b5 > 0:
                npos += 1
            if best is None or b5 / (q - 1) > best[1]:
                best = (q, b5 / (q - 1), b5)
        print(f"   {name}: #q with B5>0: {npos}/{V}; best q={best[0]} B5/q={best[1]:.4f} (B5={best[2]})")

    print("\nVERDICT: if the large-V near-dilations have measured avg-B5 > 0 and the")
    print("closed form tracks the measured deficit, (C1) is validated and the")
    print("(near-)dilated branch is a SMALL-V artifact -- (C2) rigidity unnecessary;")
    print("the aggregated theorem = [closed form + box + dispersal] at V > V0.")
    print("DONE.")


if __name__ == '__main__':
    main()
