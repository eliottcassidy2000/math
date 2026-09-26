#!/usr/bin/env python3
"""procgen_gates_20260925_expsum.py -- PART A of the cycle-gate equidistribution lane.

A1  certification of the float DP for S(h) (brute force, mpmath at 40 digits);
A2  the counting identity N = (1/M) sum_h S(h), exactly (finite-field image of Z[zeta_M]) and in floats;
A3  the census N(p,a) = #{w : |G| divides c_w} for q = 3, p <= 40, every a, by two independent
    methods (sorted meet-in-the-middle join of residues; orbit iteration over the size range),
    against the main term C(p,a)/|G|; SHEET control (the 3n-1 cycles on the side 3^a > 2^p);
A4  the DRIFT control q = 5 (gates 2^p - 5^a), p <= 40.
Every check raises on failure.  stdout = results, stderr = timing/memory.
"""
import math
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from procgen_gates_20260925_core import (S_brute, S_dp, S_dp_exact_all, S_dp_mpmath, carry, census_mitm,  # noqa: E402
                                         census_range, gate, lean_malloc, lyndon_count, report_mem, rotate_carry,
                                         words)

P_CENSUS = 40


def a1_certify():
    print("=" * 100)
    print("A1. Float DP for S(h) = sum_w e(h c_w/|G|): certification")
    print("=" * 100)
    rng = random.Random(20260925)
    worst_b = 0.0
    nb = 0
    for q in (3, 5):
        for p in range(1, 13):
            for a in range(1, p + 1):
                M = abs(gate(p, a, q))
                hs = sorted({1, 2, 3, M - 1, M // 2, rng.randrange(M)} if M > 1 else {0})
                d = S_dp(p, a, hs, q)
                for k, h in enumerate(hs):
                    err = abs(d[k] - S_brute(p, a, h, q))
                    worst_b = max(worst_b, err)
                    nb += 1
    print(f"  DP vs brute-force enumeration: {nb} (clock, h) pairs, q in {{3,5}}, p <= 12, max |diff| = {worst_b:.2e}")
    assert worst_b < 1e-9
    worst_m, nm = 0.0, 0
    for (p, a, q) in [(19, 12, 3), (27, 17, 3), (30, 19, 3), (40, 25, 3), (46, 29, 3), (65, 41, 3), (84, 53, 3),
                      (40, 17, 5), (30, 13, 5)]:
        M = abs(gate(p, a, q))
        hs = [1, 2, M // 2, M - 3] + [rng.randrange(1, M) for _ in range(3)]
        d = S_dp(p, a, hs, q)
        for k, h in enumerate(hs):
            ref = S_dp_mpmath(p, a, h, q, dps=40)
            rel = abs(d[k] - ref) / math.comb(p, a)
            worst_m = max(worst_m, rel)
            nm += 1
    print(f"  float DP vs 40-digit mpmath DP: {nm} pairs up to p = 84, max |diff|/C(p,a) = {worst_m:.2e}")
    assert worst_m < 1e-12


def a2_identity(census):
    print()
    print("=" * 100)
    print("A2. Counting identity N = (1/M) sum_{h mod M} S(h)")
    print("=" * 100)
    n_exact = 0
    maxM = 0
    for (p, a, q), (N, _) in sorted(census.items()):
        M = abs(gate(p, a, q))
        if a == 0 or M > 20000 or p > 24:
            continue
        ell, Sv = S_dp_exact_all(p, a, q)
        tot = int(Sv.sum() % ell)
        Nx = tot * pow(M, -1, ell) % ell
        assert Nx == N, (p, a, q, Nx, N)
        n_exact += 1
        maxM = max(maxM, M)
    print(f"  exact (image of Z[zeta_M] in F_l, l prime = 1 mod M): {n_exact} clocks with M <= 20000, q in {{3,5}},")
    print(f"  largest M = {maxM}: (1/M) sum_h S(h) equals the census count on every clock.")
    worst, nf = 0.0, 0
    for (p, a, q), (N, _) in sorted(census.items()):
        M = abs(gate(p, a, q))
        if a == 0 or M > 200000 or p > 26:
            continue
        hs = list(range(M))
        tot = 0.0
        for k0 in range(0, M, 4000):
            tot += S_dp(p, a, hs[k0:k0 + 4000], q).real.sum()
        err = abs(tot / M - N)
        worst = max(worst, err)
        nf += 1
    print(f"  float: {nf} clocks with M <= 200000, max |(1/M) sum_h Re S(h) - N| = {worst:.2e}")
    assert worst < 1e-6


def run_census(q, pmax):
    out = {}
    t0 = time.time()
    nm = nr = 0
    for p in range(1, pmax + 1):
        for a in range(0, p + 1):
            G = gate(p, a, q)
            r = census_range(p, a, q) if a > 0 else (1, [0])
            if abs(G) < (1 << 62):
                m = census_mitm(p, a, q, want_words=True)
                assert m == r, (p, a, q, m, r)
                nm += 1
            nr += 1
            out[(p, a, q)] = r
    print(f"  q = {q}: {nr} clocks (p <= {pmax}) by the orbit-range method; {nm} of them also by the "
          f"residue join, identical counts and points on all {nm}.   [{time.time() - t0:.1f}s]", flush=True)
    return out


def cycles_from_points(census, q):
    """group integral periodic points into cycles; return dict min-point -> (p0, a0, sorted points)."""
    cyc = {}
    for (p, a, qq), (N, pts) in census.items():
        if qq != q:
            continue
        for x in pts:
            # recompute the primitive period of x under T_q (sign-aware: x<0 uses the same formula)
            y, orbit = x, [x]
            for _ in range(p):
                y = y // 2 if y % 2 == 0 else (q * y + 1) // 2
                orbit.append(y)
            assert orbit[p] == x
            per = next(d for d in range(1, p + 1) if orbit[d] == x)
            pts0 = sorted(set(orbit[:per]))
            a0 = sum(1 for v in orbit[:per] if v % 2)
            key = min(pts0, key=abs)
            cyc[key] = (per, a0, pts0)
    return cyc


def a3_a4(census, q, label):
    print()
    print("=" * 100)
    print(f"{label}: census N(p,a) = #{{w : |2^p - {q}^a| divides c_w}} versus the main term C(p,a)/|G|")
    print("=" * 100)
    cyc = cycles_from_points(census, q)
    print(f"  distinct integral cycles found (all p <= {P_CENSUS}):")
    for key in sorted(cyc, key=lambda k: (cyc[k][0], k)):
        per, a0, pts0 = cyc[key]
        side = "dyadic 2^p > q^a (positive cycle of qx+1)" if gate(per, a0, q) > 0 else \
            "q-adic q^a > 2^p (negative for qx+1 = positive cycle of qx-1)"
        s = pts0 if len(pts0) <= 12 else pts0[:6] + ["..."] + pts0[-3:]
        print(f"    clock ({per},{a0}), G = {gate(per, a0, q)}, {side}: {s}")
    # predicted counts from the cycles: each cycle of clock (p0,a0) contributes p0 words at (k p0, k a0)
    bad = 0
    for (p, a, qq), (N, pts) in census.items():
        if qq != q:
            continue
        pred = 0
        for key, (per, a0, pts0) in cyc.items():
            if a0 == 0 and a == 0:
                pred += 1 if per == 1 else 0
                continue
            if per and p % per == 0 and a0 * (p // per) == a:
                pred += per
        if pred != N:
            bad += 1
    assert bad == 0
    print(f"  every N(p,a) equals the sum over these cycles of (period) at the multiples of their clock: "
          f"0 exceptions.")
    print()
    print("  per-period totals (a = 0 and forced clocks included in N; main term excludes a = 0):")
    print(f"  {'p':>3} | {'N dyadic':>8} {'main dyadic':>11} | {'N q-adic':>8} {'main q-adic':>11} | "
          f"sporadic necklaces (|G| > 1, primitive)")
    tot_main = {1: 0.0, -1: 0.0}
    for p in range(1, P_CENSUS + 1):
        nd = nt = 0
        md = mt = 0.0
        spor = []
        for a in range(0, p + 1):
            G = gate(p, a, q)
            N, pts = census[(p, a, q)]
            if G > 0:
                nd += N
            else:
                nt += N
            if a == 0:
                continue
            mt_ = math.comb(p, a) / abs(G)
            if G > 0:
                md += mt_
            else:
                mt += mt_
        for key, (per, a0, pts0) in cyc.items():
            if per == p and abs(gate(per, a0, q)) > 1 and a0 > 0:
                spor.append(f"({per},{a0}) min {key}")
        tot_main[1] += md
        tot_main[-1] += mt
        if p <= 30 or p % 2 == 0:
            print(f"  {p:>3} | {nd:>8} {md:>11.4f} | {nt:>8} {mt:>11.4f} | {', '.join(spor)}")
    print(f"  sum over p <= {P_CENSUS} of the word-level main term: dyadic {tot_main[1]:.3f}, "
          f"q-adic {tot_main[-1]:.3f}")
    # per-clock table for the clocks with the largest main term
    rows = []
    for (p, a, qq), (N, pts) in census.items():
        if qq != q or a == 0:
            continue
        G = gate(p, a, q)
        rows.append((math.comb(p, a) / abs(G), p, a, G, N, lyndon_count(p, a)))
    rows.sort(reverse=True)
    print()
    print("  the 16 clocks with the largest main term C/|G| (word level) and the necklace-level term L/|G|:")
    print(f"  {'(p,a)':>8} {'G':>14} {'C/|G|':>9} {'L/|G|':>8} {'N':>4}")
    for mt_, p, a, G, N, L in rows[:16]:
        print(f"  {str((p, a)):>8} {G:>14} {mt_:>9.4f} {L / abs(G):>8.4f} {N:>4}")
    return cyc


def a5_rotation():
    """exact rotation law c_{Rw} = c_w/2 or (q c_w + G)/2, hence c_{R^j w} = 2^-j q^(k_j(w)) c_w mod M,
    and the transport identity S(h) = sum_w e(2^-j q^(k_j(w)) h c_w / M) for every j."""
    print()
    print("=" * 100)
    print("A5. Rotation law and the transport identity for S(h)")
    print("=" * 100)
    import cmath
    nw = nid = 0
    worst = 0.0
    rng = random.Random(7)
    for q in (3, 5):
        for p in range(2, 15):
            for a in range(1, p):
                G = gate(p, a, q)
                M = abs(G)
                inv2 = pow(2, -1, M)
                W = list(words(p, a))
                cs = {}
                for pos in W:
                    w = [0] * p
                    for t in pos:
                        w[t] = 1
                    cs[tuple(w)] = carry(pos, a, q)
                for w, c in cs.items():
                    rw = w[1:] + w[:1]
                    assert rotate_carry(c, w[0], G, q) == cs[rw]
                    nw += 1
                hs = [1, 2, rng.randrange(M)] if M > 2 else [1]
                for h in hs:
                    S0 = sum(cmath.exp(2j * math.pi * (h * c % M) / M) for c in cs.values())
                    for j in (1, 2, p // 2, p):
                        tot = 0j
                        for w, c in cs.items():
                            k = sum(w[:j])
                            mult = pow(inv2, j, M) * pow(q, k, M) % M
                            tot += cmath.exp(2j * math.pi * (h * mult * c % M) / M)
                        worst = max(worst, abs(tot - S0))
                        nid += 1
    print(f"  c_(Rw) = c_w/2 (w_0 = 0) or (q c_w + G)/2 (w_0 = 1): exact on all {nw} words, q in {{3,5}}, p <= 14.")
    print(f"  S(h) = sum_w e(2^-j q^(k_j(w)) h c_w / M), k_j(w) = #ones among the first j letters: {nid} checks,")
    print(f"  max |difference| = {worst:.2e}.  (So S(2^-j q^k h) = sum over words with k ones in their LAST j letters of")
    print("  e(h c/M) + the rest: a large archimedean coefficient at h is transported to the frequencies 2^-j q^k h.)")
    assert worst < 1e-8


def main():
    lean_malloc()
    t = time.time()
    a1_certify()
    print(f"[A1 {time.time() - t:.1f}s]", file=sys.stderr, flush=True)
    print()
    print("=" * 100)
    print("A3/A4 census runs")
    print("=" * 100)
    cen3 = run_census(3, P_CENSUS)
    cen5 = run_census(5, P_CENSUS)
    census = dict(cen3)
    census.update(cen5)
    a2_identity(census)
    print(f"[A2 {time.time() - t:.1f}s]", file=sys.stderr, flush=True)
    a3_a4(census, 3, "A3 (q = 3; SHEET control on the side 3^a > 2^p)")
    a3_a4(census, 5, "A4 (q = 5; DRIFT control)")
    a5_rotation()
    report_mem("end")
    print(f"[A total {time.time() - t:.1f}s]", file=sys.stderr)


if __name__ == "__main__":
    main()
