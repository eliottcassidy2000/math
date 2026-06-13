#!/usr/bin/env python3
"""
two_three_tower_deep.py -- The modular tower of I(Omega, x): deeper analysis

Building on ten_eleven_deep.py results:
1. CRT tower: I(b) mod (b+1) recovers I(-1) progressively
2. 6-adic tower: H mod 6^k structure
3. Difference operator Delta^k I(b) and Stirling numbers
4. I(omega_3) cube root evaluation -- resultant connection
5. 2-adic and 3-adic valuations of H-1
6. Eisenstein integer structure: I(omega) in Z[omega]
7. Generating function G(z) = sum I(b) z^b
8. The beta basis as 3-adic Taylor expansion
9. Phase transition table across n

Author: kind-pasteur-2026-03-14-S63
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
from math import gcd, comb, factorial
import cmath


def random_tournament(n, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def get_directed_cycles(A, n):
    cycles = []
    for k in range(3, n+1, 2):
        for subset in combinations(range(n), k):
            verts = list(subset)
            nc = count_ham_cycles(A, verts)
            for _ in range(nc):
                cycles.append(frozenset(subset))
    return cycles


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp or dp[(mask, v)] == 0:
                continue
            cnt = dp[(mask, v)]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and dp[(full, v)] > 0:
            if A[verts[v]][verts[0]]:
                total += dp[(full, v)]
    return total


def compute_alpha_1_2(cycles):
    alpha_1 = len(cycles)
    alpha_2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if not (cycles[i] & cycles[j]):
                alpha_2 += 1
    return alpha_1, alpha_2


def I_omega(a1, a2, x):
    if isinstance(x, complex):
        return complex(1) + a1 * x + a2 * x * x
    return 1 + a1 * x + a2 * x * x


def crt_combine(r1, m1, r2, m2):
    g = gcd(m1, m2)
    if (r1 - r2) % g != 0:
        return None, None
    lcm = m1 * m2 // g
    def ext_gcd(a, b):
        if a == 0: return b, 0, 1
        g2, x, y = ext_gcd(b % a, a)
        return g2, y - (b // a) * x, x
    _, p, q = ext_gcd(m1, m2)
    x = (r1 + m1 * p * (r2 - r1) // g) % lcm
    return x, lcm


def main():
    n = 7
    rng = np.random.default_rng(42)

    # Collect data
    num_samples = 300
    tournaments = []
    for _ in range(num_samples):
        A = random_tournament(n, rng)
        cycles = get_directed_cycles(A, n)
        a1, a2 = compute_alpha_1_2(cycles)
        H = I_omega(a1, a2, 2)
        euler = I_omega(a1, a2, -1)
        tournaments.append({'a1': a1, 'a2': a2, 'H': H, 'euler': euler})

    # ============================================================
    print("=" * 70)
    print("PART 1: CRT Tower -- progressively recovering I(-1)")
    print("=" * 70)

    print("\n  Verification: I(b) mod (b+1) = I(-1) mod (b+1)")
    for b in range(2, 21):
        ok = all(I_omega(t['a1'], t['a2'], b) % (b+1) == t['euler'] % (b+1)
                 for t in tournaments)
        status = "OK" if ok else "FAIL"
        print(f"    b={b:2d}: {status}")

    # Build CRT tower for first tournament
    t0 = tournaments[0]
    print(f"\n  CRT tower for a1={t0['a1']}, a2={t0['a2']}, I(-1)={t0['euler']}:")
    r, m = t0['euler'] % 3, 3
    print(f"    Start: I(-1) = {r} mod {m}")
    for b in range(3, 20):
        nr = t0['euler'] % (b+1)
        cr, cm = crt_combine(r, m, nr, b+1)
        if cr is not None and cm > m:
            r, m = cr, cm
            print(f"    b={b:2d}: new info mod {b+1:2d} -> I(-1) = {r:6d} (mod {m})")

    print(f"    Actual: I(-1) = {t0['euler']}")
    print(f"    Recovered: {r} mod {m}, match: {t0['euler'] % m == r}")

    # How many evaluations to determine I(-1)?
    print("\n  Evaluations needed to uniquely determine I(-1):")
    for t in tournaments[:10]:
        r, m = t['euler'] % 3, 3
        needed = 1
        for b in range(3, 100):
            if m > abs(t['euler']) + 1:
                break
            nr = t['euler'] % (b+1)
            cr, cm = crt_combine(r, m, nr, b+1)
            if cr is not None and cm > m:
                r, m = cr, cm
                needed += 1
        print(f"    I(-1)={t['euler']:4d}: {needed} evaluations (modulus={m})")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 2: 6-adic tower -- H mod 6^k")
    print("=" * 70)

    for k in range(1, 5):
        mod = 6**k
        dist = defaultdict(int)
        for t in tournaments:
            dist[t['H'] % mod] += 1
        n_classes = len(dist)
        possible_odd = mod // 2
        print(f"\n  H mod {mod:5d}: {n_classes:3d} classes / {possible_odd} possible ({100*n_classes/possible_odd:.1f}%)")

        if mod <= 36:
            for r in sorted(dist):
                print(f"      {r:3d}: {dist[r]:3d} ({100*dist[r]/num_samples:.1f}%)")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 3: Difference operator and Stirling numbers")
    print("=" * 70)

    print("\n  I(x) = 1 + a1*x + a2*x^2")
    print("  Newton forward differences at 0:")
    print("    D^0 I(0) = 1")
    print("    D^1 I(0) = a1+a2")
    print("    D^2 I(0) = 2a2")
    print("    D^k I(0) = 0 for k>=3")
    print()
    print("  Newton: I(x) = C(x,0) + (a1+a2)*C(x,1) + 2a2*C(x,2)")
    print("  At x=2: I(2) = 1 + 2*(a1+a2) + 2a2 = 1 + 2a1 + 4a2")
    print("  At x=-1: I(-1) = 1 - (a1+a2) + 2a2 = 1 - a1 + a2")

    print("\n  Verification:")
    for t in tournaments[:5]:
        d0, d1, d2 = 1, t['a1']+t['a2'], 2*t['a2']
        H_check = d0 + 2*d1 + d2
        E_check = d0 - d1 + d2
        print(f"    a1={t['a1']:2d}, a2={t['a2']:2d}: D=[{d0},{d1},{d2}], "
              f"H={H_check}={t['H']}, I(-1)={E_check}={t['euler']}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 4: Cube root evaluation -- Eisenstein integers")
    print("=" * 70)

    w = cmath.exp(2j * cmath.pi / 3)

    print(f"\n  w = e^(2pi*i/3), w^2 + w + 1 = 0")
    print(f"  I(w) = 1 + a1*w + a2*w^2")
    print(f"  In Z[w]: I(w) = (1-a2) + (a1-a2)*w  [since w^2 = -1-w]")

    print("\n  Eisenstein norm N(I(w)) = (1-a2)^2 - (1-a2)(a1-a2) + (a1-a2)^2:")
    for t in tournaments[:10]:
        a1, a2 = t['a1'], t['a2']
        rp = 1 - a2
        wp = a1 - a2
        N_eis = rp**2 - rp*wp + wp**2
        N_dir = abs(I_omega(a1, a2, w))**2
        print(f"    a1={a1:2d}, a2={a2:2d}: I(w)=({rp})+({wp})w, "
              f"N={N_eis}, |I(w)|^2={N_dir:.1f}, match={abs(N_eis-N_dir)<0.01}")

    # N(I(w)) mod 3
    print("\n  N(I(w)) mod 3:")
    n3_dist = defaultdict(int)
    for t in tournaments:
        rp = 1 - t['a2']
        wp = t['a1'] - t['a2']
        N = rp**2 - rp*wp + wp**2
        n3_dist[N % 3] += 1
    print(f"    {dict(sorted(n3_dist.items()))}")

    # I(w) mod (1-w) = I(1) mod 3
    print("\n  I(w) mod (1-w) vs I(1) mod 3:")
    ok = all((1+t['a1']-2*t['a2']) % 3 == (1+t['a1']+t['a2']) % 3 for t in tournaments)
    print(f"    Always equal: {ok}")

    # Resultant: I(1)*|I(w)|^2
    print("\n  Resultant Res(I(x), x^3-1) = I(1)*I(w)*I(w^2) = I(1)*N(I(w)):")
    for t in tournaments[:5]:
        a1, a2 = t['a1'], t['a2']
        I1 = 1 + a1 + a2
        rp = 1 - a2
        wp = a1 - a2
        N = rp**2 - rp*wp + wp**2
        res = I1 * N
        print(f"    a1={a1:2d}, a2={a2:2d}: I(1)={I1:3d}, N(I(w))={N:5d}, Res={res:7d}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 5: 2-adic and 3-adic valuations of H-1")
    print("=" * 70)

    v2_dist = defaultdict(int)
    v3_dist = defaultdict(int)
    for t in tournaments:
        h1 = t['H'] - 1
        v2 = 0
        temp = h1
        while temp > 0 and temp % 2 == 0: v2 += 1; temp //= 2
        v2_dist[v2] += 1
        v3 = 0
        temp = abs(h1) if h1 > 0 else 0
        while temp > 0 and temp % 3 == 0: v3 += 1; temp //= 3
        v3_dist[v3] += 1

    print(f"\n  v_2(H-1) distribution:")
    for v in sorted(v2_dist):
        print(f"    v_2={v}: {v2_dist[v]:3d} ({100*v2_dist[v]/num_samples:.1f}%)")
    print(f"\n  v_3(H-1) distribution:")
    for v in sorted(v3_dist):
        print(f"    v_3={v}: {v3_dist[v]:3d} ({100*v3_dist[v]/num_samples:.1f}%)")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 6: Generating function G(z) = sum I(b) z^b")
    print("=" * 70)

    print("\n  G(z) = 1/(1-z) + a1*z/(1-z)^2 + a2*z*(1+z)/(1-z)^3")

    print("\n  G(1/2) - H = 1 + 2*a2 (disjoint-pair excess):")
    for t in tournaments[:10]:
        a1, a2 = t['a1'], t['a2']
        G_half = 2 + 2*a1 + 6*a2
        excess = G_half - t['H']
        print(f"    a1={a1:2d}, a2={a2:2d}: G(1/2)={G_half}, H={t['H']}, excess={excess}, 1+2a2={1+2*a2}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 7: The reduced formula H = 3 - 2*b0 + 6*b2")
    print("=" * 70)

    print("\n  From beta basis + constraint b0+b1+b2=1:")
    print("  H = b0 + 3*b1 + 9*b2 = b0 + 3*(1-b0-b2) + 9*b2 = 3 - 2*b0 + 6*b2")
    print()
    print("  H mod 6 depends ONLY on b0 mod 3:")
    for b0_mod3 in [0, 1, 2]:
        h_mod6 = (3 - 2*b0_mod3) % 6
        print(f"    b0 = {b0_mod3} mod 3 => H = {h_mod6} mod 6")

    print("\n  Verification:")
    for t in tournaments[:10]:
        b0 = 1 - t['a1'] + t['a2']
        b2 = t['a2']
        H_red = 3 - 2*b0 + 6*b2
        print(f"    b0={b0:4d}, b2={b2:2d}: H=3-2*{b0}+6*{b2}={H_red}={t['H']}, "
              f"b0%3={b0%3}, H%6={t['H']%6}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 8: Independence of 2-adic and 3-adic digits")
    print("=" * 70)

    for k2, k3 in [(1,1), (2,1), (1,2), (2,2)]:
        m2, m3 = 2**k2, 3**k3
        joint = defaultdict(int)
        for t in tournaments:
            joint[(t['H'] % m2, t['H'] % m3)] += 1

        marg2 = defaultdict(int)
        marg3 = defaultdict(int)
        for (r2, r3), cnt in joint.items():
            marg2[r2] += cnt
            marg3[r3] += cnt

        chi2 = 0
        for r2 in marg2:
            for r3 in marg3:
                expected = marg2[r2] * marg3[r3] / num_samples
                observed = joint.get((r2, r3), 0)
                if expected > 0:
                    chi2 += (observed - expected)**2 / expected

        df = max((len(marg2)-1) * (len(marg3)-1), 1)
        print(f"\n  mod {m2} x mod {m3}: chi2={chi2:.2f}, df={df}, chi2/df={chi2/df:.2f}")

        if m2 <= 4 and m3 <= 9:
            for r2 in sorted(marg2):
                row = " ".join(f"{joint.get((r2,r3),0):>4d}" for r3 in sorted(marg3))
                print(f"    H={r2}({m2}): {row}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 9: Phase transition across n")
    print("=" * 70)

    print("\n  n | b0<0 % | b1>0 % | <a2/a1> | C(n-3,3)/8")
    for nn in [5, 6, 7, 9]:
        rng_n = np.random.default_rng(42)
        b0_neg, b1_pos, total = 0, 0, 200
        a2a1_ratios = []
        for _ in range(total):
            A = random_tournament(nn, rng_n)
            cycles = get_directed_cycles(A, nn)
            a1, a2 = compute_alpha_1_2(cycles)
            if a1 > 0:
                a2a1_ratios.append(a2/a1)
            if 1 - a1 + a2 < 0: b0_neg += 1
            if a1 - 2*a2 > 0: b1_pos += 1

        pred = comb(nn-3, 3) / 8 if nn >= 6 else 0
        mean_ratio = np.mean(a2a1_ratios) if a2a1_ratios else 0
        print(f"  {nn:2d} | {100*b0_neg/total:5.1f}  | {100*b1_pos/total:5.1f}  | {mean_ratio:.4f}  | {pred:.4f}")

    # ============================================================
    print("\n" + "=" * 70)
    print("PART 10: Summary -- The 2-3 Architecture")
    print("=" * 70)

    print("""
  I(x) = 1 + a1*x + a2*x^2  (exact at n<=7)

  ALPHA basis (at x): I(x) = 1 + a1*x + a2*x^2
  BETA basis  (at -1): I(x) = b0 + b1*(x+1) + b2*(x+1)^2
  NEWTON basis (at 0): I(x) = C(x,0) + d1*C(x,1) + d2*C(x,2)

  Change of basis:
    b0 = 1-a1+a2,  b1 = a1-2a2,  b2 = a2
    d0 = 1,  d1 = a1+a2,  d2 = 2a2
    CONSTRAINT: b0+b1+b2 = 1

  REDUCED FORMULA: H = 3 - 2*b0 + 6*b2

  2 controls: evaluation point, Redei parity, deletion-contraction
  3 controls: mod structure, topology (I(-1)), beta step size
  6 = 2*3: Vandermonde det, natural modulus, reduced formula

  EISENSTEIN CONNECTION: I(w) in Z[w] with norm N(I(w))
  CRT TOWER: I(b) mod (b+1) = I(-1) mod (b+1) for all b
  GENERATING FUNCTION: G(1/2) = H + 1 + 2*a2
""")

    print("Done.")


if __name__ == '__main__':
    main()
