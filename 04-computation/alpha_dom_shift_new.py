#!/usr/bin/env python3
"""
alpha_dom_shift_new.py -- kind-pasteur-2026-04-16

DISCOVERY: The dominant term in H = I(Omega,2) = sum_k 2^k * alpha_k
shifts from k=1 (total odd cycles) to k=2 (disjoint pairs) between n=9 and n=11.

VERIFIED TABLE:
  n=3:  H=3,              α₁=1,           α₂=0,           dominant=α₁ (ratio=∞)
  n=5:  H=15,             α₁=7,           α₂=0,           dominant=α₁ (ratio=∞)
  n=7:  H=189,            α₁=80,          α₂=7,           dominant=α₁ (ratio=5.71)
  n=9:  H=3357,           α₁=954,         α₂=336,         dominant=α₁ (ratio=1.42)
  n=11: H=95095,          α₁=21169,       α₂=10879,       dominant=α₂ (ratio=0.97)
  n=13: H=3711175,        α₁=606027,      α₂=436748,      dominant=α₂ (ratio=0.69)
  n=15: H=198464295,      α₁=27495799,    α₂=22518662,    dominant=α₂ (ratio=0.61)
  n=17: H=13689269499,    α₁=1651334601,  α₂=1482234998,  dominant=α₂ (ratio=0.557)

KEY FORMULA:
  ratio = alpha_1 / (2 * alpha_2) crosses 1 between n=9 and n=11.
  ratio: inf, inf, 5.71, 1.42, 0.97, 0.69, 0.61, 0.557 for n=3..17.

ALPHA_4 STRUCTURE (n=17):
  [3,3,3,3]:   945,744  (2.1%)  -- four 3-cycles, 12 vertices
  [3,3,3,L≥5]: 28,778,144 (62.6%) -- three 3-cycles + one longer cycle
  [3,3,5,5]:   16,273,216 (35.4%) -- two 3-cycles + two 5-cycles (NEW at n=17)
  Total α₄:   45,997,104

ALPHA_3/ALPHA_2 RATIO TREND (α₃/α₂):
  n=9: 0.039, n=11: 0.106, n=13: 0.201, n=15: 0.260, n=17: 0.309
  Dominance threshold (8α₃ > 4α₂ i.e. α₃/α₂ > 0.5) not yet reached.
  Extrapolating: α₃ dominance over α₂ around n=21-25.

ALPHA_4 CORRECTION (n≥17): The "triple-of-3-cycles" approach misses [3,3,5,5]
quadruples when 3+3+5+5 ≤ n. Safe at n≤15 (16>15), needed at n≥17.
See n17_alpha2.py for the [3,3,5,5] fix.

CORRECTNESS NOTE: T_7 uses QR={1,2,4} NOT {1,2,3}! The standard Paley T_p uses
quadratic residues mod p. T_7 circulant({1,2,3}) gives only H=175, not 189.
n=17: 17≡1 mod 4, no Paley. Best circulants: S={1..8}, {all-odd}, {all-even}, H=13,689,269,499.
"""

from itertools import combinations
import time


def circulant(n, S):
    """Circulant tournament on Z_n: v beats (v+k) mod n for k in S."""
    adj = [0]*n
    for v in range(n):
        for k in S:
            adj[v] |= 1 << ((v+k) % n)
    return adj


def H_dp(adj):
    """Held-Karp DP for Hamiltonian path count. O(2^n * n^2)."""
    n = len(adj); N = 1 << n
    dp = [[0]*n for _ in range(N)]
    for v in range(n): dp[1 << v][v] = 1
    full = N-1
    for mask in range(1, N):
        row = dp[mask]
        for v in range(n):
            val = row[v]
            if not val: continue
            outs = adj[v] & ~mask & full
            while outs:
                ub = outs & -outs; u = ub.bit_length()-1
                dp[mask|ub][u] += val; outs ^= ub
    return sum(dp[full])


def cycle_cc(adj, n):
    """
    cc[mask] = # directed odd cycles with exactly vertex set `mask`.
    Canonical form: start at min vertex. O(n * sum_{L odd} C(n,L) * (L-1)!) in worst case,
    but DP is much faster in practice: O(n * 2^n) states.
    """
    full = (1 << n)-1
    cc = [0] * (1 << n)
    for s in range(n):
        s_bit = 1 << s
        # Only extend to vertices >= s (s = minimum vertex of cycle)
        hi_mask = full & ~(s_bit - 1)
        queue = {(s_bit, s): 1}
        while queue:
            nq = {}
            for (mask, v), cnt in queue.items():
                L = bin(mask).count('1')
                if L >= 3 and L % 2 == 1:
                    if (adj[v] >> s) & 1:
                        cc[mask] += cnt
                cands = adj[v] & ~mask & hi_mask
                while cands:
                    ub = cands & -cands; u = ub.bit_length()-1
                    key = (mask|ub, u)
                    nq[key] = nq.get(key, 0) + cnt
                    cands ^= ub
            queue = nq
    return cc


def sos(f, n):
    """Sum-over-subsets (zeta) transform: f[T] = sum_{m subset of T} f[m]."""
    f = list(f)
    for i in range(n):
        for mask in range(1 << n):
            if (mask >> i) & 1:
                f[mask] += f[mask ^ (1 << i)]
    return f


def alpha_decomp_full(adj, n):
    """
    Compute alpha_k for all k up to floor(n/3).
    Returns (alphas dict, H, cycle_counts_by_length).
    """
    full = (1 << n)-1
    kmax = n // 3

    cc = cycle_cc(adj, n)
    f = sos(cc, n)

    alpha1 = sum(cc)
    # alpha_2 = (1/2) * sum_m cc[m] * f[~m & full]
    alpha2 = sum(cc[m] * f[(~m) & full] for m in range(1 << n)) // 2
    # alpha_3 = (1/3) * sum_{m1<m2, disjoint} cc[m1]*cc[m2]*f[~(m1|m2)&full]
    nonzero = [(m, cc[m]) for m in range(1 << n) if cc[m]]
    alpha3 = sum(c1*c2*f[(~m1 & ~m2) & full]
                 for m1, c1 in nonzero
                 for m2, c2 in nonzero
                 if m2 > m1 and not (m1 & m2)) // 3

    alphas = {0: 1, 1: alpha1, 2: alpha2, 3: alpha3}

    if kmax >= 4:
        # alpha_4 at n=12-14: all quadruples must be four 3-cycles (4*3=12 ≤ n, 3*3+5 > n)
        # alpha_4 = (1/4) * sum_{disjoint triples of 3-cycles} c1*c2*c3 * f3[complement]
        cc3 = [cc[m] if bin(m).count('1') == 3 else 0 for m in range(1 << n)]
        f3 = sos(cc3, n)
        tc3 = [(m, cc3[m]) for m in range(1 << n) if cc3[m]]
        raw = 0
        for i, (m1, c1) in enumerate(tc3):
            m1_cands = [(m, c) for m, c in tc3[i+1:] if not (m1 & m)]
            for j, (m2, c2) in enumerate(m1_cands):
                m12 = m1 | m2
                for m3, c3 in m1_cands[j+1:]:
                    if not (m12 & m3):
                        raw += c1 * c2 * c3 * f3[(~(m12 | m3)) & full]
        alphas[4] = raw // 4

    H = sum(2**k * v for k, v in alphas.items())
    len_counts = {}
    for mask, cnt in enumerate(cc):
        if cnt:
            L = bin(mask).count('1')
            len_counts[L] = len_counts.get(L, 0) + cnt

    return alphas, H, len_counts


# ─── Paley quadratic residues ───────────────────────────────────────────────

def qr(p):
    """Quadratic residues mod p (for p prime, p ≡ 3 mod 4)."""
    return sorted({i*i % p for i in range(1, p)} - {0})


# ─── Main: build and verify the mechanism-shift table ───────────────────────

if __name__ == '__main__':
    print("="*72)
    print("DOMINANT MECHANISM SHIFT TABLE")
    print("H = 1 + 2α₁ + 4α₂ + 8α₃ + 16α₄")
    print("="*72)

    tournaments = [
        (3,  [1],         "Paley T_3"),
        (5,  [1, 3],      "SC reg n=5"),      # QR mod 5 invalid, use SC regular
        (7,  [1, 2, 4],   "Paley T_7"),        # QR mod 7 = {1,2,4}
        (9,  [1, 2, 3, 5], "max circulant n=9"),
        (11, [1, 3, 4, 5, 9], "Paley T_11"),   # QR mod 11
        (13, [1, 3, 5, 7, 9, 11], "max circulant n=13"),
    ]

    print(f"\n{'n':>4} {'tournament':>22} {'H':>12} {'α₁':>8} {'α₂':>8} {'α₃':>8} "
          f"{'α₄':>6} {'dom':>4} {'α₁/(2α₂)':>10}")
    print("-"*90)

    for n, S, name in tournaments:
        adj = circulant(n, S)
        H = H_dp(adj)
        alphas, H_ocf, lc = alpha_decomp_full(adj, n)
        assert H == H_ocf, f"H mismatch at n={n}: {H} vs {H_ocf}"

        terms = [2**(k+1) * alphas.get(k, 0) for k in range(1, n//3+2)]
        dom_k = 1 + max(range(len(terms)), key=lambda i: terms[i])
        ratio = (alphas[1] / (2 * alphas[2])) if alphas.get(2, 0) > 0 else float('inf')

        print(f"{n:>4} {name:>22} {H:>12,} {alphas.get(1,0):>8} {alphas.get(2,0):>8} "
              f"{alphas.get(3,0):>8} {alphas.get(4,0):>6} {'α_'+str(dom_k):>4} "
              f"{ratio:>10.4f}")

    print()
    print("CROSSOVER: α₁/(2α₂) ratio crosses 1 between n=9 (1.42) and n=11 (0.97).")
    print("MECHANISM: At n≤9, maximize total odd cycles. At n≥11, maximize disjoint pairs.")
    print("PATTERN: kmax = floor(n/3). Top-level α_kmax counts near-perfect cycle factorizations.")
    print()
    print("CORRECTNESS NOTES:")
    print("  - T_7 uses QR mod 7 = {1,2,4}, NOT {1,2,3}")
    print("  - 5≡1 mod 4: no Paley for n=5; max circulants: S={1,3} or S={2,4}, both H=15")
    print("  - 9 is not prime: no Paley; max circulant is {1,2,3,5}")
    print("  - 13≡1 mod 4: no Paley; max circulant is {1,3,5,7,9,11} (all-odd S)")


# ─── EXTENDED TABLE including n=15,17 ───────────────────────────────────────
#
# QUADRUPLE TYPE KEY:
#   [3,3,3,3]: four 3-cycles (all-3-cycle)
#   [3,3,3,5]: three 3-cycles + one 5-cycle (first appears at n=14)
#   [3,3,5,5]: two 3-cycles + two 5-cycles (first appears at n=16, needed at n=17)
#   [3,3,3,7]: three 3-cycles + one 7-cycle (first appears at n=16, in [3,3,3,L≥5])
#
# ALPHA_4 BREAKDOWN:
#   n=15: 397,720 = 83,800 [3,3,3,3] + 313,920 [3,3,3,5]  ([3,3,5,5] = 16 > 15)
#   n=17: 45,997,104 = 945,744 [3333] + 28,778,144 [333L≥5] + 16,273,216 [3355]
#
# At n=15: alpha_5 = 7,472 = 7,472 [3,3,3,3,3] (perfect triangle factorizations)
#   = directed Steiner triple system parallel classes
#   3-cycles form 2-(15,3,4) design
# At n=17: alpha_5 = 1,800,368 = 703,936 [33333] + 1,096,432 [33335]
#
# FULL MECHANISM SHIFT TABLE (all verified):
#   n     H                 α₁            α₂           α₃        α₄        α₅     ratio
#   3         3               1             0            0          0         0       ∞
#   5        15               7             0            0          0         0       ∞
#   7       189              80             7            0          0         0    5.7143
#   9      3357             954           336           13          0         0    1.4196
#  11     95095           21169         10879         1155          0         0    0.9729
#  13   3711175          606027        436748        87568       3224         0    0.6938
#  15 198464295        27495799      22518662      5849428     397720      7472    0.6105
#  17 13689269499     1651334601    1482234998    458011858  45997104   1800368    0.5570
#
# The ratio α₁/(2α₂) decreases: ∞,∞,5.71,1.42,0.97,0.69,0.61,0.557
# Crossover from α₁-dominant to α₂-dominant: n=9→11 (ratio crosses 1)
#
# α₃/α₂ trend: 0, 0.039, 0.106, 0.201, 0.260, 0.309
# Dominance of α₃ over α₂ (i.e. 8α₃ > 4α₂ ⟺ α₃/α₂ > 0.5): estimated n≈21-25
