#!/usr/bin/env python3
"""
speedups_s266.py — Finding speedups via Lie-theoretic structure
opus-2026-03-23-S266

CURRENT STATE:
  - V_n, T_n, twin_SL at any n: O(p_odd(n)) time, sub-second to n=50
  - E(G_n) exact: only known to n=9 (5 hours via nauty)
  - Residual: 0,1,8,38,222,1463,15721 — no closed form
  - V_even at odd n: closed-form Burnside (new from S262-S263)

SPEEDUP OPPORTUNITIES:

1. STABILIZER BURNSIDE for E(G_n) (kind-pasteur S20dz):
   Fix arc {0,1}, compute via S_2 × S_{n-2}. 36× faster at n=9.
   But still requires tournament enumeration.

2. PARTITION-LEVEL STABILIZER BURNSIDE:
   Instead of enumerating permutations of S_2 × S_{n-2},
   iterate over their CYCLE TYPES. This is O(p(n-2)) partitions.
   For each cycle type of S_{n-2} × identity on {0,1}:
     - Compute Fix(σ) = tournaments fixed by σ
     - Check if flipping arc {0,1} changes the class
   Problem: "check if flip changes class" is the hard part.

3. ROOT SYSTEM ORTHOGONALITY:
   From S265: 90% of arc pairs are orthogonal roots at large n.
   Orthogonal arcs produce INDEPENDENT flip effects.
   This means: the metagraph is LOCALLY a product graph.
   Speedup: local edge count ≈ product of independent sub-counts.

4. IDENTITY-DOMINATED ASYMPTOTICS:
   At n=12, identity contributes 99.96% of Burnside sum.
   So E(G_n) ≈ E_identity + small correction.
   E_identity = (labeled_E) / n! where labeled_E uses all 2^m tournaments.

5. EXACT E VIA DOUBLE BURNSIDE (kind-pasteur S20ea):
   E = (1/n!) Σ_σ |edges preserved by σ|
   Works at n=3,4 (fails at n=5 due to overcounting).
   The overcounting = MW (multi-wire). If we can compute MW...

6. EXACT E FROM THE FACTORIZATION:
   E = (T - twin_SL - residual) / 2
   We know T and twin_SL exactly. The ONLY unknown is residual.
   If we can compute residual from Burnside, E is exact.

7. RESIDUAL FROM k-CYCLE DECOMPOSITION (kind-pasteur S20ec):
   4cycle_SL = 2^n / (8*(n-4)!) — exact formula for 4-cycle contribution.
   This captures 25% of residual at n=5, 5% at n=6.
   The rest is MW (multi-wire), which is global.

8. V_EVEN FROM CLOSED-FORM (new from S262):
   V_even(n) = (1/n!) Σ_σ ccs(σ) × 2^{eo(σ) - k(σ) + 1}
   For ODD n: works for ALL cycle types (no constraint matrix!).
   O(p(n)) time — same as V_tourn.

Let me implement the fastest possible computation of everything.
"""

import sys, time
from math import comb, factorial, gcd
from functools import lru_cache
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

# ============================================================
# PARTITION GENERATION (fastest: iterative with memoization)
# ============================================================

@lru_cache(maxsize=None)
def _fact(n):
    return factorial(n)

def gen_odd_parts(n):
    """All partitions of n into odd parts."""
    result = []
    def _g(rem, mx, cur):
        if rem == 0:
            result.append(tuple(cur))
            return
        st = min(rem, mx)
        if st % 2 == 0: st -= 1
        for p in range(st, 0, -2):
            _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return result

def gen_all_parts(n):
    """All partitions of n."""
    result = []
    def _g(rem, mx, cur):
        if rem == 0:
            result.append(tuple(cur))
            return
        for p in range(min(rem, mx), 0, -1):
            _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return result

def ccs(n, ct):
    c = Counter(ct); r = _fact(n)
    for l, k in c.items(): r //= pow(l, k) * _fact(k)
    return r

def pair_orb(ct):
    s = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)):
            s += gcd(ct[i], ct[j])
    return s

def edge_orb(ct):
    s = sum(c//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)):
            s += gcd(ct[i], ct[j])
    return s

# ============================================================
# THE SPEED ENGINE v2: everything from partitions
# ============================================================

def compute_everything(n):
    """Compute ALL known quantities from partition-level Burnside."""
    nf = _fact(n)
    m = comb(n, 2)

    # === ODD PARTITIONS → tournaments, transitions, twin_SL ===
    odd_parts = gen_odd_parts(n)

    V_sum = 0; T_sum = 0; twin_sum = 0

    for ct in odd_parts:
        cc = ccs(n, ct)
        po = pair_orb(ct)
        k = len(ct)
        f = ct.count(1)
        fix = pow(2, po)

        V_sum += cc * fix
        if f >= 2:
            T_sum += cc * comb(f, 2) * fix
            twin_sum += cc * 2 * comb(f, 2) * pow(2, po - k + 1)

    V = V_sum // nf
    T = T_sum // nf
    twin_SL = twin_sum // nf
    EO = T // 2 + _fact(n-2)

    # === ALL PARTITIONS → V_even (odd n only), V_graph ===
    all_parts = gen_all_parts(n)
    V_even_sum = 0
    V_graph_sum = 0

    for ct in all_parts:
        cc = ccs(n, ct)
        eo = edge_orb(ct)
        k = len(ct)
        V_graph_sum += cc * pow(2, eo)
        # V_even closed form: works for ALL types at odd n
        V_even_sum += cc * pow(2, eo - k + 1)

    V_graph = V_graph_sum // nf
    V_even = V_even_sum // nf if n % 2 == 1 else None

    # === ANTI-AUTOMORPHISM → SC ===
    # For SC, need fix_anti which requires building permutation
    # Skip for speed (compute separately if needed)

    # === APPROXIMATE E ===
    E_master = T * (pow(2, n-1) - 2) // pow(2, n)
    E_twin = (T - twin_SL) // 2

    # === 4-CYCLE CORRECTION (kind-pasteur S20ec) ===
    if n >= 5:
        # 4cycle_SL = 2^n / (8 * (n-4)!)
        # But this is for labeled count; need to Burnside it
        # Approximate: 4cycle contribution ≈ n*(n-1)*(n-2)*(n-3) * 2^{n-5} / n!
        # = 2^{n-5} / ((n-4)!)... need to verify
        pass

    return {
        'n': n, 'V': V, 'T': T, 'twin_SL': twin_SL, 'EO': EO,
        'E_master': E_master, 'E_twin': E_twin,
        'V_graph': V_graph, 'V_even': V_even,
        'n_odd_parts': len(odd_parts), 'n_all_parts': len(all_parts),
    }

# ============================================================
# BENCHMARK: How fast can we go?
# ============================================================

if __name__ == '__main__':
    print("=" * 80)
    print("  SPEED ENGINE v2 — ALL QUANTITIES FROM PARTITIONS")
    print("  opus-2026-03-23-S266")
    print("=" * 80)

    known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

    # === Benchmark: n=3..30 ===
    print(f"\n{'n':>3} {'V':>12} {'T':>12} {'EO':>12} {'E_master':>12} "
          f"{'E_twin':>12} {'V_even':>8} {'V_graph':>10} {'#p':>5} {'time':>8}")

    total_t = time.time()
    for n in range(3, 31):
        t0 = time.time()
        r = compute_everything(n)
        dt = time.time() - t0

        ve_str = str(r['V_even']) if r['V_even'] is not None else '—'
        if r['V'] < 1e12:
            print(f"{n:>3} {r['V']:>12} {r['T']:>12} {r['EO']:>12} "
                  f"{r['E_master']:>12} {r['E_twin']:>12} {ve_str:>8} "
                  f"{r['V_graph']:>10} {r['n_all_parts']:>5} {dt:>7.3f}s")
        else:
            print(f"{n:>3} {r['V']:>12.4e} {r['T']:>12.4e} {r['EO']:>12.4e} "
                  f"{r['E_master']:>12.4e} {r['E_twin']:>12.4e} {ve_str:>8} "
                  f"{r['V_graph']:>10.4e} {r['n_all_parts']:>5} {dt:>7.3f}s")

    total_dt = time.time() - total_t
    print(f"\n  Total n=3..30: {total_dt:.3f}s")

    # === Accuracy check ===
    print(f"\n  ACCURACY:")
    print(f"  {'n':>3} {'E_exact':>10} {'E_master':>10} {'err_m%':>7} {'E_twin':>10} {'err_t%':>7}")
    for n in sorted(known_E.keys()):
        r = compute_everything(n)
        ek = known_E[n]
        em = r['E_master']
        et = r['E_twin']
        err_m = 100 * abs(em - ek) / ek
        err_t = 100 * abs(et - ek) / ek
        print(f"  {n:>3} {ek:>10} {em:>10} {err_m:>6.2f}% {et:>10} {err_t:>6.2f}%")

    # === NEW: V_even at odd n ===
    print(f"\n  V_EVEN (odd n only, closed form):")
    print(f"  {'n':>3} {'V_tourn':>12} {'V_even':>10} {'T/E ratio':>10} {'V_graph':>10}")
    for n in range(3, 22, 2):
        r = compute_everything(n)
        if r['V_even']:
            ratio = r['V'] / r['V_even']
            print(f"  {n:>3} {r['V']:>12} {r['V_even']:>10} {ratio:>10.2f} {r['V_graph']:>10}")

    # === WHAT'S STILL MISSING ===
    print(f"\n{'='*80}")
    print("  WHAT'S STILL MISSING FOR EXACT E(G_n)")
    print(f"{'='*80}")
    print("""
  KNOWN EXACTLY (from cycle index):
    V_n, T_n, twin_SL, EO, SC_n, V_merged, V_even (odd n)

  KNOWN APPROXIMATELY:
    E_master = T × (2^{n-1}-2)/2^n  (error <0.5% at n≥9)
    E_twin = (T - twin_SL)/2        (error <1.6% at n≥8)

  THE LAST UNKNOWN: residual R = T - twin_SL - 2E
    Sequence R/2: 0, 1, 8, 38, 222, 1463, 15721
    - Cannot be computed from partition-level Burnside alone
    - Requires GLOBAL information (multi-wire collisions)
    - kind-pasteur showed: 4cycle captures 25% at n=5, MW=75%
    - MW = edge multiplicity excess = genuinely global phenomenon

  SPEEDUP POSSIBILITIES:
    1. Compute 4cycle_SL exactly (another Burnside sum) → reduces residual
    2. Edge-centric Burnside at n=10 via partition of S_2 × S_8 → O(p(8)×2) time
    3. Nauty-based E at n=10 if feasible (would extend known E by 1 term)
    4. Asymptotic formula for MW from random tournament theory
""")

    # === SPEEDUP 1: How fast can we compute V_even to n=100? ===
    print(f"\n{'='*80}")
    print("  SPEEDUP: V_even to n=100 (odd n)")
    print(f"{'='*80}")

    t0 = time.time()
    for n in range(3, 102, 2):
        r_t = time.time()
        all_parts = gen_all_parts(n)
        nf = _fact(n)
        V_even_sum = 0
        for ct in all_parts:
            cc = ccs(n, ct)
            eo = edge_orb(ct)
            k = len(ct)
            V_even_sum += cc * pow(2, eo - k + 1)
        V_even = V_even_sum // nf
        dt = time.time() - r_t
        if n <= 15 or n % 10 == 1:
            print(f"  n={n:>3}: V_even = {V_even:>20} ({len(all_parts):>6} parts, {dt:.3f}s)")

    print(f"\n  Total time for V_even at odd n=3..101: {time.time()-t0:.3f}s")

    print("\nDONE.")
