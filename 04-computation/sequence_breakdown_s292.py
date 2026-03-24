#!/usr/bin/env python3
"""
sequence_breakdown_s292.py — Break ALL sequences by cycle type contribution
opus-2026-03-24-S292

The Euler product trick: each quantity (V, T, twin_SL, E_approx)
is a Burnside sum over odd-cycle-type partitions. Break each sum
into contributions from each partition type.

For each partition type λ = (c_1, ..., c_k):
  V_λ = ccs(λ) × 2^{po(λ)} / n!
  T_λ = ccs(λ) × C(f,2) × 2^{po(λ)} / n!  [f = #fixed points]
  twin_λ = ccs(λ) × 2C(f,2) × 2^{po(λ)-k+1} / n!

Then: V = Σ V_λ, T = Σ T_λ, twin = Σ twin_λ

Can we express E and the residual as sums over λ too?
"""

import sys
from math import comb, factorial, gcd
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)

def gen_odd_parts(n):
    result = []
    def _g(rem, mx, cur):
        if rem == 0: result.append(tuple(sorted(cur, reverse=True))); return
        st = min(rem, mx)
        if st % 2 == 0: st -= 1
        for p in range(st, 0, -2): _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return list(set(result))

def ccs(n, ct):
    c = Counter(ct); r = factorial(n)
    for l, k in c.items(): r //= pow(l, k) * factorial(k)
    return r

def po(ct):
    s = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): s += gcd(ct[i], ct[j])
    return s

known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}

print("=" * 80)
print("  SEQUENCE BREAKDOWN BY CYCLE TYPE")
print("  opus-2026-03-24-S292")
print("=" * 80)

for n in range(3, 11):
    nf = factorial(n)
    parts = gen_odd_parts(n)
    parts.sort(key=lambda ct: (len(ct), ct))

    print(f"\n{'='*80}")
    print(f"  n = {n}: {len(parts)} odd-cycle partitions")
    print(f"{'='*80}")

    total_V = 0; total_T = 0; total_twin = 0

    rows = []
    for ct in parts:
        cc = ccs(n, ct)
        p = po(ct)
        k = len(ct)
        f = ct.count(1)
        fix = pow(2, p)

        V_ct = cc * fix
        T_ct = cc * comb(f, 2) * fix if f >= 2 else 0
        twin_ct = cc * 2 * comb(f, 2) * pow(2, p - k + 1) if f >= 2 else 0

        total_V += V_ct
        total_T += T_ct
        total_twin += twin_ct

        # Contribution to correction R(n) = (V × n! / 2^m) - 1
        # R_ct = V_ct / (2^m / n!) for non-identity partitions
        # Actually: V = (1/n!) Σ V_ct. And V × n! / 2^m = Σ V_ct / 2^m.
        # The identity partition [1^n] contributes 2^m / n! to V,
        # i.e., V_ct = 2^m for [1^n]. So V_ct/2^m = 1 for identity.
        # Non-identity: V_ct/2^m = correction.

        m = comb(n, 2)
        R_ct = V_ct / pow(2, m)

        rows.append({
            'ct': ct, 'cc': cc, 'po': p, 'k': k, 'f': f,
            'V': V_ct, 'T': T_ct, 'twin': twin_ct,
            'R': R_ct, 'is_identity': (ct == tuple([1]*n)),
        })

    V = total_V // nf
    T = total_T // nf
    twin = total_twin // nf
    E_twin = (T - twin) // 2
    E = known_E.get(n, E_twin)

    if n <= 8:
        print(f"\n  {'partition':>20} {'ccs':>8} {'po':>4} {'k':>3} {'f':>3} "
              f"{'V_ct':>12} {'T_ct':>12} {'twin_ct':>12} {'R_ct':>10}")

        for r in rows:
            if r['V'] > 0:
                print(f"  {str(r['ct']):>20} {r['cc']:>8} {r['po']:>4} {r['k']:>3} {r['f']:>3} "
                      f"{r['V']:>12} {r['T']:>12} {r['twin']:>12} {r['R']:>10.6f}")

    print(f"\n  TOTALS (÷ {nf}):")
    print(f"    V = {V}")
    print(f"    T = {T}")
    print(f"    twin_SL = {twin}")
    print(f"    E_twin = (T-twin)/2 = {E_twin}")
    if n in known_E:
        print(f"    E_exact = {E}")
        residual = T - twin - 2*E
        print(f"    residual = T - twin - 2E = {residual}")
        print(f"    residual/2 = {residual//2}")

    # DOMINANCE: what fraction does each partition contribute?
    if n <= 8:
        identity_V = [r for r in rows if r['is_identity']][0]['V']
        print(f"\n  DOMINANCE:")
        print(f"    Identity [1^{n}] contributes {identity_V}/{total_V} = "
              f"{100*identity_V/total_V:.2f}% of V_sum")

        # Top 3 non-identity contributions
        non_id = sorted([r for r in rows if not r['is_identity'] and r['V'] > 0],
                       key=lambda r: -r['V'])
        for i, r in enumerate(non_id[:5]):
            print(f"    #{i+1} non-identity: {r['ct']} contributes {r['V']}/{total_V} = "
                  f"{100*r['V']/total_V:.2f}%")

    # The T - twin breakdown per partition type
    if n <= 8:
        print(f"\n  T - twin PER PARTITION TYPE:")
        for r in rows:
            diff = r['T'] - r['twin']
            if diff != 0:
                print(f"    {str(r['ct']):>20}: T_ct-twin_ct = {diff:>12} "
                      f"(= {r['T']:>10} - {r['twin']:>10})")

# ============================================================
# THE KEY: CAN WE EXPRESS THE RESIDUAL BY PARTITION?
# ============================================================

print(f"\n{'='*80}")
print("  THE RESIDUAL BREAKDOWN QUESTION")
print(f"{'='*80}")
print("""
  E = (T - twin_SL - residual) / 2

  We know T and twin_SL break down by partition type.
  The DIFFERENCE T - twin_SL also breaks down:

  (T - twin)_ct = ccs × C(f,2) × 2^po × [1 - 2^{1-k+1}/2^{po}]
                = ccs × C(f,2) × [2^po - 2^{po-k+2}]
                = ccs × C(f,2) × 2^{po-k+2} × [2^{k-2} - 1]

  For identity [1^n]: k=n, so [2^{n-2} - 1] × ccs × C(n,2) × 2^{m-n+2}
  For single 3-cycle [3,1^{n-3}]: k=n-2, [2^{n-4} - 1] × ...

  This means T - twin = Σ partition_contributions, each computable.
  And 2E = (T - twin) - residual.

  The residual is the PART that doesn't factor by partition type.
  It comes from MULTI-EDGE collisions: two different partitions
  contributing to the SAME metagraph edge. This is intrinsically
  non-factorizable.

  BUT: maybe we can bound the residual by bounding the collision rate.
""")

# Compute the (T-twin) per partition and check ratio with residual
print(f"  {'n':>3} {'T-twin':>12} {'residual':>10} {'ratio':>8}")
for n in range(3, 10):
    nf = factorial(n)
    T_sum = 0; twin_sum = 0
    for ct in gen_odd_parts(n):
        cc = ccs(n, ct); p = po(ct); k = len(ct); f = ct.count(1)
        T_sum += cc * comb(f, 2) * pow(2, p) if f >= 2 else 0
        twin_sum += cc * 2 * comb(f, 2) * pow(2, p - k + 1) if f >= 2 else 0
    T = T_sum // nf; twin = twin_sum // nf
    diff = T - twin
    E = known_E.get(n, 0)
    res = diff - 2*E if E else 0
    ratio = res / diff if diff > 0 else 0
    print(f"  {n:>3} {diff:>12} {res:>10} {ratio:>8.4f}")

print("\nDONE.")
