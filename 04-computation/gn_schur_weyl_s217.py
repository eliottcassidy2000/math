#!/usr/bin/env python3
"""
gn_schur_weyl_s217.py — Schur-Weyl decomposition of T_n
opus-2026-03-22-S217

KEY INSIGHT: T_n = ⟨χ_V, χ_W⟩ where:
  χ_V = tournament permutation character (Fix(σ) for odd-cycle σ, 0 otherwise)
  χ_W = arc-position permutation character (C(f(σ),2) + t(σ))

The arc-position module W decomposes as:
  W = V_{(n)} ⊕ V_{(n-1,1)} ⊕ V_{(n-2,2)}

So: T_n = m_{(n)}(V) + m_{(n-1,1)}(V) + m_{(n-2,2)}(V)

where m_λ(V) = multiplicity of irrep λ in the tournament module.

We already know: m_{(n)} = V_n (number of iso classes).

Computing m_{(n-1,1)} and m_{(n-2,2)} via Burnside gives a DECOMPOSITION
of T_n into three irreducible components.

Similarly for the anti-automorphism side:
T_anti = m_{(n)}^anti + m_{(n-1,1)}^anti + m_{(n-2,2)}^anti

This gives the complete Schur-Weyl picture of the transition orbit count.
"""

import sys
from math import comb, factorial, gcd
from collections import Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  SCHUR-WEYL DECOMPOSITION OF T_n")
print("  opus-2026-03-22-S217")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================

def gen_partitions(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0: yield []; return
    for first in range(min(n, max_part), 0, -1):
        for rest in gen_partitions(n - first, first):
            yield [first] + rest

def cycle_type_count(n, ct):
    c = Counter(ct)
    r = factorial(n)
    for l, k in c.items(): r //= (l**k) * factorial(k)
    return r

def fix_tournaments(ct):
    for c in ct:
        if c % 2 == 0: return 0
    exp = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): exp += gcd(ct[i], ct[j])
    return 2**exp

# ============================================================================
# S_n CHARACTER VALUES for key partitions
# ============================================================================

def chi_trivial(ct, n):
    """χ_{(n)}(σ) = 1 for all σ."""
    return 1

def chi_standard(ct, n):
    """χ_{(n-1,1)}(σ) = f(σ) - 1 where f = # fixed points."""
    f = ct.count(1)
    return f - 1

def chi_pair_rep(ct, n):
    """χ_{(n-2,2)}(σ) = character of action on 2-subsets minus trivial minus standard.
    χ_W(σ) = C(f,2) + t  (arc positions fixed by σ)
    χ_{(n-2,2)}(σ) = χ_W(σ) - χ_{(n)}(σ) - χ_{(n-1,1)}(σ)
                    = (C(f,2) + t) - 1 - (f-1)
                    = C(f,2) + t - f
                    = f(f-1)/2 + t - f
                    = f(f-3)/2 + t"""
    f = ct.count(1)
    t = ct.count(2)
    return f*(f-3)//2 + t

def chi_sign(ct, n):
    """χ_{(1^n)}(σ) = sign(σ) = (-1)^(n - #cycles)."""
    num_cycles = len(ct)
    return (-1)**(n - num_cycles)

def chi_sign_standard(ct, n):
    """χ_{(2,1^{n-2})}(σ) = sign(σ) × (f(σ) - 1)."""
    return chi_sign(ct, n) * chi_standard(ct, n)

# ============================================================================
# COMPUTE IRREP MULTIPLICITIES IN TOURNAMENT MODULE
# ============================================================================

print("\n  IRREP MULTIPLICITIES in tournament permutation module")
print("-" * 70)

print(f"\n  {'n':>3} {'V_n':>8} {'T_n':>10} {'m_(n)':>8} {'m_(n-1,1)':>10} "
      f"{'m_(n-2,2)':>10} {'sum':>8} {'check':>6}")

for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)

    # Compute multiplicities
    m_triv = 0      # m_{(n)}
    m_std = 0       # m_{(n-1,1)}
    m_pair = 0      # m_{(n-2,2)}
    T_sum = 0

    for ct in gen_partitions(n):
        ct_list = list(ct)
        ctc = cycle_type_count(n, ct_list)
        ft = fix_tournaments(ct_list)

        if ft == 0:
            continue

        m_triv += ctc * ft * chi_trivial(ct_list, n)
        m_std += ctc * ft * chi_standard(ct_list, n)
        m_pair += ctc * ft * chi_pair_rep(ct_list, n)

        fa = comb(ct_list.count(1), 2) + ct_list.count(2)
        T_sum += ctc * ft * fa

    m_triv //= nfact
    m_std //= nfact
    m_pair //= nfact
    T_n = T_sum // nfact

    check = m_triv + m_std + m_pair
    ok = "✓" if check == T_n else "✗"

    print(f"  {n:3d} {m_triv:8d} {T_n:10d} {m_triv:8d} {m_std:10d} "
          f"{m_pair:10d} {check:8d} {ok}")

# ============================================================================
# DEEPER ANALYSIS: What do the multiplicities mean?
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS OF IRREP MULTIPLICITIES")
print("=" * 80)

for n in range(3, 14):
    m_val = comb(n, 2)
    nfact = factorial(n)

    m_triv = 0; m_std = 0; m_pair = 0

    for ct in gen_partitions(n):
        ct_list = list(ct)
        ctc = cycle_type_count(n, ct_list)
        ft = fix_tournaments(ct_list)
        if ft == 0: continue
        m_triv += ctc * ft * chi_trivial(ct_list, n)
        m_std += ctc * ft * chi_standard(ct_list, n)
        m_pair += ctc * ft * chi_pair_rep(ct_list, n)

    m_triv //= nfact; m_std //= nfact; m_pair //= nfact
    T = m_triv + m_std + m_pair

    # Ratios
    print(f"\n  n={n}: V_n = m_(n) = {m_triv}")
    print(f"    m_(n-1,1) = {m_std:>10d}  ({m_std/m_triv:.4f} × V_n)")
    print(f"    m_(n-2,2) = {m_pair:>10d}  ({m_pair/m_triv:.4f} × V_n)")
    print(f"    T_n = {T:>10d}  = V + {m_std} + {m_pair}")
    print(f"    D_n = Vm-T = {m_triv*m_val-T}")

# ============================================================================
# THE SIGN CHARACTER AND OTHER IRREPS
# ============================================================================

print("\n" + "=" * 80)
print("  ADDITIONAL IRREP MULTIPLICITIES")
print("=" * 80)

print(f"\n  {'n':>3} {'V_n':>8} {'m_sign':>8} {'m_sign_std':>10} {'m_(n-2,1,1)':>12}")

for n in range(3, 14):
    nfact = factorial(n)

    m_triv = 0; m_sign = 0; m_sign_std = 0

    for ct in gen_partitions(n):
        ct_list = list(ct)
        ctc = cycle_type_count(n, ct_list)
        ft = fix_tournaments(ct_list)
        if ft == 0: continue
        m_triv += ctc * ft
        m_sign += ctc * ft * chi_sign(ct_list, n)
        m_sign_std += ctc * ft * chi_sign_standard(ct_list, n)

    m_triv //= nfact; m_sign //= nfact; m_sign_std //= nfact

    print(f"  {n:3d} {m_triv:8d} {m_sign:8d} {m_sign_std:10d}")

# ============================================================================
# KEY SEQUENCE: m_{(n-1,1)} and m_{(n-2,2)}
# ============================================================================

print("\n" + "=" * 80)
print("  KEY SEQUENCES")
print("=" * 80)

m_triv_seq = []; m_std_seq = []; m_pair_seq = []

for n in range(3, 14):
    nfact = factorial(n)
    mt = 0; ms = 0; mp = 0
    for ct in gen_partitions(n):
        ctl = list(ct)
        ctc = cycle_type_count(n, ctl)
        ft = fix_tournaments(ctl)
        if ft == 0: continue
        mt += ctc * ft
        ms += ctc * ft * chi_standard(ctl, n)
        mp += ctc * ft * chi_pair_rep(ctl, n)
    m_triv_seq.append(mt // nfact)
    m_std_seq.append(ms // nfact)
    m_pair_seq.append(mp // nfact)

print(f"  m_(n) = V_n:      {m_triv_seq}")
print(f"  m_(n-1,1):        {m_std_seq}")
print(f"  m_(n-2,2):        {m_pair_seq}")
print(f"  T_n = sum:        {[m_triv_seq[i]+m_std_seq[i]+m_pair_seq[i] for i in range(len(m_triv_seq))]}")

# Ratios
print(f"\n  m_(n-1,1)/V_n:    {[f'{m_std_seq[i]/m_triv_seq[i]:.6f}' for i in range(len(m_triv_seq))]}")
print(f"  m_(n-2,2)/V_n:    {[f'{m_pair_seq[i]/m_triv_seq[i]:.6f}' for i in range(len(m_triv_seq))]}")

# Check: m_(n-1,1) + m_(n-2,2) = D_n (deficit)
D_seq = [m_std_seq[i] + m_pair_seq[i] for i in range(len(m_triv_seq))]
print(f"\n  D_n = m_(n-1,1) + m_(n-2,2): {D_seq}")
print(f"  (should match: 2, 8, 32, 136, 664, 4352, 48096, 918208, 31568544, ...)")

# ============================================================================
# EDGE FORMULA CONNECTION
# ============================================================================

print("\n" + "=" * 80)
print("  CONNECTION TO EDGE COUNT")
print("=" * 80)

E_known = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161}

print(f"\n  The edge formula: |E| ≈ T/(2+ε) where ε → 0")
print(f"  T = V + m_(n-1,1) + m_(n-2,2)")
print(f"  |E| ≈ V/2 + m_(n-1,1)/2 + m_(n-2,2)/2  (if ε→0)")
print(f"  The 'correction' to V*m/2 comes from m_(n-1,1) and m_(n-2,2)")

for i, n in enumerate(range(3, 9)):
    if n not in E_known: continue
    V = m_triv_seq[i]
    ms = m_std_seq[i]
    mp = m_pair_seq[i]
    T = V + ms + mp
    E = E_known[n]
    m_val = comb(n, 2)

    print(f"\n  n={n}: V={V}, m_(n-1,1)={ms}, m_(n-2,2)={mp}, T={T}")
    print(f"    |E|={E}, T/2={T/2:.1f}, T/(2E)={T/(2*E):.4f}")
    print(f"    V*m/2={V*m_val/2}, V*m/2 - |E|={V*m_val/2-E}")
    print(f"    (V+ms+mp)/2 = {(V+ms+mp)/2:.1f}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
