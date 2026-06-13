#!/usr/bin/env python3
"""
gn_merged_burnside_s214.py — Burnside theory for the merged metagraph G_n/Z_2
opus-2026-03-22-S214

The complement-merged graph G_n/Z_2 has:
  - Vertices: SC classes (unchanged) + merged NS complement pairs
  - |V_merged| = (V_n + SC_n) / 2
  - Symmetry group: S_n × Z/2Z (2n! elements)

KEY INSIGHT: the transition orbit count for the merged graph is
  T_merged = (T_n + T_n^anti) / 2
where T_n^anti = (1/n!) Σ_σ Fix_anti(σ) × FA(σ)
      Fix_anti(σ) = #{T : σ(T) = T^op} = #{T with σ as anti-automorphism}

Computing Fix_anti(σ) requires understanding anti-automorphism orbits on arcs.
"""

import sys
from math import comb, factorial, gcd
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  BURNSIDE THEORY FOR G_n/Z_2 (MERGED METAGRAPH)")
print("  opus-2026-03-22-S214")
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
    """# tournaments fixed by σ with cycle type ct."""
    for c in ct:
        if c % 2 == 0: return 0
    exp = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): exp += gcd(ct[i], ct[j])
    return 2**exp

def fix_arcs(ct):
    """# arc positions {u,v} fixed by σ with cycle type ct.
    = C(f,2) + t where f = # fixed points, t = # 2-cycles."""
    f = ct.count(1)
    t = ct.count(2)
    return comb(f, 2) + t

# ============================================================================
# ANTI-AUTOMORPHISM FIXED POINT COUNT
# ============================================================================

def fix_anti_tournaments(ct, n):
    """# tournaments T such that σ(T) = T^op for σ with cycle type ct.

    σ is an anti-automorphism of T iff T(σi, σj) = 1 - T(i,j) for all i≠j.

    The action of σ on directed arcs (i→j) creates orbits.
    In each orbit, the T-values must alternate: v, 1-v, v, 1-v, ...

    An orbit of directed arcs of length L under σ:
    - If L is even: consistent alternation, 1 free choice
    - If L is odd: inconsistency, no valid assignment → Fix_anti = 0

    But we work with UNDIRECTED arcs {i,j} (each determines BOTH directions).
    Need to count orbits of σ on undirected arcs and check parity of each orbit.
    """

    # Build a representative permutation with the given cycle type
    # Cycle type ct = [c1, c2, ...] (sorted descending)
    perm = [0] * n
    pos = 0
    cycles = []
    for c in ct:
        cycle = list(range(pos, pos + c))
        cycles.append(cycle)
        for i in range(c):
            perm[cycle[i]] = cycle[(i+1) % c]
        pos += c

    # Find orbits of σ on directed arcs (i, j) with i ≠ j
    # Two types: (i,j) with the value T(i,j), and (j,i) with value 1-T(i,j)
    # Under the anti-aut condition: T(σi, σj) = 1 - T(i,j)
    # So (i,j) maps to (σi, σj) with value 1-T(i,j) = T(σi, σj)... wait:
    # T(σi, σj) = 1 - T(i,j). So the directed arc (i,j) with value v
    # maps to directed arc (σi, σj) with value 1-v.

    # For UNDIRECTED arcs {i,j}: T(i,j) = v implies T(j,i) = 1-v.
    # Under σ: {i,j} maps to {σi, σj}. The value T(σi, σj) = 1-T(i,j) = 1-v.
    # And T(σj, σi) = 1 - T(σi, σj) = v.
    # So the undirected arc value FLIPS under one application of σ.
    # After k applications: value = v if k even, 1-v if k odd.
    # Orbit closes after L steps. Need value after L steps = original value.
    # → L must be even for consistency.

    # So: find orbits of σ on undirected arcs {i,j}.
    # If ALL orbits have even length: Fix_anti = 2^(# orbits)
    #   (because each orbit contributes 1 free binary choice)
    # If ANY orbit has odd length: Fix_anti = 0

    # Wait, that's not quite right. Each orbit of undirected arcs of length L_undir:
    # The undirected arc value flips each step.
    # After L_undir steps, it's back to original. Need L_undir steps with alternation.
    # The value after L_undir steps = (-1)^L_undir * original value (mod 2 flip).
    # For consistency: (-1)^L_undir = 1, i.e., L_undir is even.

    # Actually for mod-2 (binary flip): after L_undir flips, value = v ⊕ L_undir.
    # Wait, the value flips each step: v, 1-v, v, 1-v, ...
    # After L_undir steps: v if L_undir even, 1-v if L_undir odd.
    # For the orbit to be consistent: need L_undir even.
    # Each even-length orbit contributes 1 free choice (the starting value v).
    # Fix_anti = 2^(# orbits of even length) if ALL orbits are even.
    # = 0 if ANY orbit is odd.

    # Compute orbits of σ on undirected arcs
    visited = set()
    n_orbits = 0
    all_even = True

    for i in range(n):
        for j in range(i+1, n):
            if (i, j) in visited:
                continue
            # Trace orbit of {i, j}
            orbit_len = 0
            a, b = i, j
            while True:
                key = (min(a, b), max(a, b))
                if key in visited and orbit_len > 0:
                    break
                visited.add(key)
                orbit_len += 1
                a, b = perm[a], perm[b]
                a, b = min(a, b), max(a, b)
                if (a, b) == (min(i,j), max(i,j)):
                    break

            n_orbits += 1
            if orbit_len % 2 != 0:
                all_even = False

    if not all_even:
        return 0
    return 2 ** n_orbits

# ============================================================================
# VERIFICATION: compute Fix_anti by brute force for small n
# ============================================================================

def fix_anti_brute(sigma, n):
    """Brute force: count tournaments T on [n] with sigma(T) = T^op."""
    m = comb(n, 2)
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    count = 0
    for bits in range(1 << m):
        # Build T
        adj = [[0]*n for _ in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        # Check: T(sigma[i], sigma[j]) = 1 - T(i, j) for all i != j
        ok = True
        for i in range(n):
            for j in range(n):
                if i == j: continue
                if adj[sigma[i]][sigma[j]] != 1 - adj[i][j]:
                    ok = False
                    break
            if not ok: break
        if ok:
            count += 1
    return count

print("\n  VERIFICATION: Fix_anti formula vs brute force")
print("-" * 60)

for n in range(3, 7):
    print(f"\n  n={n}:")
    for ct in gen_partitions(n):
        ct_list = list(ct)
        # Build representative permutation
        perm = [0] * n
        pos = 0
        for c in ct_list:
            for i in range(c):
                perm[pos + i] = pos + (i + 1) % c
            pos += c

        formula_val = fix_anti_tournaments(ct_list, n)

        # Brute force for n <= 5
        if n <= 5:
            brute_val = fix_anti_brute(perm, n)
            match = "✓" if formula_val == brute_val else "✗"
            print(f"    {str(ct_list):>20s}: formula={formula_val:6d}, brute={brute_val:6d} {match}")
        else:
            print(f"    {str(ct_list):>20s}: formula={formula_val:6d}")

# ============================================================================
# COMPUTE T_n, T_n^anti, T_merged for all n
# ============================================================================

print("\n" + "=" * 80)
print("  MAIN COMPUTATION: T_n, T_n^anti, T_merged")
print("=" * 80)

for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)

    T_sum = 0
    V_sum = 0
    T_anti_sum = 0
    SC_sum = 0

    for ct in gen_partitions(n):
        ct_list = list(ct)
        ct_count = cycle_type_count(n, ct_list)

        # Standard Burnside
        ft = fix_tournaments(ct_list)
        fa = fix_arcs(ct_list)
        V_sum += ct_count * ft
        T_sum += ct_count * ft * fa

        # Anti-automorphism
        ft_anti = fix_anti_tournaments(ct_list, n)
        T_anti_sum += ct_count * ft_anti * fa
        SC_sum += ct_count * ft_anti

    V_n = V_sum // nfact
    T_n = T_sum // nfact
    SC_n = SC_sum // nfact
    T_anti = T_anti_sum // nfact

    V_merged = (V_n + SC_n) // 2
    T_merged = (T_n + T_anti) // 2
    D_merged = V_merged * m - T_merged

    print(f"  n={n:2d}: V={V_n:>12d}, SC={SC_n:>8d}, V_merged={(V_n+SC_n)//2:>10d}")
    print(f"         T={T_n:>12d}, T_anti={T_anti:>8d}, T_merged={T_merged:>10d}")
    print(f"         m={m:3d}, D_merged={D_merged:>10d}, T_merged/V_merged={T_merged/V_merged:.4f}")

# ============================================================================
# COMPARE WITH KNOWN MERGED EDGE COUNTS
# ============================================================================

print("\n" + "=" * 80)
print("  COMPARISON WITH KNOWN DATA")
print("=" * 80)

# Known from kind-pasteur:
E_merged = {3: 1, 4: 3, 5: 21, 6: 143}
V_merged_known = {3: 2, 4: 3, 5: 10, 6: 34}
E_original = {3: 1, 4: 5, 5: 30, 6: 290, 7: 4086}

for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)

    T_sum = 0; V_sum = 0; T_anti_sum = 0; SC_sum = 0
    for ct in gen_partitions(n):
        ct_list = list(ct)
        ct_count = cycle_type_count(n, ct_list)
        ft = fix_tournaments(ct_list)
        fa = fix_arcs(ct_list)
        ft_anti = fix_anti_tournaments(ct_list, n)
        V_sum += ct_count * ft
        T_sum += ct_count * ft * fa
        T_anti_sum += ct_count * ft_anti * fa
        SC_sum += ct_count * ft_anti

    V_n = V_sum // nfact
    T_n = T_sum // nfact
    SC_n = SC_sum // nfact
    T_anti_n = T_anti_sum // nfact
    Vm = (V_n + SC_n) // 2
    Tm = (T_n + T_anti_n) // 2

    if n in E_merged:
        Em = E_merged[n]
        ratio_m = Tm / Em if Em > 0 else float('inf')
        print(f"  n={n}: V_merged={Vm}, E_merged={Em}, T_merged={Tm}, T_merged/E_merged={ratio_m:.4f}")
    elif n in E_original:
        print(f"  n={n}: V_merged={Vm}, T_merged={Tm}, E_merged=?")
        # Predict: if T_merged/E_merged ≈ some constant...
    else:
        print(f"  n={n}: V_merged={Vm}, T_merged={Tm}")

# ============================================================================
# THE MERGED T/E RATIO
# ============================================================================

print("\n" + "=" * 80)
print("  T_merged / E_merged RATIO (for formula hunting)")
print("=" * 80)

for n in sorted(E_merged):
    m = comb(n, 2)
    nfact = factorial(n)
    T_sum = 0; V_sum = 0; T_anti_sum = 0; SC_sum = 0
    for ct in gen_partitions(n):
        ct_list = list(ct)
        ct_count = cycle_type_count(n, ct_list)
        ft = fix_tournaments(ct_list)
        fa = fix_arcs(ct_list)
        ft_anti = fix_anti_tournaments(ct_list, n)
        V_sum += ct_count * ft
        T_sum += ct_count * ft * fa
        T_anti_sum += ct_count * ft_anti * fa
        SC_sum += ct_count * ft_anti

    V_n = V_sum // nfact
    T_n = T_sum // nfact
    SC_n = SC_sum // nfact
    T_anti_n = T_anti_sum // nfact
    Vm = (V_n + SC_n) // 2
    Tm = (T_n + T_anti_n) // 2
    Em = E_merged[n]

    excess = Tm - 2*Em

    print(f"  n={n}: T_m={Tm}, E_m={Em}, T_m/E_m={Tm/Em:.4f}, T_m-2E_m={excess}")

# ============================================================================
# PREDICT E_merged(7) and then derive E_original(7) for verification
# ============================================================================

print("\n" + "=" * 80)
print("  PREDICTIONS")
print("=" * 80)

# Compute T_merged for n=7
n = 7
m = comb(n, 2)
nfact = factorial(n)
T_sum = 0; V_sum = 0; T_anti_sum = 0; SC_sum = 0
for ct in gen_partitions(n):
    ct_list = list(ct)
    ct_count = cycle_type_count(n, ct_list)
    ft = fix_tournaments(ct_list)
    fa = fix_arcs(ct_list)
    ft_anti = fix_anti_tournaments(ct_list, n)
    V_sum += ct_count * ft
    T_sum += ct_count * ft * fa
    T_anti_sum += ct_count * ft_anti * fa
    SC_sum += ct_count * ft_anti

V_n = V_sum // nfact
T_n = T_sum // nfact
SC_n = SC_sum // nfact
T_anti_n = T_anti_sum // nfact
Vm7 = (V_n + SC_n) // 2
Tm7 = (T_n + T_anti_n) // 2

print(f"  n=7: V_merged={Vm7}, T_merged={Tm7}, SC={SC_n}")
print(f"  Known: E_original(7)=4086, V_original(7)=456")

# T_merged/E_merged ratios:
# Let's see if this has a cleaner pattern than T/E for original
ratios_m = [None, None, None]
for nn in [3,4,5,6]:
    mm = comb(nn, 2)
    nf = factorial(nn)
    Ts = 0; Vs = 0; Tas = 0; SCs = 0
    for ct in gen_partitions(nn):
        ctl = list(ct)
        ctc = cycle_type_count(nn, ctl)
        ft_ = fix_tournaments(ctl)
        fa_ = fix_arcs(ctl)
        fta = fix_anti_tournaments(ctl, nn)
        Vs += ctc * ft_
        Ts += ctc * ft_ * fa_
        Tas += ctc * fta * fa_
        SCs += ctc * fta
    Vn_ = Vs // nf; Tn_ = Ts // nf; Tan_ = Tas // nf; SCn_ = SCs // nf
    Vmm = (Vn_ + SCn_) // 2; Tmm = (Tn_ + Tan_) // 2
    Em_ = E_merged[nn]
    ratios_m.append(Tmm / Em_)

print(f"\n  T_merged/E_merged ratios: {[f'{r:.4f}' if r else '?' for r in ratios_m[3:7]]}")
print(f"  T_original/E_original: {[T_n_orig/E_original[nn] for nn in [3,4,5,6] for T_n_orig in [0]]}")

# Compute T_n_orig for comparison
for nn in [3,4,5,6,7]:
    mm = comb(nn, 2)
    nf = factorial(nn)
    Ts = 0
    for ct in gen_partitions(nn):
        ctl = list(ct)
        ctc = cycle_type_count(nn, ctl)
        ft_ = fix_tournaments(ctl)
        fa_ = fix_arcs(ctl)
        Ts += ctc * ft_ * fa_
    Tn_ = Ts // nf
    if nn in E_original:
        print(f"  n={nn}: T_orig/E_orig = {Tn_/E_original[nn]:.4f}")

# ============================================================================
# THE KEY SEQUENCES
# ============================================================================

print("\n" + "=" * 80)
print("  ALL SEQUENCES")
print("=" * 80)

V_seq = []; SC_seq = []; Vm_seq = []; T_seq = []; Tanti_seq = []; Tm_seq = []

for n in range(3, 14):
    m = comb(n, 2)
    nfact = factorial(n)
    T_sum = 0; V_sum = 0; T_anti_sum = 0; SC_sum = 0
    for ct in gen_partitions(n):
        ct_list = list(ct)
        ct_count = cycle_type_count(n, ct_list)
        ft = fix_tournaments(ct_list)
        fa = fix_arcs(ct_list)
        ft_anti = fix_anti_tournaments(ct_list, n)
        V_sum += ct_count * ft
        T_sum += ct_count * ft * fa
        T_anti_sum += ct_count * ft_anti * fa
        SC_sum += ct_count * ft_anti

    V_seq.append(V_sum // nfact)
    SC_seq.append(SC_sum // nfact)
    Vm_seq.append((V_sum + SC_sum) // (2 * nfact))
    T_seq.append(T_sum // nfact)
    Tanti_seq.append(T_anti_sum // nfact)
    Tm_seq.append((T_sum + T_anti_sum) // (2 * nfact))

print(f"  V_n (A000568):  {V_seq}")
print(f"  SC_n:           {SC_seq}")
print(f"  V_merged:       {Vm_seq}")
print(f"  T_n:            {T_seq}")
print(f"  T_anti:         {Tanti_seq}")
print(f"  T_merged:       {Tm_seq}")

# Deficit in merged graph
Dm_seq = [Vm_seq[i] * comb(i+3, 2) - Tm_seq[i] for i in range(len(Vm_seq))]
print(f"  D_merged:       {Dm_seq}")

print(f"\n  T_merged/V_merged: {[Tm_seq[i]/Vm_seq[i] for i in range(len(Vm_seq))]}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
