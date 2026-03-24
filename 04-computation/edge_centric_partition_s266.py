#!/usr/bin/env python3
"""
edge_centric_partition_s266.py — Partition-level edge-centric Burnside
opus-2026-03-23-S266

kind-pasteur's edge-centric formula (S20dz) shows E(G_n) = edges through
a single arc position {0,1}, computed via Burnside on Stab({0,1}) = S_2 × S_{n-2}.

Instead of enumerating permutations of Stab, iterate over CYCLE TYPES.
The cycle types of S_2 × S_{n-2} are:
  - (action on {0,1}) × (cycle type on {2,...,n-1})
  - S_2 part: either identity [1,1] or transposition [2]
  - S_{n-2} part: any partition of n-2

For each such (σ_2, σ_{n-2}), σ acts on n vertices with known cycle type.
We need: #{T fixed by σ : [T] ≠ [T flip (0,1)]}

The hard part: determining if flipping arc (0,1) changes the iso class.
This requires understanding how arc (0,1) sits within the σ-orbits.

KEY INSIGHT: For σ that fixes both 0 and 1 (identity on S_2):
  Arc (0,1) is a FIXED ARC under σ.
  Flipping it changes the class iff the class is NOT self-neutral at (0,1).
  The number of fixed T with non-neutral (0,1) = Fix_tourn(σ) - neutral(σ).

For σ that swaps 0 and 1 (transposition on S_2):
  Arc (0,1) maps to arc (1,0) under σ.
  If σ(T) = T, then T[0][1] = T[σ(0)][σ(1)] = T[1][0] = 1 - T[0][1].
  This is IMPOSSIBLE unless the tournament allows T[0][1] ≠ T[0][1].
  Wait: if σ swaps 0 and 1, and σ(T) = T, then for arc (0,1):
  T must have T[0][1] = T[σ(0)][σ(1)] = T[1][0]. But T[1][0] = 1-T[0][1].
  So T[0][1] = 1-T[0][1], impossible. So NO tournament is fixed by σ=(01).
  Actually wait: for the combined σ = (0,1) × σ_{n-2}, the fixed-point
  condition is more complex. The arc (0,1) is part of a 2-element orbit
  {(0,1), (1,0)} under σ. For a tournament fixed by σ, we need the
  orientation to be consistent: T[0][1] = T[σ(0)][σ(1)] = T[1][0] = 1-T[0][1].
  This forces T[0][1] = 1-T[0][1], which is impossible.

  So: NO tournament is fixed by any σ that swaps 0 and 1.
  Therefore: only the IDENTITY part of S_2 contributes.

This means: the stabilizer computation reduces to Burnside on S_{n-2}!

For σ ∈ S_{n-2} (fixing vertices 0 and 1):
  Fix_tourn(σ_ext) = #{T : σ_ext(T) = T}
  where σ_ext extends σ by identity on {0,1}.

  The cycle type of σ_ext has all cycles of σ plus two fixed points 0,1.
  So cycle type = σ's cycle type + [1,1].
  For this to have odd cycles only: σ's cycle type must have odd cycles only.

  Then: Fix_tourn(σ_ext) = 2^{pair_orb(σ's ct + [1,1])}

  Among these fixed tournaments:
    neutral(σ_ext) = # where flipping (0,1) preserves the class
    non_neutral = Fix_tourn - neutral

  E(G_n) = (1/|Stab|) × Σ non_neutral
         = (1/(2(n-2)!)) × Σ_{σ ∈ S_{n-2}} non_neutral(σ)

  The sum over S_{n-2} can be done via PARTITIONS of n-2!

The remaining challenge: computing neutral(σ) for each cycle type.
This is WHERE the residual lives.
"""

from math import comb, factorial, gcd
from collections import Counter
from functools import lru_cache
import time, sys

sys.stdout.reconfigure(line_buffering=True)

@lru_cache(maxsize=None)
def _f(n): return factorial(n)

def gen_parts(n):
    result = []
    def _g(rem, mx, cur):
        if rem == 0: result.append(tuple(cur)); return
        for p in range(min(rem, mx), 0, -1): _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return result

def gen_odd_parts(n):
    result = []
    def _g(rem, mx, cur):
        if rem == 0: result.append(tuple(cur)); return
        st = min(rem, mx)
        if st % 2 == 0: st -= 1
        for p in range(st, 0, -2): _g(rem - p, p, cur + [p])
    _g(n, n, [])
    return result

def ccs_val(n, ct):
    c = Counter(ct); r = _f(n)
    for l, k in c.items(): r //= pow(l, k) * _f(k)
    return r

def pair_orb(ct):
    s = sum((c-1)//2 for c in ct)
    for i in range(len(ct)):
        for j in range(i+1, len(ct)): s += gcd(ct[i], ct[j])
    return s

def compute_edge_centric(n):
    """Edge-centric E(G_n) via partition-level Burnside on S_{n-2}.

    E = (1/(2(n-2)!)) × Σ_{σ ∈ S_{n-2}, odd cycle type of ext}
        [Fix_tourn(σ_ext) - neutral(σ_ext)]

    Where σ_ext has cycle type = σ's ct + [1,1] (adding fixed pts 0,1).

    Since we DON'T know neutral(σ_ext) analytically, we compute:
    FIX_TOTAL = Σ Fix_tourn(σ_ext) = total fixed tournaments through arc (0,1)
              = T_n (transition count, scaled appropriately)
    NEUTRAL_TOTAL = twin_SL + complex_SL

    E = (FIX_TOTAL - NEUTRAL_TOTAL) / (2(n-2)!)
    """

    stab_size = 2 * _f(n - 2)

    # Iterate over cycle types of S_{n-2}
    # Extend each by adding 2 fixed points (vertices 0 and 1)
    fix_total = 0
    twin_total = 0

    for ct_inner in gen_odd_parts(n - 2):
        # Extended cycle type: ct_inner + (1, 1)
        ct_ext = tuple(sorted(list(ct_inner) + [1, 1], reverse=True))

        # Check: all odd? ct_inner has odd parts, + two 1s → still all odd
        cc_inner = ccs_val(n - 2, ct_inner)
        po_ext = pair_orb(ct_ext)
        fix = pow(2, po_ext)
        fix_total += cc_inner * fix

        # Twin contribution: twins use the pair (0,1) as the twin pair
        # When σ fixes 0 and 1, and 0,1 are twins (identical neighborhoods),
        # then arc (0,1) is neutral.
        # Number of fixed T with twin (0,1) = 2^{po_ext - 1} (arc between twins is free)
        # Wait — for twins, the arc (0,1) direction doesn't matter
        # because 0 and 1 are interchangeable. So neutral_twin = 2^{po_ext - 1}?
        # No: twins means Aut(T) contains σ_(01) transposition.
        # The twin condition: for all v ≠ 0,1: T[0][v] = T[1][v].
        # Under σ (which fixes 0,1 and has some cycle structure on 2..n-1):
        # The twin constraint reduces the free arcs by... hmm, this is complex.

        # Simple bound: twin_per_sigma ≤ fix/2 (at most half of fixed T are neutral at (0,1))
        # Actually: for identity σ (ct_inner = [1]*(n-2)):
        #   fix = 2^{C(n,2)} (all tournaments)
        #   neutral at (0,1) = #{T where flipping (0,1) preserves class}
        #   = SL at arc (0,1) = self_loops_via_e from kind-pasteur's S20dz

        # For the partition-level formula, we can use:
        k_ext = len(ct_ext)
        f_ext = ct_ext.count(1)
        # Twin SL through this arc: must have 0,1 as a twin pair
        # Condition: T[0][v] = T[1][v] for all v in 2..n-1 fixed by σ
        # Plus: arcs within σ-orbits of {2..n-1} are free
        # The twin pair (0,1) contributes 1 free arc (the direction between them)
        # The remaining arcs: same as Fix_tourn but with n-2 vertices
        # Specifically: twin_T with (0,1) twins = # tournaments on {2..n-1} fixed by σ
        # (because the arcs from 0 and arcs from 1 are identical)
        # = Fix_tourn(σ_inner_ct) IF ct_inner has all odd cycles
        po_inner = pair_orb(ct_inner)
        twin_per = pow(2, po_inner)  # tournaments on {2..n-1} × twin pair direction
        # Actually: po_inner counts pair orbits on {2..n-1}.
        # For each such tournament on {2..n-1}, we extend to a twin-pair tournament
        # on {0,1,...,n-1} by: for each v, T[0][v]=T[1][v]=value from inner tournament,
        # and T[0][1] is free. So twin extensions = 2 × (inner tournaments fixed by σ).
        # Wait: the inner tournament is on {2..n-1} with C(n-2,2) arcs.
        # But arcs from 0 and 1 to {2..n-1} are determined by the twin condition.
        # Plus the arc (0,1) itself is free.
        # So: twin_fixed(σ) = 2^{po_inner + 1} (inner + T[0][1])
        # Actually need to count arcs more carefully.
        # Arcs of the full tournament: C(n,2) = C(n-2,2) + 2(n-2) + 1
        # Under twin: arcs between {2..n-1} = C(n-2,2) free (po_inner orbits)
        # Arcs {0,v}, {1,v}: T[0][v]=T[1][v], so n-2 arcs but PAIRED, giving n-2 free bits
        # But under σ, these n-2 arcs are grouped into orbits of σ on {2..n-1}
        # Number of orbits of σ on {2..n-1} = k_inner (= number of cycles in ct_inner)
        k_inner = len(ct_inner)
        # Arc (0,1): 1 free bit (twin doesn't constrain it, but it IS the arc we're checking)
        # So twin_fixed(σ) = 2^{po_inner + k_inner + 1}
        # ... but for neutral, we need T fixed by σ AND arc (0,1) is neutral.
        # Neutral means flipping (0,1) doesn't change the class. For twins,
        # this is always true (swapping 0↔1 preserves the tournament).
        # So all twin-fixed T have neutral (0,1).
        # neutral_twin(σ) = twin_fixed(σ) = 2^{po_inner + k_inner + 1}

        twin_neutral = pow(2, po_inner + k_inner + 1)
        twin_total += cc_inner * twin_neutral

    # E = (fix_total - neutral_total) / stab_size
    # We only have twin_neutral, not full neutral
    # E_lower = (fix_total - twin_total) / stab_size

    E_upper = fix_total // stab_size  # if no neutrals
    E_twin_bound = (fix_total - twin_total) // stab_size

    return {
        'n': n,
        'fix_total': fix_total,
        'twin_total': twin_total,
        'stab_size': stab_size,
        'E_upper': E_upper,  # = T_n (all transitions / 2)
        'E_twin_bound': E_twin_bound,
    }

if __name__ == '__main__':
    print("=" * 70)
    print("  EDGE-CENTRIC PARTITION-LEVEL BURNSIDE")
    print("  opus-2026-03-23-S266")
    print("=" * 70)

    known_E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}
    known_T = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}

    print(f"\n{'n':>3} {'fix_total':>14} {'twin_total':>14} {'E_upper':>10} "
          f"{'E_twin':>10} {'T_known':>10} {'E_known':>10}")

    for n in range(3, 13):
        t0 = time.time()
        r = compute_edge_centric(n)
        dt = time.time() - t0

        # Check: E_upper should equal T_n / 2 + something?
        # Actually E_upper = fix_total / (2*(n-2)!) = T_n (transition orbits)
        # because fix_total / (n-2)! = Σ for S_{n-2} of Fix extended
        # and dividing by 2 for the S_2 factor

        tk = known_T.get(n, '?')
        ek = known_E.get(n, '?')

        print(f"{n:>3} {r['fix_total']:>14} {r['twin_total']:>14} "
              f"{r['E_upper']:>10} {r['E_twin_bound']:>10} {str(tk):>10} {str(ek):>10} "
              f"({dt:.3f}s)")

    # Analysis
    print(f"\n  KEY INSIGHT: fix_total / stab_size should relate to T_n")
    for n in range(3, 10):
        r = compute_edge_centric(n)
        print(f"  n={n}: fix/{r['stab_size']} = {r['fix_total']//r['stab_size']}, "
              f"T_known={known_T[n]}, twin/{r['stab_size']}={r['twin_total']//r['stab_size']}")

    print("\nDONE.")
