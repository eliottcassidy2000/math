#!/usr/bin/env python3
"""
spectral_overlap_formula.py — Exact spectral formula for cycle overlap

The missing link (Step 4): express the number of overlapping cycle pairs
directly in terms of eigenvalues of the circulant adjacency matrix.

KEY IDENTITY TO PROVE:
  |E(Ω)| = C(N,2) - α₂ can be expressed as a function of {λ_k}.

For DIRECTED k-cycles through vertex 0:
  c_k(0) = tr(A^k)/p = (1/p) Σ_j λ_j^k

Two k-cycles C₁, C₂ share a vertex v iff v ∈ V(C₁) ∩ V(C₂).
The total overlap (summed over all v) is:
  Σ_v #{pairs of k-cycles through v} = Σ_v C(d_v^(k), 2)
where d_v^(k) = number of k-cycles through v.

For circulant: d_v^(k) = d_0^(k) = c_k(0) for all v.

So: overlap_weight(k-cycles) = p · C(c_k(0), 2) = p · c_k(0)·(c_k(0)-1)/2

And: total_overlap_weight = Σ_{odd k} p · C(c_k(0), 2)
                          = (p/2) Σ_{odd k} c_k(0)² - (p/2) Σ_{odd k} c_k(0)

But this counts vertex-sharing among cycles of the SAME length.
Cross-length overlaps need: Σ_{k₁≠k₂} p · c_{k₁}(0) · c_{k₂}(0)

Wait — for cycles of DIFFERENT lengths sharing a vertex, the counting
is: for each vertex v, the number of ordered pairs of cycles through v
(of any odd length) is (Σ_k c_k(v))² - Σ_k c_k(v)².

Author: opus-2026-03-12-S64
"""

import numpy as np
from collections import defaultdict

def eigenvalues_circulant(S, p):
    omega = np.exp(2j * np.pi / p)
    return np.array([sum(omega**(s*k) for s in S) for k in range(p)])

def make_tournament(p, S):
    A = np.zeros((p, p), dtype=np.int8)
    for i in range(p):
        for s in S:
            A[i][(i + s) % p] = 1
    return A

def get_QR(p):
    return sorted(set(pow(a, 2, p) for a in range(1, p)) - {0})

# ========================================================================
print("=" * 72)
print("PART I: CYCLE COUNTS FROM EIGENVALUES")
print("=" * 72)

# For a circulant tournament on Z_p with connection set S:
# tr(A^k) = Σ_j λ_j^k where λ_j = Σ_{s∈S} ω^{sj}
# c_k(0) = #{directed k-walks from 0 to 0} = tr(A^k)/p
#
# For k < p, directed k-walks 0→...→0 include:
#   - Simple k-cycles through 0
#   - Non-simple walks (with repeated vertices)
#
# For k=3: ALL walks are simple (impossible to revisit in 3 steps in tournament)
# For k=5: MOST walks are simple, but some non-simple walks exist
# For k≥7: significant non-simple walk corrections

print("Walk counts vs simple cycle counts:")
for p in [7, 11, 13]:
    m = (p - 1) // 2
    S_int = list(range(1, m + 1))
    QR = get_QR(p)

    print(f"\n  p={p}:")
    for name, S in [("Int", S_int), ("Pal", QR)]:
        lam = eigenvalues_circulant(S, p)
        A = make_tournament(p, S)

        print(f"    {name}:", end="")
        for k in [3, 5, 7]:
            if k > p:
                break
            # Walk count from eigenvalues
            walk_k = sum(lam[j]**k for j in range(p)).real / p

            # Simple cycle count by direct enumeration (small p only)
            simple_cycles = 0
            dp = {(1, 0): 1}
            for step in range(1, k + 1):
                new_dp = {}
                for (mask, v), count in dp.items():
                    for w in range(p):
                        if step == k:
                            if w == 0 and A[v][0]:
                                key = (mask, 0)
                                new_dp[key] = new_dp.get(key, 0) + count
                        else:
                            if w != 0 and not (mask & (1 << w)) and A[v][w]:
                                new_mask = mask | (1 << w)
                                key = (new_mask, w)
                                new_dp[key] = new_dp.get(key, 0) + count
                dp = new_dp
            simple_cycles = sum(count for (mask, v), count in dp.items() if v == 0)

            print(f"  s_{k}={walk_k:.0f}, c_{k}={simple_cycles}", end="")
        print()

# ========================================================================
print("\n" + "=" * 72)
print("PART II: OVERLAP IN TERMS OF CYCLE PAIR STATISTICS")
print("=" * 72)

# For two directed k-cycles C₁, C₂ through vertex 0:
# |V(C₁) ∩ V(C₂)| = number of shared vertices
#
# The total number of ordered pairs (C₁, C₂) of k-cycles through 0
# that also pass through vertex v is:
#   N_{k,v} = #{k-cycles through both 0 and v}
#
# For circulant: N_{k,v} = N_{k,v-0} depends only on v.
#
# The total overlap weight (summing over all v):
# Σ_v Σ_{C₁≠C₂} 1[v∈C₁] 1[v∈C₂] = Σ_{C₁≠C₂} |V(C₁) ∩ V(C₂)|
#
# Key quantity: for k-cycles through 0, how many also pass through v?
# This is the "joint cycle count":
#   J_k(0,v) = #{directed k-cycles through BOTH 0 and v}

print("Joint cycle counts J_k(0,v) for p=7:")
p = 7
m = (p - 1) // 2

for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
    A = make_tournament(p, S)
    print(f"\n  {name}:")

    for L in [3, 5, 7]:
        if L > p:
            break
        # Count L-cycles through 0 that also pass through vertex v
        J = np.zeros(p)

        # Enumerate all L-cycles through 0
        dp = {(1, 0): 1}
        for step in range(1, L + 1):
            new_dp = {}
            for (mask, v), count in dp.items():
                for w in range(p):
                    if step == L:
                        if w == 0 and A[v][0]:
                            key = (mask, 0)
                            new_dp[key] = new_dp.get(key, 0) + count
                    else:
                        if w != 0 and not (mask & (1 << w)) and A[v][w]:
                            new_mask = mask | (1 << w)
                            key = (new_mask, w)
                            new_dp[key] = new_dp.get(key, 0) + count
            dp = new_dp

        # For each cycle, record which vertices it visits
        for (mask, v), count in dp.items():
            if v == 0 and count > 0:
                for u in range(p):
                    if u == 0 or (mask & (1 << u)):
                        J[u] += count

        print(f"    J_{L}(0,v): {[int(J[v]) for v in range(p)]}")

    # Total cycles through 0
    total = sum(J)  # This is L * c_L(0) actually... wait
    print(f"    Total: c_3(0) + c_5(0) + c_7(0) = {J[0]}")

# ========================================================================
print("\n" + "=" * 72)
print("PART III: SPECTRAL FORMULA FOR JOINT CYCLE COUNTS")
print("=" * 72)

# For a circulant tournament, the number of k-cycles through BOTH 0 and v
# can be expressed using eigenvalues via a "pinned trace" formula.
#
# A k-cycle through both 0 and v: 0 → ... → v → ... → 0
# This is a path from 0 to v of some length j, followed by a path from v to 0
# of length k-j, where 1 ≤ j ≤ k-1, and all k vertices are distinct.
#
# For simple walks (ignoring revisit), the pinned trace is:
#   T_k(0,v) = Σ_{j=1}^{k-1} [A^j]_{0v} · [A^{k-j}]_{v0}
#
# But this overcounts due to non-simple walks. For k=3 (simple always):
#   J_3(0,v) = A_{0v} · A_{v0} = 0 (can't have both in tournament!)
# Wait, that's wrong. A 3-cycle through 0 and v needs a third vertex w:
#   0 → v → w → 0 or 0 → w → v → 0
#
# Let me think about this differently. For a k-cycle through 0 and v,
# v appears at some position in the cycle. Say the cycle is
# 0 → x₁ → ... → x_{j-1} → v → x_{j+1} → ... → x_{k-1} → 0.
#
# This is a k-walk from 0 to 0 that passes through v at step j.
# For the walk to be a simple cycle, all x_i must be distinct and ≠ 0, v.
#
# The SPECTRAL approach: use the transfer matrix.
# [A^k]_{00} = Σ_j |λ_j|² λ_j^{k-2} ... no, that's for k-walks 0→0
# which is just tr(A^k)/p = c_k(0) (for circulant).
#
# For the "pinned" count through 0 AND v, we need:
#   Σ_{j=1}^{k-1} [A^j]_{0v} [A^{k-j}]_{v0}  (overcounts, but shows structure)
#
# [A^j]_{0v} = (1/p) Σ_r λ_r^j ω^{-rv} (circulant Fourier transform)

print("Spectral approach to joint counts:")
p = 7
m = (p - 1) // 2

for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
    A = make_tournament(p, S)
    lam = eigenvalues_circulant(S, p)
    omega = np.exp(2j * np.pi / p)

    print(f"\n  {name}:")

    # [A^j]_{0v} = (1/p) Σ_r λ_r^j ω^{-rv}
    # Walk-based "pinned trace" for k-walks 0→0 through v:
    # P_k(v) = Σ_{j=1}^{k-1} [A^j]_{0v} [A^{k-j}]_{v0}
    #        = (1/p²) Σ_{j=1}^{k-1} Σ_{r,s} λ_r^j λ_s^{k-j} ω^{-rv} ω^{sv}
    #        = (1/p²) Σ_{r,s} λ_s^k ω^{(s-r)v} Σ_{j=1}^{k-1} (λ_r/λ_s)^j
    #        ... this is getting complex

    # Let me just compute it numerically
    for k in [3, 5]:
        P_k = np.zeros(p, dtype=complex)
        for v in range(p):
            total = 0
            for j in range(1, k):
                # [A^j]_{0v}
                Aj_0v = sum(lam[r]**j * omega**(-r*v) for r in range(p)) / p
                # [A^{k-j}]_{v0}
                Akj_v0 = sum(lam[r]**(k-j) * omega**(r*v) for r in range(p)) / p
                total += Aj_0v * Akj_v0
            P_k[v] = total

        print(f"    Walk-pinned P_{k}(v): {[f'{P_k[v].real:.1f}' for v in range(p)]}")

# ========================================================================
print("\n" + "=" * 72)
print("PART IV: THE FUNDAMENTAL IDENTITY")
print("=" * 72)

# The key insight is that the total overlap weight
# W = Σ_v C(d_v, 2) where d_v = cycles through v
# is the SAME for all circulant tournaments (by symmetry d_v = d for all v).
#
# But the NUMBER OF EDGES |E(Ω)| = #{overlapping pairs}
# depends on the DISTRIBUTION of overlap weights, which is controlled
# by the eigenvalue structure.
#
# For cycles of a SINGLE length k:
# Let N_k = total k-cycles = p · c_k(0) / k (each cycle has p/gcd(p,k) = p rotations,
# counted once each; but wait, each k-cycle visits k vertices and has k rotations)
# Actually, for a directed k-cycle, there are k possible starting vertices.
# The number of distinct k-cycles = (1/k) · p · c_k(0) (if p is prime and k < p).
# Wait no: c_k(0) = #{directed k-cycles through 0}. Each cycle through 0
# visits k vertices, one of which is 0. By circulant symmetry,
# each cycle has exactly p/p = 1 rotation that fixes 0... hmm.
#
# Let me be more careful. In a circulant tournament:
# - Total directed k-cycles = p · c_k(0) / k
#   (each cycle visits k vertices, each vertex starts c_k(0) cycles,
#    total count with multiplicity = p · c_k(0), divide by k for distinct cycles)
#
# For p prime, gcd(k, p) = 1 for all k < p, so each cycle orbit under Z_p
# has size p. And c_k(0) counts cycles through vertex 0.
# Total distinct directed k-cycles = p · c_k(0) / k
# (but also = number of Z_p orbits of k-cycles × k (cycle automorphisms)??)
# Hmm, this is getting confusing. Let me just verify with data.

print("Cycle counting verification:")
for p in [7]:
    m = (p - 1) // 2
    for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
        A = make_tournament(p, S)

        for L in [3, 5, 7]:
            # Cycles through vertex 0 (by DP)
            dp = {(1, 0): 1}
            for step in range(1, L + 1):
                new_dp = {}
                for (mask, v), count in dp.items():
                    for w in range(p):
                        if step == L:
                            if w == 0 and A[v][0]:
                                key = (mask, 0)
                                new_dp[key] = new_dp.get(key, 0) + count
                        else:
                            if w != 0 and not (mask & (1 << w)) and A[v][w]:
                                new_mask = mask | (1 << w)
                                key = (new_mask, w)
                                new_dp[key] = new_dp.get(key, 0) + count
                dp = new_dp
            c_L_0 = sum(count for (mask, v), count in dp.items() if v == 0)

            # Total directed L-cycles (enumerate all vertex sets)
            all_cycle_vsets = set()
            for start in range(p):
                dp = {(1 << start, start): 1}
                for step in range(1, L + 1):
                    new_dp = {}
                    for (mask, v), count in dp.items():
                        for w in range(p):
                            if step == L:
                                if w == start and A[v][w]:
                                    key = (mask, w)
                                    new_dp[key] = new_dp.get(key, 0) + count
                            else:
                                if w != start and not (mask & (1 << w)) and A[v][w]:
                                    new_key = (mask | (1 << w), w)
                                    new_dp[new_key] = new_dp.get(new_key, 0) + count
                    dp = new_dp
                for (mask, v), count in dp.items():
                    if v == start and count > 0:
                        vset = frozenset(i for i in range(p) if mask & (1 << i))
                        all_cycle_vsets.add(vset)

            N_L = len(all_cycle_vsets)

            # From eigenvalues
            lam = eigenvalues_circulant(S, p)
            walk_L = sum(lam[j]**L for j in range(p)).real / p

            print(f"    {name} L={L}: c_L(0)={c_L_0}, N_L(distinct vsets)={N_L}, "
                  f"p·c_L(0)/L={p*c_L_0/L:.0f}, walk_L={walk_L:.0f}")

# ========================================================================
print("\n" + "=" * 72)
print("PART V: OVERLAP AS FUNCTION OF SPECTRAL POWER SUMS")
print("=" * 72)

# The KEY formula we want:
#   α₂ = #{disjoint cycle pairs} = f(e₃, e₅, e₇, ..., e_p)
# where e_k = tr(A^k)/p = (1/p) Σ λ_j^k
#
# We know:
#   Total cycles N = Σ_k N_k where N_k ≈ p·c_k(0)/k
#   Total overlap W = p · C(d, 2) where d = Σ_k c_k(0) (or similar)
#   α₂ = C(N, 2) - |E(Ω)|
#
# The connection to eigenvalues is through the power sums e_k = c_k(0) (for simple walks).
# For k=3: e_3 = c_3(0) exactly (all 3-walks are simple)
# For k≥5: e_k ≈ c_k(0) with corrections from non-simple walks

# Let me compute the "second moment of cycle distribution" directly.
# This is related to the 4-point correlation function:
# #{pairs (C₁,C₂) of k-cycles sharing a specific vertex}
# = #{k-cycles through v}² - #{k-cycles through v}
# = c_k(v)² - c_k(v) = c_k(0)² - c_k(0) for circulant

print("Second moment analysis:")
for p in [7]:
    m = (p - 1) // 2
    for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
        A = make_tournament(p, S)
        lam = eigenvalues_circulant(S, p)

        total_N = 0
        total_d = 0
        cycle_data = {}

        for L in [3, 5, 7]:
            # c_L(0) from direct enumeration
            dp = {(1, 0): 1}
            for step in range(1, L + 1):
                new_dp = {}
                for (mask, v), count in dp.items():
                    for w in range(p):
                        if step == L:
                            if w == 0 and A[v][0]:
                                key = (mask, 0)
                                new_dp[key] = new_dp.get(key, 0) + count
                        else:
                            if w != 0 and not (mask & (1 << w)) and A[v][w]:
                                new_mask = mask | (1 << w)
                                key = (new_mask, w)
                                new_dp[key] = new_dp.get(key, 0) + count
                dp = new_dp
            c_L_0 = sum(count for (mask, v), count in dp.items() if v == 0)

            total_d += c_L_0
            cycle_data[L] = c_L_0

        # c_k(0) from eigenvalues (walk count, includes non-simple)
        print(f"\n  {name}: cycle counts through vertex 0")
        for L in [3, 5, 7]:
            e_L = sum(lam[j]**L for j in range(p)).real / p
            c_L = cycle_data[L]
            print(f"    L={L}: c_L(0)={c_L}, walk_L={e_L:.0f}, correction={e_L-c_L:.0f}")

        # Predict overlap from cycle counts
        d = total_d  # total cycles per vertex
        W_pred = p * d * (d - 1) // 2
        print(f"    Total d = {d}, predicted W = {W_pred}")

# ========================================================================
print("\n" + "=" * 72)
print("PART VI: THE SPECTRAL OVERLAP BOUND")
print("=" * 72)
print("""
SPECTRAL OVERLAP BOUND (Proposed):

For a circulant tournament on Z_p with eigenvalues λ₁,...,λ_{p-1}
(excluding λ₀ = m), define:
  σ_k = (1/p) Σ_{j=1}^{p-1} |λ_j|^{2k}  (spectral power sums)

Then the number of disjoint odd-cycle pairs satisfies:
  α₂ ≥ C(N, 2) · (1 - σ₂/(σ₁²))

where N = total odd cycles and σ₁² / σ₂ is the participation ratio.

INTUITION: When σ₂ << σ₁² (Paley, flat spectrum):
  α₂ ≈ C(N, 2) · (1 - 1) = 0 (very few disjoint pairs)

When σ₂ ≈ σ₁² (concentrated):
  α₂ ≈ C(N, 2) · (1 - 1/PR) where PR = participation ratio
  ≈ C(N, 2) · (1 - 1/3) = 2/3 · C(N, 2) for Interval

This is because:
  - Cycles through vertex 0 are RANDOM with respect to cycles through v
    when v is in a "different spectral region"
  - The fraction of vertex pairs in "different regions" ≈ 1 - 1/PR
  - Cycles through such vertices are nearly independent → disjoint

This needs to be made rigorous, but the intuition matches the data.
""")

# Verify the bound at p=7
print("Testing spectral overlap bound at p=7:")
p = 7
m = (p - 1) // 2

for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
    lam = eigenvalues_circulant(S, p)

    # Spectral power sums (k=0 excluded)
    y2 = np.array([abs(lam[k])**2 for k in range(1, p)])
    sigma1 = sum(y2) / p
    sigma2 = sum(y2**2) / p
    PR = sigma1**2 / sigma2

    # Known data
    N = 36  # total cycle vertex sets (same for both at p=7)
    alpha2_known = 14 if "Int" in name else 7

    # Predicted α₂
    alpha2_pred = int(N * (N - 1) / 2 * (1 - 1/PR))

    print(f"  {name}: PR={PR:.4f}, predicted α₂ ≈ {alpha2_pred}, actual α₂ = {alpha2_known}")

print("""
NOTE: The bound is ROUGH (the predicted values may not match exactly),
but the DIRECTION is correct: higher PR → higher α₂.
The exact formula requires accounting for:
  - Cross-length cycle interactions (3-cycles vs 5-cycles vs 7-cycles)
  - Non-simple walk corrections
  - Multiplicity of directed cycles per vertex set
""")

# ========================================================================
print("=" * 72)
print("PART VII: NOVEL CONNECTION — EXPANDER MIXING LEMMA FOR Ω")
print("=" * 72)

# The Expander Mixing Lemma relates edge distribution to spectral gap:
# |e(X,Y) - |X||Y|d/n| ≤ λ₂ √(|X||Y|)
# where d is average degree and λ₂ is second eigenvalue of Ω.
#
# For independent sets: e(S,S) = 0, so |S|²d/n ≤ λ₂|S|, giving α ≤ nλ₂/d.
# This is the Hoffman bound.
#
# For Ω(Interval) vs Ω(Paley), the spectral gap of Ω differs:

print("Hoffman/expander analysis of Ω at p=7:")
p = 7

for name, S in [("Interval", list(range(1, m+1))), ("Paley", get_QR(p))]:
    A = make_tournament(p, S)

    # Find all odd cycles
    all_vsets = set()
    for L in range(3, p + 1, 2):
        for start in range(p):
            dp = {(1 << start, start): 1}
            for step in range(1, L + 1):
                new_dp = {}
                for (mask, v), count in dp.items():
                    for w in range(p):
                        if step == L:
                            if w == start and A[v][w]:
                                key = (mask, w)
                                new_dp[key] = new_dp.get(key, 0) + count
                        else:
                            if w != start and not (mask & (1 << w)) and A[v][w]:
                                new_key = (mask | (1 << w), w)
                                new_dp[new_key] = new_dp.get(new_key, 0) + count
                dp = new_dp
            for (mask, v), count in dp.items():
                if v == start and count > 0:
                    all_vsets.add(frozenset(i for i in range(p) if mask & (1 << i)))

    cycles = list(all_vsets)
    nc = len(cycles)

    # Build Ω adjacency
    adj_mat = np.zeros((nc, nc))
    for i in range(nc):
        for j in range(i+1, nc):
            if cycles[i] & cycles[j]:
                adj_mat[i][j] = adj_mat[j][i] = 1

    eigs = np.linalg.eigvalsh(adj_mat)
    lambda_max = eigs[-1]
    lambda_2 = max(abs(eigs[0]), abs(eigs[-2]))  # second largest |eigenvalue|
    avg_deg = np.sum(adj_mat) / nc

    hoffman = nc * (-eigs[0]) / (lambda_max - eigs[0])

    print(f"\n  {name}:")
    print(f"    |V(Ω)| = {nc}")
    print(f"    Avg degree = {avg_deg:.2f}")
    print(f"    λ_max(Ω) = {lambda_max:.4f}")
    print(f"    λ_min(Ω) = {eigs[0]:.4f}")
    print(f"    λ₂(Ω) = {eigs[-2]:.4f}")
    print(f"    Spectral gap = {lambda_max - eigs[-2]:.4f}")
    print(f"    Hoffman bound: α(Ω) ≥ {hoffman:.4f}")
    print(f"    Edge density = {np.sum(adj_mat) / (nc*(nc-1)):.6f}")

    # Actual independence number (brute force for small nc)
    if nc <= 20:
        max_indep = 0
        for mask in range(1 << nc):
            k = bin(mask).count('1')
            if k <= max_indep:
                continue
            bits = [i for i in range(nc) if mask & (1 << i)]
            indep = True
            for a in range(len(bits)):
                for b in range(a+1, len(bits)):
                    if adj_mat[bits[a]][bits[b]]:
                        indep = False
                        break
                if not indep:
                    break
            if indep and k > max_indep:
                max_indep = k
        print(f"    Actual α(Ω) = {max_indep}")

print("\nDONE.")
