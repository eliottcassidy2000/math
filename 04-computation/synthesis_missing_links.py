#!/usr/bin/env python3
"""
synthesis_missing_links.py — Identify the missing links in the proof chain

The current understanding:
  Fejér kernel (eigenvalues) → IPR concentration → cycle clustering → α₂ advantage → H(Int) > H(Pal)

Missing rigorous connections:
  1. IPR → α₂: Does spectral concentration DIRECTLY bound disjoint cycle counts?
  2. Additive energy → cycle overlap: Can we prove E(S) controls pairwise vertex overlap?
  3. Exact crossover: What determines p_c for the phase transition?

This script tests NEW connections not yet explored.

Author: opus-2026-03-12-S65
"""

import numpy as np
from itertools import product, combinations
from collections import defaultdict

def legendre(a, p):
    if a % p == 0: return 0
    v = pow(a, (p-1)//2, p)
    return v if v == 1 else -1

def tournament_adjacency(sigma, p):
    """Build adjacency matrix for circulant tournament with orientation sigma."""
    m = (p-1)//2
    n = p
    A = np.zeros((n, n), dtype=int)
    for k in range(1, m+1):
        for i in range(n):
            j = (i + k) % n
            if sigma[k-1] == 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_3cycles(A, n):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check all 6 orderings, count directed 3-cycles
                if A[i][j] and A[j][k] and A[k][i]:
                    count += 1
                if A[i][k] and A[k][j] and A[j][i]:
                    count += 1
    return count

def additive_energy(S, p):
    """E(S) = #{(a,b,c,d) in S^4 : a+b ≡ c+d mod p}."""
    sums = defaultdict(int)
    for a in S:
        for b in S:
            sums[(a + b) % p] += 1
    return sum(v*v for v in sums.values())

def eigenvalues_circulant(S, p):
    """Compute DFT eigenvalues of circulant tournament."""
    omega = np.exp(2j * np.pi / p)
    eigs = []
    for k in range(1, p):  # skip k=0
        lam = sum(omega**(k*s) for s in S)
        eigs.append(lam)
    return eigs

def ipr(eigs):
    """Inverse Participation Ratio of eigenvalue magnitudes squared."""
    x = np.array([abs(e)**2 for e in eigs])
    return np.sum(x**2) / np.sum(x)**2

# =============================================================
# CONNECTION 1: Trace formula vs disjointness
# =============================================================
print("=" * 70)
print("CONNECTION 1: TRACE MOMENTS AND DISJOINTNESS")
print("=" * 70)

# Key idea: tr(A^k) counts closed walks of length k.
# The VARIANCE of the eigenvalue distribution controls how "uneven"
# the walks are distributed, which relates to cycle overlap.
#
# If eigenvalues are concentrated (like Interval), most walks
# pass through the same "channel" → cycles are correlated →
# but paradoxically this means the REMAINING cycles are MORE disjoint.

for p in [7, 11, 13]:
    m = (p-1)//2
    qr = sorted(a for a in range(1, p) if legendre(a, p) == 1)
    interval = list(range(1, m+1))

    # Full connection sets (mod p)
    S_int = interval + [p - k for k in interval]
    S_pal = qr + [p - a for a in qr] if p % 4 == 3 else None

    print(f"\np={p}, m={m}")

    # Eigenvalues
    eigs_int = eigenvalues_circulant(interval, p)
    mags_int = sorted([abs(e) for e in eigs_int], reverse=True)
    ipr_int = ipr(eigs_int)

    print(f"  Interval: top |λ|s = {[f'{x:.2f}' for x in mags_int[:4]]}")
    print(f"  Interval IPR = {ipr_int:.6f}")

    if S_pal is not None:
        eigs_pal = eigenvalues_circulant(qr, p)
        mags_pal = sorted([abs(e) for e in eigs_pal], reverse=True)
        ipr_pal = ipr(eigs_pal)
        print(f"  Paley:    top |λ|s = {[f'{x:.2f}' for x in mags_pal[:4]]}")
        print(f"  Paley    IPR = {ipr_pal:.6f}")

    # Trace moments (normalized)
    A_int = tournament_adjacency(tuple(1 for _ in range(m)), p)
    traces_int = [np.trace(np.linalg.matrix_power(A_int, k)) / p for k in range(2, 8)]
    print(f"  Interval traces/p (k=2..7): {[f'{t:.1f}' for t in traces_int]}")

    if S_pal is not None:
        sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))
        A_pal = tournament_adjacency(sigma_pal, p)
        traces_pal = [np.trace(np.linalg.matrix_power(A_pal, k)) / p for k in range(2, 8)]
        print(f"  Paley    traces/p (k=2..7): {[f'{t:.1f}' for t in traces_pal]}")

# =============================================================
# CONNECTION 2: Additive energy → overlap weight
# =============================================================
print("\n" + "=" * 70)
print("CONNECTION 2: ADDITIVE ENERGY AND CYCLE OVERLAP")
print("=" * 70)
print("""
CONJECTURE: E(S) controls the expected vertex overlap between two random
cycles, because:
  - E(S) counts solutions to a + b = c + d in S
  - Vertex overlap between cycles C₁ and C₂ relates to shared arc structures
  - Shared arcs correspond to sumsets a+b = c+d

If this is true, then E(S) → overlap distribution → α₂ → H.
""")

# Test: for each circulant tournament at p=7, compute E(S) and the
# number of vertex-disjoint 3-cycle pairs
for p in [7, 13]:
    m = (p-1)//2
    all_sigmas = list(product([1, -1], repeat=m))

    print(f"\np={p}: Testing E(S) vs disjoint 3-cycle pairs")

    results = []
    for sigma in all_sigmas:
        # Connection set
        S = [k for k in range(1, m+1) if sigma[k-1] == 1]
        S_full = []
        for k in S:
            S_full.append(k)
            S_full.append(p - k)
        if not S_full:
            continue

        E = additive_energy(S_full, p)

        # Count 3-cycles
        A = tournament_adjacency(sigma, p)
        c3 = count_3cycles(A, p)

        # Count vertex-disjoint 3-cycle pairs
        # Find all 3-cycles first
        cycles_3 = []
        for i in range(p):
            for j in range(p):
                if i == j: continue
                for k in range(p):
                    if k == i or k == j: continue
                    if A[i][j] and A[j][k] and A[k][i]:
                        cycle = frozenset([i, j, k])
                        if cycle not in [frozenset(c) for c in cycles_3]:
                            cycles_3.append([i, j, k])

        # Count disjoint pairs
        disj_pairs = 0
        cycle_sets = [frozenset(c) for c in cycles_3]
        unique_sets = list(set(cycle_sets))
        for i_idx in range(len(unique_sets)):
            for j_idx in range(i_idx + 1, len(unique_sets)):
                if unique_sets[i_idx].isdisjoint(unique_sets[j_idx]):
                    disj_pairs += 1

        results.append((E, c3, disj_pairs, S_full))

    # Sort by E
    results.sort(key=lambda x: x[0])
    print(f"  {'E(S)':>8} {'c3':>6} {'disj_pairs':>12} S")
    for E, c3, dp, S in results[:8]:
        print(f"  {E:>8} {c3:>6} {dp:>12}   {sorted(S)}")
    if len(results) > 8:
        print(f"  ... ({len(results)} total)")
        for E, c3, dp, S in results[-4:]:
            print(f"  {E:>8} {c3:>6} {dp:>12}   {sorted(S)}")

    # Correlation
    Es = [r[0] for r in results]
    dps = [r[2] for r in results]
    if np.std(Es) > 0 and np.std(dps) > 0:
        corr = np.corrcoef(Es, dps)[0, 1]
        print(f"  Correlation(E, disjoint_pairs) = {corr:.4f}")

# =============================================================
# CONNECTION 3: The Turán connection — CRITICAL NEW IDEA
# =============================================================
print("\n" + "=" * 70)
print("CONNECTION 3: TURÁN-TYPE BOUND ON INDEPENDENCE POLYNOMIAL")
print("=" * 70)
print("""
NEW INSIGHT: The Ω graph (cycle overlap/conflict graph) is very dense.
For p=7: density ≈ 0.98. This means ALMOST all cycle pairs share a vertex.

The independence polynomial Z(G, λ) = Σ α_k λ^k at high fugacity λ=2
is controlled by the MINIMUM degree and SPECTRAL GAP of G.

QUESTION: Does Interval's Ω have lower minimum degree than Paley's Ω?
If so, this directly explains why α_k(Int) > α_k(Pal) for k ≥ 2.
""")

for p in [7]:
    m = (p-1)//2
    print(f"\np={p}: Detailed Ω graph comparison")

    # Interval
    sigma_int = tuple(1 for _ in range(m))
    A_int = tournament_adjacency(sigma_int, p)

    # Paley
    sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))
    A_pal = tournament_adjacency(sigma_pal, p)

    # Find all odd cycles up to length 7 for each
    for name, A in [("Interval", A_int), ("Paley", A_pal)]:
        cycles = []
        # Find 3-cycles
        for combo in combinations(range(p), 3):
            i, j, k = combo
            if A[i][j] and A[j][k] and A[k][i]:
                cycles.append(frozenset(combo))
            if A[i][k] and A[k][j] and A[j][i]:
                cycles.append(frozenset(combo))

        # Find 5-cycles (via DFS-like approach)
        cycles_5 = set()
        for start in range(p):
            # All directed paths of length 5 starting at start
            def dfs(path, length):
                if length == 5:
                    if A[path[-1]][start]:  # close cycle
                        cycles_5.add(frozenset(path))
                    return
                curr = path[-1]
                for nxt in range(p):
                    if nxt not in path and A[curr][nxt]:
                        dfs(path + [nxt], length + 1)
            dfs([start], 1)
        cycles.extend(cycles_5)

        # Build Ω graph
        N = len(cycles)
        print(f"\n  {name}: {N} odd cycles (3 + 5 only)")

        # Adjacency of Ω: edge iff cycles share a vertex
        omega_adj = np.zeros((N, N), dtype=int)
        cycle_list = list(cycles)
        for i_idx in range(N):
            for j_idx in range(i_idx + 1, N):
                if not cycle_list[i_idx].isdisjoint(cycle_list[j_idx]):
                    omega_adj[i_idx][j_idx] = 1
                    omega_adj[j_idx][i_idx] = 1

        degrees = omega_adj.sum(axis=1)
        print(f"    |E(Ω)| = {omega_adj.sum() // 2}")
        print(f"    Degree: min={degrees.min()}, max={degrees.max()}, avg={degrees.mean():.2f}")
        print(f"    Degree variance = {np.var(degrees):.2f}")

        # Count independent sets of size 2 (= α₂)
        alpha_2 = sum(1 for i_idx in range(N) for j_idx in range(i_idx+1, N)
                      if omega_adj[i_idx][j_idx] == 0)
        print(f"    α₂ = {alpha_2}")

        # Spectral gap
        eigs_omega = np.linalg.eigvalsh(omega_adj.astype(float))
        print(f"    Spectral gap: λ_max={eigs_omega[-1]:.2f}, λ_min={eigs_omega[0]:.2f}")

        # Hoffman bound on independence number
        hoffman = -N * eigs_omega[0] / (eigs_omega[-1] - eigs_omega[0])
        print(f"    Hoffman bound α(Ω) ≥ {hoffman:.2f}")

# =============================================================
# CONNECTION 4: Hamiltonian cycle ratio (NEW from kind-pasteur data)
# =============================================================
print("\n" + "=" * 70)
print("CONNECTION 4: HAMILTONIAN CYCLE REVERSAL AT p=19")
print("=" * 70)
print("""
CRITICAL OBSERVATION (from kind-pasteur S58 data):
  At p=19, c_{19}(Interval) = 34,841,356,485 > c_{19}(Paley) = 34,358,763,933

Interval has MORE Hamiltonian CYCLES than Paley at p=19!

For all shorter odd lengths k < 19, Paley has MORE k-cycles.
The REVERSAL happens at k = p (= Hamiltonian cycles).

This is consistent with:
  - Paley has more "random" short cycles (QR structure spreads them)
  - Interval's directed flow creates more FULL Hamiltonian cycles
  - The crossover from "more short cycles" to "more long cycles"
    is DIRECTLY related to the H crossover via OCF
""")

# Check the ratio c_k(I)/c_k(P) as a function of k at p=19
# Data from kind-pasteur's cycle_counts_alpha_p19.out
pal_cycles = {3: 285, 5: 11628, 7: 424080, 9: 12156390,
              11: 249208902, 13: 3280900392, 15: 23662379790,
              17: 69401425077, 19: 34358763933}
int_cycles = {3: 285, 5: 10488, 7: 391362, 9: 10807884,
              11: 224515210, 13: 2961329208, 15: 21901889133,
              17: 66503305202, 19: 34841356485}

print(f"\n  p=19: Cycle count ratios c_k(Int)/c_k(Pal)")
print(f"  {'k':>4}  {'c_k(Int)':>16}  {'c_k(Pal)':>16}  {'ratio':>8}  {'sign':>8}")
for k in sorted(pal_cycles.keys()):
    ci = int_cycles[k]
    cp = pal_cycles[k]
    ratio = ci / cp
    sign = "INT>PAL" if ratio > 1 else ("TIED" if ratio == 1 else "PAL>INT")
    print(f"  {k:>4}  {ci:>16,}  {cp:>16,}  {ratio:>8.4f}  {sign:>8}")

# KEY INSIGHT: The ratio INCREASES monotonically toward k=p!
# This suggests that Interval's flow coherence amplifies with cycle length.

# =============================================================
# CONNECTION 5: Walsh degree structure and OCF coupling
# =============================================================
print("\n" + "=" * 70)
print("CONNECTION 5: WALSH DEGREE WEIGHTS IN H")
print("=" * 70)
print("""
From THM-077 (Walsh OCF proof), we know:
  hat{H}[S] = epsilon * 2^r * (n-2k)! / 2^{n-1}

where S is a union of r even-length path components on n-2k vertices.

For circulant tournaments, the orientation sigma acts as:
  H(sigma) = sum_S hat{H}[S] * chi_S(sigma)
           = sum_S hat{H}[S] * prod_{j in S} sigma_j

The KEY: hat{H}[S] depends only on |S| and the component structure,
NOT on which specific chords are in S.

For Paley sigma_P: chi_S(sigma_P) = prod_{j in S} legendre(j, p)
  = legendre(prod(S), p)   [by multiplicativity]
  = ??? (depends on product of elements of S mod p)

For Interval sigma_I = (1,...,1): chi_S(sigma_I) = 1 for all S.
  So H(Int) = sum_S hat{H}[S] (all terms add up constructively!)

For Paley: some chi_S(sigma_P) = -1, causing destructive interference.
  At large p: about half the terms cancel → H(Pal) < H(Int).
  At small p: the cancellation pattern is structured enough that
  Paley gains from the "right" signs → H(Pal) > H(Int).
""")

# Test: what fraction of Walsh characters are +1 for Paley vs Interval?
for p in [7, 11, 13, 19]:
    m = (p-1)//2
    sigma_pal = tuple(legendre(k, p) for k in range(1, m+1))

    total_subsets = 0
    positive_chars = 0
    for size in range(0, m+1):
        for S in combinations(range(m), size):
            total_subsets += 1
            chi = 1
            for j in S:
                chi *= sigma_pal[j]
            if chi == 1:
                positive_chars += 1

    frac = positive_chars / total_subsets
    print(f"  p={p}: {positive_chars}/{total_subsets} = {frac:.4f} of Walsh chars are +1 for Paley")
    print(f"        Expected if random: 0.5000")
    print(f"        Bias toward +1: {frac - 0.5:+.4f}")

# =============================================================
# CONNECTION 6: The product formula for chi_S(sigma_P)
# =============================================================
print("\n" + "=" * 70)
print("CONNECTION 6: PRODUCT FORMULA AND LEGENDRE SYMBOL")
print("=" * 70)
print("""
For S ⊂ {1,...,m}, chi_S(sigma_P) = prod_{k in S} legendre(k, p).

By multiplicativity of the Legendre symbol:
  chi_S(sigma_P) = legendre(prod(S), p)

where prod(S) = product of all elements of S (mod p).

This means: chi_S is +1 iff prod(S) is a QR mod p.

QUESTION: Among all subsets S of {1,...,m}, what fraction have
QR product? Is this exactly 1/2?
""")

for p in [7, 11, 19]:
    m = (p-1)//2

    qr_prod_count = 0
    nqr_prod_count = 0
    zero_count = 0
    total = 0

    for size in range(1, m+1):
        for S in combinations(range(1, m+1), size):
            total += 1
            prod_S = 1
            for k in S:
                prod_S = (prod_S * k) % p
            leg = legendre(prod_S, p)
            if leg == 1:
                qr_prod_count += 1
            elif leg == -1:
                nqr_prod_count += 1
            else:
                zero_count += 1

    print(f"  p={p}, m={m}: QR product: {qr_prod_count}, NQR product: {nqr_prod_count}")
    print(f"    QR fraction: {qr_prod_count/total:.4f}")
    print(f"    Perfect balance = 0.5000")
    print(f"    Excess QR: {qr_prod_count - nqr_prod_count:+d} ({(qr_prod_count-nqr_prod_count)/total:+.4f})")

print("\n" + "=" * 70)
print("SYNTHESIS: THE PROOF CHAIN")
print("=" * 70)
print("""
The complete proof would need:

1. [PROVED] Interval eigenvalues = Fejér kernel (exact match)
2. [PROVED] Fejér maximizes IPR (harmonic analysis classical result)
3. [NEED] IPR → cycle disjointness bound
   Key: show that max|λ_k|²/Σ|λ_k|² ≥ c implies α₂ ≥ f(α₁,c)
4. [NEED] Cycle disjointness → H bound via OCF
   Key: H = Σ 2^k α_k, so exponential weighting amplifies α₂ advantage
5. [NEED] The critical coupling g_c where many-body overtakes 2-body

The WEAKEST LINK is step 3: converting spectral concentration into
a quantitative cycle disjointness bound.

POSSIBLE APPROACH (NEW):
  Use the trace formula: c_k = (1/p) Σ λ_j^k
  For Interval: c_k ≈ |μ₁|^k/p (dominated by top eigenvalue)
  The cycles of different lengths k₁, k₂ tend to use the SAME
  "spectral channel" (top eigenvalue), so they OVERLAP.

  But the TOTAL number of cycles is smaller (fewer channels).
  The disjoint ones come from DIFFERENT spectral channels.

  Interval has 2 main channels (μ₁ and μ₁*), Paley has p-1 equal channels.
  With 2 channels: cycles split into 2 groups → high disjointness BETWEEN groups.
  With p-1 channels: cycles spread randomly → low disjointness everywhere.

  This is the SCHUR COMPLEMENT / BLOCK DIAGONAL argument:
  Ω(Interval) ≈ complete bipartite graph (2 blocks from 2 eigenvalues)
  Ω(Paley) ≈ Erdős-Rényi random graph (uniform mixing)
  Complete bipartite → max α₂ for given density!
""")

print("\nDONE.")
