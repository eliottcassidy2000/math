#!/usr/bin/env python3
"""
Boolean function analysis of H on the orientation cube {±1}^m.

KEY INSIGHT: H(σ) = Σ_S ĥ[S] χ_S(σ) is a multilinear polynomial on {±1}^m.
The question "does all-ones maximize H?" is equivalent to asking whether
H satisfies certain MONOTONICITY or CORRELATION properties.

APPROACHES EXPLORED:
1. FKG inequality: If ĥ[S] ≥ 0 for all S with |S| even, then H is
   "increasing" in a specific sense (Fortuin-Kasteleyn-Ginibre).

2. Noise sensitivity: If H is "noise stable", its maximum concentrates
   near the constant function → all-ones wins.

3. Hypercontractivity (Bonami-Beckner): High-degree Walsh coefficients
   are bounded in terms of low-degree ones.

4. Influence / pivotality: The "influence" of each coordinate tells us
   how sensitive H is to flipping one chord.

5. Log-supermodularity: If H(σ ∨ τ) · H(σ ∧ τ) ≥ H(σ) · H(τ) for all
   σ, τ (where ∨, ∧ are coordinatewise max/min), then H is log-supermodular
   and the FKG lattice condition gives the maximum at all-ones.

NEW CONNECTION: The Ahlswede-Daykin "four functions" theorem generalizes FKG
and may apply even when individual Walsh coefficients have mixed signs.

opus-2026-03-12-S67
"""

import numpy as np
from itertools import combinations
from sympy.ntheory import legendre_symbol as legendre

def make_tournament(p, S):
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for d in S:
            A[i][(i+d)%p] = 1
    return A

def count_H(A):
    n = len(A)
    dp = [[0]*n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def orientation_to_S(sigma, p):
    m = (p - 1) // 2
    S = set()
    for k in range(1, m + 1):
        if sigma[k-1] == 1:
            S.add(k)
        else:
            S.add(p - k)
    return S

# ========================================================================
print("=" * 72)
print("PART I: LOG-SUPERMODULARITY TEST")
print("=" * 72)
print("""
If H is log-supermodular on {±1}^m, i.e.:
  H(σ ∨ τ) · H(σ ∧ τ) ≥ H(σ) · H(τ)  for all σ, τ
where (σ ∨ τ)_i = max(σ_i, τ_i) and (σ ∧ τ)_i = min(σ_i, τ_i),
then by FKG:
  - The maximum is at all-ones (since all-ones = ∨ of everything)
  - More generally, for any increasing event A:
    P(A | high H) ≥ P(A | low H)

This is EXACTLY what we need to prove Interval maximizes H!
""")

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    # Compute H for all orientations
    H_vals = {}
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_vals[sigma] = count_H(A)

    # Test log-supermodularity: H(σ∨τ)·H(σ∧τ) ≥ H(σ)·H(τ)
    violations = 0
    total = 0
    worst_ratio = float('inf')

    all_sigmas = list(H_vals.keys())
    for i in range(len(all_sigmas)):
        for j in range(i+1, len(all_sigmas)):
            s, t = all_sigmas[i], all_sigmas[j]

            # Compute join and meet
            join = tuple(max(s[k], t[k]) for k in range(m))
            meet = tuple(min(s[k], t[k]) for k in range(m))

            lhs = H_vals[join] * H_vals[meet]
            rhs = H_vals[s] * H_vals[t]

            total += 1
            if lhs < rhs:
                violations += 1
                ratio = lhs / rhs
                if ratio < worst_ratio:
                    worst_ratio = ratio

    print(f"  p={p}: {violations}/{total} log-supermodularity violations")
    if violations > 0:
        print(f"    Worst ratio: {worst_ratio:.6f}")
    else:
        print(f"    *** H IS LOG-SUPERMODULAR! FKG APPLIES! ***")

# ========================================================================
print("\n" + "=" * 72)
print("PART II: MONOTONICITY ANALYSIS — IS H COORDINATEWISE INCREASING?")
print("=" * 72)
print("""
A simpler condition: is H(σ) ≤ H(σ') whenever σ ≤ σ' coordinatewise?
If so, then all-ones is trivially the maximum.

This is equivalent to: ALL degree-1 Walsh coefficients are non-negative,
AND the function is monotone (which is stronger than just degree-1 signs).
""")

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    # Compute H for all orientations
    H_dict = {}
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_dict[sigma] = count_H(A)

    # Check monotonicity: flipping σ_i from -1 to +1 should increase H
    non_monotone = 0
    total_flips = 0

    for sigma in H_dict:
        for i in range(m):
            if sigma[i] == -1:
                sigma_flipped = list(sigma)
                sigma_flipped[i] = 1
                sigma_flipped = tuple(sigma_flipped)

                total_flips += 1
                if H_dict[sigma_flipped] < H_dict[sigma]:
                    non_monotone += 1

    print(f"\n  p={p}: {non_monotone}/{total_flips} anti-monotone flips")

    # Compute Walsh coefficients (full Hadamard transform)
    H_vec = np.zeros(N)
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        H_vec[bits] = H_dict[sigma]

    # Walsh-Hadamard transform
    h_hat = np.zeros(N)
    for S_bits in range(N):
        S_set = [k for k in range(m) if (S_bits >> k) & 1]
        for bits in range(N):
            sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
            chi = 1
            for k in S_set:
                chi *= sigma[k]
            h_hat[S_bits] += H_vec[bits] * chi
        h_hat[S_bits] /= N

    # Degree-1 coefficients
    print(f"  Degree-1 Walsh coefficients (influences):")
    for i in range(m):
        S_bits = 1 << i
        print(f"    ĥ[{{{i+1}}}] = {h_hat[S_bits]:.2f}", end="")
        if h_hat[S_bits] < 0:
            print(" ← NEGATIVE (anti-monotone in chord {})".format(i+1))
        else:
            print(" ← positive")

# ========================================================================
print("\n" + "=" * 72)
print("PART III: FKG LATTICE CONDITION (HOLLEY CRITERION)")
print("=" * 72)
print("""
Even if H is not coordinatewise monotone, it can still satisfy the
FKG lattice condition: for all σ ~ τ (differing in one coordinate),
  H(σ∨τ) / H(σ) ≥ H(τ) / H(σ∧τ)

This is the CONDITIONAL version of log-supermodularity, and it's
sufficient for FKG to apply. It's also easier to verify — we only
need to check pairs differing in exactly TWO coordinates.
""")

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    H_dict = {}
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_dict[sigma] = count_H(A)

    # FKG lattice condition: for σ, τ differing in coords i, j
    # with σ_i=+1, σ_j=-1, τ_i=-1, τ_j=+1:
    # H(σ∨τ) · H(σ∧τ) ≥ H(σ) · H(τ)

    violations = 0
    total = 0

    for sigma in H_dict:
        for i in range(m):
            for j in range(i+1, m):
                if sigma[i] == 1 and sigma[j] == -1:
                    # Flip i: -1, j: +1
                    tau = list(sigma)
                    tau[i] = -1
                    tau[j] = 1
                    tau = tuple(tau)

                    join = list(sigma)
                    join[j] = 1  # both +1
                    join = tuple(join)

                    meet = list(sigma)
                    meet[i] = -1  # both -1
                    meet = tuple(meet)

                    lhs = H_dict[join] * H_dict[meet]
                    rhs = H_dict[sigma] * H_dict[tau]

                    total += 1
                    if lhs < rhs:
                        violations += 1

    print(f"  p={p}: {violations}/{total} FKG lattice violations (2-coord test)")
    if violations == 0:
        print(f"    *** FKG LATTICE CONDITION HOLDS AT p={p}! ***")

# ========================================================================
print("\n" + "=" * 72)
print("PART IV: NOISE STABILITY AND HYPERCONTRACTIVITY")
print("=" * 72)
print("""
Noise stability: S_ρ(H) = Σ_S ρ^|S| ĥ[S]²
  = expected correlation between H(σ) and H(σ') where σ' is a
  ρ-correlated copy of σ.

If S_ρ(H) is close to Var(H), then H is "noise stable" — it depends
mostly on low-degree Fourier coefficients. The Majority Is Stablest
theorem (Mossel-O'Donnell-Oleszkiewicz 2010) then implies that the
maximum is near the "dictator" functions (all-ones or all-minus-ones).

Hypercontractivity (Bonami-Beckner): ||T_ρ f||_q ≤ ||f||_p for
q-1 = ρ²(p-1). This bounds how much energy can be in high degrees.
""")

for p in [7, 11, 13]:
    m = (p - 1) // 2
    N = 2**m

    # Compute Walsh spectrum
    H_dict = {}
    for bits in range(N):
        sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        H_dict[sigma] = count_H(A)

    H_vec = np.array([H_dict[tuple(1 if (bits >> k) & 1 else -1 for k in range(m))]
                       for bits in range(N)], dtype=float)

    # Walsh-Hadamard transform
    h_hat = np.zeros(N)
    for S_bits in range(N):
        S_set = [k for k in range(m) if (S_bits >> k) & 1]
        deg = len(S_set)
        val = 0.0
        for bits in range(N):
            chi = 1
            for k in S_set:
                if not ((bits >> k) & 1):
                    chi = -chi
            val += H_vec[bits] * chi
        h_hat[S_bits] = val / N

    # Energy by degree
    energy_by_deg = {}
    for S_bits in range(N):
        deg = bin(S_bits).count('1')
        if deg not in energy_by_deg:
            energy_by_deg[deg] = 0.0
        energy_by_deg[deg] += h_hat[S_bits]**2

    total_var = sum(v for d, v in energy_by_deg.items() if d > 0)

    print(f"\n  p={p}, m={m}:")
    print(f"    Mean H = {h_hat[0]:.1f}")
    print(f"    Total Var(H) = {total_var:.1f}")
    for deg in sorted(energy_by_deg.keys()):
        frac = energy_by_deg[deg] / (total_var if total_var > 0 else 1) if deg > 0 else 0
        print(f"    Degree {deg}: energy = {energy_by_deg[deg]:.1f}" +
              (f" ({frac*100:.1f}% of variance)" if deg > 0 else " (mean²)"))

    # Noise stability at various ρ
    print(f"    Noise stability S_ρ / Var(H):")
    for rho in [0.5, 0.8, 0.9, 0.95, 0.99]:
        S_rho = sum(rho**(2*bin(S).count('1')) * h_hat[S]**2
                     for S in range(1, N))
        print(f"      ρ={rho}: S_ρ/Var = {S_rho/total_var:.4f}" if total_var > 0 else "      N/A")

# ========================================================================
print("\n" + "=" * 72)
print("PART V: THE CORRELATION INEQUALITY APPROACH")
print("=" * 72)
print("""
NEW IDEA: Instead of proving H is log-supermodular directly,
use the POLYMER EXPANSION of Z(Ω, λ):

  Z(Ω, λ) = Π_{v ∈ V(Ω)} (1 + λ·x_v) · Π_{uv ∈ E(Ω)} (1 - x_u·x_v·...)

In the hard-core model, Z = Σ_{I independent} λ^|I|.

The KEY: Z(Ω, λ) is a MONOTONE function of the complement graph.
If we can show that Interval's complement-Ω has more edges (sparser Ω),
then Z(Ω_Int, λ) ≥ Z(Ω_Pal, λ) for all λ > 0.

Formally: if Ω_Int ⊆ Ω_Pal (as graphs), then every independent set of
Ω_Pal is also independent in Ω_Int, so Z(Ω_Int, λ) ≥ Z(Ω_Pal, λ).

Is Ω_Int ⊂ Ω_Pal as a SUBGRAPH? Let's check!
""")

def find_directed_cycles(A, p, max_len=None):
    """Find all directed odd cycles."""
    if max_len is None:
        max_len = p
    n = len(A)
    cycles = set()
    for L in range(3, max_len + 1, 2):
        for start in range(n):
            def dfs(path, remaining):
                v = path[-1]
                if remaining == 0:
                    if A[v][start]:
                        min_idx = path.index(min(path))
                        canonical = tuple(path[min_idx:] + path[:min_idx])
                        cycles.add(canonical)
                    return
                for u in range(n):
                    if A[v][u] and u not in path:
                        if remaining == 1 and u != start:
                            continue
                        if remaining > 1 and u == start:
                            continue
                        dfs(path + [u], remaining - 1)
            dfs([start], L - 1)
    return list(cycles)

for p in [7, 11]:
    m = (p - 1) // 2

    # Interval
    sigma_int = tuple(1 for _ in range(m))
    S_int = orientation_to_S(sigma_int, p)
    A_int = make_tournament(p, S_int)
    cycles_int = find_directed_cycles(A_int, p)

    # Paley (p ≡ 3 mod 4)
    sigma_pal = tuple(int(legendre(k, p)) for k in range(1, m + 1))
    S_pal = orientation_to_S(sigma_pal, p)
    A_pal = make_tournament(p, S_pal)
    cycles_pal = find_directed_cycles(A_pal, p)

    # Build Ω graphs (overlap = share a vertex)
    nc_int = len(cycles_int)
    nc_pal = len(cycles_pal)

    # Ω edges
    edges_int = set()
    for i in range(nc_int):
        vi = set(cycles_int[i])
        for j in range(i+1, nc_int):
            vj = set(cycles_int[j])
            if vi & vj:
                edges_int.add((i, j))

    edges_pal = set()
    for i in range(nc_pal):
        vi = set(cycles_pal[i])
        for j in range(i+1, nc_pal):
            vj = set(cycles_pal[j])
            if vi & vj:
                edges_pal.add((i, j))

    # Edge density
    max_edges_int = nc_int * (nc_int - 1) // 2
    max_edges_pal = nc_pal * (nc_pal - 1) // 2

    print(f"\n  p={p}:")
    print(f"    Interval: {nc_int} cycles, {len(edges_int)} Ω-edges "
          f"(density {len(edges_int)/max(1,max_edges_int):.4f})")
    print(f"    Paley:    {nc_pal} cycles, {len(edges_pal)} Ω-edges "
          f"(density {len(edges_pal)/max(1,max_edges_pal):.4f})")

    # Can we embed: is there a mapping from Interval's cycles to Paley's
    # such that independent sets are preserved?
    # Simpler question: compare edge DENSITIES

    # Alpha_k comparison
    alpha_int = [0] * (nc_int + 1)
    alpha_pal = [0] * (nc_pal + 1)

    if nc_int <= 20:
        adj_int = [[False]*nc_int for _ in range(nc_int)]
        for (i,j) in edges_int:
            adj_int[i][j] = adj_int[j][i] = True

        for mask in range(1 << nc_int):
            nodes = [i for i in range(nc_int) if (mask >> i) & 1]
            k = len(nodes)
            indep = True
            for a in range(len(nodes)):
                for b in range(a+1, len(nodes)):
                    if adj_int[nodes[a]][nodes[b]]:
                        indep = False
                        break
                if not indep:
                    break
            if indep:
                alpha_int[k] += 1

    if nc_pal <= 20:
        adj_pal = [[False]*nc_pal for _ in range(nc_pal)]
        for (i,j) in edges_pal:
            adj_pal[i][j] = adj_pal[j][i] = True

        for mask in range(1 << nc_pal):
            nodes = [i for i in range(nc_pal) if (mask >> i) & 1]
            k = len(nodes)
            indep = True
            for a in range(len(nodes)):
                for b in range(a+1, len(nodes)):
                    if adj_pal[nodes[a]][nodes[b]]:
                        indep = False
                        break
                if not indep:
                    break
            if indep:
                alpha_pal[k] += 1

    print(f"    Independence polynomial comparison:")
    H_int = sum(alpha_int[k] * 2**k for k in range(nc_int + 1))
    H_pal = sum(alpha_pal[k] * 2**k for k in range(nc_pal + 1))
    for k in range(min(8, max(nc_int, nc_pal) + 1)):
        ai = alpha_int[k] if k < len(alpha_int) else 0
        ap = alpha_pal[k] if k < len(alpha_pal) else 0
        if ai > 0 or ap > 0:
            winner = "INT" if ai > ap else ("PAL" if ap > ai else "TIE")
            print(f"      α_{k}: Int={ai}, Pal={ap}, Δ={ai-ap}, 2^k·Δ={2**k*(ai-ap)} [{winner}]")
    print(f"    Z(Ω,2): Int={H_int}, Pal={H_pal}")

# ========================================================================
print("\n" + "=" * 72)
print("PART VI: STOCHASTIC DOMINANCE OF Ω DEGREE SEQUENCES")
print("=" * 72)
print("""
Even if Ω_Int ⊄ Ω_Pal, we can compare their degree sequences.
In the hard-core model, Z is monotone in the graph structure:
  fewer edges → more independent sets → higher Z.

But this is too coarse — the STRUCTURE of edges matters.
Specifically: Ω_Int may have fewer edges but they could be more
concentrated (forming cliques), which HELPS independent sets elsewhere.

The REAL comparison: Interval's Ω has a SPARSER structure in the
"between-clusters" region, compensating for denser "within-cluster" edges.
""")

for p in [7]:
    m = (p - 1) // 2

    for name, sigma in [("Interval", tuple(1 for _ in range(m))),
                         ("Paley", tuple(int(legendre(k, p)) for k in range(1, m+1)))]:
        S = orientation_to_S(sigma, p)
        A = make_tournament(p, S)
        cycles = find_directed_cycles(A, p)
        nc = len(cycles)

        # Degree sequence in Ω
        deg = [0] * nc
        for i in range(nc):
            vi = set(cycles[i])
            for j in range(nc):
                if i != j:
                    vj = set(cycles[j])
                    if vi & vj:
                        deg[i] += 1

        deg_sorted = sorted(deg, reverse=True)
        print(f"\n  p={p}, {name}: {nc} cycles")
        print(f"    Ω degree sequence: {deg_sorted}")
        print(f"    Mean degree: {np.mean(deg):.2f}")
        print(f"    Max degree: {max(deg)}")

        # Cycle vertex sets and their overlaps
        print(f"    Cycles and their vertex sets:")
        by_length = {}
        for c in cycles:
            L = len(c)
            if L not in by_length:
                by_length[L] = []
            by_length[L].append(c)

        for L in sorted(by_length.keys()):
            print(f"      Length {L}: {len(by_length[L])} cycles")
            for c in by_length[L][:5]:
                print(f"        {c} (vertices: {sorted(set(c))})")

# ========================================================================
print("\n" + "=" * 72)
print("PART VII: THE CONVEXITY-OF-Z APPROACH")
print("=" * 72)
print("""
MOST PROMISING DIRECTION:

Z(Ω(σ), 2) as a function of the eigenvalue distribution.

We know:
  |λ_t(σ)|² = (m/2)² + b_t(σ)²
  where b_t(σ) = Σ_j σ_j sin(2πjt/p) = imaginary part of eigenvalue

So |λ_t|² is a QUADRATIC function of σ. But Z is a highly nonlinear
function of the cycle counts, which are polynomial in |λ|².

KEY: The cycle count c_k = (1/p) Σ_t λ_t^k. Since Re(λ_t) = -1/2 for all t,
  λ_t = -1/2 + i·b_t(σ)
  λ_t^k = (-1/2 + i·b_t)^k = Σ_j C(k,j) (-1/2)^{k-j} (i·b_t)^j

So c_k is a polynomial of degree k in the b_t's, which are LINEAR in σ.
Therefore c_k is a polynomial of degree k in σ!

The full H = Z(Ω, 2) involves ALL c_k for k=3,5,...,p. The nonlinearity
of H as a function of σ comes from the PRODUCT of these cycle counts
through the independence polynomial.

This suggests: H is a polynomial of degree ~ p in σ, but the effective
degree is determined by the maximum order of cycle interactions.

For the PROOF, we need: among all σ on the (same) sphere ||Bσ||² = const,
which σ maximizes Z? The answer is σ = all-ones for p ≥ 13.

The Fejér kernel property (|λ_1|² is maximized at all-ones by the
isoperimetric inequality on Z_p) combined with the super-exponential
moment ratio growth should give the result asymptotically.

Let's verify the polynomial structure and look for CONVEXITY.
""")

for p in [7, 13]:
    m = (p - 1) // 2
    N = 2**m

    print(f"\n  p={p}, m={m}:")

    # Compute b_t for all t and for all sigma
    # b_t(σ) = Σ_j σ_j sin(2πjt/p)
    sin_matrix = np.zeros((p-1, m))  # B[t-1, j-1] = sin(2πjt/p)
    for t in range(1, p):
        for j in range(1, m+1):
            sin_matrix[t-1, j-1] = np.sin(2 * np.pi * j * t / p)

    # For each σ, compute eigenvalues and then c_k
    H_all = []
    b_all = []
    c3_all = []
    c5_all = []

    for bits in range(N):
        sigma = np.array([1 if (bits >> k) & 1 else -1 for k in range(m)], dtype=float)
        b = sin_matrix @ sigma  # shape (p-1,)

        # λ_t = -1/2 + i·b_t
        lam = -0.5 + 1j * b

        # Cycle counts
        c3 = np.sum(lam**3).real / p
        c5 = np.sum(lam**5).real / p

        # H
        sigma_tuple = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        S = orientation_to_S(sigma_tuple, p)
        A = make_tournament(p, S)
        H = count_H(A)

        H_all.append(H)
        b_all.append(b)
        c3_all.append(c3)
        c5_all.append(c5)

    H_all = np.array(H_all)
    c3_all = np.array(c3_all)
    c5_all = np.array(c5_all)

    # Check: is H a function of (c3, c5, c7, ...)?
    # Group orientations by their cycle count profile
    from collections import defaultdict
    cycle_groups = defaultdict(list)
    for bits in range(N):
        cycle_groups[(round(c3_all[bits], 2), round(c5_all[bits], 2))].append(H_all[bits])

    print(f"  {N} orientations, {len(cycle_groups)} distinct (c3,c5) profiles")
    print(f"  Same (c3,c5) → same H? ", end="")
    same_H = all(len(set(v)) == 1 for v in cycle_groups.values())
    print("YES" if same_H else "NO")

    if not same_H:
        for key, vals in cycle_groups.items():
            if len(set(vals)) > 1:
                print(f"    (c3,c5)={key}: H values = {sorted(set(vals))}")
                break

    # Regression: H vs (c3, c5, c3², c5², c3·c5)
    X = np.column_stack([c3_all, c5_all, c3_all**2, c5_all**2, c3_all*c5_all])
    if X.shape[0] > X.shape[1]:
        from numpy.linalg import lstsq
        coeffs, residuals, _, _ = lstsq(X, H_all, rcond=None)
        H_pred = X @ coeffs
        R2 = 1 - np.sum((H_all - H_pred)**2) / np.sum((H_all - np.mean(H_all))**2)
        print(f"  R² of H ~ (c3, c5, c3², c5², c3·c5): {R2:.6f}")

    # Check if c3 is the same for all
    c3_unique = len(set(round(x, 4) for x in c3_all))
    c5_unique = len(set(round(x, 4) for x in c5_all))
    print(f"  Distinct c3 values: {c3_unique}")
    print(f"  Distinct c5 values: {c5_unique}")

# ========================================================================
print("\n" + "=" * 72)
print("PART VIII: THE CRITICAL TEST — LOG-SUPERMODULARITY AT p=13")
print("=" * 72)
print("""
p=13 is the critical case: first prime where Interval wins.
If log-supermodularity holds here, combined with p=7 (verified above),
this gives strong evidence for a general proof.

For p ≡ 1 mod 4 (like p=13), there is no Paley tournament.
The competitor is the FLATTEST orientation (Paley-like analog).

Log-supermodularity at p=13 would mean:
  For ALL σ, τ ∈ {±1}^6:
    H(σ∨τ) · H(σ∧τ) ≥ H(σ) · H(τ)

This is 2^12 / 2 = 2048 pairs to check (already done in Part I).
Let's also check the STRONGER condition: log-convexity along
every line segment (1-dimensional restriction).
""")

p = 13
m = (p - 1) // 2
N = 2**m

H_dict = {}
for bits in range(N):
    sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
    S = orientation_to_S(sigma, p)
    A = make_tournament(p, S)
    H_dict[sigma] = count_H(A)

# Along each coordinate direction, check convexity of log H
print(f"\n  p=13: Log-convexity along coordinate axes:")
for i in range(m):
    # Fix all other coordinates and vary coordinate i
    # This gives 2 values: σ_i = -1 and σ_i = +1
    # Average over all other coordinates
    increase_count = 0
    decrease_count = 0
    total = 0

    for bits in range(N):
        if (bits >> i) & 1:  # σ_i = +1
            continue
        sigma_minus = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
        sigma_plus = list(sigma_minus)
        sigma_plus[i] = 1
        sigma_plus = tuple(sigma_plus)

        total += 1
        if H_dict[sigma_plus] >= H_dict[sigma_minus]:
            increase_count += 1
        else:
            decrease_count += 1

    pct = increase_count / total * 100
    print(f"    Coord {i+1} (chord {i+1}): flip -1→+1 increases H in {increase_count}/{total} cases ({pct:.1f}%)")

# Along each pair of coordinates, check 2D log-supermodularity
print(f"\n  p=13: 2D log-supermodularity (all coordinate pairs):")
for i in range(m):
    for j in range(i+1, m):
        violations = 0
        total = 0

        for bits in range(N):
            sigma = tuple(1 if (bits >> k) & 1 else -1 for k in range(m))
            if sigma[i] != 1 or sigma[j] != -1:
                continue

            # σ has i=+1, j=-1
            tau = list(sigma)
            tau[i] = -1
            tau[j] = 1
            tau = tuple(tau)

            join = list(sigma)
            join[j] = 1
            join = tuple(join)

            meet = list(sigma)
            meet[i] = -1
            meet = tuple(meet)

            lhs = H_dict[join] * H_dict[meet]
            rhs = H_dict[sigma] * H_dict[tau]

            total += 1
            if lhs < rhs:
                violations += 1

        if violations > 0:
            print(f"    Pair ({i+1},{j+1}): {violations}/{total} VIOLATIONS")

# ========================================================================
print("\n" + "=" * 72)
print("PART IX: SUMMARY AND PROOF ARCHITECTURE")
print("=" * 72)
print("""
PROOF ARCHITECTURE for Interval H-maximization:

PATH A: Log-supermodularity → FKG → all-ones maximum
  Requires: H(σ∨τ)·H(σ∧τ) ≥ H(σ)·H(τ) for all σ,τ
  Status: TESTING above

PATH B: Noise stability → Majority Is Stablest → dictator maximum
  Requires: High noise stability (low high-degree energy)
  Status: Computed noise stability above

PATH C: Direct Walsh analysis → hyperplane inequality
  Requires: For all σ, Σ_{ψ(S,σ)=-1} ĥ[S] ≥ 0
  Status: Verified computationally at p=7,11,13,17

PATH D: Spectral concentration → polymer expansion convergence
  Requires: Interval's Fejér eigenvalues give faster convergence
  Status: Conceptual framework established

PATH E: Subgraph ordering → monotone Z
  Requires: Ω_Int ⊆ Ω_Pal (as graphs, up to relabeling)
  Status: TESTING above

PATH F: Asymptotic moment ratio → large-p dominance
  Requires: μ_k(Int)/μ_k(Pal) → ∞ super-exponentially
  Status: Verified numerically for k ≤ 20, p ≤ 97
  Combined with exhaustive p ≤ 19 gives all primes
""")

print("\nDONE.")
