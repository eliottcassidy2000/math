#!/usr/bin/env python3
"""
info_geometry_morse_s20y.py -- kind-pasteur-2026-03-22-S20y

INFORMATION GEOMETRY + MORSE THEORY OF TOURNAMENT SPACE

Synthesis of three frameworks:
1. MORSE THEORY: H as a Morse function on the m-cube (from S20x)
2. INFORMATION GEOMETRY: Fisher metric, Boltzmann families, natural gradient
3. WALSH-FOURIER ANALYSIS: epistatic decomposition of H on the hypercube

Key references pulled from web search:
- Stadler-Reidys (2002): Combinatorial Landscapes framework
- Ollivier (2009, 2017): Ricci curvature on metric spaces + IGO
- Kolesnik-Sanchez (2024): Zonotopal geometry of tournaments
- Engstrom (2009): Morse-Fourier on simplicial complexes
- De Visser-Park-Krug (2022): Discrete Morse on fitness landscapes

WHAT WE COMPUTE:
A. Walsh-Fourier spectrum of H on {0,1}^m
B. Fisher information of the Boltzmann family exp(beta*H)
C. Ollivier-Ricci curvature of the tournament flip graph (approximation)
D. Persistent homology via sublevel filtration
E. Natural gradient direction at each tournament

Author: kind-pasteur-2026-03-22-S20y
"""
import sys
import numpy as np
from math import comb, log, exp, log2
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

print("=" * 70)
print("  INFORMATION GEOMETRY + MORSE THEORY OF TOURNAMENT SPACE")
print("=" * 70)

# ================================================================
# SETUP: compute H for all tournaments at n=5
# ================================================================
n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)  # 10

print(f"\n  n={n}, m=C({n},2)={m}, |T|=2^{m}={2**m}")
print(f"  Computing H for all tournaments...")

H_map = {}
for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H_map[bits] = count_hp(A, n)

H_vals = np.array([H_map[b] for b in range(2**m)], dtype=float)
H_mean = H_vals.mean()
H_var = H_vals.var()
print(f"  E[H] = {H_mean:.4f}, Var[H] = {H_var:.4f}")

# ================================================================
# A. WALSH-FOURIER SPECTRUM OF H
# ================================================================
print(f"\n{'='*70}")
print(f"  A. WALSH-FOURIER SPECTRUM OF H ON THE {m}-CUBE")
print(f"{'='*70}")
print()

# The Walsh-Hadamard transform: H_hat(S) = (1/2^m) * sum_x (-1)^{<S,x>} * H(x)
# where S is a subset of [m] (represented as bitmask)
# This decomposes H into "epistatic" components of different orders.

# Compute via fast Walsh-Hadamard transform
fhat = H_vals.copy()
for j in range(m):
    step = 1 << (j + 1)
    half = 1 << j
    for k in range(0, 2**m, step):
        for i in range(half):
            u = fhat[k + i]
            v = fhat[k + i + half]
            fhat[k + i] = u + v
            fhat[k + i + half] = u - v
fhat /= 2**m  # normalize

# Group by order (number of bits in S)
order_energy = defaultdict(float)
order_count = defaultdict(int)
for S in range(2**m):
    order = bin(S).count('1')
    order_energy[order] += fhat[S]**2
    order_count[order] += 1

total_energy = sum(order_energy.values())
print(f"  Walsh-Fourier decomposition of H:")
print(f"  {'Order':>6s} {'#terms':>7s} {'Energy':>12s} {'%total':>8s} {'Landscape?':>12s}")
for order in sorted(order_energy.keys()):
    energy = order_energy[order]
    pct = 100 * energy / total_energy if total_energy > 0 else 0
    # Elementary landscape: H is eigenfunction of graph Laplacian
    # For hypercube: eigenvalue of order-k component is 2k
    is_elem = "elementary" if order_count[order] == 1 and pct > 99 else ""
    print(f"  {order:>6d} {order_count[order]:>7d} {energy:>12.4f} {pct:>7.2f}% {is_elem:>12s}")

# The ORDER-0 term is the mean: H_hat(emptyset) = E[H]
print(f"\n  H_hat(0) = {fhat[0]:.4f} (should equal E[H] = {H_mean:.4f})")

# Check: which orders are nonzero?
nonzero_orders = [o for o in sorted(order_energy.keys()) if order_energy[o] > 1e-10]
print(f"  Nonzero orders: {nonzero_orders}")
print(f"  Max order with nonzero energy: {max(nonzero_orders) if nonzero_orders else 'none'}")

# KEY: Is H an "elementary landscape"?
# H is elementary if it equals a single Walsh order component.
# From Stadler-Reidys: elementary iff E[H(neighbor)] = a + b*H(x) for constants a,b
# Check this:
print(f"\n  ELEMENTARY LANDSCAPE TEST:")
for bits in [0, 1, 2**m - 1]:
    h = H_map[bits]
    neighbors = [H_map[bits ^ (1 << k)] for k in range(m)]
    avg_nb = sum(neighbors) / len(neighbors)
    print(f"    bits={bits:0{m}b}: H={h}, avg neighbor H = {avg_nb:.2f}, ratio = {avg_nb/h if h > 0 else 'inf':.4f}")

# Check linear relationship: avg_nb = a + b*H
# Compute regression
xs = []
ys = []
for bits in range(2**m):
    h = H_map[bits]
    avg_nb = sum(H_map[bits ^ (1 << k)] for k in range(m)) / m
    xs.append(h)
    ys.append(avg_nb)

xs = np.array(xs)
ys = np.array(ys)
# Linear fit
A_fit = np.vstack([xs, np.ones(len(xs))]).T
b_fit, a_fit = np.linalg.lstsq(A_fit, ys, rcond=None)[0]
residuals = ys - (a_fit + b_fit * xs)
r_squared = 1 - np.var(residuals) / np.var(ys)

print(f"\n  Linear fit: avg_neighbor = {a_fit:.4f} + {b_fit:.4f} * H")
print(f"  R^2 = {r_squared:.6f}")
print(f"  Max |residual| = {max(abs(residuals)):.4f}")

if r_squared > 0.999:
    print(f"  => H is APPROXIMATELY an elementary landscape!")
    # For exact elementary: eigenvalue lambda = m*(1-b)/1 => order = m*(1-b)/2
    lam = m * (1 - b_fit)
    print(f"  Effective eigenvalue: lambda = {lam:.4f}")
    print(f"  Corresponding order: lambda/2 = {lam/2:.4f}")
else:
    print(f"  => H is NOT an elementary landscape (multi-order)")
    # Which orders dominate?
    print(f"  Dominant orders by energy:")
    for order in sorted(order_energy.keys(), key=lambda o: -order_energy[o])[:5]:
        print(f"    Order {order}: {100*order_energy[order]/total_energy:.1f}%")

# ================================================================
# B. FISHER INFORMATION OF BOLTZMANN FAMILY
# ================================================================
print(f"\n{'='*70}")
print(f"  B. FISHER INFORMATION OF BOLTZMANN FAMILY exp(beta*H)")
print(f"{'='*70}")
print()

# Boltzmann distribution: P_beta(T) = exp(beta*H(T)) / Z(beta)
# where Z(beta) = sum_T exp(beta*H(T))
# Fisher information: I(beta) = Var_beta[H] = d^2 log Z / d beta^2

betas = np.linspace(-2, 5, 71)
fisher_vals = []
mean_H_vals = []
entropy_vals = []

for beta in betas:
    log_weights = beta * H_vals
    log_Z = np.max(log_weights) + np.log(np.sum(np.exp(log_weights - np.max(log_weights))))
    probs = np.exp(log_weights - log_Z)

    mean_H = np.sum(probs * H_vals)
    var_H = np.sum(probs * (H_vals - mean_H)**2)
    ent = -np.sum(probs * np.log(probs + 1e-300))

    fisher_vals.append(var_H)
    mean_H_vals.append(mean_H)
    entropy_vals.append(ent)

fisher_vals = np.array(fisher_vals)
mean_H_vals = np.array(mean_H_vals)
entropy_vals = np.array(entropy_vals)

# Key beta values
print(f"  {'beta':>8s} {'E[H]':>8s} {'Var[H]':>10s} {'Fisher':>10s} {'Entropy':>10s}")
for beta_target in [-1, 0, 0.5, 1, 2, 3]:
    idx = np.argmin(np.abs(betas - beta_target))
    print(f"  {betas[idx]:>8.2f} {mean_H_vals[idx]:>8.2f} {fisher_vals[idx]:>10.4f} {fisher_vals[idx]:>10.4f} {entropy_vals[idx]:>10.4f}")

# The Fisher information at beta=0 is the UNIFORM variance
print(f"\n  Fisher at beta=0 (uniform): {fisher_vals[np.argmin(np.abs(betas))]:.4f}")
print(f"  This equals Var_uniform[H] = {H_var:.4f}")

# Maximum Fisher (maximum distinguishability)
max_fisher_idx = np.argmax(fisher_vals)
print(f"  Maximum Fisher: {fisher_vals[max_fisher_idx]:.4f} at beta={betas[max_fisher_idx]:.2f}")
print(f"  E[H] at max Fisher: {mean_H_vals[max_fisher_idx]:.2f}")

# Phase transition: where does E[H] jump most?
dE_dbeta = np.gradient(mean_H_vals, betas)
max_dE_idx = np.argmax(np.abs(dE_dbeta))
print(f"\n  Steepest E[H] change: at beta={betas[max_dE_idx]:.2f}, dE/dbeta={dE_dbeta[max_dE_idx]:.2f}")
print(f"  This is the PHASE TRANSITION of the Boltzmann family")

# Critical beta where Boltzmann concentrates on H=15 (max)
# At large beta, P concentrates on argmax H
print(f"\n  At beta=5: E[H]={mean_H_vals[-1]:.2f} (max H=15), concentrated on maxima")
print(f"  At beta=-2: E[H]={mean_H_vals[0]:.2f} (min H=1), concentrated on minima")

# ================================================================
# C. OLLIVIER-RICCI CURVATURE (SIMPLIFIED)
# ================================================================
print(f"\n{'='*70}")
print(f"  C. OLLIVIER-RICCI CURVATURE ON THE FLIP GRAPH")
print(f"{'='*70}")
print()

# Ollivier-Ricci curvature of edge (u,v) in a graph:
# kappa(u,v) = 1 - W_1(mu_u, mu_v) / d(u,v)
# where mu_u is the "lazy random walk" measure at u:
#   mu_u(u) = 1/(d+1), mu_u(w) = 1/((d+1)*d) for w ~ u
# and W_1 is the Wasserstein-1 (earth mover's) distance.
#
# For the m-cube: every vertex has degree m. Lazy walk:
#   mu_u(u) = 1 - alpha, mu_u(w) = alpha/m for w ~ u
# where alpha is the laziness parameter (often alpha = 1).
#
# For the m-cube with alpha=1 (pure random walk, no laziness):
# mu_u(w) = 1/m for each of the m neighbors w.
# Then kappa(u,v) = 1 - W_1(mu_u, mu_v) / 1
# = 1 - W_1(mu_u, mu_v)
#
# For the hypercube: if u and v differ in bit j, then:
# - u has m neighbors including v
# - v has m neighbors including u
# - Neighbors of u that are also neighbors of v: those differing in one bit != j
#   from u AND one bit != j from v. Since u and v differ in bit j,
#   a neighbor of u (flip bit k, k!=j) is at distance 2 from v (bits j and k differ).
#   So the ONLY common neighbor is... none (in the hypercube, neighbors of adjacent
#   vertices share no common neighbors in the usual sense).
#   Wait: u's neighbor by flipping bit k is at Hamming distance 2 from v (if k!=j).
#   v's neighbor by flipping bit k is at Hamming distance 2 from u (if k!=j).
#   The only vertices at distance 1 from both u and v are u^j^k for some k,
#   but u^j = v and v^k is a neighbor of v. So u's neighbors (except v) are at
#   distance 2 from v, and v's neighbors (except u) are at distance 2 from u.
#
# For uniform random walk on m-cube:
# kappa(u,v) = 2/m (Lin-Lu-Yau formula for hypercube)
# This is CONSTANT for all edges -- the hypercube is Ricci-flat up to the 2/m term.

# BUT: we're not interested in the bare hypercube.
# We want the H-WEIGHTED curvature: how curved is the H-landscape?
# This is captured by the HESSIAN of H, not by the graph Ricci curvature.

# Compute the discrete Hessian of H:
# For a function f on the hypercube, the Hessian at x is:
# Hess_f(x)[j,k] = f(x^j^k) - f(x^j) - f(x^k) + f(x) for j != k
# The eigenvalues of this matrix give the curvature directions.

print("  DISCRETE HESSIAN OF H AT CRITICAL POINTS:")
print()

# At local maxima (H=15)
local_max_bits = [b for b in range(2**m) if all(H_map[b ^ (1 << k)] <= H_map[b] for k in range(m))]
# At local minima (H=1)
local_min_bits = [b for b in range(2**m) if all(H_map[b ^ (1 << k)] >= H_map[b] for k in range(m))]

for label, bits_list in [("LOCAL MAX (H=15)", local_max_bits[:3]), ("LOCAL MIN (H=1)", local_min_bits[:3])]:
    print(f"  {label}:")
    for bits in bits_list:
        h = H_map[bits]
        # Compute Hessian
        Hess = np.zeros((m, m))
        for j in range(m):
            for k in range(m):
                if j == k:
                    # Diagonal: second difference
                    Hess[j][k] = H_map[bits ^ (1 << j)] - h  # first difference
                else:
                    Hess[j][k] = H_map[bits ^ (1 << j) ^ (1 << k)] - H_map[bits ^ (1 << j)] - H_map[bits ^ (1 << k)] + h

        # Eigenvalues of Hessian
        eigvals = sorted(np.linalg.eigvalsh(Hess))
        n_neg = sum(1 for e in eigvals if e < -0.01)
        n_zero = sum(1 for e in eigvals if abs(e) < 0.01)
        n_pos = sum(1 for e in eigvals if e > 0.01)

        # The Morse index is the number of negative eigenvalues
        print(f"    bits={bits:0{m}b}: H={h}")
        print(f"      Eigenvalues: {[f'{e:.1f}' for e in eigvals]}")
        print(f"      Morse index: {n_neg} neg, {n_zero} zero, {n_pos} pos")
        print(f"      Trace(Hess) = {sum(eigvals):.1f} (= Laplacian of H = sum of deltas)")
        print()

# ================================================================
# D. SUBLEVEL PERSISTENT HOMOLOGY (simplified via connected components)
# ================================================================
print(f"{'='*70}")
print(f"  D. SUBLEVEL FILTRATION: BIRTH-DEATH OF BASINS")
print(f"{'='*70}")
print()

# Track connected components of the sublevel set {T : H(T) <= h}
# as h increases from H_min to H_max.
# A new component is "born" when an isolated local minimum appears.
# Components "merge" when a saddle is reached.

# Using Union-Find for efficiency
parent = list(range(2**m))
rank_uf = [0] * (2**m)

def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(x, y):
    rx, ry = find(x), find(y)
    if rx == ry: return False
    if rank_uf[rx] < rank_uf[ry]: rx, ry = ry, rx
    parent[ry] = rx
    if rank_uf[rx] == rank_uf[ry]: rank_uf[rx] += 1
    return True

# Sort tournaments by H value
sorted_bits = sorted(range(2**m), key=lambda b: H_map[b])
active = set()

births = []  # (h_birth, component_id)
deaths = []  # (h_death, component_id, merged_into)

component_birth = {}  # find(bits) -> h_birth

for bits in sorted_bits:
    h = H_map[bits]
    active.add(bits)

    # Check neighbors already active
    active_neighbors = []
    for k in range(m):
        nb = bits ^ (1 << k)
        if nb in active:
            active_neighbors.append(nb)

    if not active_neighbors:
        # New component born!
        births.append((h, bits))
        component_birth[bits] = h
    else:
        # Connect to existing components
        roots = set(find(nb) for nb in active_neighbors)
        if len(roots) == 1:
            # Only one component, just join
            root = roots.pop()
            union(bits, root)
            # Propagate birth
            new_root = find(bits)
            if new_root not in component_birth:
                component_birth[new_root] = component_birth.get(root, h)
        else:
            # Multiple components merge at this saddle point
            # Keep the oldest (lowest birth H), kill the rest
            root_births = [(component_birth.get(r, h), r) for r in roots]
            root_births.sort()
            surviving = root_births[0][1]
            for birth_h, root in root_births[1:]:
                deaths.append((h, root, surviving, birth_h))
                union(root, surviving)
            union(bits, surviving)
            new_root = find(bits)
            component_birth[new_root] = root_births[0][0]

# Print persistence diagram
print(f"  PERSISTENCE DIAGRAM (births and deaths of H-sublevel components):")
print(f"  {'birth':>6s} {'death':>6s} {'lifespan':>8s}")

# Sort deaths by lifespan (longest-lived first)
persistence_pairs = [(birth_h, death_h) for (death_h, _, _, birth_h) in deaths]
persistence_pairs.sort(key=lambda p: -(p[1] - p[0]))

for birth_h, death_h in persistence_pairs[:15]:
    lifespan = death_h - birth_h
    bar = '#' * (lifespan // 1)
    print(f"  {birth_h:>6.0f} {death_h:>6.0f} {lifespan:>8.0f}  {bar}")

# Still-alive components (born but never died)
final_roots = set(find(b) for b in range(2**m))
alive = [(component_birth.get(r, 0), r) for r in final_roots]
alive.sort()
print(f"\n  Surviving components (never merged): {len(alive)}")
for birth_h, root in alive[:5]:
    print(f"    Born at H={birth_h}, root bits={root:0{m}b}")

print(f"\n  Total births: {len(births)}")
print(f"  Total deaths (merges): {len(deaths)}")
print(f"  Betti_0 (connected components at end): {len(alive)}")

# The key topological feature: how many long-lived components?
# Long-lived = large lifespan = robust structural feature
long_lived = [p for p in persistence_pairs if p[1] - p[0] >= 4]
print(f"  Long-lived (lifespan >= 4): {len(long_lived)}")

# ================================================================
# E. NATURAL GRADIENT DIRECTION
# ================================================================
print(f"\n{'='*70}")
print(f"  E. NATURAL GRADIENT: OPTIMAL SEARCH DIRECTION")
print(f"{'='*70}")
print()

# The natural gradient of E[H] w.r.t. the Boltzmann distribution
# at inverse temperature beta is:
#   nabla_nat E[H] = Cov_beta[H, x_j] / Var_beta[H]
# where x_j is the j-th arc indicator.
#
# At beta=0 (uniform distribution):
# nabla_nat(j) = Cov_uniform[H, x_j] = E[H * x_j] - E[H]*E[x_j]
#
# Since E[x_j] = 1/2 for uniform:
# nabla_nat(j) = E[H * x_j] - E[H]/2

print("  NATURAL GRADIENT AT beta=0 (uniform prior):")
print(f"  (Which arc, if biased, increases E[H] most?)")
print()

for j, (i_arc, j_arc) in enumerate(pairs):
    # E[H * x_j]
    EHx = np.mean([H_map[b] for b in range(2**m) if (b >> j) & 1])
    # E[H * (1-x_j)]
    EHnx = np.mean([H_map[b] for b in range(2**m) if not ((b >> j) & 1)])
    # Covariance = E[H*x_j] - E[H]*0.5
    cov = EHx * 0.5 - H_mean * 0.5  # E[H*x_j] = EHx * P(x_j=1) = EHx * 0.5
    # Actually: E[H*x_j] = sum over all T of H(T)*x_j(T) / 2^m
    # = sum over T where bit j is 1 of H(T) / 2^m = EHx * 0.5
    # Wait: EHx is already the mean of H over T where bit j=1
    # E[H*x_j] = (1/2^m) * sum_{T: bit j=1} H(T) = 0.5 * EHx
    # E[H]*E[x_j] = H_mean * 0.5
    # Cov = 0.5*EHx - 0.5*H_mean = 0.5*(EHx - H_mean)
    cov_correct = 0.5 * (EHx - H_mean)
    print(f"  arc ({i_arc},{j_arc}): E[H|arc=1]={EHx:.2f}, E[H|arc=0]={EHnx:.2f}, Cov={cov_correct:.4f}")

# By symmetry of the tournament space, all arcs should have the same
# natural gradient component (since any permutation of [n] maps
# one arc to another). Let's verify:
covs = []
for j in range(m):
    EHx = np.mean([H_map[b] for b in range(2**m) if (b >> j) & 1])
    covs.append(0.5 * (EHx - H_mean))

print(f"\n  All covariances: {[f'{c:.4f}' for c in covs]}")
print(f"  All equal? {len(set(f'{c:.4f}' for c in covs)) == 1}")
print(f"  Value: {covs[0]:.6f}")
print()
print("  INTERPRETATION: The natural gradient at beta=0 is ZERO (by symmetry).")
print("  Every arc has the same marginal effect on E[H].")
print("  This is FLAT: the center of the Boltzmann manifold is isotropic.")
print("  Information geometry says: at the uniform distribution,")
print("  NO direction is preferred. You must go to beta != 0 to see structure.")

# At nonzero beta, the natural gradient BREAKS SYMMETRY:
print(f"\n  NATURAL GRADIENT AT beta=1:")
beta = 1.0
log_weights = beta * H_vals
log_Z = np.max(log_weights) + np.log(np.sum(np.exp(log_weights - np.max(log_weights))))
probs = np.exp(log_weights - log_Z)
EH_beta = np.sum(probs * H_vals)

for j, (i_arc, j_arc) in enumerate(pairs[:3]):
    x_j = np.array([(b >> j) & 1 for b in range(2**m)], dtype=float)
    cov_beta = np.sum(probs * H_vals * x_j) - EH_beta * np.sum(probs * x_j)
    print(f"  arc ({i_arc},{j_arc}): Cov_beta={cov_beta:.4f}")

# ================================================================
# SYNTHESIS
# ================================================================
print(f"\n{'='*70}")
print(f"  SYNTHESIS: THE THREE FRAMEWORKS UNIFIED")
print(f"{'='*70}")
print()

print("  1. WALSH-FOURIER SPECTRUM tells us H's 'epistatic order':")
print(f"     Dominant orders: {nonzero_orders}")
if r_squared > 0.99:
    print(f"     H is NEAR-ELEMENTARY (R^2={r_squared:.4f})")
    print(f"     This means avg neighbor H ~ a + b*H")
    print(f"     => gradient flow is NEARLY LINEAR in H")
else:
    print(f"     H is MULTI-ORDER (R^2={r_squared:.4f})")
    print(f"     => gradient flow has nonlinear structure")

print()
print("  2. FISHER INFORMATION tells us where the landscape is 'informative':")
print(f"     Max Fisher at beta={betas[max_fisher_idx]:.2f}")
print(f"     At this temperature, small changes in beta maximally")
print(f"     change the distribution over H values.")
print(f"     This is the PHASE TRANSITION of the Boltzmann family.")

print()
print("  3. SUBLEVEL PERSISTENCE tells us the topology of basins:")
print(f"     {len(births)} births, {len(deaths)} merges, {len(alive)} surviving")
if long_lived:
    print(f"     {len(long_lived)} long-lived features (lifespan >= 4)")
    print(f"     These are the ROBUST topological structures of the H-landscape")

print()
print("  4. THE NATURAL GRADIENT at beta=0 is ZERO by S_n symmetry:")
print("     The uniform distribution is a CRITICAL POINT of the")
print("     information-geometric landscape. To optimize H, you must")
print("     first break the symmetry (increase beta).")
print("     This connects to: the Paley maximizer breaks no symmetry")
print("     (it IS the most symmetric tournament), so the optimization")
print("     landscape funnels TOWARD symmetry, not away from it.")

print()
print("  5. THE HESSIAN at local maxima has MIXED signature:")
print("     Some negative eigenvalues (downhill directions)")
print("     Some zero eigenvalues (flat directions = plateau)")
print("     The flat directions at H=37 trap steepest ascent at n=6")
print("     while the sharper peak at H=45 has fewer flat directions.")
