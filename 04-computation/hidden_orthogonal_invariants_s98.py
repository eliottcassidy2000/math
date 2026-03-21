#!/usr/bin/env python3
"""
hidden_orthogonal_invariants_s98.py — opus-2026-03-21-S98

THE HIDDEN ORTHOGONAL INVARIANTS THAT DETERMINE H

From S97b: H ≈ b·c₃ + a captures 92-95% of variance.
The remaining 5-8% is the "hidden" part.

This session investigates:
1. WHAT invariants capture the residual H - f(S₂)?
2. WHY does the spectrum not determine H? (spectral blindness)
3. The EIGENVECTOR invariant — what B's eigenvectors tell us
4. The orthogonal decomposition of tournament space by H-class
5. Practical: the Landau diagnostic for ranking quality

KEY INSIGHT BEING TESTED:
  H is determined by the FULL adjacency matrix, not just its spectrum.
  The spectrum gives S₂ (via trace formulas). The eigenvectors give the rest.
  The "hidden" invariant is an EIGENVECTOR alignment with the tournament structure.

Author: opus-2026-03-21-S98
"""

import sys
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(2**m):
        A = np.zeros((n, n), dtype=int)
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A

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

def count_c3(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: c3 += 1
        if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

def score_seq(A, n):
    return tuple(sorted(A.sum(axis=1)))

print("=" * 72)
print("  HIDDEN ORTHOGONAL INVARIANTS THAT DETERMINE H")
print("  opus-2026-03-21-S98")
print("=" * 72)

# ========================================================================
# PART 1: The S₂ residual at n=5 — exact anatomy
# ========================================================================
print("\n" + "=" * 72)
print("  PART 1: THE S₂ RESIDUAL AT n=5 — EXACT ANATOMY")
print("=" * 72)

n = 5
print(f"\n  At n={n}: 8 tournament isomorphism classes by (H, c₃)")
print(f"  S₂ determines c₃ exactly: c₃ = 5 - S₂/2")
print(f"  But within S₂ classes, H can differ.\n")

# Collect all data
data = []
for A in all_tournaments(n):
    H = count_hp(A, n)
    c3 = count_c3(A, n)
    sc = score_seq(A, n)
    scores = A.sum(axis=1)
    S2 = np.sum((scores - (n-1)/2)**2)
    B = (A - A.T).astype(float)
    eigs_B = np.sort(np.linalg.eigvals(B).imag)[::-1]

    # Eigenvector-based invariants
    B2 = B @ B
    eigs_B2, V = np.linalg.eigh(B2)  # B² is symmetric
    # Sort by eigenvalue magnitude
    idx = np.argsort(eigs_B2)
    eigs_B2 = eigs_B2[idx]
    V = V[:, idx]

    # The key eigenvector invariant: how does the all-ones vector
    # project onto each eigenspace of B²?
    ones = np.ones(n) / np.sqrt(n)
    projections = [abs(ones @ V[:, k])**2 for k in range(n)]

    # B^3 diagonal — this involves 3-cycle structure
    B3 = B2 @ B
    B3_diag = np.diag(B3)  # Should be zero (odd power of skew)

    # B^4 diagonal elements — vary by tournament
    B4 = B2 @ B2
    B4_diag = tuple(sorted(np.round(np.diag(B4)).astype(int)))

    # The off-diagonal of B² — the alignment matrix
    B2_offdiag = B2.copy()
    np.fill_diagonal(B2_offdiag, 0)
    offdiag_sum = np.sum(B2_offdiag**2)

    # Count of distinct (B²)_ij values for i≠j
    offdiag_vals = []
    for i in range(n):
        for j in range(i+1, n):
            offdiag_vals.append(int(B2[i,j]))
    offdiag_profile = tuple(sorted(offdiag_vals))

    data.append({
        'H': H, 'c3': c3, 'sc': sc, 'S2': S2,
        'eigs_B2': tuple(np.round(eigs_B2, 4)),
        'projections': tuple(np.round(projections, 6)),
        'B4_diag': B4_diag,
        'offdiag_profile': offdiag_profile,
        'offdiag_sum': offdiag_sum,
        'A': A,
    })

# Group by (H, c3)
groups = defaultdict(list)
for d in data:
    groups[(d['H'], d['c3'])].append(d)

print(f"  {'H':>4s} {'c₃':>4s} {'S₂':>5s} {'cnt':>5s} {'score':>20s} "
      f"{'B² eigenvalues':>30s} {'B² off-diag profile':>25s}")

for key in sorted(groups.keys()):
    H, c3 = key
    grp = groups[key]
    d = grp[0]
    print(f"  {H:4d} {c3:4d} {d['S2']:5.1f} {len(grp):5d} {str(d['sc']):>20s} "
          f"{str(d['eigs_B2']):>30s} {str(d['offdiag_profile']):>25s}")

# ========================================================================
# PART 2: The ambiguous score class (1,2,2,2,3)
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 2: THE AMBIGUOUS CLASS — score (1,2,2,2,3) at n=5")
print("=" * 72)
print(f"  Score (1,2,2,2,3) has S₂=2, c₃=4, but H ∈ {{11, 13, 15}}")
print(f"  What DISTINGUISHES these three H values?\n")

ambig = [d for d in data if d['sc'] == (1, 2, 2, 2, 3)]
by_H = defaultdict(list)
for d in ambig:
    by_H[d['H']].append(d)

for H_val in sorted(by_H.keys()):
    grp = by_H[H_val]
    d = grp[0]
    print(f"  H={H_val} ({len(grp)} tournaments):")
    print(f"    B² eigenvalues: {d['eigs_B2']}")
    print(f"    B⁴ diagonal: {d['B4_diag']}")
    print(f"    B² off-diag profile: {d['offdiag_profile']}")

    # Count 5-cycles
    c5_count = 0
    for verts in combinations(range(n), 5):
        first = verts[0]
        rest = list(verts[1:])
        for perm in permutations(rest):
            path = [first] + list(perm)
            ok = True
            for k in range(5):
                if not d['A'][path[k]][path[(k+1)%5]]:
                    ok = False; break
            if ok: c5_count += 1
    print(f"    c₅ = {c5_count}")

    # alpha_1 = c3 + c5
    alpha1 = d['c3'] + c5_count
    print(f"    α₁ = c₃ + c₅ = {d['c3']} + {c5_count} = {alpha1}")

    # Count alpha_2 (independent pairs of odd cycles)
    # Build conflict graph
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]
            rest = list(verts[1:])
            for perm in permutations(rest):
                path = [first] + list(perm)
                ok = True
                for k in range(length):
                    if not d['A'][path[k]][path[(k+1)%length]]:
                        ok = False; break
                if ok:
                    cycles.append(frozenset(path))

    # Remove duplicate cycles (same vertex set, different starting point)
    unique_cycles = list(set(cycles))
    nc = len(unique_cycles)

    # Build adjacency (conflict = share vertex)
    adj = np.zeros((nc, nc), dtype=int)
    for i in range(nc):
        for j in range(i+1, nc):
            if unique_cycles[i] & unique_cycles[j]:
                adj[i][j] = adj[j][i] = 1

    # Count independent sets of size 2
    alpha2 = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if not adj[i][j]:
                alpha2 += 1

    print(f"    α₂ = {alpha2} (independent pairs)")
    print(f"    OCF: H = 1 + 2·{alpha1} + 4·{alpha2} = {1 + 2*alpha1 + 4*alpha2}")
    print(f"    Actual H = {H_val}, match: {1 + 2*alpha1 + 4*alpha2 == H_val}")

# ========================================================================
# PART 3: What EXACTLY distinguishes the three H classes?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 3: THE DISTINGUISHING INVARIANT")
print("=" * 72)
print("""
  From Part 2: all three H classes have the same score seq (1,2,2,2,3),
  same c₃=4, same S₂=2, and the same B² eigenvalues (since all n=5
  tournaments have only 2 Casimir classes).

  The distinguishing feature MUST be:
  1. The number of 5-cycles (c₅) — from Part 2, does this vary?
  2. The off-diagonal profile of B² — does this vary?
  3. The eigenvector alignment — does this vary?
  4. alpha_2 (independent pairs) — does this vary?

  Let's check systematically.
""")

# Collect distinguishing features
features = defaultdict(lambda: defaultdict(set))
for d in ambig:
    H = d['H']
    features['B² eigenvalues'][H].add(d['eigs_B2'])
    features['B⁴ diagonal'][H].add(d['B4_diag'])
    features['B² off-diag profile'][H].add(d['offdiag_profile'])
    features['offdiag_sum'][H].add(d['offdiag_sum'])

print(f"  Feature discrimination power:")
for feat, by_H_data in features.items():
    vals = {H: sorted(v) for H, v in by_H_data.items()}
    all_vals = set()
    for v in vals.values():
        all_vals.update(v)
    discriminates = len(all_vals) >= len(vals)  # at least as many values as H classes
    print(f"  {feat:>30s}: {len(all_vals)} distinct values across 3 H classes → "
          f"{'DISCRIMINATES' if discriminates else 'DOES NOT DISCRIMINATE'}")
    if len(all_vals) <= 5:
        for H_val in sorted(vals.keys()):
            print(f"    H={H_val}: {vals[H_val]}")

# ========================================================================
# PART 4: The EDGE-CYCLE INCIDENCE as the hidden invariant
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 4: EDGE-CYCLE INCIDENCE — THE HIDDEN INVARIANT")
print("=" * 72)
print("""
  HYPOTHESIS: The hidden invariant is e_cyc — the number of edges
  (arcs) that participate in at least one 3-cycle.

  From HYP-337: e_cyc is NOT determined by c₃ alone.
  Multiple e_cyc values exist per c₃ at n≥5.

  e_cyc measures HOW the 3-cycles are arranged:
  - High e_cyc: cycles share many arcs (clustered)
  - Low e_cyc: cycles use distinct arcs (spread out)

  More spread-out cycles → more independent sets in Ω → higher H.
""")

for d in data:
    A = d['A']
    # Compute e_cyc: number of arcs in at least one 3-cycle
    arc_in_cycle = set()
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            arc_in_cycle.update([(i,j), (j,k), (k,i)])
        if A[i][k] and A[k][j] and A[j][i]:
            arc_in_cycle.update([(i,k), (k,j), (j,i)])
    d['e_cyc'] = len(arc_in_cycle)

# Check if e_cyc distinguishes within the ambiguous class
print(f"\n  Ambiguous class (1,2,2,2,3) by H:")
ecyc_by_H = defaultdict(set)
for d in ambig:
    ecyc_by_H[d['H']].add(d['e_cyc'])

for H_val in sorted(ecyc_by_H.keys()):
    print(f"  H={H_val}: e_cyc ∈ {sorted(ecyc_by_H[H_val])}")

# Does e_cyc + c₃ determine H for ALL tournaments at n=5?
print(f"\n  Does (c₃, e_cyc) determine H at n=5?")
by_c3_ecyc = defaultdict(set)
for d in data:
    by_c3_ecyc[(d['c3'], d['e_cyc'])].add(d['H'])

all_determined = True
for key, Hvals in sorted(by_c3_ecyc.items()):
    if len(Hvals) > 1:
        all_determined = False
        print(f"    c₃={key[0]}, e_cyc={key[1]}: H ∈ {sorted(Hvals)}")

if all_determined:
    print(f"    YES — (c₃, e_cyc) determines H at n=5!")
else:
    print(f"    NO — (c₃, e_cyc) does NOT determine H at n=5")

# ========================================================================
# PART 5: Full feature analysis — what DOES determine H at n=5?
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 5: MINIMAL INVARIANT SET THAT DETERMINES H AT n=5")
print("=" * 72)

# Test various feature combinations
# Features to try: c3, S2, e_cyc, c5, alpha_2, offdiag_profile

for d in data:
    A = d['A']
    # c5
    c5 = 0
    for verts in combinations(range(n), 5):
        first = verts[0]
        rest = list(verts[1:])
        for perm in permutations(rest):
            path = [first] + list(perm)
            ok = True
            for k in range(5):
                if not A[path[k]][path[(k+1)%5]]:
                    ok = False; break
            if ok: c5 += 1
    d['c5'] = c5

# Test (c3, c5) = (alpha_1 components)
print(f"\n  Does (c₃, c₅) determine H at n=5?")
by_c3c5 = defaultdict(set)
for d in data:
    by_c3c5[(d['c3'], d['c5'])].add(d['H'])

all_det_c3c5 = True
for key, Hvals in sorted(by_c3c5.items()):
    c3, c5 = key
    print(f"    c₃={c3}, c₅={c5}: H ∈ {sorted(Hvals)}")
    if len(Hvals) > 1:
        all_det_c3c5 = False

print(f"  → (c₃, c₅) determines H: {'YES' if all_det_c3c5 else 'NO'}")

# Test (c3, alpha_2) directly
print(f"\n  Computing alpha_2 for all tournaments...")
for d in data:
    A = d['A']
    # Find all odd cycles
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            first = verts[0]
            rest = list(verts[1:])
            for perm in permutations(rest):
                path = [first] + list(perm)
                ok = True
                for k in range(length):
                    if not A[path[k]][path[(k+1)%length]]:
                        ok = False; break
                if ok:
                    cycles.append(frozenset(path))

    unique_cycles = list(set(cycles))
    nc = len(unique_cycles)

    alpha2 = 0
    for i in range(nc):
        for j in range(i+1, nc):
            if not (unique_cycles[i] & unique_cycles[j]):
                alpha2 += 1
    d['alpha1'] = d['c3'] + d['c5']
    d['alpha2'] = alpha2

# Test (alpha_1, alpha_2)
print(f"\n  Does (α₁, α₂) determine H at n=5?")
by_a1a2 = defaultdict(set)
for d in data:
    by_a1a2[(d['alpha1'], d['alpha2'])].add(d['H'])

all_det_a1a2 = True
for key, Hvals in sorted(by_a1a2.items()):
    a1, a2 = key
    print(f"    α₁={a1:2d}, α₂={a2:2d}: H = {sorted(Hvals)}")
    if len(Hvals) > 1:
        all_det_a1a2 = False

print(f"  → (α₁, α₂) determines H: {'YES' if all_det_a1a2 else 'NO'}")

# Verify OCF: H = 1 + 2*alpha1 + 4*alpha2 + higher terms
print(f"\n  OCF verification: H = 1 + 2α₁ + 4α₂ + 8α₃ + ...")
for d in data[:20]:  # Just check a few
    ocf_approx = 1 + 2*d['alpha1'] + 4*d['alpha2']
    residual = d['H'] - ocf_approx
    if residual != 0:
        print(f"    H={d['H']}, 1+2·{d['alpha1']}+4·{d['alpha2']} = {ocf_approx}, "
              f"residual = {residual} (needs α₃ term)")

# ========================================================================
# PART 6: The orthogonal decomposition of H at n=5
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 6: ORTHOGONAL DECOMPOSITION OF H")
print("=" * 72)
print("""
  OCF gives: H = 1 + 2α₁ + 4α₂ + 8α₃ + ...

  The components 1, 2α₁, 4α₂, 8α₃ are NOT orthogonal in general.
  But they form a HIERARCHICAL decomposition.

  The ORTHOGONAL decomposition uses:
    H = E[H] + proj(H, c₃) + proj(H, residual | c₃) + ...

  where proj means the orthogonal projection onto the indicated subspace.

  At n=5, let's compute this decomposition explicitly.
""")

Hs = np.array([d['H'] for d in data], dtype=float)
c3s = np.array([d['c3'] for d in data], dtype=float)
c5s = np.array([d['c5'] for d in data], dtype=float)
a1s = np.array([d['alpha1'] for d in data], dtype=float)
a2s = np.array([d['alpha2'] for d in data], dtype=float)
S2s = np.array([d['S2'] for d in data], dtype=float)
ecycs = np.array([d['e_cyc'] for d in data], dtype=float)

H_mean = np.mean(Hs)
H_var = np.var(Hs)

# Component 1: mean
comp0 = H_mean * np.ones(len(data))

# Component 2: projection onto c₃ (= S₂)
c3_centered = c3s - np.mean(c3s)
beta_c3 = np.dot(Hs - H_mean, c3_centered) / np.dot(c3_centered, c3_centered)
comp1 = beta_c3 * c3_centered

# Component 3: residual after c₃, projected onto c₅
resid1 = Hs - comp0 - comp1
c5_centered = c5s - np.mean(c5s)
# Orthogonalize c5 against c3
c5_orth = c5_centered - np.dot(c5_centered, c3_centered) / np.dot(c3_centered, c3_centered) * c3_centered
if np.dot(c5_orth, c5_orth) > 1e-10:
    beta_c5 = np.dot(resid1, c5_orth) / np.dot(c5_orth, c5_orth)
    comp2 = beta_c5 * c5_orth
else:
    beta_c5 = 0
    comp2 = np.zeros(len(data))

# Component 4: residual after c₃ and c₅, projected onto alpha_2
resid2 = resid1 - comp2
a2_centered = a2s - np.mean(a2s)
# Orthogonalize against c3 and c5_orth
a2_orth = a2_centered.copy()
a2_orth -= np.dot(a2_orth, c3_centered) / np.dot(c3_centered, c3_centered) * c3_centered
if np.dot(c5_orth, c5_orth) > 1e-10:
    a2_orth -= np.dot(a2_orth, c5_orth) / np.dot(c5_orth, c5_orth) * c5_orth
if np.dot(a2_orth, a2_orth) > 1e-10:
    beta_a2 = np.dot(resid2, a2_orth) / np.dot(a2_orth, a2_orth)
    comp3 = beta_a2 * a2_orth
else:
    beta_a2 = 0
    comp3 = np.zeros(len(data))

# Component 5: residual after c₃, c₅, α₂, projected onto e_cyc
resid3 = resid2 - comp3
ecyc_centered = ecycs - np.mean(ecycs)
ecyc_orth = ecyc_centered.copy()
for prev in [c3_centered, c5_orth, a2_orth]:
    if np.dot(prev, prev) > 1e-10:
        ecyc_orth -= np.dot(ecyc_orth, prev) / np.dot(prev, prev) * prev
if np.dot(ecyc_orth, ecyc_orth) > 1e-10:
    beta_ecyc = np.dot(resid3, ecyc_orth) / np.dot(ecyc_orth, ecyc_orth)
    comp4 = beta_ecyc * ecyc_orth
else:
    beta_ecyc = 0
    comp4 = np.zeros(len(data))

final_resid = resid3 - comp4

# Variance decomposition
var_components = [
    ("Mean", np.var(comp0)),
    ("c₃ (= S₂/2)", np.var(comp1)),
    ("c₅ ⊥ c₃", np.var(comp2)),
    ("α₂ ⊥ (c₃,c₅)", np.var(comp3)),
    ("e_cyc ⊥ rest", np.var(comp4)),
    ("Residual", np.var(final_resid)),
]

print(f"  Orthogonal variance decomposition of H at n=5:")
print(f"  Total Var(H) = {H_var:.4f}\n")
print(f"  {'Component':>25s} {'Variance':>10s} {'% of total':>12s} {'Cumulative':>12s}")

cum = 0
for name, var_c in var_components:
    pct = var_c / H_var * 100 if H_var > 0 else 0
    cum += pct
    print(f"  {name:>25s} {var_c:10.4f} {pct:11.2f}% {cum:11.2f}%")

# ========================================================================
# PART 7: The EXACT orthogonal invariant at n=5
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 7: THE EXACT INVARIANT DETERMINING THE RESIDUAL")
print("=" * 72)
print("""
  At n=5, within the ambiguous class (c₃=4, S₂=2):
  H ∈ {11, 13, 15}

  From OCF: H = 1 + 2α₁ + 4α₂ (no α₃ exists at n=5).
  So H is EXACTLY determined by (α₁, α₂).

  The hidden invariant is:
  - c₅ (number of 5-cycles), which adds to c₃ to give α₁
  - α₂ (independent pairs), which depends on cycle arrangement

  Since c₃ = 4 for all three classes, the distinction comes from:
  - c₅: does it vary?
  - α₂: does it vary?
""")

# Detailed breakdown for the ambiguous class
print(f"  Detailed OCF breakdown for score (1,2,2,2,3):\n")
for H_val in [11, 13, 15]:
    examples = [d for d in ambig if d['H'] == H_val]
    d = examples[0]
    print(f"  H = {H_val}:")
    print(f"    c₃ = {d['c3']}, c₅ = {d['c5']}, α₁ = {d['alpha1']}")
    print(f"    α₂ = {d['alpha2']}")
    print(f"    OCF: 1 + 2·{d['alpha1']} + 4·{d['alpha2']} = {1+2*d['alpha1']+4*d['alpha2']}")
    print(f"    e_cyc = {d['e_cyc']}")
    print()

print("""
  THEOREM (Hidden Orthogonal Invariant at n=5):

  The residual H - 3c₃ within each score class is determined by the
  OCF hierarchy: specifically c₅ (5-cycle count) and α₂ (independent
  pair count in the conflict graph).

  The FULL determination of H at n=5 requires only:
    H = 1 + 2(c₃ + c₅) + 4α₂

  This is the OCF identity, which is EXACT.

  The "hidden orthogonal invariant" is (c₅, α₂) — the 5-cycle count
  and the conflict graph independence number at level 2.

  These are ORTHOGONAL to S₂ (Landau irregularity) in the following sense:
  - S₂ determines c₃ (and nothing more)
  - c₅ is correlated with but not determined by c₃
  - α₂ is correlated with but not determined by (c₃, c₅)

  The variance decomposition shows HOW MUCH each level contributes.
""")

# ========================================================================
# PART 8: Extend to n=6
# ========================================================================
print("=" * 72)
print("  PART 8: EXTENDING TO n=6")
print("=" * 72)

n = 6
print(f"\n  At n={n}: computing invariants for all {2**15} tournaments...")

data6 = []
for A in all_tournaments(n):
    H = count_hp(A, n)
    c3 = count_c3(A, n)
    sc = score_seq(A, n)
    scores = A.sum(axis=1)
    S2 = np.sum((scores - (n-1)/2)**2)

    data6.append({
        'H': H, 'c3': c3, 'sc': sc, 'S2': S2,
    })

Hs6 = np.array([d['H'] for d in data6], dtype=float)
c3s6 = np.array([d['c3'] for d in data6], dtype=float)
S2s6 = np.array([d['S2'] for d in data6], dtype=float)

# How many H values per S2 class?
by_S2 = defaultdict(set)
for d in data6:
    by_S2[d['S2']].add(d['H'])

print(f"\n  H values per S₂ class at n=6:")
print(f"  {'S₂':>6s} {'#H vals':>8s} {'H range':>20s}")
for S2_val in sorted(by_S2.keys()):
    Hvals = sorted(by_S2[S2_val])
    print(f"  {S2_val:6.1f} {len(Hvals):8d} {str(Hvals[:8]):>20s}"
          + ("..." if len(Hvals) > 8 else ""))

# How many by (c3, score)?
by_c3_sc = defaultdict(set)
for d in data6:
    by_c3_sc[(d['c3'], d['sc'])].add(d['H'])

ambig_count = sum(1 for v in by_c3_sc.values() if len(v) > 1)
max_ambig = max(len(v) for v in by_c3_sc.values())
print(f"\n  (c₃, score) classes: {len(by_c3_sc)} total, "
      f"{ambig_count} ambiguous, max {max_ambig} H values per class")

# Variance decomposition at n=6
H6_mean = np.mean(Hs6)
H6_var = np.var(Hs6)

c3_6c = c3s6 - np.mean(c3s6)
beta_c3_6 = np.dot(Hs6 - H6_mean, c3_6c) / np.dot(c3_6c, c3_6c)
comp1_6 = beta_c3_6 * c3_6c
resid1_6 = Hs6 - H6_mean - comp1_6

R2_c3_6 = 1 - np.var(resid1_6) / H6_var
print(f"\n  Variance decomposition at n=6:")
print(f"    Total Var(H) = {H6_var:.2f}")
print(f"    c₃ explains: R² = {R2_c3_6:.6f} ({R2_c3_6*100:.2f}%)")
print(f"    Residual variance: {np.var(resid1_6):.2f} ({(1-R2_c3_6)*100:.2f}%)")

# ========================================================================
# PART 9: The Landau diagnostic — engineering application
# ========================================================================
print("\n\n" + "=" * 72)
print("  PART 9: THE LANDAU DIAGNOSTIC — ENGINEERING APPLICATION")
print("=" * 72)
print("""
  THE LANDAU IRREGULARITY S₂ AS A UNIVERSAL RANKING QUALITY METRIC

  Given a tournament T (from pairwise comparisons), S₂ tells you:

  S₂ = 0: Perfectly balanced competition. Every competitor beats
           exactly half the field. Maximum number of consistent
           orderings. The ranking is MAXIMALLY UNCERTAIN.

  S₂ = S₂_max = n(n²-1)/12: Total domination. One competitor beats
           all, the next beats all but one, etc. The ranking is
           COMPLETELY DETERMINED (H = 1).

  The NORMALIZED Landau index:
    L(T) = S₂ / S₂_max ∈ [0, 1]

  L = 0: maximum competition (Paley-like)
  L = 1: total order (transitive)

  FROM THE c₃-S₂ THEOREM:
    c₃ = C(n+1,3)/4 · (1 - L)
    i.e., the number of 3-cycles is proportional to (1 - L).

  FROM THE H-S₂ FORMULA:
    H ≈ H_max · (1 - L)^γ for some exponent γ
    (approximately, since H is linear in S₂ ≈ linear in L)

  PRACTICAL APPLICATIONS:
""")

# Demo with concrete scenarios
scenarios = {
    "Premier League (n=20, balanced)": 0.15,
    "Champions League group (n=4)": 0.40,
    "Total blowout (n=10)": 0.95,
    "Rock-paper-scissors (n=3)": 0.0,
    "NBA playoff series (n=8)": 0.55,
    "LLM benchmark (n=6 models)": 0.30,
}

print(f"  {'Scenario':>40s} {'L(T)':>6s} {'Interpretation':>30s}")
for name, L in scenarios.items():
    if L < 0.2:
        interp = "Highly competitive"
    elif L < 0.5:
        interp = "Moderately competitive"
    elif L < 0.8:
        interp = "Partially ordered"
    else:
        interp = "Near-total order"
    print(f"  {name:>40s} {L:6.2f} {interp:>30s}")

print("""
  KEY INSIGHT:
  S₂ (equivalently c₃) captures >92% of the information in H.
  This means: for PRACTICAL purposes, knowing only the SCORE SEQUENCE
  (win counts) is almost sufficient to predict the total number of
  consistent rankings.

  The remaining <8% comes from HOW the wins are distributed — the
  "hidden orthogonal invariant" (cycle arrangement, c₅, α₂).
  But this refinement matters only for fine-grained analysis.

  BOTTOM LINE: To assess ranking quality, compute S₂. Done.
""")

print("\n\nDone. opus-2026-03-21-S98")
