#!/usr/bin/env python3
"""
N=6 DISCRIMINANT SEARCH — What invariant completes H = f(t3, t5, ???) at n=6?

opus-2026-03-14-S89b

At n≤4: H = 1 + 2t₃ (exact)
At n=5: H = 1 + 2t₃ + 2t₅ (exact, R²=1)
At n=6: H = 1 + 2(t₃+t₅) only 53% correct — need MORE.

Candidates for the missing invariant:
1. Edge-overlap count: how many pairs of 3-cycles share an edge
2. Triangle graph density: the clique number or chromatic number of the triangle overlap graph
3. Score sequence entropy
4. Number of directed 4-paths (not cycles — path structure)
5. Hypercube distance from transitive tournament (Hamming)
6. The NUMBER of Hamiltonian 5-cycles (not just any 5-cycles)
7. π-related: the binary expansion of H values

Key insight from THM-208: the coefficient 2 = Φ₂(1) at n≤4, and the slope
becomes 3 = Φ₃(1) at n=5 and 6 = |S₃| at n=6. But these are REGRESSION
slopes, not exact formula coefficients. We need the EXACT formula.

APPROACH: Enumerate all 32768 tournaments at n=6. For each, compute:
- H (via DP bitmask)
- t₃ (3-cycle count)
- t₅ (5-cycle count)
- e₃₃ (number of edge-sharing 3-cycle pairs)
- v₃ (number of vertices in at least one 3-cycle)
- t₃₃_share (number of 3-cycle pairs sharing exactly 2 vertices = sharing an edge)
- score sequence as tuple
- Hamming distance to nearest transitive tournament

Then search for f such that H = 1 + 2t₃ + 2t₅ + f(extra) exactly.
"""

from math import factorial, comb
from collections import Counter, defaultdict
from itertools import combinations
import sys

def compute_H(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj.get((v, u), 0) == 1:
                    new_mask = mask | (1 << u)
                    dp[(new_mask, u)] = dp.get((new_mask, u), 0) + dp[(mask, v)]
    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))

def tournament_from_bits(n, bits):
    adj = {}
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
            idx += 1
    return adj

def count_3cycles(n, adj):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check if i->j->k->i or i->k->j->i
                a = adj.get((i,j),0)
                b = adj.get((j,k),0)
                c = adj.get((k,i),0)
                if a == 1 and b == 1 and c == 1:
                    count += 1
                elif a == 0 and b == 0 and c == 0:
                    count += 1
    return count

def count_5cycles(n, adj):
    """Count directed 5-cycles in a tournament."""
    count = 0
    for combo in combinations(range(n), 5):
        verts = list(combo)
        # Check all possible 5-cycle orderings
        from itertools import permutations as perms
        for perm in perms(verts):
            is_cycle = True
            for idx in range(5):
                if adj.get((perm[idx], perm[(idx+1)%5]), 0) != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    # Each 5-cycle is counted 5 times (cyclic rotations) * 2 (directions... no, tournament direction is fixed)
    # Actually each DIRECTED 5-cycle appears 5 times from cyclic rotation
    return count // 5

def get_3cycle_triples(n, adj):
    """Return list of (i,j,k) triples that form 3-cycles."""
    triples = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                a = adj.get((i,j),0)
                b = adj.get((j,k),0)
                c = adj.get((k,i),0)
                if (a == 1 and b == 1 and c == 1) or (a == 0 and b == 0 and c == 0):
                    triples.append((i,j,k))
    return triples

def edge_sharing_pairs(triples):
    """Count pairs of 3-cycles sharing at least one edge (= sharing exactly 2 vertices)."""
    count = 0
    for idx1 in range(len(triples)):
        for idx2 in range(idx1+1, len(triples)):
            s1 = set(triples[idx1])
            s2 = set(triples[idx2])
            if len(s1 & s2) == 2:
                count += 1
    return count

def vertex_sharing_pairs(triples):
    """Count pairs of 3-cycles sharing exactly 1 vertex."""
    count = 0
    for idx1 in range(len(triples)):
        for idx2 in range(idx1+1, len(triples)):
            s1 = set(triples[idx1])
            s2 = set(triples[idx2])
            if len(s1 & s2) == 1:
                count += 1
    return count

def disjoint_pairs(triples):
    """Count pairs of 3-cycles sharing no vertex."""
    count = 0
    for idx1 in range(len(triples)):
        for idx2 in range(idx1+1, len(triples)):
            s1 = set(triples[idx1])
            s2 = set(triples[idx2])
            if len(s1 & s2) == 0:
                count += 1
    return count

def score_sequence(n, adj):
    scores = []
    for i in range(n):
        s = sum(adj.get((i,j), 0) for j in range(n) if j != i)
        scores.append(s)
    return tuple(sorted(scores))

def hamming_to_transitive(n, bits):
    """Min Hamming distance to any transitive tournament."""
    m = n*(n-1)//2
    # Generate all transitive tournaments
    from itertools import permutations
    min_dist = m
    for perm in permutations(range(n)):
        # Transitive: perm[0] beats all, perm[1] beats all except perm[0], etc.
        trans_bits = 0
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                # i->j iff perm^{-1}(i) < perm^{-1}(j), i.e., i is "higher ranked"
                inv_perm = [0]*n
                for pos, v in enumerate(perm):
                    inv_perm[v] = pos
                if inv_perm[i] < inv_perm[j]:
                    trans_bits |= (1 << idx)
                idx += 1
        dist = bin(bits ^ trans_bits).count('1')
        min_dist = min(min_dist, dist)
    return min_dist

print("="*70)
print("N=6 DISCRIMINANT SEARCH")
print("What invariant completes H = 1 + 2t₃ + 2t₅ + ??? at n=6?")
print("="*70)

n = 6
m = n*(n-1)//2  # 15
total = 1 << m   # 32768

print(f"\nEnumerating all {total} tournaments on n={n} vertices (m={m} arcs)...")

# Collect all data
data = []  # list of (bits, H, t3, t5, e33, v33, d33, score, ham_dist)

# Pre-compute transitive tournament bit patterns for Hamming distance
from itertools import permutations
trans_patterns = set()
for perm in permutations(range(n)):
    inv_perm = [0]*n
    for pos, v in enumerate(perm):
        inv_perm[v] = pos
    bits_val = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if inv_perm[i] < inv_perm[j]:
                bits_val |= (1 << idx)
            idx += 1
    trans_patterns.add(bits_val)

print(f"Number of transitive tournaments: {len(trans_patterns)} (should be {factorial(n)}={factorial(n)})")
# Hmm, n=6 has 720 transitive tournaments but we only need the bit patterns
# Actually 720 permutations might give duplicates... no, each permutation gives a distinct tournament
# because the ranking fully determines all arcs.

progress_step = total // 10
for bits in range(total):
    if bits % progress_step == 0:
        pct = 100 * bits // total
        print(f"  {pct}% done...", file=sys.stderr)

    adj = tournament_from_bits(n, bits)
    H = compute_H(n, adj)
    t3 = count_3cycles(n, adj)

    # 5-cycles: use optimized counting
    t5 = count_5cycles(n, adj)

    # 3-cycle overlap structure
    triples = get_3cycle_triples(n, adj)
    e33 = edge_sharing_pairs(triples)
    v33 = vertex_sharing_pairs(triples)
    d33 = disjoint_pairs(triples)

    sc = score_sequence(n, adj)

    # Hamming to nearest transitive
    ham = min(bin(bits ^ tp).count('1') for tp in trans_patterns)

    data.append((bits, H, t3, t5, e33, v33, d33, sc, ham))

print(f"\nDone! Computed {len(data)} tournaments.")

# ===== ANALYSIS =====

print("\n" + "="*70)
print("PART 1: Residual ε = H - 1 - 2(t₃ + t₅)")
print("="*70)

residuals = defaultdict(list)
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    eps = H - 1 - 2*(t3 + t5)
    residuals[(t3, t5)].append((eps, e33, v33, d33, sc, ham, H, bits))

# Show groups where residual is not constant
non_const = 0
for key in sorted(residuals.keys()):
    vals = residuals[key]
    eps_vals = [v[0] for v in vals]
    if len(set(eps_vals)) > 1:
        non_const += 1
        if non_const <= 10:
            print(f"\n  (t₃={key[0]}, t₅={key[1]}): ε ∈ {sorted(set(eps_vals))}, count={len(vals)}")
            # Check if e33 discriminates
            e33_groups = defaultdict(list)
            for eps_v, e33_v, v33_v, d33_v, sc_v, ham_v, H_v, bits_v in vals:
                e33_groups[e33_v].append(eps_v)
            e33_unique = all(len(set(v)) == 1 for v in e33_groups.values())
            print(f"    e₃₃ (edge-sharing pairs) discriminates? {e33_unique}")
            if e33_unique:
                for e33_val in sorted(e33_groups.keys()):
                    print(f"      e₃₃={e33_val} → ε={e33_groups[e33_val][0]}")

            # Check if (e33, v33) discriminates
            ev_groups = defaultdict(list)
            for eps_v, e33_v, v33_v, d33_v, sc_v, ham_v, H_v, bits_v in vals:
                ev_groups[(e33_v, v33_v)].append(eps_v)
            ev_unique = all(len(set(v)) == 1 for v in ev_groups.values())
            print(f"    (e₃₃, v₃₃) discriminates? {ev_unique}")

            # Check if Hamming distance discriminates
            ham_groups = defaultdict(list)
            for eps_v, e33_v, v33_v, d33_v, sc_v, ham_v, H_v, bits_v in vals:
                ham_groups[ham_v].append(eps_v)
            ham_unique = all(len(set(v)) == 1 for v in ham_groups.values())
            print(f"    Hamming dist discriminates? {ham_unique}")

            # Check if score sequence discriminates
            sc_groups = defaultdict(list)
            for eps_v, e33_v, v33_v, d33_v, sc_v, ham_v, H_v, bits_v in vals:
                sc_groups[sc_v].append(eps_v)
            sc_unique = all(len(set(v)) == 1 for v in sc_groups.values())
            print(f"    Score seq discriminates? {sc_unique}")

print(f"\nTotal non-constant (t₃,t₅) groups: {non_const} / {len(residuals)}")

print("\n" + "="*70)
print("PART 2: Does (t₃, t₅, e₃₃) determine H?")
print("="*70)

groups_3 = defaultdict(set)
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    groups_3[(t3, t5, e33)].add(H)

non_det_3 = sum(1 for v in groups_3.values() if len(v) > 1)
print(f"(t₃, t₅, e₃₃) groups: {len(groups_3)}")
print(f"Non-determined: {non_det_3}")

if non_det_3 == 0:
    print("*** (t₃, t₅, e₃₃) DETERMINES H at n=6! ***")
else:
    print(f"\n(t₃, t₅, e₃₃) leaves {non_det_3} groups non-determined")
    # Show first few
    cnt = 0
    for key in sorted(groups_3.keys()):
        if len(groups_3[key]) > 1:
            cnt += 1
            if cnt <= 5:
                print(f"  {key} → H ∈ {sorted(groups_3[key])}")

print("\n" + "="*70)
print("PART 3: Does (t₃, t₅, e₃₃, v₃₃) determine H?")
print("="*70)

groups_4 = defaultdict(set)
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    groups_4[(t3, t5, e33, v33)].add(H)

non_det_4 = sum(1 for v in groups_4.values() if len(v) > 1)
print(f"(t₃, t₅, e₃₃, v₃₃) groups: {len(groups_4)}")
print(f"Non-determined: {non_det_4}")

if non_det_4 == 0:
    print("*** (t₃, t₅, e₃₃, v₃₃) DETERMINES H at n=6! ***")

print("\n" + "="*70)
print("PART 4: Does (t₃, t₅, e₃₃, d₃₃) determine H?")
print("="*70)

groups_4b = defaultdict(set)
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    groups_4b[(t3, t5, e33, d33)].add(H)

non_det_4b = sum(1 for v in groups_4b.values() if len(v) > 1)
print(f"(t₃, t₅, e₃₃, d₃₃) groups: {len(groups_4b)}")
print(f"Non-determined: {non_det_4b}")

print("\n" + "="*70)
print("PART 5: Does (t₃, t₅, score_seq) determine H?")
print("="*70)

groups_sc = defaultdict(set)
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    groups_sc[(t3, t5, sc)].add(H)

non_det_sc = sum(1 for v in groups_sc.values() if len(v) > 1)
print(f"(t₃, t₅, score_seq) groups: {len(groups_sc)}")
print(f"Non-determined: {non_det_sc}")

if non_det_sc == 0:
    print("*** (t₃, t₅, score_seq) DETERMINES H at n=6! ***")
else:
    print(f"\n(t₃, t₅, score_seq) leaves {non_det_sc} groups non-determined")
    cnt = 0
    for key in sorted(groups_sc.keys()):
        if len(groups_sc[key]) > 1:
            cnt += 1
            if cnt <= 5:
                print(f"  t₃={key[0]}, t₅={key[1]}, score={key[2]} → H ∈ {sorted(groups_sc[key])}")

print("\n" + "="*70)
print("PART 6: Does (t₃, t₅, ham_dist) determine H?")
print("="*70)

groups_ham = defaultdict(set)
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    groups_ham[(t3, t5, ham)].add(H)

non_det_ham = sum(1 for v in groups_ham.values() if len(v) > 1)
print(f"(t₃, t₅, ham_dist) groups: {len(groups_ham)}")
print(f"Non-determined: {non_det_ham}")

print("\n" + "="*70)
print("PART 7: Linear regression H ~ a + b·t₃ + c·t₅ + d·e₃₃")
print("="*70)

# Collect arrays
import numpy as np
Hs = np.array([d[1] for d in data], dtype=float)
t3s = np.array([d[2] for d in data], dtype=float)
t5s = np.array([d[3] for d in data], dtype=float)
e33s = np.array([d[4] for d in data], dtype=float)
v33s = np.array([d[5] for d in data], dtype=float)
d33s = np.array([d[6] for d in data], dtype=float)
hams = np.array([d[8] for d in data], dtype=float)

# OLS: H ~ 1 + t3 + t5 + e33
X = np.column_stack([np.ones(len(data)), t3s, t5s, e33s])
beta = np.linalg.lstsq(X, Hs, rcond=None)[0]
pred = X @ beta
resid = Hs - pred
r2 = 1 - np.var(resid)/np.var(Hs)
print(f"H = {beta[0]:.4f} + {beta[1]:.4f}·t₃ + {beta[2]:.4f}·t₅ + {beta[3]:.4f}·e₃₃")
print(f"R² = {r2:.6f}")
print(f"Max |residual| = {np.max(np.abs(resid)):.4f}")
exact = np.sum(np.abs(resid) < 0.001)
print(f"Exact matches: {exact}/{len(data)}")

# Try H ~ 1 + t3 + t5 + e33 + v33
X2 = np.column_stack([np.ones(len(data)), t3s, t5s, e33s, v33s])
beta2 = np.linalg.lstsq(X2, Hs, rcond=None)[0]
pred2 = X2 @ beta2
resid2 = Hs - pred2
r2_2 = 1 - np.var(resid2)/np.var(Hs)
print(f"\nH = {beta2[0]:.4f} + {beta2[1]:.4f}·t₃ + {beta2[2]:.4f}·t₅ + {beta2[3]:.4f}·e₃₃ + {beta2[4]:.4f}·v₃₃")
print(f"R² = {r2_2:.6f}")
print(f"Max |residual| = {np.max(np.abs(resid2)):.4f}")

# Try H ~ 1 + t3 + t5 + e33 + d33 + v33
X3 = np.column_stack([np.ones(len(data)), t3s, t5s, e33s, d33s, v33s])
beta3 = np.linalg.lstsq(X3, Hs, rcond=None)[0]
pred3 = X3 @ beta3
resid3 = Hs - pred3
r2_3 = 1 - np.var(resid3)/np.var(Hs)
print(f"\nH = {beta3[0]:.4f} + {beta3[1]:.4f}·t₃ + {beta3[2]:.4f}·t₅ + {beta3[3]:.4f}·e₃₃ + {beta3[4]:.4f}·d₃₃ + {beta3[5]:.4f}·v₃₃")
print(f"R² = {r2_3:.6f}")
print(f"Max |residual| = {np.max(np.abs(resid3)):.4f}")

# Try with ham distance too
X4 = np.column_stack([np.ones(len(data)), t3s, t5s, e33s, d33s, v33s, hams])
beta4 = np.linalg.lstsq(X4, Hs, rcond=None)[0]
pred4 = X4 @ beta4
resid4 = Hs - pred4
r2_4 = 1 - np.var(resid4)/np.var(Hs)
print(f"\nWith ham_dist: R² = {r2_4:.6f}, Max |resid| = {np.max(np.abs(resid4)):.4f}")

print("\n" + "="*70)
print("PART 8: EXACT residual ε = H - 1 - 2(t₃+t₅) vs e₃₃")
print("="*70)

# Check if ε is a linear function of e33 within each (t3,t5) group
for key in sorted(residuals.keys()):
    vals = residuals[key]
    eps_vals = [v[0] for v in vals]
    if len(set(eps_vals)) > 1:
        # Try linear fit ε ~ a + b·e33 within this group
        eps_arr = np.array([v[0] for v in vals], dtype=float)
        e33_arr = np.array([v[1] for v in vals], dtype=float)
        if np.var(e33_arr) > 0:
            slope = np.cov(eps_arr, e33_arr)[0,1] / np.var(e33_arr)
            intercept = np.mean(eps_arr) - slope * np.mean(e33_arr)
            pred_eps = intercept + slope * e33_arr
            resid_eps = eps_arr - pred_eps
            if np.max(np.abs(resid_eps)) < 0.01:
                print(f"  (t₃={key[0]}, t₅={key[1]}): ε = {intercept:.1f} + {slope:.1f}·e₃₃  [EXACT!]")
            else:
                print(f"  (t₃={key[0]}, t₅={key[1]}): ε ~ {intercept:.2f} + {slope:.2f}·e₃₃  max_resid={np.max(np.abs(resid_eps)):.2f}")

print("\n" + "="*70)
print("PART 9: CROWN JEWEL TEST — does H = 1 + 2(t₃+t₅) + 2·e₃₃ work?")
print("="*70)

exact_count = 0
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    pred_h = 1 + 2*(t3 + t5) + 2*e33
    if pred_h == H:
        exact_count += 1

print(f"H = 1 + 2(t₃+t₅+e₃₃): {exact_count}/{len(data)} exact")

# Try other coefficients for e33
for coeff in [1, 2, 3, 4, 6]:
    exact_c = sum(1 for bits, H, t3, t5, e33, v33, d33, sc, ham in data
                  if H == 1 + 2*(t3+t5) + coeff*e33)
    print(f"H = 1 + 2(t₃+t₅) + {coeff}·e₃₃: {exact_c}/{len(data)} exact")

print("\n" + "="*70)
print("PART 10: The ε distribution — does ε depend on e₃₃ alone?")
print("="*70)

eps_by_e33 = defaultdict(list)
for bits, H, t3, t5, e33, v33, d33, sc, ham in data:
    eps = H - 1 - 2*(t3 + t5)
    eps_by_e33[e33].append(eps)

for e33_val in sorted(eps_by_e33.keys()):
    vals = eps_by_e33[e33_val]
    print(f"  e₃₃={e33_val}: ε range [{min(vals)}, {max(vals)}], mean={sum(vals)/len(vals):.2f}, unique={len(set(vals))}")

print("\n" + "="*70)
print("PART 11: Fibonaccian / π connections in H-spectrum")
print("="*70)

H_counter = Counter(d[1] for d in data)
H_vals = sorted(H_counter.keys())
print(f"Distinct H values at n=6: {len(H_vals)}")
print(f"H range: [{min(H_vals)}, {max(H_vals)}]")
print(f"All H values: {H_vals}")

# Check: are H values always odd?
all_odd = all(h % 2 == 1 for h in H_vals)
print(f"All H values odd? {all_odd}")

# Fibonacci connection: which H values are Fibonacci numbers?
fibs = set()
a, b = 1, 1
while b <= max(H_vals) + 100:
    fibs.add(b)
    a, b = b, a+b
fib_hits = [h for h in H_vals if h in fibs]
print(f"H values that are Fibonacci: {fib_hits}")

# Check H mod 6 distribution (Pisano period)
H_mod6 = Counter(h % 6 for h in [d[1] for d in data])
print(f"H mod 6 distribution: {dict(sorted(H_mod6.items()))}")

# π/6 = arcsin(1/2) — tournament connection?
import math
print(f"\nπ connections:")
print(f"  C(6,2) = 15 = m, and 15/π ≈ {15/math.pi:.4f}")
print(f"  32768 / π² ≈ {32768/math.pi**2:.4f}")
print(f"  Mean H = {sum(d[1] for d in data)/len(data):.4f}")
mean_H = sum(d[1] for d in data)/len(data)
print(f"  Mean H / 6! = {mean_H/720:.6f}")
print(f"  6!/π³ = {720/math.pi**3:.4f}")
print(f"  Sum of all H = {sum(d[1] for d in data)}")
total_H = sum(d[1] for d in data)
print(f"  Sum H / 2^15 = {total_H/32768:.4f}")
print(f"  Sum H / (6! · 2^15) = {total_H/(720*32768):.6f}")

print("\n" + "="*70)
print("DONE")
print("="*70)
