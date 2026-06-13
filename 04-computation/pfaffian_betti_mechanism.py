#!/usr/bin/env python3
"""
pfaffian_betti_mechanism.py — WHY does Pfaffian determine topology?

At n=6 (exhaustive):
  C-phase (beta_1>0): |Pf(S)| in {1,3} => Pf^2 in {1,9}
  S-phase (beta_3>0): |Pf(S)| in {7,9} => Pf^2 in {49,81}

Key observation: |Pf|^2 = |det(S)| = product of eigenvalue magnitudes squared.
For n=6 skew S, eigenvalues are +-i*a, +-i*b, +-i*c with a,b,c > 0.
|det(S)| = a^2 * b^2 * c^2.

So |Pf| = a*b*c = product of positive skew eigenvalues.

The spectral gap is max(a,b,c) - min(a,b,c).
The spectral product is a*b*c = |Pf|.

Hypothesis: beta_1>0 (C-phase) requires ONE dominant eigenvalue and two small ones.
This gives large gap AND small product (dominant * small * small).

beta_3>0 (S-phase) requires near-degenerate eigenvalues (a ≈ b ≈ c).
This gives small gap AND large product (all similar).

Let's verify this and find the exact constraints.

Author: kind-pasteur-2026-03-08-S40
"""
import sys, os, random, time
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w')
from path_homology_v2 import path_betti_numbers
sys.stdout = _saved

def skew_adj(A, n):
    S = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            S[i][j] = A[i][j] - A[j][i]
    return S

def spectral_triple(S):
    """For 6x6 skew S, return the 3 positive imaginary eigenvalues sorted."""
    eigs = np.linalg.eigvals(S)
    pos = sorted([abs(e.imag) for e in eigs if e.imag > 0.01])
    return pos

def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                if A[i][k] and A[k][j] and A[j][i]: c3 += 1
    return c3

def H_tournament(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full])

# ===== EXHAUSTIVE n=6 with full eigenvalue data =====
print("=" * 70)
print("n=6 EXHAUSTIVE: Eigenvalue triples by topological phase")
print("=" * 70)

n = 6
m = n * (n-1) // 2
total = 1 << m
t0 = time.time()

phase_eigen = {'P': [], 'C': [], 'S': []}
phase_data = {'P': [], 'C': [], 'S': []}

for bits in range(total):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    try:
        beta = path_betti_numbers(A, n, max_dim=4)
    except:
        continue

    b1 = int(beta[1]) if len(beta) > 1 else 0
    b3 = int(beta[3]) if len(beta) > 3 else 0

    S = skew_adj(A, n)
    eig3 = spectral_triple(S)

    phase = 'S' if b3 > 0 else ('C' if b1 > 0 else 'P')
    phase_eigen[phase].append(tuple(round(e, 6) for e in eig3))

    c3 = count_3cycles(A, n)
    pf = round(np.prod(eig3))  # |Pf| = product of positive skew eigenvalues
    phase_data[phase].append({'eig3': eig3, 'c3': c3, 'pf': pf})

elapsed = time.time() - t0
print(f"Done in {elapsed:.1f}s")

# Unique eigenvalue triples per phase
print("\n" + "=" * 70)
print("UNIQUE EIGENVALUE TRIPLES PER PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    eig_set = Counter(phase_eigen[phase])
    print(f"\n{phase} ({len(phase_eigen[phase])} tournaments, {len(eig_set)} unique triples):")
    for triple, cnt in eig_set.most_common(10):
        a, b, c = triple
        product = round(a*b*c, 2)
        gap = round(c - a, 3)
        ratio = round(a/c, 3) if c > 0.01 else 0
        print(f"  ({a:.3f}, {b:.3f}, {c:.3f}): cnt={cnt}, |Pf|={product}, gap={gap}, ratio={ratio}")

# Eigenvalue product = |Pf|
print("\n" + "=" * 70)
print("EIGENVALUE PRODUCT (= |Pf|) BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    products = Counter(d['pf'] for d in phase_data[phase])
    print(f"\n{phase}: |Pf| = {dict(sorted(products.items()))}")

# Eigenvalue ratio min/max
print("\n" + "=" * 70)
print("EIGENVALUE RATIO (min/max) BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    ratios = [d['eig3'][0]/d['eig3'][2] if d['eig3'][2] > 0.01 else 0 for d in phase_data[phase]]
    if ratios:
        print(f"\n{phase}: ratio mean={np.mean(ratios):.4f}, min={min(ratios):.4f}, max={max(ratios):.4f}")

# Key question: what eigenvalue constraint forces beta_1>0?
print("\n" + "=" * 70)
print("EIGENVALUE CONSTRAINTS FOR EACH PHASE")
print("=" * 70)

# For C-phase: check if the SMALLEST eigenvalue is always 1
c_min_eigs = [d['eig3'][0] for d in phase_data['C']]
c_max_eigs = [d['eig3'][2] for d in phase_data['C']]
if c_min_eigs:
    print(f"\nC-phase: min_eig in {sorted(set(round(e,3) for e in c_min_eigs))}")
    print(f"C-phase: max_eig in {sorted(set(round(e,3) for e in c_max_eigs))}")

# For S-phase: check if all eigenvalues are the same
s_triples = set(phase_eigen['S'])
if s_triples:
    print(f"\nS-phase: all unique triples:")
    for t in sorted(s_triples):
        print(f"  {t}")

# c3 correlation with phase
print("\n" + "=" * 70)
print("c3 DISTRIBUTION BY PHASE")
print("=" * 70)

for phase in ['P', 'C', 'S']:
    c3_counts = Counter(d['c3'] for d in phase_data[phase])
    print(f"\n{phase}: c3 = {dict(sorted(c3_counts.items()))}")

# Score sequence analysis
print("\n" + "=" * 70)
print("tr(S^2) BY PHASE (= -2 * sum_ij S[i][j]^2 / 2 = -n(n-1))")
print("=" * 70)
# Actually tr(S^2) = -2*sum_{i<j} s_{ij}^2 = -2 * C(n,2) for tournament (all s_{ij} = ±1)
# So tr(S^2) is CONSTANT for all n-tournaments: -n(n-1)
# tr(S^4) = sum eigenvalues^4 = 2*(a^4 + b^4 + c^4) — varies
# Newton's: a^2+b^2+c^2 = -tr(S^2)/2 = n(n-1)/2 = 15 for n=6

for phase in ['P', 'C', 'S']:
    if not phase_data[phase]:
        continue
    sums_sq = [sum(e**2 for e in d['eig3']) for d in phase_data[phase]]
    sums_4th = [sum(e**4 for e in d['eig3']) for d in phase_data[phase]]
    print(f"\n{phase}:")
    print(f"  a^2+b^2+c^2: {sorted(set(round(s,2) for s in sums_sq))}")
    print(f"  a^4+b^4+c^4: {sorted(set(round(s,2) for s in sums_4th))[:10]}")

# KEY INSIGHT: a^2 + b^2 + c^2 = 15 is FIXED for all n=6 tournaments.
# So the eigenvalue triple lies on the sphere a^2+b^2+c^2 = 15 with 0 < a <= b <= c.
# The product abc = |Pf| is maximized when a=b=c=sqrt(5) (AM-GM), giving |Pf| = 5*sqrt(5) ≈ 11.18
# But |Pf| must be integer (Pfaffian of integer matrix), so max |Pf| = 9.
# The minimum |Pf| with all eigs > 0 approaches 0 as the smallest eig approaches 0.

print("\n" + "=" * 70)
print("ALGEBRAIC EXPLANATION")
print("=" * 70)
print("""
For n=6 tournament with skew eigenvalues +-ia, +-ib, +-ic (a<=b<=c):
  CONSTRAINT: a^2 + b^2 + c^2 = C(n,2) = 15 (from tr(S^2) = -2*C(6,2) = -30)

  |Pf(S)| = a*b*c (integer, from Pfaffian of integer matrix)
  spectral gap = c - a

  AM-GM: abc <= (a^2+b^2+c^2)/3)^{3/2}/sqrt(27) = 5^{3/2}/sqrt(27) = 11.18/5.20 ≈ 2.15
  Wait, AM-GM gives abc <= ((a^2+b^2+c^2)/3)^{3/2} = 5*sqrt(5) ≈ 11.18 when a=b=c=sqrt(5)
  But Pfaffian must be ODD integer. Achievable: 1,3,5,7,9.

  C-phase (beta_1>0): |Pf| in {1,3} => abc small => ONE eigenvalue near 0 or very spread
  S-phase (beta_3>0): |Pf| in {7,9} => abc large => eigenvalues near-degenerate

  The topology "sees" the eigenvalue distribution shape on the constraint sphere!
""")
