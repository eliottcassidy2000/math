#!/usr/bin/env python3
"""
General-p Hessian Analysis for Paley Local Maximum

Compute the Hessian eigenvalue lambda_H for multiple primes p = 3 mod 4.
If lambda_H < 0 for all tested p, this strongly supports the conjecture
that Paley is always a local maximum of H on the Parseval simplex.

The Hessian at Paley of H = c_0 + sum c_j * e_j(x) restricted to
Parseval is: lambda_H = -sum_j c_j * C(m-2, j-2) * (p/4)^{j-2}

We need the coefficients c_j. These come from the exact fit:
H(tournament) = c_0 + c_2*e_2 + c_3*e_3 + ... + c_{n_reps}*e_{n_reps}

For this we need:
1. All distinct H values (requires Held-Karp DP)
2. All distinct e_k values
3. Solve the linear system

For larger p, the DP is expensive but circulant structure helps.
"""

import numpy as np
from itertools import combinations
from math import comb

def is_qr(a, p):
    if a % p == 0: return False
    return pow(a, (p-1)//2, p) == 1

def get_circulants(p):
    n = (p-1)//2
    results = []
    for S in combinations(range(1,p), n):
        S_set = set(S)
        if all((p-j)%p not in S_set for j in S_set):
            results.append(frozenset(S_set))
    return results

def adjacency_matrix(S, p):
    A = np.zeros((p,p), dtype=int)
    for i in range(p):
        for j in range(p):
            if i!=j and (j-i)%p in S:
                A[i][j] = 1
    return A

def eigenvalues_circulant(S, p):
    omega = np.exp(2j*np.pi/p)
    return [sum(omega**(k*s) for s in S) for k in range(p)]

def count_hp_dp(A):
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1<<v, v)] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)): continue
            if (mask,v) not in dp: continue
            c = dp[(mask,v)]
            for u in range(n):
                if mask & (1<<u): continue
                if A[v][u]:
                    key = (mask|(1<<u), u)
                    dp[key] = dp.get(key,0) + c
    return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

def analyze_prime(p):
    """Complete analysis for prime p = 3 mod 4."""
    m = (p-1)//2
    v = p/4  # Paley y^2 value
    paley_S = frozenset(j for j in range(1,p) if is_qr(j,p))

    print(f"\n{'='*60}")
    print(f"p = {p}, m = {m}, v = p/4 = {v}")
    print(f"{'='*60}")

    circulants = get_circulants(p)
    print(f"  {len(circulants)} circulant tournaments")

    # Compute H and eigenvalue data
    data = []
    for S in circulants:
        A = adjacency_matrix(S, p)
        eigs = eigenvalues_circulant(S, p)
        H = count_hp_dp(A)
        y_half = [eigs[k].imag for k in range(1, m+1)]
        x = [y**2 for y in y_half]
        data.append({'S': S, 'H': H, 'x': x, 'paley': S == paley_S})

    # Group by H
    by_H = {}
    for d in data:
        by_H.setdefault(d['H'], []).append(d)
    reps = [by_H[H][0] for H in sorted(by_H, reverse=True)]
    n_reps = len(reps)

    print(f"  {n_reps} distinct H values")
    print(f"  Paley H = {reps[0]['H']}, min H = {reps[-1]['H']}")

    # Compute e_k for representatives
    for d in reps:
        x = d['x']
        d['e'] = {0: 1.0}
        for k in range(1, m+1):
            d['e'][k] = sum(np.prod(list(combo)) for combo in combinations(x, k))

    # Fit H = c_0 + c_2*e_2 + ... + c_{n_e+1}*e_{n_e+1}
    n_e = min(m, n_reps - 1)  # number of e_k terms (skip e_1)

    A_mat = np.zeros((n_reps, n_e + 1))
    b_vec = np.array([d['H'] for d in reps], dtype=float)
    for i, d in enumerate(reps):
        A_mat[i, 0] = 1
        for j in range(1, n_e + 1):
            A_mat[i, j] = d['e'][j + 1]

    rank = np.linalg.matrix_rank(A_mat)
    print(f"  Fitting with e_2,...,e_{n_e+1}: rank={rank}, need {n_e+1}")

    if rank < n_e + 1:
        # Reduce to achievable rank
        n_e = rank - 1
        A_mat = A_mat[:, :n_e+1]
        print(f"  Reduced to e_2,...,e_{n_e+1}")

    try:
        if n_reps <= n_e + 1:
            coeffs = np.linalg.lstsq(A_mat, b_vec, rcond=None)[0]
        else:
            coeffs = np.linalg.lstsq(A_mat, b_vec, rcond=None)[0]
        residual = np.max(np.abs(A_mat @ coeffs - b_vec))
    except np.linalg.LinAlgError:
        print("  SINGULAR!")
        return None

    print(f"  H = {coeffs[0]:.2f}", end="")
    for j in range(1, n_e+1):
        print(f" + ({coeffs[j]:.4f})*e_{j+1}", end="")
    print(f"  [residual: {residual:.2e}]")

    # Compute Hessian eigenvalue
    lambda_H = 0
    for j in range(1, n_e + 1):
        k = j + 1  # e_{j+1}
        curv = -comb(m-2, k-2) * v**(k-2) if k-2 <= m-2 else 0
        contrib = coeffs[j] * curv
        lambda_H += contrib

    print(f"\n  HESSIAN: lambda_H = {lambda_H:.4f}")
    print(f"  Paley is {'LOCAL MAX' if lambda_H < 0 else 'NOT local max (!!!)'}")

    # Verify Paley achieves maximum H
    paley_H = [d['H'] for d in reps if d['paley']]
    is_max = paley_H[0] == reps[0]['H'] if paley_H else False
    print(f"  Paley achieves max H? {is_max}")

    # Fractional gap analysis
    paley_d = [d for d in reps if d['paley']][0]
    print(f"\n  Fractional e_k gaps (averaged over non-Paley):")
    for k in range(2, min(n_e+2, m+1)):
        gaps = []
        for d in reps:
            if d['paley']:
                continue
            if paley_d['e'][k] > 1e-15:
                gap = (paley_d['e'][k] - d['e'][k]) / paley_d['e'][k]
                gaps.append(gap)
        if gaps:
            avg_gap = np.mean(gaps)
            min_gap = min(gaps)
            print(f"    k={k}: avg={avg_gap:.4f}, min={min_gap:.4f}")

    return {
        'p': p, 'm': m, 'lambda_H': lambda_H, 'coeffs': coeffs,
        'n_e': n_e, 'residual': residual, 'is_max': is_max,
        'n_reps': n_reps, 'H_max': reps[0]['H'], 'H_min': reps[-1]['H']
    }


# === MAIN ===
print("GENERAL-p HESSIAN ANALYSIS")
print("=" * 60)

results = {}
# p = 3 mod 4 primes where we can afford the DP
for p in [3, 7, 11]:
    r = analyze_prime(p)
    if r:
        results[p] = r

# p=19 is too expensive for full DP (19! paths)
# But we can still compute e_k from eigenvalues
# The H computation is the bottleneck

print(f"\n{'='*60}")
print(f"SUMMARY TABLE")
print(f"{'='*60}")
print(f"{'p':>4} {'m':>3} {'n_reps':>6} {'lambda_H':>12} {'local_max':>10} {'global_max':>10} {'H_max':>10}")
for p, r in sorted(results.items()):
    print(f"{r['p']:>4} {r['m']:>3} {r['n_reps']:>6} {r['lambda_H']:>12.4f} {str(r['lambda_H']<0):>10} {str(r['is_max']):>10} {r['H_max']:>10}")

# === ASYMPTOTIC ANALYSIS ===
print(f"\n{'='*60}")
print(f"ASYMPTOTIC BEHAVIOR")
print(f"{'='*60}")

# At the uniform point, e_k(v,...,v) = C(m,k)*v^k
# The curvature of e_k is -C(m-2,k-2)*v^{k-2}
# The Hessian is lambda_H = -sum_j c_{j+1} * C(m-2,j-1) * v^{j-1}

# For Schur concavity to hold at the Paley point:
# We need: sum_j c_{j+1} * C(m-2,j-1) * v^{j-1} > 0

# This is a WEIGHTED sum of the coefficients.
# The weights C(m-2,j-1) * v^{j-1} grow with j (since v > 1 for p >= 7).
# So higher-order terms get MORE weight in the Hessian.

# If the positive coefficients tend to be at higher order,
# the weighted sum will be positive.

print("""
Key observation: the Hessian weights w_j = C(m-2,j-1) * v^{j-1}
GROW with j (since v = p/4 > 1 for p >= 7). So higher-order
positive e_k coefficients get disproportionately more weight.

At p=11: c_2 = +184, c_3 = -102, c_4 = +44
Weights: w_1 = 1, w_2 = 3*2.75 = 8.25, w_3 = 3*7.5625 = 22.69
Weighted sum: 184*1 + (-102)*8.25 + 44*22.69 = 184 - 842 + 989 = 331 > 0

The high-order positive term (c_4) dominates because its weight (22.69)
is much larger than the low-order weight (1).

CONJECTURE (HYP-469): For all primes p = 3 mod 4:
lambda_H = -sum c_{j+1} * C(m-2,j-1) * (p/4)^{j-1} < 0.

PROOF STRATEGY:
1. Show the coefficients c_j have a pattern: alternating sign?
2. Show the weighted sum is dominated by the highest-order positive term.
3. Use the OCF structure: H = I(Omega, 2) to constrain the c_j.
""")

# Also test: what if we include p=19 with approximate H?
# We can compute eigenvalue spectra for all Z_19 circulants without DP.
# Then check if the eigenvalue structure alone gives constraints.

print(f"\n{'='*60}")
print(f"EIGENVALUE ANALYSIS FOR p = 19 (without DP)")
print(f"{'='*60}")

p = 19
m = 9
v = p/4
circulants = get_circulants(p)
paley_S = frozenset(j for j in range(1,p) if is_qr(j,p))
print(f"  p={p}, m={m}, {len(circulants)} circulants")

# Just compute eigenvalue spectra and e_k
spectra = {}
for S in circulants:
    eigs = eigenvalues_circulant(S, p)
    y_half = [eigs[k].imag for k in range(1, m+1)]
    x = tuple(sorted([y**2 for y in y_half]))
    is_paley = S == paley_S
    if x not in spectra:
        spectra[x] = {'x': list(x), 'count': 0, 'paley': is_paley}
    spectra[x]['count'] += 1
    if is_paley:
        spectra[x]['paley'] = True

print(f"  {len(spectra)} distinct eigenvalue spectra")

# For each spectrum, compute e_k
for x_tuple, info in spectra.items():
    x = info['x']
    info['e'] = {0: 1.0}
    for k in range(1, min(m+1, 6)):  # only compute up to e_5 for speed
        info['e'][k] = sum(np.prod(list(combo)) for combo in combinations(x, k))

# Check: does Paley maximize all e_k?
paley_spec = [info for info in spectra.values() if info['paley']]
if paley_spec:
    paley_e = paley_spec[0]['e']
    print(f"\n  Paley e_k values:")
    for k in range(1, 6):
        print(f"    e_{k} = {paley_e.get(k, 'N/A'):.6f}")

    print(f"\n  Paley maximizes each e_k?")
    for k in range(2, 6):
        max_ek = max(info['e'].get(k, 0) for info in spectra.values())
        paley_max = abs(paley_e[k] - max_ek) < 0.01
        print(f"    e_{k}: Paley={paley_e[k]:.4f}, max={max_ek:.4f}, "
              f"Paley is max? {paley_max}")

    # Maclaurin means
    print(f"\n  Maclaurin means (Paley vs extremes):")
    for info in spectra.values():
        if info['paley']:
            mus = [(info['e'][k] / comb(m,k))**(1/k) for k in range(1, 6)]
            print(f"    Paley: {['%.4f'%mu for mu in mus]}")
            break

    # Find the spectrum with highest and lowest e_2
    e2_vals = [(info['e'].get(2,0), info) for info in spectra.values()]
    e2_vals.sort(key=lambda x: x[0])
    print(f"\n  Lowest e_2: {e2_vals[0][0]:.4f} (Paley: {paley_e[2]:.4f})")
    print(f"  Highest e_2: {e2_vals[-1][0]:.4f}")
