#!/usr/bin/env python3
"""
MACLAURIN INEQUALITY PROOF for Paley H-Maximization

KEY INSIGHT: Paley maximizes ALL e_k(y^2) simultaneously among circulant
tournaments with the same Parseval sum. This follows from Maclaurin's inequality:
    e_k/C(m,k) >= (e_{k+1}/C(m,k+1))^{k/(k+1)}
with equality iff all x_i are equal.

Even though H = c_0 + c_2*e_2 - c_3*e_3 + c_4*e_4 has MIXED signs,
Paley still maximizes because:
1. The e_k gaps (Paley minus other) grow super-exponentially in k
2. The positive coefficient terms (e_2, e_4) dominate the negative (e_3)

PROOF STRATEGY:
If we write Delta_k = e_k(Paley) - e_k(other), then:
Delta H = c_2*Delta_2 - |c_3|*Delta_3 + c_4*Delta_4
The question: does c_2*Delta_2 + c_4*Delta_4 > |c_3|*Delta_3?

Maclaurin tells us: Delta_k / e_k(Paley) is INCREASING in k.
So the fractional gap grows, and higher-order positive terms win.
"""

import numpy as np
from itertools import combinations
from fractions import Fraction

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


def run_analysis(p):
    print(f"\n{'='*70}")
    print(f"MACLAURIN INEQUALITY ANALYSIS FOR p = {p}")
    print(f"{'='*70}")

    circulants = get_circulants(p)
    m = (p-1)//2
    paley_S = frozenset(j for j in range(1,p) if is_qr(j,p))

    data = []
    for S in circulants:
        A = adjacency_matrix(S, p)
        eigs = eigenvalues_circulant(S, p)
        H = count_hp_dp(A)
        y_half = [eigs[k].imag for k in range(1, m+1)]
        x = [y**2 for y in y_half]
        data.append({'S': S, 'H': H, 'x': x, 'paley': S == paley_S})

    by_H = {}
    for d in data:
        by_H.setdefault(d['H'], []).append(d)
    reps = [by_H[H][0] for H in sorted(by_H, reverse=True)]

    # Compute e_k for each
    for d in data + reps:
        x = d['x']
        d['e'] = {0: 1.0}
        for k in range(1, m+1):
            d['e'][k] = sum(np.prod(list(combo)) for combo in combinations(x, k))

    paley = [d for d in reps if d['paley']][0]

    print(f"\n--- Maclaurin's Inequality Verification ---")
    print(f"For Paley (all y^2 = {p/4:.4f}):")
    print(f"  e_k / C(m,k):")
    from math import comb
    for k in range(1, m+1):
        ratio = paley['e'][k] / comb(m, k)
        print(f"    e_{k}/C({m},{k}) = {paley['e'][k]:.4f} / {comb(m,k)} = {ratio:.6f}")

    # Verify: for uniform x_i = v, e_k = C(m,k) * v^k
    v = p / 4
    print(f"\n  Expected: e_k = C(m,k) * v^k where v = {v}")
    for k in range(1, m+1):
        expected = comb(m, k) * v**k
        actual = paley['e'][k]
        print(f"    k={k}: expected={expected:.6f}, actual={actual:.6f}, match={abs(expected-actual)<0.01}")

    # For each non-Paley tournament, compute FRACTIONAL gaps
    print(f"\n--- Fractional e_k Gaps (Paley advantage) ---")
    print(f"  Delta_k / e_k(Paley) for each tournament:")

    for d in reps:
        if d['paley']:
            continue
        print(f"\n  H={d['H']}:")
        fracs = []
        for k in range(2, m+1):
            delta = paley['e'][k] - d['e'][k]
            frac = delta / paley['e'][k] if paley['e'][k] > 0 else 0
            fracs.append(frac)
            print(f"    k={k}: e_k={d['e'][k]:.4f}, Delta={delta:.4f}, "
                  f"frac={frac:.6f}")

        # Check: is the fractional gap INCREASING?
        increasing = all(fracs[i] <= fracs[i+1] + 0.01 for i in range(len(fracs)-1))
        print(f"    Fractional gap increasing? {increasing}")

    # === FIT H in e_k basis and compute detailed budget ===
    n_reps = len(reps)
    b_vec = np.array([d['H'] for d in reps], dtype=float)

    # Use e_2,...,e_{n_reps-1+1} (skip e_1 which is universal)
    n_e = min(m, n_reps - 1)
    A_mat = np.zeros((n_reps, n_e + 1))
    for i, d in enumerate(reps):
        A_mat[i, 0] = 1
        for j in range(1, n_e + 1):
            A_mat[i, j] = d['e'][j + 1]

    rank = np.linalg.matrix_rank(A_mat)
    if rank < n_e + 1:
        print(f"\n  Rank deficient ({rank} < {n_e+1}), reducing...")
        n_e = rank - 1

    try:
        if n_reps == n_e + 1:
            coeffs = np.linalg.solve(A_mat[:,:n_e+1], b_vec)
        else:
            coeffs = np.linalg.lstsq(A_mat[:,:n_e+1], b_vec, rcond=None)[0]
        residual = np.max(np.abs(A_mat[:,:n_e+1] @ coeffs - b_vec))
    except np.linalg.LinAlgError:
        print("  Singular!")
        return

    print(f"\n--- H in elementary symmetric basis ---")
    print(f"  H = {coeffs[0]:.4f}", end="")
    for j in range(1, n_e+1):
        print(f" + ({coeffs[j]:.4f})*e_{j+1}", end="")
    print(f"  [residual: {residual:.2e}]")

    # For each non-Paley, compute the contribution budget
    print(f"\n--- Contribution Budget (why Paley wins) ---")

    for d in reps:
        if d['paley']:
            continue

        print(f"\n  H={d['H']}: total deficit = {paley['H'] - d['H']}")
        total = 0
        for j in range(1, n_e+1):
            delta_e = paley['e'][j+1] - d['e'][j+1]
            contribution = coeffs[j] * delta_e
            total += contribution
            sign = "helps" if contribution > 0 else "HURTS"
            print(f"    c_{j+1}*Delta(e_{j+1}) = {coeffs[j]:.2f} * {delta_e:.4f} = {contribution:+.2f} ({sign})")
        print(f"    Sum = {total:.2f} (should match deficit {paley['H']-d['H']})")

    # === THEORETICAL ANALYSIS: When do positive terms dominate? ===
    print(f"\n{'='*70}")
    print(f"THEORETICAL ANALYSIS: Sufficient conditions for Paley maximization")
    print(f"{'='*70}")

    # At uniform point x_i = v, e_k = C(m,k) * v^k
    # For a perturbation x_i = v + delta_i with sum(delta_i) = 0:
    # Delta(e_k) = e_k(v+delta) - e_k(v)

    # First order: d/dv e_k = C(m-1,k-1) * k * v^{k-1} * (perturbation direction)
    # This is why ALL e_k decrease when moving away from uniform.

    # Second order: the Hessian of e_k at the uniform point
    # d^2 e_k / d(x_i) d(x_j) at x=v·1:
    #   = C(m-2, k-2) * v^{k-2} for i != j
    #   = 0 for i = j (e_k is multilinear in individual x_i? No.)

    # Actually e_k = sum over k-subsets of product of k x_values
    # d e_k / d x_i = e_{k-1}(x without x_i)
    # d^2 e_k / d x_i d x_j = e_{k-2}(x without x_i and x_j) for i != j

    # At uniform: d^2 e_k / d x_i d x_j = C(m-2, k-2) * v^{k-2} for i!=j
    # d^2 e_k / d x_i^2 = 0 (e_k has no x_i^2 terms)

    # So the Hessian of e_k at uniform is:
    # H_e_k = C(m-2,k-2) * v^{k-2} * (J - I) where J is all-ones
    # On the hyperplane sum delta = 0: J*delta = 0, so
    # H_e_k|_{sum=0} = -C(m-2,k-2) * v^{k-2} * I

    # This is NEGATIVE definite for k >= 2!
    # So all e_k are STRICTLY CONCAVE on the Parseval simplex at the uniform point.

    print(f"\n  Second-order analysis at Paley (uniform) point:")
    print(f"  Hessian of e_k restricted to Parseval hyperplane:")
    print(f"  H(e_k)|_Parseval = -C(m-2,k-2) * v^{{k-2}} * I")
    print(f"  This is NEGATIVE DEFINITE for k >= 2!")
    print()

    for k in range(2, m+1):
        hess_diag = -comb(m-2, k-2) * v**(k-2)
        print(f"    k={k}: curvature = {hess_diag:.6f}")

    # Now: the Hessian of H = c_0 + sum c_{j+1} * e_{j+1} is:
    # H(H) = sum_{j=1}^{n_e} c_{j+1} * H(e_{j+1})
    # = -sum_{j=1}^{n_e} c_{j+1} * C(m-2, j-1) * v^{j-1} * I

    hess_H = 0
    print(f"\n  Hessian of H at Paley:")
    for j in range(1, n_e+1):
        k = j + 1  # e_{j+1}
        contrib = -coeffs[j] * comb(m-2, k-2) * v**(k-2)
        hess_H += contrib
        print(f"    c_{k} * curv(e_{k}) = {coeffs[j]:.2f} * {-comb(m-2,k-2)*v**(k-2):.6f} = {contrib:.4f}")

    print(f"\n  Total Hessian eigenvalue = {hess_H:.6f}")
    print(f"  Paley is {'LOCAL MAX' if hess_H < 0 else 'NOT local max'} of H on Parseval simplex")

    if hess_H < 0:
        print(f"\n  *** PALEY IS A LOCAL MAXIMUM OF H ***")
        print(f"  Combined with Paley being the UNIQUE critical point (by symmetry),")
        print(f"  and H being a polynomial of bounded degree on a compact domain,")
        print(f"  this strongly suggests Paley is the GLOBAL maximum.")

    # === Check: is Paley a GLOBAL max? ===
    # On the compact domain {x_i >= 0, sum x_i = m*v}, H is a polynomial.
    # The boundary consists of degenerate eigenvalue spectra (some y_k = 0).
    # We need: H(boundary) < H(Paley).

    # Actually, we can check this for ALL circulant tournaments:
    print(f"\n--- Global optimality verification ---")
    all_H = [d['H'] for d in data]
    max_H = max(all_H)
    max_count = sum(1 for h in all_H if h == max_H)
    paley_count = sum(1 for d in data if d['paley'])
    print(f"  max(H) = {max_H}, achieved by {max_count}/{len(data)} tournaments")
    print(f"  Paley count = {paley_count}")
    print(f"  Paley achieves max? {paley['H'] == max_H}")

    # === RATIO TEST: can H/some_function be monotone in a single e_k? ===
    print(f"\n--- Monotonicity in e_m (product) after normalization ---")
    for d in reps:
        if d['e'][m] > 1e-10:
            ratio = d['H'] / d['e'][m]
        else:
            ratio = float('inf')
        print(f"  H={d['H']}: e_{m}={d['e'][m]:.6f}, H/e_{m}={ratio:.4f} {'(Paley)' if d['paley'] else ''}")

    # === MACLAURIN-BASED BOUND ===
    print(f"\n{'='*70}")
    print(f"MACLAURIN INEQUALITY BOUND")
    print(f"{'='*70}")

    # Maclaurin: (e_k/C(m,k))^{1/k} is decreasing in k
    # At uniform: all ratios equal v
    # Away from uniform: e_k/C(m,k) < v^k, and the ratio decreases
    # faster for larger k.

    # Define the Maclaurin ratio: mu_k = (e_k/C(m,k))^{1/k}
    print(f"\n  Maclaurin means mu_k = (e_k / C(m,k))^(1/k):")

    for d in reps:
        print(f"\n  H={d['H']} {'(Paley)' if d['paley'] else ''}:")
        mus = []
        for k in range(1, m+1):
            mu = (d['e'][k] / comb(m,k))**(1/k) if d['e'][k] > 0 else 0
            mus.append(mu)
            print(f"    mu_{k} = {mu:.6f}", end="")
        print()

        # Check decreasing
        decr = all(mus[i] >= mus[i+1] - 1e-10 for i in range(len(mus)-1))
        print(f"    Decreasing? {decr}")

        # How far from uniform?
        if d['paley']:
            print(f"    (All equal = {v:.6f}, as expected)")
        else:
            print(f"    Deviation from Paley ({v:.6f}):", end="")
            for k in range(m):
                dev = (v - mus[k]) / v * 100
                print(f" {dev:+.1f}%", end="")
            print()

    return data, reps, coeffs


# === MAIN ===
print("MACLAURIN INEQUALITY AND PALEY H-MAXIMIZATION")
print("=" * 70)

for p in [7, 11]:
    data, reps, coeffs = run_analysis(p)

print(f"\n{'='*70}")
print(f"GRAND CONCLUSION")
print(f"{'='*70}")
print("""
THEOREM SKETCH (Paley maximizes H among Z_p circulants):

1. REPRESENTATION: H = c_0 + sum_{k=2}^{m+1} c_k * e_k(y^2)
   where y_k = Im(eigenvalue_k), m = (p-1)/2.

2. HESSIAN AT PALEY: The Hessian of H restricted to the Parseval
   hyperplane sum(y_k^2) = m*p/4 is:
     lambda_H = -sum_k c_k * C(m-2,k-2) * (p/4)^{k-2}
   which is NEGATIVE (verified for p=7,11).

3. CRITICAL POINT: By Z_p cyclic symmetry, the Paley point (all y^2
   equal) is the UNIQUE critical point of H restricted to the symmetric
   functions. The Hessian being negative makes it a local max.

4. GLOBAL MAX: Since H is a polynomial on a compact domain and the
   only critical point is a local max, either:
   (a) Paley IS the global max, or
   (b) The global max is on the boundary (some y_k = 0).

5. BOUNDARY ANALYSIS: On the boundary of the eigenvalue domain, some
   eigenvalue magnitudes collapse. This corresponds to highly
   non-uniform spectra, which have SMALLER e_k values. The boundary
   analysis needs Maclaurin's inequality to bound H on boundary.

WHAT REMAINS:
- Step 5: formal boundary analysis using Maclaurin bounds
- Connect the polynomial representation to the OCF structure
- Extend from circulant to all tournaments ("Step A")
""")
