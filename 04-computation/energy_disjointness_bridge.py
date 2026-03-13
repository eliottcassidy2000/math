#!/usr/bin/env python3
"""
energy_disjointness_bridge.py -- kind-pasteur-2026-03-13-S60

DISCOVERY from overlap-gauss bridge:
  Corr(disj3, additive_energy) = 1.000 at p=7,11

This script investigates:
1. Is disj3 a LINEAR function of additive_energy? What's the exact formula?
2. Does the perfect correlation extend to p=13, 17, 19, 23?
3. What algebraic identity connects them?
4. Connection to Fourier L4 norm (energy = sum |S_hat(t)|^4 by Parseval)
5. Can we express the relationship purely in terms of Gauss sums?

Key identity to test:
  additive_energy(S) = sum_t |S_hat(t)|^4
  disj3 = (c3^2 - overlap) / 2  where overlap = sum_{i<j} 1_{V(C_i) cap V(C_j) != empty}

Since disj3 = -c5/2 + const (THM-155 Walsh restatement),
and c5 is related to the 5th moment of eigenvalues,
and energy is the 4th moment of Fourier...
there should be a DIRECT algebraic bridge.
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict


def legendre(a, p):
    a = a % p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s)%p] = 1
    return A


def compute_c3_c5_disj3(A, p):
    """Compute c3, c5, disj3 for a tournament on Z_p."""
    # 3-cycle vertex sets
    c3_sets = []
    for a, b, c in combinations(range(p), 3):
        if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
            c3_sets.append(frozenset([a, b, c]))
    c3_sets = list(set(c3_sets))
    c3 = len(c3_sets)

    # c5 via trace
    A_np = np.array(A, dtype=np.float64)
    A5 = np.linalg.matrix_power(A_np, 5)
    c5 = int(round(np.trace(A5))) // 5

    # disj3 (disjoint 3-cycle pairs)
    disj3 = 0
    for i in range(c3):
        for j in range(i+1, c3):
            if not (c3_sets[i] & c3_sets[j]):
                disj3 += 1

    return c3, c5, disj3, c3_sets


def additive_energy(S, p):
    """Compute additive energy E(S) = |{(a,b,c,d) in S^4 : a+b=c+d mod p}|."""
    S_set = set(S)
    energy = 0
    for a in S:
        for b in S:
            for c in S:
                d = (a + b - c) % p
                if d in S_set:
                    energy += 1
    return energy


def additive_energy_fourier(S, p):
    """Compute additive energy via Fourier: E(S) = (1/p) * sum_t |S_hat(t)|^4."""
    omega = cmath.exp(2j * cmath.pi / p)
    total = 0.0
    for t in range(p):
        val = sum(omega ** (t * s) for s in S)
        total += abs(val)**4
    return total / p  # Should equal integer energy


def fourier_L4(S, p):
    """Compute sum |S_hat(t)|^4 and individual |S_hat(t)|^2."""
    omega = cmath.exp(2j * cmath.pi / p)
    spectrum = []
    for t in range(p):
        val = sum(omega ** (t * s) for s in S)
        spectrum.append(abs(val)**2)
    L4 = sum(x**2 for x in spectrum)
    return L4, spectrum


def eigenvalue_D_sum(p, S):
    """Compute sum_t D_t^2 and sum_t D_t^4 where D_t = Im(lambda_t).

    lambda_t = sum_{s in S} omega^{st}, Re(lambda_t) = -1/2 for t!=0.
    """
    omega = cmath.exp(2j * cmath.pi / p)
    D2_sum = 0.0
    D4_sum = 0.0
    D_vals = []
    for t in range(1, p):
        lam = sum(omega ** (s * t) for s in S)
        D = lam.imag
        D_vals.append(D)
        D2_sum += D**2
        D4_sum += D**4
    return D2_sum, D4_sum, D_vals


def all_orientations_analysis(p):
    """Analyze disj3 vs additive energy across ALL 2^m orientations at prime p."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"ALL-ORIENTATION ANALYSIS at p={p}, m={m}, N={N}")
    print(f"{'='*70}")

    disj3_vals = []
    energy_vals = []
    c5_vals = []
    D2_vals = []
    D4_vals = []
    L4_vals = []

    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        c3, c5, disj3, _ = compute_c3_c5_disj3(A, p)
        E = additive_energy(S, p)
        D2, D4, _ = eigenvalue_D_sum(p, S)
        L4, _ = fourier_L4(S, p)

        disj3_vals.append(disj3)
        energy_vals.append(E)
        c5_vals.append(c5)
        D2_vals.append(D2)
        D4_vals.append(D4)
        L4_vals.append(L4)

    # Correlation analysis
    d3 = np.array(disj3_vals, dtype=float)
    en = np.array(energy_vals, dtype=float)
    c5a = np.array(c5_vals, dtype=float)
    d2a = np.array(D2_vals, dtype=float)
    d4a = np.array(D4_vals, dtype=float)
    l4a = np.array(L4_vals, dtype=float)

    def corr(x, y):
        if np.std(x) < 1e-12 or np.std(y) < 1e-12:
            return float('nan')
        return np.corrcoef(x, y)[0, 1]

    print(f"\n  Correlation matrix:")
    print(f"    Corr(disj3, energy)    = {corr(d3, en):.6f}")
    print(f"    Corr(disj3, c5)        = {corr(d3, c5a):.6f}")
    print(f"    Corr(disj3, sum D^2)   = {corr(d3, d2a):.6f}")
    print(f"    Corr(disj3, sum D^4)   = {corr(d3, d4a):.6f}")
    print(f"    Corr(disj3, L4)        = {corr(d3, l4a):.6f}")
    print(f"    Corr(energy, c5)       = {corr(en, c5a):.6f}")
    print(f"    Corr(energy, sum D^2)  = {corr(en, d2a):.6f}")
    print(f"    Corr(energy, sum D^4)  = {corr(en, d4a):.6f}")
    print(f"    Corr(energy, L4)       = {corr(en, l4a):.6f}")
    print(f"    Corr(c5, sum D^4)      = {corr(c5a, d4a):.6f}")

    # Linear regression: disj3 = a * energy + b
    if np.std(en) > 1e-12:
        slope = np.cov(d3, en)[0, 1] / np.var(en)
        intercept = np.mean(d3) - slope * np.mean(en)
        residual = d3 - (slope * en + intercept)
        max_err = np.max(np.abs(residual))
        print(f"\n  Linear fit: disj3 = {slope:.6f} * energy + {intercept:.2f}")
        print(f"    Max residual = {max_err:.6f}")
        if max_err < 0.01:
            print(f"    *** EXACT LINEAR RELATION ***")

    # Check disj3 = a * c5 + b (we know this from THM-155)
    if np.std(c5a) > 1e-12:
        slope_c5 = np.cov(d3, c5a)[0, 1] / np.var(c5a)
        int_c5 = np.mean(d3) - slope_c5 * np.mean(c5a)
        res_c5 = d3 - (slope_c5 * c5a + int_c5)
        print(f"\n  Linear fit: disj3 = {slope_c5:.6f} * c5 + {int_c5:.2f}")
        print(f"    Max residual = {np.max(np.abs(res_c5)):.6f}")

    # Check energy = a * c5 + b
    if np.std(c5a) > 1e-12 and np.std(en) > 1e-12:
        slope_ec = np.cov(en, c5a)[0, 1] / np.var(c5a)
        int_ec = np.mean(en) - slope_ec * np.mean(c5a)
        res_ec = en - (slope_ec * c5a + int_ec)
        print(f"\n  Linear fit: energy = {slope_ec:.6f} * c5 + {int_ec:.2f}")
        print(f"    Max residual = {np.max(np.abs(res_ec)):.6f}")

    # Check energy = a * sum_D^4 + b
    if np.std(d4a) > 1e-12 and np.std(en) > 1e-12:
        slope_ed = np.cov(en, d4a)[0, 1] / np.var(d4a)
        int_ed = np.mean(en) - slope_ed * np.mean(d4a)
        res_ed = en - (slope_ed * d4a + int_ed)
        print(f"\n  Linear fit: energy = {slope_ed:.6f} * sum_D^4 + {int_ed:.2f}")
        print(f"    Max residual = {np.max(np.abs(res_ed)):.6f}")

    # DEEPER: Check energy = a*D4 + b*D2 + c (quadratic in D2?)
    X = np.column_stack([d4a, d2a, np.ones(N)])
    coeffs, res, _, _ = np.linalg.lstsq(X, en, rcond=None)
    pred = X @ coeffs
    max_err_3 = np.max(np.abs(en - pred))
    print(f"\n  Multivariate: energy = {coeffs[0]:.6f}*D4 + {coeffs[1]:.6f}*D2 + {coeffs[2]:.2f}")
    print(f"    Max residual = {max_err_3:.6f}")

    # Print unique (disj3, energy, c5, D2, D4) tuples
    print(f"\n  Unique orientation data (first 20):")
    data = list(zip(disj3_vals, energy_vals, c5_vals,
                     [round(x, 2) for x in D2_vals],
                     [round(x, 2) for x in D4_vals]))
    seen = set()
    count = 0
    for d in sorted(set(data)):
        if count >= 20:
            break
        mult = data.count(d)
        print(f"    disj3={d[0]:>6d}  energy={d[1]:>6d}  c5={d[2]:>6d}  D2={d[3]:>10.2f}  D4={d[4]:>10.2f}  (×{mult})")
        count += 1

    return disj3_vals, energy_vals, c5_vals, D2_vals, D4_vals


def algebraic_identity_test(p):
    """Test whether energy = f(c5) has a CLOSED FORM depending only on p.

    We know:
      c5 = tr(A^5)/5  where A is circulant with eigenvalues lambda_t = -1/2 + i*D_t
      tr(A^5) = m^5 + 2*sum_t Re((-1/2 + i*D_t)^5)

    Re(z^5) where z = -1/2 + iD:
      = (-1/2)^5 + C(5,2)(-1/2)^3*(-D^2) + C(5,4)(-1/2)*D^4
      = -1/32 + 10*(1/8)*D^2 + 5*(-1/2)*D^4
      = -1/32 - 10/8*D^2 - 5/2*D^4

    Wait, let me redo. (−1/2+iD)^5 by binomial:
    j=0: C(5,0)(−1/2)^5 (iD)^0 = −1/32                    [real, j even]
    j=1: C(5,1)(−1/2)^4 (iD)^1 = 5/16 * iD                [imag]
    j=2: C(5,2)(−1/2)^3 (iD)^2 = 10*(−1/8)*(−D^2) = 10D^2/8 = 5D^2/4  [real, i^2=−1]
    j=3: C(5,3)(−1/2)^2 (iD)^3 = 10*(1/4)*(−iD^3) = −5iD^3/2  [imag]
    j=4: C(5,4)(−1/2)^1 (iD)^4 = 5*(−1/2)*D^4 = −5D^4/2   [real, i^4=1]
    j=5: C(5,5)(−1/2)^0 (iD)^5 = iD^5                     [imag]

    Re(z^5) = −1/32 + 5D^2/4 − 5D^4/2

    And energy(S) = sum_t |S_hat(t)|^4 / p (but we need to be more careful).

    Actually, |S_hat(t)|^2 = |lambda_t|^2 = 1/4 + D_t^2
    So |S_hat(t)|^4 = (1/4 + D_t^2)^2 = 1/16 + D_t^2/2 + D_t^4

    Therefore:
    sum_{t=1}^{p-1} |lambda_t|^4 = (p-1)/16 + sum D_t^2/2 + sum D_t^4

    We know sum D_t^2 = mp/2 (from THM-155 eigenvalue identity).
    So: sum |lambda_t|^4 = (p-1)/16 + mp/4 + sum D_t^4

    And: c5 = (m^5 + 2*sum Re(z^5))/5
         = (m^5 - (p-1)/16 + 5*sum D^2/2 - 5*sum D^4)/5
         = (m^5 - (p-1)/16 + 5mp/4 - 5*sum D^4)/5

    So: 5*c5 = m^5 - (p-1)/16 + 5mp/4 - 5*sum D^4
    => sum D^4 = (m^5 - (p-1)/16 + 5mp/4 - 5*c5) / 5

    And the additive energy of S is related to |S_hat|^4 via Parseval.
    But S_hat(t) != lambda_t in general...

    Actually S_hat(t) = sum_{s in S} omega^{st} = lambda_t for the circulant!
    So: additive_energy = (1/p) * sum_t |S_hat(t)|^4 = (1/p)(m^4 + sum_{t!=0} |lambda_t|^4)

    This gives a DIRECT connection between energy and sum D^4, hence to c5.
    """
    m = (p - 1) // 2

    print(f"\n{'='*70}")
    print(f"ALGEBRAIC IDENTITY TEST at p={p}")
    print(f"{'='*70}")

    # Predicted: sum D_t^4 = (m^5 - (p-1)/16 + 5mp/4 - 5c5) / 5
    # And: energy = (1/p)(m^4 + (p-1)/16 + mp/4 + sum D^4)

    # Let's verify on a few orientations
    N = 1 << m

    for bits in [0, 1, N-1, N//3]:
        if bits >= N:
            continue
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        A_np = np.array(A, dtype=np.float64)

        # Direct c5
        A5 = np.linalg.matrix_power(A_np, 5)
        c5 = int(round(np.trace(A5))) // 5

        # Direct energy
        E = additive_energy(S, p)

        # Eigenvalue computation
        omega = cmath.exp(2j * cmath.pi / p)
        D_vals = []
        for t in range(1, p):
            lam = sum(omega ** (s * t) for s in S)
            D_vals.append(lam.imag)

        sum_D2 = sum(d**2 for d in D_vals)
        sum_D4 = sum(d**4 for d in D_vals)

        # Predicted sum_D4 from c5
        # 5c5 = m^5 + 2*sum Re(z^5) where sum Re = -(p-1)/32 + 5*sum_D2/4 - 5*sum_D4/2
        # Actually 2*sum Re(z^5) = 2*(-(p-1)/32 + 5*sum_D2/4 - 5*sum_D4/2)
        # So: 5c5 = m^5 - (p-1)/16 + 5*sum_D2/2 - 5*sum_D4
        # => sum_D4 = (m^5 - (p-1)/16 + 5*sum_D2/2 - 5*c5) / 5
        pred_D4 = (m**5 - (p-1)/16 + 5*sum_D2/2 - 5*c5) / 5

        # Predicted energy from eigenvalues
        # E*p = |S|^4 + sum_{t!=0} |lambda_t|^4
        #     = m^4 + sum (1/4 + D_t^2)^2
        #     = m^4 + (p-1)/16 + sum_D2/2 + sum_D4
        pred_E = (m**4 + (p-1)/16 + sum_D2/2 + sum_D4) / p

        print(f"\n  bits={bits:0{m}b} S={S}")
        print(f"    c5 = {c5}")
        print(f"    energy = {E}")
        print(f"    sum_D2 = {sum_D2:.4f} (expected mp/2 = {m*p/2})")
        print(f"    sum_D4 = {sum_D4:.4f}")
        print(f"    pred_D4 from c5 = {pred_D4:.4f} (diff = {abs(sum_D4 - pred_D4):.6f})")
        print(f"    pred_E from eigenvals = {pred_E:.4f} (diff = {abs(E - pred_E):.6f})")

        # Now derive exact formula: E = f(c5, p)
        # E*p = m^4 + (p-1)/16 + mp/4 + sum_D4
        # sum_D4 = (m^5 - (p-1)/16 + 5mp/4 - 5c5) / 5
        # So: E*p = m^4 + (p-1)/16 + mp/4 + (m^5 - (p-1)/16 + 5mp/4 - 5c5)/5
        #        = m^4 + (p-1)/16 + mp/4 + m^5/5 - (p-1)/80 + mp/4 - c5
        #        = m^4 + m^5/5 + (p-1)(1/16 - 1/80) + mp/2 - c5
        #        = m^4 + m^5/5 + (p-1)*4/80 + mp/2 - c5
        #        = m^4 + m^5/5 + (p-1)/20 + mp/2 - c5

        pred_Ep = m**4 + m**5/5 + (p-1)/20 + m*p/2 - c5
        print(f"    PREDICTED E*p = m^4 + m^5/5 + (p-1)/20 + mp/2 - c5")
        print(f"    = {pred_Ep:.4f}")
        print(f"    Actual E*p = {E*p}")
        print(f"    Match: {abs(E*p - pred_Ep) < 0.01}")


def derive_closed_form():
    """Derive and display the closed-form identity connecting energy, c5, and disj3."""
    print(f"\n{'='*70}")
    print("CLOSED-FORM IDENTITY DERIVATION")
    print("="*70)

    print("""
  Starting from eigenvalue structure of circulant tournament on Z_p:

  lambda_0 = m = (p-1)/2
  lambda_t = -1/2 + i*D_t  for t = 1, ..., p-1

  Known: sum_{t!=0} D_t^2 = mp/2   (from sum |lambda_t|^2 = m)

  IDENTITY 1: c5 from eigenvalues
    c5 = tr(A^5)/5 = [m^5 + 2*sum Re((-1/2+iD_t)^5)] / 5
    Re(z^5) = -1/32 + 5D^2/4 - 5D^4/2

    5c5 = m^5 - (p-1)/16 + 5mp/4 - 5*sum D^4
    => sum D^4 = (m^5 - (p-1)/16 + 5mp/4 - 5c5) / 5   ...(*)

  IDENTITY 2: energy from eigenvalues
    E(S) = (1/p)*sum_t |S_hat(t)|^4 = (1/p)*[m^4 + sum_{t!=0}(1/4+D_t^2)^2]
    = (1/p)*[m^4 + (p-1)/16 + mp/4 + sum D^4]          ...(**)

  COMBINING (*) and (**):
    E*p = m^4 + (p-1)/16 + mp/4 + [m^5/5 - (p-1)/80 + mp/4 - c5]
    E*p = m^4 + m^5/5 + (p-1)/20 + mp/2 - c5

    E = [m^4 + m^5/5 + (p-1)/20 + mp/2 - c5] / p

  IDENTITY 3: disj3 from c5 (THM-155 Walsh version)
    disj3 = -c5/2 + const(p)

  COMBINING:
    disj3 = -c5/2 + const
    E = const'(p) - c5/p

    => disj3 = p*E/2 + const''(p)

  This is a LINEAR relationship with slope p/2!
""")


def verify_slope(p):
    """Verify that disj3 = (p/2)*energy + const across orientations."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n  Testing disj3 = (p/2)*energy + const at p={p}:")

    pairs = []
    for bits in range(min(N, 64)):  # limit for large p
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        c3, c5, disj3, _ = compute_c3_c5_disj3(A, p)
        E = additive_energy(S, p)

        pairs.append((disj3, E, c5))

    d3 = np.array([x[0] for x in pairs], dtype=float)
    en = np.array([x[1] for x in pairs], dtype=float)
    c5a = np.array([x[2] for x in pairs], dtype=float)

    if np.std(en) > 1e-12:
        slope = np.cov(d3, en)[0, 1] / np.var(en)
        intercept = np.mean(d3) - slope * np.mean(en)
        residual = d3 - (slope * en + intercept)
        max_err = np.max(np.abs(residual))

        # Expected slope from derivation
        # disj3 = -c5/2 + A  and  E = B - c5/p  =>  c5 = p*(B-E)
        # disj3 = -p*(B-E)/2 + A = p*E/2 + (A - pB/2)
        expected_slope = p / 2

        print(f"    Measured slope     = {slope:.6f}")
        print(f"    Expected slope p/2 = {expected_slope:.6f}")
        print(f"    Slope match: {abs(slope - expected_slope) < 0.01}")
        print(f"    Intercept = {intercept:.2f}")
        print(f"    Max residual = {max_err:.6f}")

        # Verify: c5 = const - p*E (linear in E with slope -p)
        if np.std(en) > 1e-12:
            slope_ce = np.cov(c5a, en)[0, 1] / np.var(en)
            print(f"    c5 vs energy slope = {slope_ce:.6f} (expected {-p})")


def higher_moment_investigation(p):
    """Investigate whether higher cycle counts (c7, c9) also relate to energy moments."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"HIGHER MOMENT INVESTIGATION at p={p}")
    print(f"{'='*70}")

    data = []
    for bits in range(min(N, 64)):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        A_np = np.array(A, dtype=np.float64)

        ck_vals = {}
        for k in range(3, min(p+1, 14), 2):
            Ak = np.linalg.matrix_power(A_np, k)
            ck_vals[k] = int(round(np.trace(Ak))) // k

        E = additive_energy(S, p)
        _, D4, D_vals = eigenvalue_D_sum(p, S)
        D6 = sum(d**6 for d in D_vals)

        data.append({**ck_vals, 'energy': E, 'D4': D4, 'D6': D6})

    # For each c_k, check if it's an affine function of energy
    for k in range(3, min(p+1, 14), 2):
        ck_arr = np.array([d[k] for d in data], dtype=float)
        en_arr = np.array([d['energy'] for d in data], dtype=float)

        if np.std(en_arr) < 1e-12 or np.std(ck_arr) < 1e-12:
            print(f"\n  c_{k}: constant (std=0)")
            continue

        r = np.corrcoef(ck_arr, en_arr)[0, 1]
        print(f"\n  c_{k} vs energy: corr = {r:.6f}", end='')

        if abs(r) > 0.999:
            slope = np.cov(ck_arr, en_arr)[0, 1] / np.var(en_arr)
            resid = ck_arr - (slope * en_arr + np.mean(ck_arr) - slope*np.mean(en_arr))
            print(f"  EXACT? slope={slope:.4f}, max_resid={np.max(np.abs(resid)):.4f}")
        else:
            # Check if c_k is function of (energy, D6) or similar
            D4_arr = np.array([d['D4'] for d in data], dtype=float)
            D6_arr = np.array([d['D6'] for d in data], dtype=float)

            X = np.column_stack([en_arr, D4_arr, np.ones(len(data))])
            try:
                coeffs, _, _, _ = np.linalg.lstsq(X, ck_arr, rcond=None)
                pred = X @ coeffs
                max_err = np.max(np.abs(ck_arr - pred))
                print(f"  2-var fit: {coeffs[0]:.4f}*E + {coeffs[1]:.4f}*D4 + {coeffs[2]:.2f}, max_err={max_err:.4f}")
            except:
                print()

            # Try with D6
            X = np.column_stack([en_arr, D6_arr, np.ones(len(data))])
            try:
                coeffs, _, _, _ = np.linalg.lstsq(X, ck_arr, rcond=None)
                pred = X @ coeffs
                max_err = np.max(np.abs(ck_arr - pred))
                print(f"    + D6 fit: {coeffs[0]:.4f}*E + {coeffs[1]:.4f}*D6 + {coeffs[2]:.2f}, max_err={max_err:.4f}")
            except:
                pass


# ========================================================
# MAIN
# ========================================================

print("=" * 70)
print("ENERGY-DISJOINTNESS BRIDGE")
print("Testing: disj3 = (p/2)*energy + const")
print("=" * 70)

derive_closed_form()

for p in [7, 11, 13]:
    all_orientations_analysis(p)

print("\n\n" + "=" * 70)
print("SLOPE VERIFICATION: disj3 = (p/2)*energy + const")
print("=" * 70)

for p in [7, 11, 13, 17, 19]:
    verify_slope(p)

print("\n\n" + "=" * 70)
print("ALGEBRAIC IDENTITY VERIFICATION")
print("=" * 70)

for p in [7, 11, 13]:
    algebraic_identity_test(p)

print("\n")
higher_moment_investigation(7)
higher_moment_investigation(11)

print("\nDONE.")
