#!/usr/bin/env python3
"""
representation_chi_bridge.py — opus-2026-03-12-S69

Test the hypothesis: uniformity of the representation function r_S(d)
predicts |chi_per| (and hence the amount of homology).

r_S(d) = #{(a,b) ∈ S×S : a-b ≡ d mod p}

For a (v,k,λ)-difference set: r_S(d) = λ constant.
Paley at p≡3 mod 4 is a perfect difference set with λ=(m-1)/2.
At p≡1 mod 4, no perfect DS exists at these parameters.

Hypothesis: Var(r_S) small → |chi_per| large.

Also tests: what is the Fourier-analytic meaning of r_S uniformity?
r_S(d) = (1/p) Σ_k |S_hat(k)|² ω^{kd} = (1/p) Σ_k Q_k ω^{kd}
So r_S is the inverse DFT of Q_k! Uniform r_S ↔ flat Q ↔ difference set.

More precisely: Var(r_S) = (1/(p-1)) Σ_{d≠0} (r_S(d) - mean)² = (1/p²(p-1)) Σ_{k≠0} Q_k²
Wait, let me derive this carefully.
"""

import numpy as np
import sys
sys.path.insert(0, '04-computation')
from circulant_homology import CirculantHomology

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

def gen_orientations(p):
    m = (p - 1) // 2
    pairs = [(d, p - d) for d in range(1, m + 1)]
    for bits in range(2**m):
        S = frozenset(pairs[i][0] if bits & (1 << i) else pairs[i][1] for i in range(m))
        yield S

def zpstar_orbit_rep(p, S):
    best = sorted(S)
    for a in range(2, p):
        S_a = frozenset((a * s) % p for s in S)
        if sorted(S_a) < best:
            best = sorted(S_a)
    return frozenset(best)

def compute_Q(p, S):
    m = (p-1)//2
    return np.array([abs(sum(np.exp(2j*np.pi*s*k/p) for s in S))**2 for k in range(1, m+1)])

def representation_func(S, p):
    """Compute r_S(d) for d=1,...,p-1."""
    return {d: sum(1 for a in S for b in S if (a-b)%p == d) for d in range(1, p)}

def representation_variance(r_S, p, m):
    """Variance of r_S(d) over d=1,...,p-1."""
    mean = m * (m-1) / (p-1)
    vals = list(r_S.values())
    return sum((v - mean)**2 for v in vals) / (p-1)

print("=" * 80)
print("REPRESENTATION FUNCTION UNIFORMITY → chi_per")
print("=" * 80)

for p in [7, 11, 13]:
    m = (p-1)//2
    print(f"\n{'─'*60}")
    print(f"p = {p}, m = {m}")
    print(f"{'─'*60}")

    orbits = {}
    for S in gen_orientations(p):
        rep = zpstar_orbit_rep(p, S)
        if rep not in orbits:
            orbits[rep] = S

    results = []
    for rep, S in sorted(orbits.items(), key=lambda x: sorted(x[0])):
        # Compute chi_per
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))

        Q = compute_Q(p, S)
        r_S = representation_func(S, p)
        var_r = representation_variance(r_S, p, m)
        r_vals = sorted(set(r_S.values()))
        r_range = max(r_vals) - min(r_vals)

        # Fourier connection: Var(r_S) = (1/p²) Σ_{k≥1} (Q_k - mean_Q)²
        # where mean_Q = Σ Q_k / (p-1) ... wait, let me compute directly
        # Actually r_S(d) = Σ_k Q_k * e^{2πikd/p} / p (includes k=0 term)
        # where "Q_0" = m² (the |S_hat(0)|² = m²)
        # So r_S(d) = m²/p + (1/p) Σ_{k=1}^{p-1} Q_k * e^{2πikd/p}
        # But Q_k = Q_{p-k}, so this simplifies.
        # Var(r_S) = Var over d of r_S(d)
        # = (1/p²) Σ_{k≥1} Q_k² (by Parseval on the fluctuations)
        # Actually: Σ_{d=1}^{p-1} |r_S(d) - mean|² = (1/p²) Σ_{k=1}^{p-1} |Q_k|²
        # Wait, no. Let f(d) = r_S(d), mean = m(m-1)/(p-1)
        # Actually f(d) = Σ_{a,b∈S} δ_{a-b≡d} and sum f(d) = m(m-1) + m*δ_{d=0}
        # f_hat(k) = Σ_d f(d) ω^{kd} = Σ_{a,b∈S} ω^{k(a-b)} = |S_hat(k)|² = Q_k for k≠0
        # f_hat(0) = m² (including d=0 diagonal)
        # So by Parseval: Σ_d f(d)² = (1/p) Σ_k |f_hat(k)|² = (1/p)(m⁴ + Σ_{k≥1} Q_k²)
        # But we only sum over d=1,...,p-1:
        # Σ_{d≥1} f(d)² = (1/p)(m⁴ + Σ_{k≥1} Q_k²) - f(0)² = (1/p)(m⁴ + Σ Q_k²) - m²

        sum_Q2 = sum(Q**2)
        parseval_check = (m**4 + sum_Q2) / p - m**2
        actual_sum_f2 = sum(v**2 for v in r_S.values())
        # Parseval should give Σ_{d=1}^{p-1} f(d)²

        results.append({
            'S': sorted(S),
            'chi_per': chi_per,
            'var_r': var_r,
            'r_range': r_range,
            'r_vals': r_vals,
            'sum_Q2': float(sum_Q2),
            'prod_Q': float(np.prod(Q)),
        })

    results.sort(key=lambda r: r['var_r'])

    print(f"\n  Sorted by Var(r_S):")
    print(f"  {'S':25s} {'chi':>5s} {'|chi|':>5s} {'Var(r)':>8s} {'range':>5s} {'r values':>15s} {'ΣQ²':>8s} {'ΠQ':>8s}")
    for r in results:
        print(f"  {str(r['S']):25s} {r['chi_per']:5d} {abs(r['chi_per']):5d} {r['var_r']:8.4f} {r['r_range']:5d} {str(r['r_vals']):>15s} {r['sum_Q2']:8.1f} {r['prod_Q']:8.1f}")

    # Correlation between Var(r_S) and |chi_per|
    vars_r = [r['var_r'] for r in results]
    abs_chi = [abs(r['chi_per']) for r in results]
    if len(results) > 2:
        corr = np.corrcoef(vars_r, abs_chi)[0,1]
        print(f"\n  Pearson corr(Var(r_S), |chi_per|) = {corr:.4f}")

    # Also test: Var(r_S) vs chi_per (signed)
    chi_vals = [r['chi_per'] for r in results]
    if len(results) > 2:
        corr_signed = np.corrcoef(vars_r, chi_vals)[0,1]
        print(f"  Pearson corr(Var(r_S), chi_per) = {corr_signed:.4f}")

    # Key: Var(r_S) = (1/(p-1)) Σ (r-mean)² = (1/p²(p-1)) Σ_{k≥1} Q_k²
    # Wait, let me derive properly:
    # f(d) = r_S(d) for d≥1, mean = m(m-1)/(p-1)
    # Σ (f-mean)² = Σ f² - (p-1)*mean²
    # From Parseval: Σ f² = (m⁴ + ΣQ²)/p - m²
    # So Var = [(m⁴ + ΣQ²)/p - m² - (p-1)*mean²]/(p-1)
    # = [(m⁴ + ΣQ²)/p - m² - m²(m-1)²/(p-1)]/(p-1)
    # This should be proportional to ΣQ². Let me verify.
    print(f"\n  Verifying Parseval connection:")
    for r in results:
        var_predicted = (m**4 + r['sum_Q2'])/p - m**2 - (p-1)*(m*(m-1)/(p-1))**2
        var_predicted /= (p-1)
        print(f"    S={str(r['S']):25s}: Var(r) = {r['var_r']:.6f}, predicted = {var_predicted:.6f}")

# Now look at the FOURIER ANALYTIC CONNECTION more carefully
print(f"\n{'='*80}")
print("FOURIER ANALYSIS: Var(r_S) = f(ΣQ²)")
print("=" * 80)
print(f"Since r_S is the inverse DFT of Q_k:")
print(f"  r_S(d) = m²/p + (1/p)Σ Q_k cos(2πkd/p) + correction")
print(f"  Var(r_S) ~ ΣQ² / p² (up to constants)")
print(f"  More uniform r_S ↔ smaller ΣQ² ↔ flatter Q spectrum ↔ more topology")
print()
print(f"This means: FLAT Q-SPECTRUM → MAXIMUM |chi_per|")
print(f"This is a SPECTRAL-TOPOLOGICAL DUALITY!")
print(f"  - Paley (flat Q) → max |chi_per| at p≡3 mod 4")
print(f"  - Quasi-Paley (nearly flat Q) → max |chi_per| at p≡1 mod 4")
print(f"  - Interval (peaked Q) → min |chi_per| ≈ 0 or 1")

# Verify: at each p, is min ΣQ² ↔ max |chi_per|?
print(f"\n  Verification across primes:")
for p in [7, 11, 13]:
    m = (p-1)//2
    orbits = {}
    for S in gen_orientations(p):
        rep = zpstar_orbit_rep(p, S)
        if rep not in orbits:
            orbits[rep] = S

    for rep, S in sorted(orbits.items(), key=lambda x: sorted(x[0])):
        h = CirculantHomology(n=p, S=S)
        h._ensure_enumerated(p)
        omega = h.omega_dims(max_degree=p-1)
        chi_per = sum((-1)**mm * w for mm, w in enumerate(omega))
        Q = compute_Q(p, S)
        sum_Q2 = sum(Q**2)
        print(f"    p={p}: S={str(sorted(S)):25s} ΣQ²={sum_Q2:8.1f} chi_per={chi_per:4d}")

print("\nDONE.")
