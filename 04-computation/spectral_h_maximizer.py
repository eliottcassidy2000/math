# ⚠️ WARNING: The connection set {1,2,3,5,8} does NOT give a tournament on Z_11
# because S ∩ (-S mod 11) = {3,8} ≠ ∅. All results for this "DRT" are INVALID.
# The ONLY valid circulant DRT at n=11 is the Paley tournament (QR={1,3,4,5,9}).
# See MISTAKE-017.

import os; os.environ['PYTHONIOENCODING'] = 'utf-8'
"""
spectral_h_maximizer.py
kind-pasteur-2026-03-07-S39b

Deep spectral analysis of H-maximizing tournaments.

INV-055 lead: Linial-Morgenstern uses spectral methods to bound cycle densities.
Can spectral properties EXPLAIN why Paley maximizes H?

Key questions:
1. Does the skew spectrum (eigenvalues of A-A^T) characterize Paley among
   regular tournaments?
2. Is there a spectral formula for alpha_1 (total directed odd cycles)?
3. Does the "flatness" of the spectrum correlate with H-maximization?
4. At n=11 (two DRT classes): what spectral property distinguishes them?

Using the fact that for tournament A:
  tr(A^k) = k*c_k for k=3,4,5 (THM-118)
  tr((A-A^T)^k) = tr(S^k) where S is the skew-adjacency matrix
  For odd k: tr(S^k) = 2*tr(A^k) - ... (relates to cycle counts)
"""

import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from tournament_fast import c3_from_score, c4_fast, c5_fast, alpha2_from_trace
from math import comb
import random

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    print("numpy not available. Exiting.")
    sys.exit(0)


def paley_tournament(p):
    """Construct Paley tournament T_p for prime p = 3 mod 4."""
    QR = set()
    for x in range(1, p):
        QR.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in QR:
                T[i][j] = 1
    return T


def skew_spectrum(T):
    """Compute eigenvalues of S = A - A^T (purely imaginary for tournaments)."""
    A = np.array(T, dtype=float)
    S = A - A.T
    evals = np.linalg.eigvals(S)
    # Should be purely imaginary (real parts ~ 0)
    imag = sorted([e.imag for e in evals], reverse=True)
    return imag


def adjacency_spectrum(T):
    """Eigenvalues of the adjacency matrix A."""
    A = np.array(T, dtype=float)
    return sorted(np.linalg.eigvals(A), key=lambda x: -abs(x))


print("=" * 70)
print("SPECTRAL H-MAXIMIZER ANALYSIS")
print("=" * 70)

# ============================================================
# 1. Paley vs Regular tournaments at n=7
# ============================================================
print("\n--- n=7: Paley vs all regular tournaments ---")

n = 7
m = n*(n-1)//2
T7_paley = paley_tournament(7)
h_paley = hamiltonian_path_count(T7_paley)

# Skew spectrum of Paley
sk_paley = skew_spectrum(T7_paley)
print(f"Paley T_7: H={h_paley}")
print(f"  Skew spectrum (imag): {[f'{x:.4f}' for x in sk_paley]}")
print(f"  |Skew eigenvalues|: {sorted([abs(x) for x in sk_paley], reverse=True)[:4]}")

# Regular tournaments at n=7 have score (3,3,3,3,3,3,3)
random.seed(42)
all_bits = range(1 << m)

# Sample and find regular tournaments
regular_data = []
random_sample = [random.randint(0, (1 << m) - 1) for _ in range(5000)]

for bits in random_sample:
    T = tournament_from_bits(n, bits)
    scores = sorted(sum(row) for row in T)
    if scores == [3]*7:  # Regular
        h = hamiltonian_path_count(T)
        sk = skew_spectrum(T)
        c3 = c3_from_score(T)
        c5 = c5_fast(T)
        a2 = alpha2_from_trace(T)
        regular_data.append({
            'bits': bits, 'H': h, 'skew': sk,
            'c3': c3, 'c5': c5, 'alpha2': a2,
            'skew_flat': max(abs(x) for x in sk) - min(abs(x) for x in sk if abs(x) > 0.01)
        })

print(f"\nFound {len(regular_data)} regular tournaments in sample")

if regular_data:
    # Group by H value (should be 3 classes: H=189, 175, 171)
    h_groups = {}
    for d in regular_data:
        h_groups.setdefault(d['H'], []).append(d)

    print(f"\nRegular tournament classes by H:")
    for h_val in sorted(h_groups.keys(), reverse=True):
        group = h_groups[h_val]
        sk_example = group[0]['skew']
        sk_max = max(abs(x) for x in sk_example)
        sk_nonzero = sorted([abs(x) for x in sk_example if abs(x) > 0.01], reverse=True)
        print(f"  H={h_val}: count={len(group)}, c3={group[0]['c3']}, c5={group[0]['c5']}, "
              f"alpha2={group[0]['alpha2']}")
        print(f"    |Skew eigenvalues|: {[f'{x:.4f}' for x in sk_nonzero]}")

    # Check: do all tours with same H have same skew spectrum?
    print(f"\n  Spectral invariance within H-class:")
    for h_val in sorted(h_groups.keys(), reverse=True):
        group = h_groups[h_val]
        spectra = set()
        for d in group:
            sk_rounded = tuple(round(abs(x), 3) for x in d['skew'] if abs(x) > 0.01)
            spectra.add(sk_rounded)
        print(f"    H={h_val}: {len(spectra)} distinct spectra among {len(group)} tournaments")


# ============================================================
# 2. Spectral formula for cycle counts
# ============================================================
print("\n" + "=" * 70)
print("SPECTRAL FORMULAS FOR CYCLE COUNTS")
print("=" * 70)

# For skew matrix S = A - A^T:
# S^k for odd k: tr(S^k) relates to cycle counts
# S = A - A^T, so S^3 = A^3 - A^2*A^T - A*A^T*A - A^T*A^2 + A*A^T*A^T + A^T*A*A^T + A^T*A^2 - A^T^3
# Actually simpler: for tournament, A + A^T = J - I
# So A^T = J - I - A
# And S = A - (J-I-A) = 2A - J + I

# tr(S^k) for odd k:
# At k=1: tr(S) = tr(2A-J+I) = 0 - n + n = 0
# At k=3: tr(S^3) = tr((2A-J+I)^3)

# Let's verify: tr(S^3) = ? * c3
print("\nChecking tr(S^k) vs cycle counts:")

random.seed(42)
sample = [random.randint(0, (1 << m) - 1) for _ in range(200)]

s3_formula = []
s5_formula = []
for bits in sample[:100]:
    T = tournament_from_bits(n, bits)
    A = np.array(T, dtype=float)
    S = A - A.T
    S3 = S @ S @ S
    S5 = S3 @ S @ S
    tr_S3 = np.trace(S3)
    tr_S5 = np.trace(S5)

    c3 = c3_from_score(T)
    c5 = c5_fast(T)

    # tr(A^3) = 3*c3, tr(A^5) = 5*c5
    tr_A3 = 3 * c3
    tr_A5 = 5 * c5

    s3_formula.append((tr_S3, c3, tr_A3))
    s5_formula.append((tr_S5, c5, tr_A5))

# Check if tr(S^3) / c3 is constant
ratios_3 = [x[0]/x[1] if x[1] != 0 else None for x in s3_formula]
ratios_3 = [r for r in ratios_3 if r is not None]
if ratios_3:
    r3_min, r3_max = min(ratios_3), max(ratios_3)
    print(f"  tr(S^3)/c3: min={r3_min:.4f}, max={r3_max:.4f} "
          f"({'CONSTANT' if abs(r3_max - r3_min) < 0.001 else 'VARIES'})")

# Check tr(S^3) vs tr(A^3)
ratios_SA3 = [x[0]/x[2] if x[2] != 0 else None for x in s3_formula]
ratios_SA3 = [r for r in ratios_SA3 if r is not None]
if ratios_SA3:
    r_min, r_max = min(ratios_SA3), max(ratios_SA3)
    print(f"  tr(S^3)/tr(A^3): min={r_min:.4f}, max={r_max:.4f} "
          f"({'CONSTANT' if abs(r_max - r_min) < 0.001 else 'VARIES'})")
    if abs(r_max - r_min) < 0.001:
        print(f"    => tr(S^3) = {r_min:.0f} * tr(A^3) = {r_min:.0f} * 3 * c3")

# Same for k=5
ratios_5 = [x[0]/x[1] if x[1] != 0 else None for x in s5_formula]
ratios_5 = [r for r in ratios_5 if r is not None]
if ratios_5:
    r5_min, r5_max = min(ratios_5), max(ratios_5)
    print(f"  tr(S^5)/c5: min={r5_min:.4f}, max={r5_max:.4f} "
          f"({'CONSTANT' if abs(r5_max - r5_min) < 0.001 else 'VARIES'})")

ratios_SA5 = [x[0]/x[2] if x[2] != 0 else None for x in s5_formula]
ratios_SA5 = [r for r in ratios_SA5 if r is not None]
if ratios_SA5:
    r_min, r_max = min(ratios_SA5), max(ratios_SA5)
    print(f"  tr(S^5)/tr(A^5): min={r_min:.4f}, max={r_max:.4f} "
          f"({'CONSTANT' if abs(r_max - r_min) < 0.001 else 'VARIES'})")

# ============================================================
# 3. Skew spectrum flatness and H
# ============================================================
print("\n" + "=" * 70)
print("SPECTRAL FLATNESS vs H")
print("=" * 70)

data_spec = []
for bits in sample:
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    sk = skew_spectrum(T)
    sk_abs = sorted([abs(x) for x in sk], reverse=True)

    # Various spectral measures
    spec_max = sk_abs[0] if sk_abs else 0
    spec_gap = sk_abs[0] - sk_abs[1] if len(sk_abs) >= 2 else 0
    spec_var = np.var(sk_abs)
    # "Flatness" = 1 - normalized variance
    spec_mean = np.mean(sk_abs) if sk_abs else 0
    flatness = 1 - (np.std(sk_abs) / spec_mean) if spec_mean > 0 else 0

    data_spec.append({
        'H': h, 'spec_max': spec_max, 'spec_gap': spec_gap,
        'spec_var': spec_var, 'flatness': flatness, 'sk_abs': sk_abs
    })

hs = [d['H'] for d in data_spec]
spec_maxs = [d['spec_max'] for d in data_spec]
spec_gaps = [d['spec_gap'] for d in data_spec]
flatnesses = [d['flatness'] for d in data_spec]

corr_max = np.corrcoef(spec_maxs, hs)[0,1]
corr_gap = np.corrcoef(spec_gaps, hs)[0,1]
corr_flat = np.corrcoef(flatnesses, hs)[0,1]

print(f"  Correlation(|skew_max|, H): {corr_max:.4f}")
print(f"  Correlation(spectral_gap, H): {corr_gap:.4f}")
print(f"  Correlation(flatness, H): {corr_flat:.4f}")

# Top 5 by H
top5 = sorted(data_spec, key=lambda d: -d['H'])[:5]
print(f"\n  Top 5 by H:")
for d in top5:
    print(f"    H={d['H']:>5}, |sk_max|={d['spec_max']:.3f}, "
          f"gap={d['spec_gap']:.3f}, flat={d['flatness']:.3f}")

# Bottom 5 by H
bot5 = sorted(data_spec, key=lambda d: d['H'])[:5]
print(f"\n  Bottom 5 by H:")
for d in bot5:
    print(f"    H={d['H']:>5}, |sk_max|={d['spec_max']:.3f}, "
          f"gap={d['spec_gap']:.3f}, flat={d['flatness']:.3f}")


# ============================================================
# 4. n=11 DRT spectral comparison
# ============================================================
print("\n" + "=" * 70)
print("n=11 DRT SPECTRAL COMPARISON")
print("=" * 70)

n11 = 11
T11_paley = paley_tournament(11)

# Non-Paley DRT at n=11: difference set {1,2,3,5,8}
T11_np = [[0]*11 for _ in range(11)]
DS_np = {1, 2, 3, 5, 8}
for i in range(11):
    for j in range(11):
        if i != j and (j - i) % 11 in DS_np:
            T11_np[i][j] = 1

# Verify they're different
h_paley11 = hamiltonian_path_count(T11_paley)
h_np11 = hamiltonian_path_count(T11_np)
print(f"  Paley T_11: H = {h_paley11}")
print(f"  Non-Paley DRT: H = {h_np11}")

# Skew spectra
sk_p = skew_spectrum(T11_paley)
sk_np = skew_spectrum(T11_np)
sk_p_abs = sorted([abs(x) for x in sk_p], reverse=True)
sk_np_abs = sorted([abs(x) for x in sk_np], reverse=True)

print(f"\n  Paley skew |eigenvalues|: {[f'{x:.4f}' for x in sk_p_abs[:6]]}")
print(f"  Non-Paley skew |eigenvalues|: {[f'{x:.4f}' for x in sk_np_abs[:6]]}")

# Adjacency spectra
adj_p = adjacency_spectrum(T11_paley)
adj_np = adjacency_spectrum(T11_np)
adj_p_abs = sorted([abs(x) for x in adj_p], reverse=True)
adj_np_abs = sorted([abs(x) for x in adj_np], reverse=True)

print(f"\n  Paley adj |eigenvalues|: {[f'{abs(x):.4f}' for x in adj_p[:6]]}")
print(f"  Non-Paley adj |eigenvalues|: {[f'{abs(x):.4f}' for x in adj_np[:6]]}")

# Cycle counts
c3_p = c3_from_score(T11_paley)
c5_p = c5_fast(T11_paley)
c3_np = c3_from_score(T11_np)
c5_np = c5_fast(T11_np)

print(f"\n  Paley: c3={c3_p}, c5={c5_p}")
print(f"  Non-Paley: c3={c3_np}, c5={c5_np}")

# Spectral formula: for DRTs, eigenvalues of S should relate to Gauss sums
# For Paley at prime p: S eigenvalues should be ±i*sqrt(p) (multiplicity (p-1)/2 each) + 0
# Check:
print(f"\n  Paley T_11 expected skew eigenvalue: ±i*sqrt(11) = ±{np.sqrt(11):.4f}i")
print(f"  Actual unique |values|: {sorted(set(round(x, 3) for x in sk_p_abs), reverse=True)}")

print(f"\n  Non-Paley T_11 unique |skew values|: {sorted(set(round(x, 3) for x in sk_np_abs), reverse=True)}")


# ============================================================
# 5. Spectral characterization of H-maximizer
# ============================================================
print("\n" + "=" * 70)
print("CONJECTURE: Paley = minimal spectral gap among regular tournaments")
print("=" * 70)

# For Paley T_p, all nonzero skew eigenvalues have the same |value| = sqrt(p)
# This means Paley has ZERO spectral gap in the skew spectrum
# (all nonzero eigenvalues equal in absolute value)

# Hypothesis: This "flat spectrum" property characterizes the H-maximizer

if regular_data:
    print(f"\n  Regular n=7 tournaments: spectral gap analysis")
    for h_val in sorted(h_groups.keys(), reverse=True):
        group = h_groups[h_val]
        gaps = []
        for d in group:
            sk = d['skew']
            sk_nz = sorted([abs(x) for x in sk if abs(x) > 0.01], reverse=True)
            if len(sk_nz) >= 2:
                gap = sk_nz[0] - sk_nz[-1]
                gaps.append(gap)
        if gaps:
            print(f"    H={h_val}: spectral gap = {np.mean(gaps):.4f} "
                  f"(min={min(gaps):.4f}, max={max(gaps):.4f})")

# Check at n=11
sk_p_nz = [x for x in sk_p_abs if x > 0.01]
sk_np_nz = [x for x in sk_np_abs if x > 0.01]
gap_p = max(sk_p_nz) - min(sk_p_nz) if len(sk_p_nz) >= 2 else 0
gap_np = max(sk_np_nz) - min(sk_np_nz) if len(sk_np_nz) >= 2 else 0
print(f"\n  n=11 Paley: skew spectral gap = {gap_p:.6f}")
print(f"  n=11 Non-Paley: skew spectral gap = {gap_np:.6f}")
print(f"  Paley has SMALLER gap: {gap_p < gap_np}")

# The Paley conference matrix property:
# For Paley T_p: S = A - A^T satisfies S^2 = -pI + J
# This means ALL nonzero eigenvalues of S have |value| = sqrt(p)
A_p = np.array(T11_paley, dtype=float)
S_p = A_p - A_p.T
S2 = S_p @ S_p
expected = -11 * np.eye(11) + np.ones((11, 11))
error = np.max(np.abs(S2 - expected))
print(f"\n  Conference matrix check for Paley T_11:")
print(f"    S^2 = -p*I + J: max error = {error:.10f} ({'PASS' if error < 1e-8 else 'FAIL'})")

# Same for non-Paley
A_np = np.array(T11_np, dtype=float)
S_np = A_np - A_np.T
S2_np = S_np @ S_np
error_np = np.max(np.abs(S2_np - expected))
print(f"    Non-Paley S^2 = -p*I + J: max error = {error_np:.10f} "
      f"({'PASS' if error_np < 1e-8 else 'FAIL'})")

if error_np > 1e-8:
    print(f"    => Non-Paley DRT does NOT have conference matrix property!")
    # What does S^2 look like for non-Paley?
    S2_np_diag = [S2_np[i][i] for i in range(11)]
    print(f"    Non-Paley S^2 diagonal: {[int(round(x)) for x in S2_np_diag]}")
    # Off-diagonal
    off_diag = set()
    for i in range(11):
        for j in range(11):
            if i != j:
                off_diag.add(int(round(S2_np[i][j])))
    print(f"    Non-Paley S^2 off-diagonal values: {sorted(off_diag)}")


print("\nDone.")
