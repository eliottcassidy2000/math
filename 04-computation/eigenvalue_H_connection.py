#!/usr/bin/env python3
"""
eigenvalue_H_connection.py -- kind-pasteur-2026-03-13-S60

For circulant tournaments on Z_p, the adjacency matrix has eigenvalues:
  lambda_j = sum_{s in S} omega^{js}  (j = 0, 1, ..., p-1)
where omega = exp(2*pi*i/p).

lambda_0 = |S| = (p-1)/2 (always).
For Paley: lambda_j = Gauss sum / character values.

Question: does the H-class structure at p=11 correspond to
some eigenvalue invariant?

Also: H(T) is related to the permanent of (I + A(T)), and
permanents of circulant matrices have nice eigenvalue formulas
(Minc's theorem, Schur's inequality, etc.).

In fact: perm(C) for circulant C with eigenvalues lambda_j is
  perm(C) = sum over permutations sigma of prod lambda_sigma(j)
which doesn't simplify directly, but the ABSOLUTE values
|lambda_j| are related to H.
"""

import cmath
import math
from collections import defaultdict

p = 11
m = (p - 1) // 2
omega = cmath.exp(2j * cmath.pi / p)

print(f"EIGENVALUE ANALYSIS AT p={p}")
print(f"{'='*70}")

# Known H-classes at p=11
H_classes = {
    95095: [2, 29],     # Class A (Paley/anti-Paley)
    93467: [4, 5, 7, 11, 15, 16, 20, 24, 26, 27],  # Class B
    93027: [0, 3, 6, 10, 13, 18, 21, 25, 28, 31],   # Class C
    92411: [1, 8, 9, 12, 14, 17, 19, 22, 23, 30],   # Class D
}

def bits_to_S(bits):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))
    return S

def eigenvalues(S, p):
    """Compute eigenvalues of circulant adjacency matrix."""
    eigs = []
    for j in range(p):
        lam = sum(omega**(j * s) for s in S)
        eigs.append(lam)
    return eigs

# Compute eigenvalue spectra for representatives of each class
print(f"\nEigenvalue spectra (|lambda_j| for j=0,...,{p-1}):")
for H_val, bits_list in sorted(H_classes.items(), reverse=True):
    bits = bits_list[0]
    S = bits_to_S(bits)
    eigs = eigenvalues(S, p)

    abs_eigs = sorted([abs(e) for e in eigs], reverse=True)
    print(f"\n  H={H_val} (bits={bits}), S={sorted(S)}:")
    print(f"    |lambda|: {', '.join(f'{x:.4f}' for x in abs_eigs)}")

    # Phase structure
    phases = sorted([(j, eigs[j].real, eigs[j].imag, abs(eigs[j]),
                       cmath.phase(eigs[j]) * 180 / cmath.pi)
                      for j in range(p)], key=lambda x: x[0])
    print(f"    j:  re(lam)  im(lam)  |lam|  phase(deg)")
    for j, re, im, ab, ph in phases:
        print(f"    {j}: {re:>8.4f} {im:>8.4f} {ab:>7.4f} {ph:>8.2f}")

# Key invariants from eigenvalues
print(f"\n{'='*70}")
print(f"  EIGENVALUE INVARIANTS")
print(f"{'='*70}")

for H_val, bits_list in sorted(H_classes.items(), reverse=True):
    bits = bits_list[0]
    S = bits_to_S(bits)
    eigs = eigenvalues(S, p)

    # Product of eigenvalues = det(A)
    det_A = 1
    for e in eigs:
        det_A *= e

    # Sum of |lambda_j|^2 = tr(A A^T) = number of edges p*(p-1)/2
    sum_sq = sum(abs(e)**2 for e in eigs)

    # Sum of |lambda_j|^4 = tr((AA^T)^2) = sum_j (d_j)^2 for circulant
    sum_4th = sum(abs(e)**4 for e in eigs)

    # Product of |lambda_j| (related to permanent via vdW bound)
    prod_abs = 1
    for e in eigs:
        prod_abs *= abs(e)

    # Multiset of |lambda| values (up to rounding)
    abs_multiset = sorted([round(abs(e), 6) for e in eigs], reverse=True)

    print(f"\n  H={H_val}:")
    print(f"    det(A) = {det_A:.4f}")
    print(f"    sum |lam|^2 = {sum_sq:.4f} (expected: p*(p-1)/2 = {p*(p-1)/2})")
    print(f"    sum |lam|^4 = {sum_4th:.4f}")
    print(f"    prod |lam| = {prod_abs:.4f}")
    print(f"    |lam| multiset = {abs_multiset}")

# For Paley: eigenvalues are Gauss sums
# lambda_j = sum_{s QR} omega^{js} = (chi(j)*g - 1)/2 where g is the Gauss sum
# g = sum_{t=1}^{p-1} chi(t) omega^t = i^((p-1)/2) * sqrt(p)  (for p=3 mod 4)
# At p=11: p = 3 mod 4, so g = i * sqrt(11) (Gauss sum)
# lambda_j = (chi(j) * i * sqrt(11) - 1) / 2  for j != 0
# lambda_0 = (p-1)/2 = 5

print(f"\n{'='*70}")
print(f"  PALEY GAUSS SUM STRUCTURE")
print(f"{'='*70}")

# Gauss sum for p=11
g = sum((1 if pow(t, (p-1)//2, p) == 1 else -1) * omega**t for t in range(1, p))
print(f"  Gauss sum g = {g:.6f}")
print(f"  |g| = {abs(g):.6f}, sqrt(p) = {math.sqrt(p):.6f}")
print(f"  g/sqrt(p) = {g/math.sqrt(p):.6f}")
print(f"  Expected: i^{(p-1)//2} = i^{(p-1)//2}")
i_power = 1j ** ((p-1)//2)
print(f"  i^{(p-1)//2} = {i_power:.6f}")
print(f"  g - i^{(p-1)//2} * sqrt(p) = {g - i_power * math.sqrt(p):.10f}")

# Check all Paley eigenvalues against the Gauss sum formula
S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
eigs_paley = eigenvalues(S_qr, p)

print(f"\n  Paley eigenvalue verification:")
for j in range(p):
    computed = eigs_paley[j]
    if j == 0:
        expected = (p-1)/2
    else:
        chi_j = 1 if pow(j, (p-1)//2, p) == 1 else -1
        expected = (chi_j * g - 1) / 2
    diff = abs(computed - expected)
    print(f"    j={j}: computed={computed:.6f}, expected={expected:.6f}, diff={diff:.2e}")

# For non-Paley orientations, the eigenvalues are NOT Gauss sums.
# They are character sums of the specific orientation set.
# The ABSOLUTE VALUES determine the tournament's spectral properties.

# Key question: does sum |lambda_j|^k for some k correlate with H?
print(f"\n{'='*70}")
print(f"  SPECTRAL INVARIANTS vs H")
print(f"{'='*70}")

spectral_data = []
for bits in range(1 << m):
    S = bits_to_S(bits)
    eigs = eigenvalues(S, p)

    # Various spectral invariants
    sum_2 = sum(abs(e)**2 for e in eigs).real
    sum_4 = sum(abs(e)**4 for e in eigs).real
    sum_6 = sum(abs(e)**6 for e in eigs).real
    prod_abs = 1
    for e in eigs:
        prod_abs *= abs(e)
    prod_abs = prod_abs.real

    # Sum of real parts
    sum_re = sum(e.real for e in eigs).real

    # Sum of |lambda_j|^2 * lambda_j (cubic invariant)
    sum_cubic = sum(abs(e)**2 * e for e in eigs).real

    spectral_data.append({
        'bits': bits, 'sum_2': sum_2, 'sum_4': sum_4, 'sum_6': sum_6,
        'prod': prod_abs, 'sum_re': sum_re, 'sum_cubic': sum_cubic
    })

# Group by H-class and check which spectral invariants distinguish classes
for inv_name in ['sum_2', 'sum_4', 'sum_6', 'prod', 'sum_re', 'sum_cubic']:
    vals_per_class = {}
    for H_val, bits_list in H_classes.items():
        class_vals = set()
        for bits in bits_list:
            v = next(s[inv_name] for s in spectral_data if s['bits'] == bits)
            class_vals.add(round(v, 6))
        vals_per_class[H_val] = sorted(class_vals)

    # Check if this invariant distinguishes all 4 classes
    all_vals = []
    for H_val in sorted(vals_per_class, reverse=True):
        all_vals.extend(vals_per_class[H_val])

    distinct = len(set(round(v, 4) for v in all_vals))
    separates = distinct == 4

    print(f"\n  {inv_name}: ", end="")
    for H_val in sorted(vals_per_class, reverse=True):
        print(f"H={H_val}: {vals_per_class[H_val]}, ", end="")
    print(f"  {'SEPARATES' if separates else 'does not separate'}")

print("\nDONE.")
