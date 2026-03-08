#!/usr/bin/env python3
"""
gaussian_integer_investigation.py — Why is I(Ω(T), i) always a Gaussian integer?

DISCOVERY: For all tournaments T tested (n≤5 exhaustive, n=6 sampled),
I(Ω(T), i) ∈ Z[i] (Gaussian integers).

Moreover at n≤5: Re(I(Ω, i)) ∈ {0, 1} always!

This script investigates:
1. The pattern in Re and Im parts
2. Whether |I(Ω, i)|² has a combinatorial interpretation
3. Connection to the independence polynomial's algebraic structure
4. Whether I(Ω, i) · I(Ω, -i) = |I(Ω, i)|² relates to H(T) or t₃

Author: opus-2026-03-07-S46f
"""

import sys
import os
import cmath
from itertools import combinations, permutations
from collections import Counter, defaultdict
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '03-artifacts', 'code'))
from tournament_lib import all_tournaments, random_tournament, hamiltonian_path_count

def find_odd_cycle_vertex_sets(T):
    n = len(T)
    cycle_vsets = set()
    for length in range(3, n + 1, 2):
        for combo in combinations(range(n), length):
            vset = frozenset(combo)
            if vset in cycle_vsets:
                continue
            first = combo[0]
            for perm in permutations(combo[1:]):
                path = (first,) + perm
                if all(T[path[i]][path[(i + 1) % length]] for i in range(length)):
                    cycle_vsets.add(vset)
                    break
    return list(cycle_vsets)

def indep_poly_coefficients(cycles):
    m = len(cycles)
    if m == 0:
        return [1]
    vsets = [frozenset(c) for c in cycles]
    nbr = [0] * m
    for i in range(m):
        for j in range(i + 1, m):
            if vsets[i] & vsets[j]:
                nbr[i] |= 1 << j
                nbr[j] |= 1 << i
    counts = [0] * (m + 1)
    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
        if ok:
            counts[bin(mask).count('1')] += 1
    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return counts

def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

def count_3cycles(T, n):
    return sum(1 for i,j,k in combinations(range(n), 3)
               if (T[i][j] and T[j][k] and T[k][i]) or
                  (T[i][k] and T[k][j] and T[j][i]))

def score_seq(T, n):
    return tuple(sorted(sum(T[i][j] for j in range(n)) for i in range(n)))

# =====================================================================
# PART 1: Exhaustive Gaussian integer analysis
# =====================================================================
print("=" * 70)
print("GAUSSIAN INTEGER ANALYSIS OF I(Ω(T), i)")
print("=" * 70)

for n in range(3, 7):
    print(f"\n--- n = {n} ---")

    gen = all_tournaments(n) if n <= 5 else (random_tournament(n) for _ in range(5000))
    label = "exhaustive" if n <= 5 else "5000 samples"

    real_parts = Counter()
    imag_parts = Counter()
    gi_values = Counter()  # (Re, Im) pairs
    norm_values = Counter()  # |I|²

    # Track correlation with other invariants
    data = []

    for T in gen:
        nn = len(T)
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        val = eval_poly(coeffs, 1j)
        H = eval_poly(coeffs, 2)  # = H(T) by OCF
        t3 = count_3cycles(T, nn)

        re = int(round(val.real))
        im = int(round(val.imag))
        norm = re*re + im*im

        real_parts[re] += 1
        imag_parts[im] += 1
        gi_values[(re, im)] += 1
        norm_values[norm] += 1

        data.append({
            'coeffs': coeffs, 'H': H, 't3': t3,
            're': re, 'im': im, 'norm': norm,
            'score': score_seq(T, nn), 'n_cycles': len(cycles)
        })

    total = len(data)
    print(f"  {total} tournaments ({label})")

    print(f"\n  Real parts: {dict(sorted(real_parts.items()))}")
    print(f"  Imag parts: {dict(sorted(imag_parts.items()))}")

    print(f"\n  (Re, Im) distribution:")
    for key in sorted(gi_values.keys()):
        if gi_values[key] >= max(1, total // 100):
            print(f"    ({key[0]:+d}, {key[1]:+d}i): {gi_values[key]}")

    print(f"\n  |I(Ω, i)|² = Re² + Im² distribution:")
    for norm in sorted(norm_values.keys()):
        print(f"    |I|² = {norm}: {norm_values[norm]}")

    # KEY: does |I(Ω, i)|² relate to H(T)?
    print(f"\n  Checking |I(Ω, i)|² vs H(T):")
    norm_by_H = defaultdict(Counter)
    for d in data:
        norm_by_H[d['H']][d['norm']] += 1

    for H in sorted(norm_by_H.keys()):
        norms = dict(sorted(norm_by_H[H].items()))
        if len(norms) == 1:
            norm_val = list(norms.keys())[0]
            print(f"    H={H}: |I|² = {norm_val} always")
        else:
            print(f"    H={H}: |I|² = {norms}")

    # Check: is |I(Ω, i)|² determined by t₃?
    norm_by_t3 = defaultdict(Counter)
    for d in data:
        norm_by_t3[d['t3']][d['norm']] += 1

    determined = all(len(v) == 1 for v in norm_by_t3.values())
    print(f"\n  |I(Ω, i)|² determined by t₃ alone? {determined}")
    if not determined:
        for t3 in sorted(norm_by_t3.keys()):
            if len(norm_by_t3[t3]) > 1:
                print(f"    t₃={t3}: |I|² = {dict(sorted(norm_by_t3[t3].items()))}")

    # Check: Re(I(Ω, i)) pattern
    re_vals = sorted(set(d['re'] for d in data))
    print(f"\n  Re(I(Ω, i)) values: {re_vals}")

    # IMPORTANT: Is Re always 0 or 1?
    always_01 = all(d['re'] in [0, 1] for d in data)
    print(f"  Re ∈ {{0, 1}} always? {always_01}")

    if not always_01:
        re_by_coeffs_degree = defaultdict(Counter)
        for d in data:
            deg = len(d['coeffs']) - 1
            re_by_coeffs_degree[deg][d['re']] += 1
        print(f"  Re by polynomial degree:")
        for deg in sorted(re_by_coeffs_degree.keys()):
            print(f"    degree {deg}: Re values = {dict(sorted(re_by_coeffs_degree[deg].items()))}")

    # Check: Im(I(Ω, i)) pattern
    # At x=i: I(Ω, i) = α₀ + α₁·i + α₂·i² + α₃·i³ + ...
    # = α₀ + α₁·i - α₂ - α₃·i + α₄ + α₅·i - ...
    # Re = α₀ - α₂ + α₄ - ...
    # Im = α₁ - α₃ + α₅ - ...
    print(f"\n  Algebraic structure:")
    print(f"    I(Ω, i) = Σ αₖ iᵏ")
    print(f"    Re(I) = α₀ - α₂ + α₄ - ...")
    print(f"    Im(I) = α₁ - α₃ + α₅ - ...")

    # Verify
    for d in data[:5]:
        c = d['coeffs']
        re_manual = sum(c[k] * ((-1)**(k//2)) for k in range(0, len(c), 2))
        im_manual = sum(c[k] * ((-1)**((k-1)//2)) for k in range(1, len(c), 2))
        print(f"    coeffs={c}: Re={re_manual} (check: {d['re']}), Im={im_manual} (check: {d['im']})")

# =====================================================================
# PART 2: WHY is Re always 0 or 1?
# =====================================================================
print("\n\n" + "=" * 70)
print("PART 2: WHY IS Re(I(Ω, i)) ∈ {0, 1}?")
print("=" * 70)

print("""
Re(I(Ω, i)) = α₀ - α₂ + α₄ - α₆ + ...

Since α₀ = 1 always (empty set), we need:
  α₂ - α₄ + α₆ - ... ∈ {0, 1}

This is the ALTERNATING sum of even-indexed independence numbers.
""")

for n in range(3, 7):
    print(f"\n--- n = {n} ---")
    gen = all_tournaments(n) if n <= 5 else (random_tournament(n) for _ in range(3000))

    even_alt_sum = Counter()  # α₂ - α₄ + α₆ - ...

    for T in gen:
        nn = len(T)
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)

        # α₂ - α₄ + α₆ - ... = 1 - Re(I(Ω, i))
        val = sum((-1)**(k//2) * coeffs[k] for k in range(2, len(coeffs), 2))
        even_alt_sum[val] += 1

    print(f"  α₂ - α₄ + α₆ - ... distribution: {dict(sorted(even_alt_sum.items()))}")
    print(f"  This equals 1 - Re(I(Ω, i))")

# =====================================================================
# PART 3: Is Im(I(Ω, i)) related to t₃ or H?
# =====================================================================
print("\n\n" + "=" * 70)
print("PART 3: Im(I(Ω, i)) PATTERNS")
print("=" * 70)

for n in range(3, 7):
    print(f"\n--- n = {n} ---")
    gen = all_tournaments(n) if n <= 5 else (random_tournament(n) for _ in range(3000))

    im_by_t3 = defaultdict(Counter)
    im_by_H = defaultdict(Counter)
    im_by_ncycles = defaultdict(Counter)

    for T in gen:
        nn = len(T)
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        val = eval_poly(coeffs, 1j)
        im = int(round(val.imag))
        H = eval_poly(coeffs, 2)
        t3 = count_3cycles(T, nn)

        im_by_t3[t3][im] += 1
        im_by_H[H][im] += 1
        im_by_ncycles[len(cycles)][im] += 1

    print(f"  Im by t₃:")
    for t3 in sorted(im_by_t3.keys()):
        vals = dict(sorted(im_by_t3[t3].items()))
        if len(vals) == 1:
            print(f"    t₃={t3}: Im = {list(vals.keys())[0]} always")
        else:
            print(f"    t₃={t3}: Im = {vals}")

    print(f"\n  Im by |Ω| (number of cycles):")
    for nc in sorted(im_by_ncycles.keys()):
        vals = dict(sorted(im_by_ncycles[nc].items()))
        print(f"    |Ω|={nc}: Im = {vals}")

# =====================================================================
# PART 4: |I(Ω, i)|² and the connection to Rédei parity
# =====================================================================
print("\n\n" + "=" * 70)
print("PART 4: |I(Ω, i)|² AND PARITY")
print("=" * 70)

print("""
Since I(Ω, i) ∈ Z[i], we have |I(Ω, i)|² ∈ Z.
|I(Ω, i)|² = Re² + Im² = N(I(Ω, i)) (Gaussian norm).

H(T) = I(Ω, 2) is always odd (Rédei).
I(Ω, -1) = Re - Im (from i→-1: need to re-derive).
Actually I(Ω, -1) = Σ αₖ(-1)ᵏ.

Key question: What is |I(Ω, i)|² mod 2?
""")

for n in range(3, 7):
    print(f"\n--- n = {n} ---")
    gen = all_tournaments(n) if n <= 5 else (random_tournament(n) for _ in range(3000))

    norm_mod2 = Counter()
    norm_mod4 = Counter()
    H_mod2_vs_norm = defaultdict(Counter)

    for T in gen:
        nn = len(T)
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        val = eval_poly(coeffs, 1j)
        re = int(round(val.real))
        im = int(round(val.imag))
        norm = re*re + im*im
        H = eval_poly(coeffs, 2)

        norm_mod2[norm % 2] += 1
        norm_mod4[norm % 4] += 1
        H_mod2_vs_norm[H % 2][norm % 2] += 1

    print(f"  |I(Ω, i)|² mod 2: {dict(sorted(norm_mod2.items()))}")
    print(f"  |I(Ω, i)|² mod 4: {dict(sorted(norm_mod4.items()))}")

    # Is norm always odd?
    always_odd = all(k == 1 for k in norm_mod2.keys())
    print(f"  |I(Ω, i)|² always odd: {always_odd}")

# =====================================================================
# PART 5: Connection to I(Ω, x) mod (x²+1)
# =====================================================================
print("\n\n" + "=" * 70)
print("PART 5: I(Ω, x) IN Z[x]/(x²+1) ≅ Z[i]")
print("=" * 70)

print("""
The fact that I(Ω, i) ∈ Z[i] is AUTOMATIC: since all αₖ are integers,
evaluating at x = i gives an element of Z[i].

But the INTERESTING question is: what patterns do Re and Im follow?

Re(I) = α₀ - α₂ + α₄ - ... = I_even(-1) where I_even(y) = Σ α_{2k} y^k
Im(I) = α₁ - α₃ + α₅ - ... = I_odd(-1) where I_odd(y) = Σ α_{2k+1} y^k

So Re(I(Ω,i)) = I_even(Ω, -1) = alternating sum of EVEN independence numbers
   Im(I(Ω,i)) = I_odd(Ω, -1) = alternating sum of ODD independence numbers

And I(Ω, -1) = α₀ - α₁ + α₂ - ... = Re - Im.
""")

# Verify
print("Verification: I(Ω,-1) = Re(I(Ω,i)) - Im(I(Ω,i))")
for n in [3, 4, 5]:
    all_match = True
    for T in all_tournaments(n):
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        val_i = eval_poly(coeffs, 1j)
        val_m1 = eval_poly(coeffs, -1)
        re = int(round(val_i.real))
        im = int(round(val_i.imag))
        if re - im != val_m1:
            all_match = False
    print(f"  n={n}: {all_match}")

print("""
STRUCTURAL INSIGHT:
  I(Ω, i) decomposes I(Ω, x) into "even" and "odd" parts.
  Re = alternating sum of even-sized independent sets: α₀ - α₂ + α₄ - ...
  Im = alternating sum of odd-sized independent sets: α₁ - α₃ + α₅ - ...

  The constraint Re ∈ {0, 1} (verified for n≤5) means:
    α₂ - α₄ + α₆ - ... ∈ {0, 1}

  This is a NON-TRIVIAL constraint on the independence polynomial!
  It says the even-alternating partial sum is ALMOST exactly α₀ = 1.
""")

# =====================================================================
# PART 6: The norm |I(Ω, i)|² and OCF
# =====================================================================
print("=" * 70)
print("PART 6: |I(Ω, i)|² = NORM AND ITS MEANING")
print("=" * 70)

print("""
|I(Ω, i)|² = (α₀ - α₂ + α₄ - ...)² + (α₁ - α₃ + α₅ - ...)²

At n=3: norms = {1, 2} → I(Ω,i) = 1 or 1+i
At n=4: norms = {1, 2, 5} → I(Ω,i) = 1, 1+i, or 1+2i
At n=5: norms = {1, 2, 5, 17, 26, 37}

Are these norms related to primes of the form a²+b²?
  1 = 0²+1², 2 = 1²+1², 5 = 1²+2², 17 = 1²+4², 26 = 1²+5², 37 = 1²+6²

Note: 26 = 2 × 13. And 13 = 2² + 3² is a Gaussian prime norm.
So the factorization structure in Z[i] matters!
""")

for n in range(3, 7):
    print(f"\n--- n = {n} ---")
    gen = all_tournaments(n) if n <= 5 else (random_tournament(n) for _ in range(3000))

    gi_vals = Counter()
    for T in gen:
        cycles = find_odd_cycle_vertex_sets(T)
        coeffs = indep_poly_coefficients(cycles)
        val = eval_poly(coeffs, 1j)
        re = int(round(val.real))
        im = int(round(val.imag))
        gi_vals[(re, im)] += 1

    print(f"  Gaussian integer values (Re, Im·i) and counts:")
    for (re, im), cnt in sorted(gi_vals.items()):
        norm = re*re + im*im
        print(f"    {re:+d} + {im:+d}i  (norm={norm}): {cnt}")

print("\nDone.")
