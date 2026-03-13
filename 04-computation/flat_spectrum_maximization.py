"""
FLAT SPECTRUM MAXIMIZATION — opus-2026-03-12-S67e

KEY INSIGHT FROM OTHER AGENT'S DATA:
  The PALEY tournament (QR set) BEATS the Interval for H!
  - p=7: Interval H=175, Paley H=189
  - p=11: Interval H=93027, Paley H=95095

  The Paley tournament has FLAT spectrum: Q_k = (p+1)/4 for all k.
  The Interval has PEAKED spectrum: Q_k = Fejér kernel samples.

THIS SCRIPT EXPLORES:
  1. Why flat spectrum maximizes H (convexity / Schur-concavity)
  2. Difference set characterization (flat ⟺ difference set)
  3. Automorphism group structure (AGL(1,p) for Paley)
  4. Dihedral subgroups and their role
  5. Entropy maximization interpretation
  6. The resonance cascade for FLAT spectra (contrast with Fibonacci)
"""

import numpy as np
from itertools import combinations, permutations
from math import gcd, factorial, comb
from fractions import Fraction

def is_qr(a, p):
    """Check if a is a quadratic residue mod p."""
    if a % p == 0:
        return False
    return pow(a, (p-1)//2, p) == 1

def qr_set(p):
    """Return the set of quadratic residues mod p."""
    return sorted([a for a in range(1, p) if is_qr(a, p)])

def fourier_magnitudes(S, p):
    """Compute Q_k = |S_hat(k)|^2 for k=1,...,(p-1)/2."""
    m = (p - 1) // 2
    omega = np.exp(2j * np.pi / p)
    Q = []
    for k in range(1, m + 1):
        S_hat = sum(omega ** (s * k) for s in S)
        Q.append(abs(S_hat) ** 2)
    return np.array(Q)

def count_hamiltonian_paths(S, p):
    """Count Hamiltonian paths in circulant tournament C_p^S starting from 0."""
    # DP over subsets
    n = p
    dp = {}
    dp[(1 << 0, 0)] = 1
    for mask_size in range(1, n):
        for (mask, v), count in list(dp.items()):
            if bin(mask).count('1') != mask_size:
                continue
            for s in S:
                w = (v + s) % n
                if not (mask & (1 << w)):
                    new_mask = mask | (1 << w)
                    key = (new_mask, w)
                    dp[key] = dp.get(key, 0) + count
    full_mask = (1 << n) - 1
    H = sum(dp.get((full_mask, v), 0) for v in range(n))
    return H // n, H  # H_from_0, H_total

def difference_multiset(S, p):
    """Compute the difference multiset {a-b mod p : a,b in S, a!=b}."""
    diffs = {}
    for a in S:
        for b in S:
            if a != b:
                d = (a - b) % p
                diffs[d] = diffs.get(d, 0) + 1
    return diffs

def elementary_symmetric(Q):
    """Compute elementary symmetric functions of Q."""
    m = len(Q)
    esf = [Fraction(1)]
    for j in range(1, m + 1):
        val = Fraction(0)
        for combo in combinations(range(m), j):
            prod = Fraction(1)
            for i in combo:
                prod *= Fraction(Q[i]).limit_denominator(10**12)
            val += prod
        esf.append(val)
    return esf

print("=" * 72)
print("PART 1: FLAT vs PEAKED SPECTRA — EMPIRICAL COMPARISON")
print("=" * 72)

for p in [7, 11, 13]:
    m = (p - 1) // 2
    print(f"\n--- p = {p}, m = {m} ---")

    # Interval
    S_int = list(range(1, m + 1))
    # Also equivalent to S_int_alt = [p-m, ..., p-1]
    Q_int = fourier_magnitudes(S_int, p)

    # Paley (QR set)
    S_qr = qr_set(p)
    Q_qr = fourier_magnitudes(S_qr, p)

    print(f"  Interval S = {S_int}")
    print(f"    Q_k = [{', '.join(f'{q:.4f}' for q in Q_int)}]")
    print(f"    prod(1+Q) = {np.prod(1 + Q_int):.4f}")
    print(f"    max/min Q = {max(Q_int)/min(Q_int):.4f}")
    print(f"    σ²(Q)/mean² = {np.var(Q_int)/np.mean(Q_int)**2:.4f} (coefficient of variation squared)")

    print(f"  Paley S = {S_qr}")
    print(f"    Q_k = [{', '.join(f'{q:.4f}' for q in Q_qr)}]")
    print(f"    prod(1+Q) = {np.prod(1 + Q_qr):.4f}")
    print(f"    max/min Q = {max(Q_qr)/min(Q_qr):.4f}")
    print(f"    σ²(Q)/mean² = {np.var(Q_qr)/np.mean(Q_qr)**2:.4f}")

    # Difference multiset analysis
    diffs_int = difference_multiset(S_int, p)
    diffs_qr = difference_multiset(S_qr, p)

    print(f"\n  Difference multiset (how many ways d = a-b, a,b ∈ S):")
    print(f"    Interval: {dict(sorted(diffs_int.items()))}")
    print(f"    Paley:    {dict(sorted(diffs_qr.items()))}")

    # λ = (m-1)/2 for a perfect difference set
    lambda_ideal = m * (m - 1) // (p - 1)
    print(f"    Perfect difference set λ = m(m-1)/(p-1) = {lambda_ideal}")

    # Check if QR is a difference set
    qr_vals = set(diffs_qr.values())
    print(f"    QR difference values: {qr_vals} ({'FLAT = difference set!' if len(qr_vals) <= 2 else 'NOT flat'})")


print("\n" + "=" * 72)
print("PART 2: WHY FLAT SPECTRUM → H MAXIMIZATION (CONVEXITY)")
print("=" * 72)

print("""
SCHUR-CONCAVITY ARGUMENT:

For the partition function Z = I(Ω, z) at fixed graph Ω, the quantity
  prod(z + Q_k) = z^m + e_1 z^{m-1} + ... + e_m
is a symmetric function of Q_k with fixed e_1 = m(p-m)/2.

By AM-GM: prod(z+Q_k) ≤ (z + mean(Q))^m  iff Q_k all equal.

Wait — that's the WRONG direction! AM-GM says:
  (geometric mean) ≤ (arithmetic mean)
  prod(z+Q_k)^{1/m} ≤ (1/m) Σ(z+Q_k) = z + mean(Q)

So prod(z+Q_k) ≤ (z + mean(Q))^m with equality iff all Q_k equal.

This means FLAT spectrum MAXIMIZES prod(1+Q_k)!

But H = I(Ω, 2) ≠ prod(1+Q_k) in general. The relationship is:
  H = p · amplification · prod(1+Q_k)

For the Paley tournament:
  prod(1+Q_k) is maximized (flat Q_k)
  amplification factor varies

Let's check if the H advantage comes entirely from the prod term.
""")

for p in [7, 11]:
    m = (p - 1) // 2
    S_int = list(range(1, m + 1))
    S_qr = qr_set(p)
    Q_int = fourier_magnitudes(S_int, p)
    Q_qr = fourier_magnitudes(S_qr, p)

    if p <= 11:
        H0_int, H_int = count_hamiltonian_paths(S_int, p)
        H0_qr, H_qr = count_hamiltonian_paths(S_qr, p)

        prod_int = np.prod(1 + Q_int)
        prod_qr = np.prod(1 + Q_qr)

        amp_int = H_int / (p * prod_int)
        amp_qr = H_qr / (p * prod_qr)

        print(f"\n  p = {p}:")
        print(f"    Interval: H={H_int}, prod(1+Q)={prod_int:.1f}, amp = H/(p·prod) = {amp_int:.6f}")
        print(f"    Paley:    H={H_qr}, prod(1+Q)={prod_qr:.1f}, amp = H/(p·prod) = {amp_qr:.6f}")
        print(f"    H ratio: Paley/Interval = {H_qr/H_int:.6f}")
        print(f"    prod ratio: {prod_qr/prod_int:.6f}")
        print(f"    amp ratio: {amp_qr/amp_int:.6f}")


print("\n" + "=" * 72)
print("PART 3: DIFFERENCE SETS AND GAUSS SUMS")
print("=" * 72)

print("""
THEOREM (classical): For p ≡ 3 mod 4, the quadratic residues mod p
form a (p, (p-1)/2, (p-3)/4)-difference set.

Equivalently: |Ŝ(k)|² = constant for all k ≠ 0.

PROOF SKETCH:
  Ŝ(k) = Σ_{a∈QR} ω^{ak} = (1/2)(Σ_{a=1}^{p-1} ω^{ak} + Σ χ(a)ω^{ak})
        = (1/2)(-1 + χ(k)^{-1} · g)  where g = Gauss sum

  |Ŝ(k)|² = (1/4)|χ(k)g - 1|² = (1/4)(p + 1 - 2·Re(χ(k)·g̅))

  For p ≡ 3 mod 4: g² = χ(-1)·p = (-1)·p = -p
  So g = i√p (times a root of unity). Then g̅ = -g.

  Re(χ(k)·g̅) = Re(-χ(k)·g) = -Re(χ(k)·g)

  But Σ_k Re(χ(k)·g) = Re(g · Σ χ(k)) = 0, and each |Re(χ(k)·g)| ≤ √p.

  Actually for p ≡ 3 mod 4:
  g = i^{(p-1)/2} · √p · (some 4th root of unity)

  The exact computation gives |Ŝ(k)|² = (p+1)/4 for ALL k.
""")

# Verify the Gauss sum computation
for p in [7, 11, 19, 23, 31, 43]:
    if p % 4 != 3:
        continue
    m = (p - 1) // 2
    S = qr_set(p)
    Q = fourier_magnitudes(S, p)

    expected = (p + 1) / 4
    is_flat = np.allclose(Q, expected, atol=1e-6)
    print(f"  p = {p}: Q_k = {Q[0]:.6f}, (p+1)/4 = {expected:.6f}, flat = {is_flat}")
    print(f"    prod(1+Q) = {np.prod(1+Q):.1f} = ((p+5)/4)^m = {((p+5)/4)**m:.1f}")

print("\n  For p ≡ 1 mod 4:")
for p in [5, 13, 17, 29, 37, 41]:
    m = (p - 1) // 2
    S = qr_set(p)
    Q = fourier_magnitudes(S, p)
    print(f"  p = {p}: Q range = [{min(Q):.4f}, {max(Q):.4f}], ratio = {max(Q)/min(Q):.4f}")
    print(f"    Note: QR does NOT form a tournament for p ≡ 1 mod 4 (symmetric)")


print("\n" + "=" * 72)
print("PART 4: AUTOMORPHISM GROUPS AND DIHEDRAL STRUCTURE")
print("=" * 72)

print("""
The Paley tournament C_p^{QR} has automorphism group:
  Aut(C_p^QR) = AGL(1, p)_QR = {x ↦ ax + b : a ∈ QR(p), b ∈ F_p}

This has order p · (p-1)/2 = p · m.

Group structure:
  AGL(1,p)_QR ≅ F_p ⋊ QR(p)* ≅ Z_p ⋊ Z_m

  Translations: x ↦ x + b (form Z_p, act on vertices)
  Dilations: x ↦ ax, a ∈ QR (form Z_m ≅ Z_{(p-1)/2})

DIHEDRAL SUBGROUP:
  For p ≡ 3 mod 4, -1 is NOT a QR (since (-1)^{(p-1)/2} = -1).
  So the map x ↦ -x is NOT in the automorphism group.

  However, if we consider the UNDIRECTED Paley graph (p ≡ 1 mod 4),
  then x ↦ -x IS an automorphism, giving a dihedral structure.

  For the tournament (p ≡ 3 mod 4):
  The "reflections" are x ↦ ax where a has order 2 in QR*.
  If (p-1)/2 is even, there exist involutions in QR*.

  Example: p = 7, m = 3. QR = {1, 2, 4}.
  Powers: 2^1=2, 2^2=4, 2^3=1. So QR* = <2> ≅ Z_3.
  No involution in Z_3 → no dihedral subgroup from dilations alone.

  BUT: Consider the extended group including x ↦ -x/a for a ∈ QR.
  These are "orientation-reversing" maps.
""")

# Compute automorphism group structure for small primes
for p in [7, 11, 23]:
    m = (p - 1) // 2
    qr = qr_set(p)

    # Find the generator of QR* (multiplicative group)
    for g in qr:
        if g == 1:
            continue
        powers = set()
        x = 1
        for _ in range(m):
            x = (x * g) % p
            powers.add(x)
        if powers == set(qr):
            print(f"\n  p = {p}: QR* = <{g}> ≅ Z_{m}")
            print(f"    |Aut| = p · m = {p * m}")

            # Check for elements of order 2 (involutions)
            involutions = [a for a in qr if pow(a, 2, p) == 1 and a != 1]
            if involutions:
                print(f"    Involutions in QR*: {involutions}")
            else:
                print(f"    No involutions in QR* (m = {m} is odd)")

            # Orbits on pairs under the automorphism group
            # For flat spectrum: all non-trivial Fourier modes are equivalent
            print(f"    Orbits on frequency space: single orbit (flat Q_k)")
            break


print("\n" + "=" * 72)
print("PART 5: ENTROPY AND INFORMATION THEORY")
print("=" * 72)

print("""
INFORMATION-THEORETIC VIEW:

Define the normalized Q-distribution: q_k = Q_k / Σ Q_k

Shannon entropy: H_ent = -Σ q_k log(q_k)

Maximum entropy (uniform): H_max = log(m)
Interval (peaked): H_ent << H_max
Paley (flat): H_ent = H_max exactly
""")

for p in [7, 11, 13, 17, 23, 29]:
    m = (p - 1) // 2

    # Interval
    S_int = list(range(1, m + 1))
    Q_int = fourier_magnitudes(S_int, p)
    q_int = Q_int / Q_int.sum()
    H_ent_int = -np.sum(q_int * np.log(q_int))

    # Paley (for p ≡ 3 mod 4 only)
    S_qr = qr_set(p)
    Q_qr = fourier_magnitudes(S_qr, p)
    q_qr = Q_qr / Q_qr.sum()
    H_ent_qr = -np.sum(q_qr * np.log(q_qr))

    H_max = np.log(m)

    print(f"  p = {p}, m = {m}:")
    print(f"    Interval: H_ent = {H_ent_int:.4f} / {H_max:.4f} = {H_ent_int/H_max:.4f}")
    print(f"    Paley:    H_ent = {H_ent_qr:.4f} / {H_max:.4f} = {H_ent_qr/H_max:.4f}")

    # Rényi entropy of order 2 (log of 1/IPR)
    IPR_int = np.sum(q_int**2)
    IPR_qr = np.sum(q_qr**2)
    print(f"    Rényi-2: Interval IPR = {IPR_int:.4f}, Paley IPR = {IPR_qr:.4f} (= 1/m = {1/m:.4f})")


print("\n" + "=" * 72)
print("PART 6: AM-GM BOUND AND THE FLAT SPECTRUM THEOREM")
print("=" * 72)

print("""
FLAT SPECTRUM THEOREM (to prove):

Among all connection sets S ⊂ {1,...,p-1} with |S| = m,
the one maximizing prod(1+Q_k) is the one with Q_k = constant.

PROOF: By AM-GM inequality applied to (1+Q_k):
  [Π(1+Q_k)]^{1/m} ≤ (1/m) Σ(1+Q_k) = 1 + (1/m)Σ Q_k = 1 + m/2

  (since Σ Q_k = m²/2 = m(p-m)/2 at m=(p-1)/2... wait, Σ Q_k = m·mean(Q))

Actually: Σ Q_k = Σ_{k=1}^m |Ŝ(k)|² = Σ_{a,b∈S} Σ_{k=1}^m ω^{(a-b)k}
For a=b: contributes m·|S| = m²
For a≠b: Σ_{k=1}^m ω^{(a-b)k} depends on whether sum over HALF the roots

Parseval: Σ_{k=0}^{p-1} |Ŝ(k)|² = p·|S| = pm
But Ŝ(0) = m, so Σ_{k=1}^{p-1} |Ŝ(k)|² = pm - m² = m(p-m)
By symmetry |Ŝ(k)|² = |Ŝ(p-k)|², so Σ_{k=1}^m |Ŝ(k)|² = m(p-m)/2

So Σ Q_k = m(p-m)/2 for ALL S.
Mean(Q) = (p-m)/2 = (p+1)/4 (at m = (p-1)/2).
""")

# Verify Parseval
for p in [7, 11, 13, 17]:
    m = (p - 1) // 2
    S_int = list(range(1, m + 1))
    S_qr = qr_set(p)
    Q_int = fourier_magnitudes(S_int, p)
    Q_qr = fourier_magnitudes(S_qr, p)

    expected = m * (p - m) / 2
    print(f"  p={p}: Σ Q_int = {Q_int.sum():.4f}, Σ Q_qr = {Q_qr.sum():.4f}, m(p-m)/2 = {expected}")

print("""
So: prod(1+Q_k) ≤ (1 + (p+1)/4)^m = ((p+5)/4)^m

Equality iff Q_k = (p+1)/4 for all k iff S is a DIFFERENCE SET.

For p ≡ 3 mod 4: QR(p) IS a difference set → achieves the bound!
For p ≡ 1 mod 4: QR(p) is NOT a connection set (since -1 ∈ QR,
  the set is symmetric and doesn't give a tournament).

QUESTION: For p ≡ 1 mod 4, can any S achieve Q_k = constant?
""")

# For p ≡ 1 mod 4, search for flat-spectrum sets
for p in [5, 13]:
    m = (p - 1) // 2
    print(f"\n  p = {p}, searching for flattest S among all {comb(p-1, m)} sets:")

    best_flatness = float('inf')
    best_S = None
    best_prod = 0

    all_elements = list(range(1, p))
    count = 0
    for S in combinations(all_elements, m):
        S = list(S)
        # Check if S is a valid connection set (if s ∈ S then p-s ∉ S)
        valid = True
        for s in S:
            if (p - s) in S:
                valid = False
                break
        if not valid:
            continue
        count += 1

        Q = fourier_magnitudes(S, p)
        flatness = np.var(Q)
        prod_val = np.prod(1 + Q)

        if prod_val > best_prod:
            best_prod = prod_val
            best_S_prod = S
            best_Q_prod = Q

        if flatness < best_flatness:
            best_flatness = flatness
            best_S = S
            best_Q = Q

    print(f"    Total valid connection sets: {count}")
    print(f"    Flattest: S = {best_S}, Q = [{', '.join(f'{q:.4f}' for q in best_Q)}]")
    print(f"      var(Q) = {np.var(best_Q):.6f}, prod(1+Q) = {np.prod(1+best_Q):.4f}")
    print(f"    Max prod(1+Q): S = {best_S_prod}, Q = [{', '.join(f'{q:.4f}' for q in best_Q_prod)}]")
    print(f"      var(Q) = {np.var(best_Q_prod):.6f}, prod(1+Q) = {np.prod(1+best_Q_prod):.4f}")

    # Also check H for the best prod set if small enough
    if p <= 11:
        H0, H = count_hamiltonian_paths(best_S_prod, p)
        print(f"      H = {H}, H/p = {H0}")


print("\n" + "=" * 72)
print("PART 7: H vs prod(1+Q) — IS H MONOTONE IN prod?")
print("=" * 72)

# For p = 7 and 11, compute H and prod(1+Q) for ALL valid connection sets
for p in [7, 11]:
    m = (p - 1) // 2
    print(f"\n  p = {p}:")

    results = []
    all_elements = list(range(1, p))
    for S in combinations(all_elements, m):
        S = list(S)
        valid = True
        for s in S:
            if (p - s) in S:
                valid = False
                break
        if not valid:
            continue

        Q = fourier_magnitudes(S, p)
        prod_val = np.prod(1 + Q)

        if p <= 11:
            H0, H = count_hamiltonian_paths(S, p)
            results.append((S, Q, prod_val, H, H0))

    # Sort by prod(1+Q)
    results.sort(key=lambda x: x[2])

    print(f"  {'S':<25} {'prod(1+Q)':>12} {'H':>10} {'H/p':>8} {'H/(p·prod)':>12}")
    for S, Q, prod_val, H, H0 in results:
        print(f"  {str(S):<25} {prod_val:>12.1f} {H:>10} {H0:>8} {H/(p*prod_val):>12.4f}")

    # Correlation between prod and H
    prods = [r[2] for r in results]
    Hs = [r[3] for r in results]
    corr = np.corrcoef(prods, Hs)[0, 1]
    print(f"\n  Correlation(prod(1+Q), H) = {corr:.6f}")


print("\n" + "=" * 72)
print("PART 8: THE DIHEDRAL CONNECTION — SYMMETRY & MAXIMIZATION")
print("=" * 72)

print("""
DIHEDRAL GROUPS AND TOURNAMENT AUTOMORPHISMS:

The dihedral group D_p = <r, s | r^p = s^2 = (rs)^2 = 1> has order 2p.
  r = rotation (vertex relabeling v ↦ v+1 mod p)
  s = reflection (vertex relabeling v ↦ -v mod p)

For a CIRCULANT tournament C_p^S:
  - r is ALWAYS an automorphism (by definition of circulant)
  - s is an automorphism iff S = p - S (the "palindrome" condition)

For the Paley tournament (p ≡ 3 mod 4):
  - S = QR(p). Is QR(p) = p - QR(p)?
  - (p - a) mod p ∈ QR iff (-1)·a ∈ QR iff χ(-1)·χ(a) = 1
  - χ(-1) = (-1)^{(p-1)/2} = -1 for p ≡ 3 mod 4
  - So p - a ∈ QR iff a ∈ QNR. So QR ≠ p - QR.
  - Therefore D_p is NOT in Aut(Paley) for p ≡ 3 mod 4.

For the Interval S = {1,...,m}:
  - p - S = {m+1,...,p-1}. Clearly S ≠ p - S.
  - So D_p is NOT in Aut(Interval) either.

HOWEVER: D_p acts on the UNDIRECTED Paley graph (p ≡ 1 mod 4).
And for p ≡ 3 mod 4, the dihedral group D_{(p-1)/2} acts on the
FREQUENCY SPACE {1,...,m} via the Galois group of Q(ζ_p).

Let's check which circulant tournaments have D_p-symmetry.
""")

for p in [7, 11]:
    m = (p - 1) // 2
    print(f"\n  p = {p}: Connection sets with palindrome S = p - S:")
    all_elements = list(range(1, p))
    palindrome_sets = []
    for S in combinations(all_elements, m):
        S_set = set(S)
        complement = set((p - s) % p for s in S)
        # For valid tournament: S ∩ (p-S) = ∅ and S ∪ (p-S) = {1,...,p-1}
        # Palindrome means S = p - S, which contradicts S ∩ (p-S) = ∅
        # So NO circulant tournament has D_p symmetry!
        if S_set == complement:
            palindrome_sets.append(S)

    if palindrome_sets:
        print(f"    Found: {palindrome_sets}")
    else:
        print(f"    NONE — D_p cannot be an automorphism of any circulant tournament on p vertices")
        print(f"    (because S = p-S contradicts the tournament property S ∩ (p-S) = ∅)")


print("\n" + "=" * 72)
print("PART 9: GENERALIZED RESONANCE FOR FLAT SPECTRA")
print("=" * 72)

print("""
For the Paley tournament with flat Q_k = c = (p+1)/4:

  prod(1+Q_k) = (1+c)^m = ((p+5)/4)^m

  prod(a+Q_k) = (a+c)^m   for all a

  This is a MONOMIAL in (a+c), not a polynomial in a with
  interesting structure like the Interval's Morgan-Voyce.

The Interval's resonance cascade arises because the RECURRENCE
  B_m = (2+x)B_{m-1} - B_{m-2}
has a non-trivial interaction between levels (the -B_{m-2} term).

For flat spectrum: prod(a+c) = (a+c)^m has NO recurrence — it's
just exponentiation. The "resonance" is trivial (pure amplification).

COMPARISON:
  Interval: prod(1+Q) = F_p ~ φ^p/√5 (Fibonacci — needs resonance)
  Paley: prod(1+Q) = ((p+5)/4)^m (exponential — no resonance needed)

  Which grows faster?
""")

# Compare growth rates
print("  Growth comparison:")
print(f"  {'p':>4} {'F_p':>15} {'((p+5)/4)^m':>15} {'ratio Paley/Fib':>18}")
for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    m = (p - 1) // 2

    # Fibonacci numbers
    a, b = 1, 1
    for _ in range(p - 1):
        a, b = b, a + b
    F_p = b

    # Paley product (only for p ≡ 3 mod 4)
    if p % 4 == 3:
        paley_prod = ((p + 5) / 4) ** m
        ratio = paley_prod / F_p
        print(f"  {p:>4} {F_p:>15} {paley_prod:>15.1f} {ratio:>18.6f}")
    else:
        print(f"  {p:>4} {F_p:>15} {'(p≡1 mod 4)':>15} {'N/A':>18}")


print("\n" + "=" * 72)
print("PART 10: THE DEEP CONNECTION — LATTICE GAS AT FLAT POTENTIAL")
print("=" * 72)

print("""
LATTICE GAS INTERPRETATION:

H(T) = I(Ω(T), 2) = independence polynomial of Ω at z=2.

For the Paley tournament with Q_k = c = (p+1)/4 for all k:
  - The spectral input is UNIFORM across all channels
  - Each Fourier mode contributes equally to path counting
  - The amplification factor H/(p · prod(1+Q)) measures how
    the odd-cycle structure of Ω adds to the spectral base

KEY QUESTION: Is the amplification factor H/(p·(1+c)^m) smaller
for Paley than for Interval? If so, the Paley tournament wins by
brute-force spectral superiority (larger prod(1+Q)) despite having
LESS graph-level amplification.

From the data:
  p=7: Interval amp = 1.923, Paley amp = 1.000
  p=11: Interval amp = 95.02, Paley amp = 8.442

So YES: Paley has LOWER amplification but HIGHER base.
As p grows, does the base advantage overcome the amplification gap?

ANALOGY: Two investment strategies:
  - Interval: Low principal (F_p) × high multiplier (grows factorially)
  - Paley: High principal ((1+c)^m) × low multiplier

  For large p, which dominates?
""")

# Asymptotic comparison
print("  Asymptotic growth rates:")
print(f"  {'p':>4} {'log(F_p)':>10} {'m·log(1+c)':>12} {'Fib growth rate':>16} {'Paley growth rate':>18}")
for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
    m = (p - 1) // 2
    log_phi = np.log((1 + np.sqrt(5)) / 2)
    log_Fp = p * log_phi - np.log(np.sqrt(5))  # asymptotic

    if p % 4 == 3:
        c = (p + 1) / 4
        log_paley = m * np.log(1 + c)
        fib_rate = log_Fp / p
        paley_rate = log_paley / p
        print(f"  {p:>4} {log_Fp:>10.2f} {log_paley:>12.2f} {fib_rate:>16.6f} {paley_rate:>18.6f}")


print("\n" + "=" * 72)
print("SYNTHESIS: FLAT SPECTRUM MAXIMIZATION")
print("=" * 72)

print("""
MAIN FINDINGS:

1. FLAT SPECTRUM THEOREM: prod(1+Q_k) is maximized when Q_k = constant
   = (p+1)/4. This follows from AM-GM with constraint Σ Q_k = m(p-m)/2.

2. DIFFERENCE SET CONNECTION: Q_k = constant ⟺ S is a difference set.
   For p ≡ 3 mod 4, the quadratic residues QR(p) form a (p,m,(m-1)/2)-
   difference set, achieving the bound.

3. PALEY BEATS INTERVAL: At p=7 and p=11, the Paley tournament has
   higher H than the Interval, despite the Interval's Fibonacci structure.
   The advantage comes from the spectral base prod(1+Q), not from
   the graph amplification factor.

4. GROWTH RATE CROSSOVER: Paley's growth rate m·log((p+5)/4) ~ (p/2)·log(p/4)
   dominates Fibonacci's p·log(φ) for large p. The crossover is:
   (1/2)·log(p/4) > log(φ) ⟹ p > 4φ² ≈ 10.5
   So for ALL p ≥ 11, the Paley spectral base grows faster than Fibonacci.

5. TWO RESONANCE MECHANISMS:
   - INTERVAL: Fibonacci resonance cascade (recurrence → φ-locking → F_p)
     Beautiful mathematics but SUBOPTIMAL for H.
   - PALEY: Flat resonance (AM-GM optimality → (1+c)^m)
     Less romantic but DOMINANT for H.

6. DIHEDRAL CONNECTION: The automorphism group AGL(1,p) for Paley
   is NOT dihedral, but is a semidirect product Z_p ⋊ Z_m.
   Dihedral D_p cannot be an automorphism of ANY circulant tournament
   (since the reflection v↦-v would require S = p-S, contradicting
   the tournament property).

7. OPEN QUESTION: For p ≡ 1 mod 4, there is no Paley tournament.
   What connection set achieves the flattest spectrum and highest H?
   Does it involve quadratic residues in a modified way?
""")
