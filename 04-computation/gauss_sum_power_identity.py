#!/usr/bin/env python3
"""
Gauss Sum Power Identity: Why Paley tournaments maximize 5-cycles among circulants.

For Z_p circulant tournaments, eigenvalues of adjacency matrix A are:
  λ_k = Σ_{s∈S} ω^{ks}   where ω = e^{2πi/p}, S = connection set with |S|=(p-1)/2

The trace formula: c_m = (1/m) Σ_{k=0}^{p-1} λ_k^m  counts m-cycles.

Key question: Why does S = QR_p (quadratic residues) maximize c_5?

Author: opus research agent
"""

import numpy as np
from itertools import combinations
from fractions import Fraction
import cmath

def omega(p):
    """Primitive p-th root of unity."""
    return np.exp(2j * np.pi / p)

def eigenvalues_circulant(p, S):
    """Compute eigenvalues λ_k = Σ_{s∈S} ω^{ks} for k=0,...,p-1."""
    w = omega(p)
    eigs = []
    for k in range(p):
        lam = sum(w**(k*s) for s in S)
        eigs.append(lam)
    return eigs

def cycle_count_from_eigenvalues(eigs, m):
    """c_m = (1/m) * Σ λ_k^m (should be a real integer for tournaments)."""
    total = sum(lam**m for lam in eigs)
    return total / m

def quadratic_residues(p):
    """QR mod p."""
    qr = set()
    for x in range(1, p):
        qr.add((x*x) % p)
    return sorted(qr)

def all_circulant_tournaments(p):
    """
    Generate all connection sets S ⊂ Z_p^* with |S| = (p-1)/2
    such that S ∩ (-S) = ∅ (tournament condition: exactly one of s, p-s in S).
    """
    half = (p - 1) // 2
    # Pairs: {s, p-s} for s = 1,...,(p-1)/2
    pairs = []
    for s in range(1, half + 1):
        pairs.append((s, p - s))

    # Each tournament picks one from each pair
    tournaments = []
    for bits in range(2**half):
        S = []
        for i in range(half):
            if bits & (1 << i):
                S.append(pairs[i][1])
            else:
                S.append(pairs[i][0])
        tournaments.append(sorted(S))
    return tournaments

def is_isomorphic_circulant(p, S1, S2):
    """Check if S1 and S2 give isomorphic circulants (via multiplication by units)."""
    S1_set = set(S1)
    for a in range(1, p):
        if set((a * s) % p for s in S2) == S1_set:
            return True
    return False

def gauss_sum(p, S, k):
    """G_k(S) = Σ_{s∈S} ω^{ks}."""
    w = omega(p)
    return sum(w**(k*s) for s in S)

def print_separator(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}\n")

# ============================================================
#  SECTION 1: All Z_7 circulant tournaments
# ============================================================
print_separator("SECTION 1: Z_7 Circulant Tournaments — Eigenvalue Analysis")

p = 7
qr7 = quadratic_residues(p)
print(f"QR_7 = {qr7}")

all_tours = all_circulant_tournaments(p)
print(f"Total circulant tournaments on Z_7: {len(all_tours)}")

# Group by isomorphism class
classes = []
used = [False] * len(all_tours)
for i, S in enumerate(all_tours):
    if used[i]:
        continue
    cls = [S]
    for j in range(i+1, len(all_tours)):
        if not used[j] and is_isomorphic_circulant(p, S, all_tours[j]):
            cls.append(all_tours[j])
            used[j] = True
    used[i] = True
    classes.append(cls)

print(f"Isomorphism classes: {len(classes)}\n")

for ci, cls in enumerate(classes):
    S = cls[0]
    is_paley = (set(S) == set(qr7))
    label = " ← PALEY" if is_paley else ""
    eigs = eigenvalues_circulant(p, S)

    print(f"Class {ci+1}: S = {S}{label}  (size {len(cls)})")

    # Print eigenvalues
    for k in range(p):
        lam = eigs[k]
        mag = abs(lam)
        phase = cmath.phase(lam)
        print(f"  λ_{k} = {lam.real:+8.4f} {lam.imag:+8.4f}i   |λ|={mag:.4f}  arg={phase:.4f}")

    # Cycle counts
    for m in [3, 5, 7]:
        cm = cycle_count_from_eigenvalues(eigs, m)
        print(f"  c_{m} = {cm.real:.1f}  (imag: {cm.imag:.2e})")

    # Power sums
    for m in [3, 5]:
        ps = sum(lam**m for lam in eigs)
        print(f"  Σλ^{m} = {ps.real:.1f}")
    print()

# ============================================================
#  SECTION 2: Verify c_3 universality and compute c_5
# ============================================================
print_separator("SECTION 2: c_3 and c_5 for ALL Z_7 circulants")

for S in all_tours:
    eigs = eigenvalues_circulant(p, S)
    c3 = cycle_count_from_eigenvalues(eigs, 3).real
    c5 = cycle_count_from_eigenvalues(eigs, 5).real
    is_paley = (set(S) == set(qr7))
    tag = " ← PALEY" if is_paley else ""
    print(f"S={S}:  c_3={c3:5.0f}  c_5={c5:5.0f}{tag}")

# ============================================================
#  SECTION 3: Eigenvalue magnitude analysis
# ============================================================
print_separator("SECTION 3: Why Paley Maximizes c_5 — Eigenvalue Magnitudes")

print("Key insight: c_5 = (1/5) Σ λ_k^5")
print("For circulants, λ_0 = |S| = (p-1)/2 always.")
print("The k≠0 eigenvalues are Gauss-type sums.\n")

for ci, cls in enumerate(classes):
    S = cls[0]
    is_paley = (set(S) == set(qr7))
    label = "PALEY" if is_paley else "NON-PALEY"
    eigs = eigenvalues_circulant(p, S)

    print(f"--- {label}: S = {S} ---")
    mags = [abs(e) for e in eigs]
    mag5 = [abs(e)**5 for e in eigs]
    print(f"  |λ_k|:   {['%.4f' % m for m in mags]}")
    print(f"  |λ_k|^5: {['%.2f' % m for m in mag5]}")
    print(f"  Σ|λ_k|^5 = {sum(mag5):.2f}")

    # Actual 5th power sum (complex)
    ps5 = sum(e**5 for e in eigs)
    print(f"  Σλ_k^5   = {ps5.real:.2f} (actual, determines c_5)")
    print(f"  Ratio actual/max_possible = {ps5.real / sum(mag5):.4f}")
    print()

print("For Paley (QR): ALL |λ_k| = √p for k≠0 (Gauss sum property)")
print("  → |λ_k|^5 = p^{5/2} for each k≠0")
print("  → Maximum possible 5th power contribution from k≠0 terms")
print()
print("For non-Paley: |λ_k| vary, and by power-mean inequality,")
print("  Σ|λ_k|^5 < (p-1)*p^{5/2} when magnitudes are unequal.")

# ============================================================
#  SECTION 4: Gauss sum identity for QR
# ============================================================
print_separator("SECTION 4: Gauss Sum Identity for QR Sets")

print("For S = QR_p, the eigenvalue is a Gauss sum:")
print("  G_k = Σ_{s∈QR} ω^{ks}")
print()
print("Classical result: G_k = (χ(k)*g(χ) - 1) / 2")
print("  where g(χ) = Σ_{t=1}^{p-1} χ(t)ω^t is the quadratic Gauss sum")
print("  and χ is the Legendre symbol mod p")
print()
print("For p ≡ 3 (mod 4): g(χ) = i√p")
print("For p ≡ 1 (mod 4): g(χ) = √p")
print()

for pp in [5, 7, 11, 13]:
    qr = quadratic_residues(pp)
    print(f"\np = {pp}, QR = {qr}")
    w = omega(pp)

    # Compute quadratic Gauss sum
    g_chi = sum(pow(t, (pp-1)//2, pp) * (1 if pow(t, (pp-1)//2, pp) == 1 else -1) * w**t
                for t in range(1, pp))
    # Actually let me compute it properly
    def legendre(a, p):
        if a % p == 0: return 0
        return 1 if pow(a, (p-1)//2, p) == 1 else -1

    g_chi = sum(legendre(t, pp) * w**t for t in range(1, pp))
    print(f"  g(χ) = {g_chi.real:.4f} + {g_chi.imag:.4f}i")
    print(f"  |g(χ)|² = {abs(g_chi)**2:.2f}  (should be {pp})")

    # Verify eigenvalue formula
    eigs_direct = eigenvalues_circulant(pp, qr)
    for k in range(pp):
        chi_k = legendre(k, pp) if k > 0 else 0
        if k > 0:
            lam_formula = (chi_k * g_chi - 1) / 2
        else:
            lam_formula = len(qr)
        lam_direct = eigs_direct[k]
        err = abs(lam_formula - lam_direct)
        if k <= 3 or err > 1e-10:
            print(f"  k={k}: λ_direct={lam_direct.real:+.4f}{lam_direct.imag:+.4f}i  "
                  f"λ_formula={lam_formula.real:+.4f}{lam_formula.imag:+.4f}i  err={err:.2e}")

# ============================================================
#  SECTION 5: The 5th power sum via Gauss sums
# ============================================================
print_separator("SECTION 5: Fifth Power Sum Identity")

print("For S = QR_p, λ_k = (χ(k)·g - 1)/2 where g = g(χ).")
print("So λ_k^5 = (1/32)(χ(k)g - 1)^5")
print()
print("Expanding (χ(k)g - 1)^5 by binomial theorem:")
print("  = Σ_{j=0}^{5} C(5,j) (χ(k)g)^j (-1)^{5-j}")
print("  = -1 + 5χg - 10g² + 10χg³ - 5g⁴ + χg⁵")
print("  (using χ(k)^j = χ(k) for odd j, 1 for even j, when k≠0)")
print()

def legendre(a, p):
    if a % p == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

for pp in [5, 7, 11, 13, 17, 19, 23]:
    qr = quadratic_residues(pp)
    w = omega(pp)
    g = sum(legendre(t, pp) * w**t for t in range(1, pp))

    # Sum of χ(k) over k=1..p-1
    sum_chi = sum(legendre(k, pp) for k in range(1, pp))  # = 0 always

    # Compute Σ_{k=1}^{p-1} λ_k^5
    eigs = eigenvalues_circulant(pp, qr)
    sum5_direct = sum(eigs[k]**5 for k in range(1, pp))

    # Formula: Σ_{k≠0} (χg-1)^5 / 32
    # = (1/32) Σ_{k≠0} [-1 + 5χ(k)g - 10g² + 10χ(k)g³ - 5g⁴ + χ(k)g⁵]
    # Note: Σχ(k) = 0, so terms with odd power of χ vanish!
    # = (1/32)(p-1)(-1 - 10g² - 5g⁴) + (1/32)·0·(5g + 10g³ + g⁵)
    # = (p-1)/32 * (-1 - 10g² - 5g⁴)

    # Wait, let me be more careful. χ(k)^2 = 1 for k≠0.
    # χ(k)^{odd} = χ(k). So:
    # (χg)^0 = 1
    # (χg)^1 = χ(k)·g
    # (χg)^2 = χ(k)²·g² = g²
    # (χg)^3 = χ(k)³·g³ = χ(k)·g³
    # (χg)^4 = g⁴
    # (χg)^5 = χ(k)·g⁵

    # (χg-1)^5 = Σ C(5,j)(χg)^j(-1)^{5-j}
    # j=0: (-1)^5 = -1
    # j=1: 5·χg·(-1)^4 = 5χg
    # j=2: 10·g²·(-1)^3 = -10g²
    # j=3: 10·χg³·(-1)^2 = 10χg³
    # j=4: 5·g⁴·(-1)^1 = -5g⁴
    # j=5: χg⁵·1 = χg⁵

    # Sum over k=1..p-1:
    # Σ(χg-1)^5 = (p-1)(-1) + 0·5g + (p-1)(-10g²) + 0·10g³ + (p-1)(-5g⁴) + 0·g⁵
    #            = (p-1)(-1 - 10g² - 5g⁴)

    formula_sum = (pp - 1) * (-1 - 10*g**2 - 5*g**4) / 32

    # Also add λ_0^5 = ((p-1)/2)^5
    lam0 = (pp - 1) / 2
    total_formula = formula_sum + lam0**5
    total_direct = sum(e**5 for e in eigs)

    c5_direct = total_direct.real / 5

    # For g² when p ≡ 3 mod 4: g = i√p, g² = -p
    # For g² when p ≡ 1 mod 4: g = √p, g² = p
    # But actually g² = χ(-1)·p always (standard result)
    chi_neg1 = legendre(-1, pp)  # = (-1)^{(p-1)/2}
    g2_expected = chi_neg1 * pp

    print(f"p={pp:2d}: g²={g**2:.2f} (expect χ(-1)·p = {g2_expected}), "
          f"c_5(Paley)={c5_direct.real:.0f}, "
          f"formula check: {abs(total_formula - total_direct):.2e}")

# ============================================================
#  SECTION 6: Clean formula for c_5(Paley)
# ============================================================
print_separator("SECTION 6: Closed Form for c_5(Paley_p)")

print("Using g² = χ(-1)·p:")
print("  For p ≡ 3 mod 4: g² = -p, g⁴ = p²")
print("  For p ≡ 1 mod 4: g² = p,  g⁴ = p²")
print()
print("Σ_{k≠0} λ_k^5 = (p-1)/32 · (-1 - 10·χ(-1)·p - 5p²)")
print()
print("Total: 5·c_5 = ((p-1)/2)^5 + (p-1)/32·(-1 - 10·χ(-1)·p - 5p²)")
print()

for pp in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
    chi_neg1 = 1 if (pp - 1) // 2 % 2 == 0 else -1
    lam0_5 = ((pp-1)//2)**5
    rest = (pp - 1) * (-1 - 10*chi_neg1*pp - 5*pp**2) // 32
    five_c5 = lam0_5 + rest
    c5 = five_c5 // 5 if five_c5 % 5 == 0 else five_c5 / 5

    # Also compute directly to verify
    qr = quadratic_residues(pp)
    eigs = eigenvalues_circulant(pp, qr)
    c5_direct = round(sum(e**5 for e in eigs).real / 5)

    match = "✓" if c5 == c5_direct else "✗"
    print(f"  p={pp:2d}: χ(-1)={chi_neg1:+d}, c_5 = {c5:8}, direct={c5_direct:8}  {match}")

# ============================================================
#  SECTION 6b: Simplify the formula
# ============================================================
print_separator("SECTION 6b: Simplified c_5(Paley_p) Formula")

print("5·c_5(Paley_p) = ((p-1)/2)^5 + (p-1)/32 · (-1 - 10χ(-1)p - 5p²)")
print()
print("Let's separate by residue class of p mod 4:")
print()

# p ≡ 1 mod 4: χ(-1) = 1
print("Case p ≡ 1 (mod 4):  χ(-1) = +1")
print("  5c_5 = (p-1)^5/32 + (p-1)(-1 - 10p - 5p²)/32")
print("       = (p-1)/32 · [(p-1)^4 - 1 - 10p - 5p²]")
print()
# p ≡ 3 mod 4: χ(-1) = -1
print("Case p ≡ 3 (mod 4):  χ(-1) = -1")
print("  5c_5 = (p-1)^5/32 + (p-1)(-1 + 10p - 5p²)/32")
print("       = (p-1)/32 · [(p-1)^4 - 1 + 10p - 5p²]")

print()
print("Expanding (p-1)^4 = p^4 - 4p^3 + 6p^2 - 4p + 1:")
print()
print("Case p ≡ 1 (mod 4):")
print("  bracket = p^4 - 4p^3 + 6p^2 - 4p + 1 - 1 - 10p - 5p^2")
print("         = p^4 - 4p^3 + p^2 - 14p")
print("         = p(p^3 - 4p^2 + p - 14)")
print()
print("Case p ≡ 3 (mod 4):")
print("  bracket = p^4 - 4p^3 + 6p^2 - 4p + 1 - 1 + 10p - 5p^2")
print("         = p^4 - 4p^3 + p^2 + 6p")
print("         = p(p^3 - 4p^2 + p + 6)")

print()
print("So: 5c_5(Paley_p) = p(p-1)/32 × f(p)")
print("  where f(p) = p³ - 4p² + p - 14  if p ≡ 1 (mod 4)")
print("        f(p) = p³ - 4p² + p + 6   if p ≡ 3 (mod 4)")
print()

for pp in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
    chi_neg1 = 1 if (pp - 1) // 2 % 2 == 0 else -1
    if chi_neg1 == 1:
        f_p = pp**3 - 4*pp**2 + pp - 14
    else:
        f_p = pp**3 - 4*pp**2 + pp + 6
    c5_formula = pp * (pp - 1) * f_p // (32 * 5)
    remainder = (pp * (pp - 1) * f_p) % (32 * 5)

    qr = quadratic_residues(pp)
    eigs = eigenvalues_circulant(pp, qr)
    c5_direct = round(sum(e**5 for e in eigs).real / 5)

    match = "✓" if c5_formula == c5_direct and remainder == 0 else f"✗ (rem={remainder})"
    print(f"  p={pp:2d} ≡ {pp%4} mod 4: f(p)={f_p:6d}, "
          f"c_5 = {c5_formula:8d}, direct = {c5_direct:8d}  {match}")

# ============================================================
#  SECTION 7: Compare Paley vs ALL other circulants
# ============================================================
print_separator("SECTION 7: Paley vs Non-Paley c_5 Ratios at Various p")

for pp in [5, 7, 11, 13]:
    qr = quadratic_residues(pp)
    tours = all_circulant_tournaments(pp)

    # Get Paley c_5
    eigs_paley = eigenvalues_circulant(pp, qr)
    c5_paley = round(cycle_count_from_eigenvalues(eigs_paley, 5).real)
    c3_paley = round(cycle_count_from_eigenvalues(eigs_paley, 3).real)

    print(f"\np = {pp}, QR = {qr}")
    print(f"  Paley: c_3 = {c3_paley}, c_5 = {c5_paley}")

    # Check all others (by isomorphism class)
    seen = set()
    c5_values = {}
    for S in tours:
        S_tuple = tuple(S)
        # Check if we've seen an isomorphic one
        found = False
        for S2 in seen:
            if is_isomorphic_circulant(pp, S, list(S2)):
                found = True
                break
        if found:
            continue
        seen.add(S_tuple)

        eigs = eigenvalues_circulant(pp, S)
        c3 = round(cycle_count_from_eigenvalues(eigs, 3).real)
        c5 = round(cycle_count_from_eigenvalues(eigs, 5).real)
        is_paley = (set(S) == set(qr))
        label = " PALEY" if is_paley else ""

        ratio = c5 / c5_paley if c5_paley != 0 else float('inf')
        print(f"  S={S}: c_3={c3}, c_5={c5}, c_5/c_5(Paley) = {ratio:.4f}{label}")

        if c5 not in c5_values:
            c5_values[c5] = []
        c5_values[c5].append(S)

    if c5_paley > 0:
        min_c5 = min(c5_values.keys())
        print(f"  → c_5(Paley)/c_5(min) = {c5_paley/min_c5:.4f}" if min_c5 > 0 else "  → min c_5 = 0")

# ============================================================
#  SECTION 8: Deeper analysis — why equal magnitudes maximize
# ============================================================
print_separator("SECTION 8: Power-Mean Inequality Explanation")

print("""
WHY Paley maximizes c_5 among circulants:

For ANY circulant tournament on Z_p:
  λ_0 = (p-1)/2  (always)
  Σ_{k=1}^{p-1} |λ_k|² = Σ_{k=1}^{p-1} Σ_{s,t∈S} ω^{k(s-t)}
                         = (p-1)·|S| - |S|²/...

Actually, let's compute Σ|λ_k|² directly via Parseval:
  Σ_{k=0}^{p-1} |λ_k|² = p · |S| = p(p-1)/2

So Σ_{k=1}^{p-1} |λ_k|² = p(p-1)/2 - ((p-1)/2)² = (p-1)(p+1)/4

For Paley: each |λ_k| = √((p+1)/4·(p-1)/(p-1)) — wait let me just verify:
""")

for pp in [5, 7, 11, 13]:
    qr = quadratic_residues(pp)
    eigs = eigenvalues_circulant(pp, qr)
    sum_sq = sum(abs(e)**2 for e in eigs[1:])
    expected = pp * (pp-1) / 2 - ((pp-1)/2)**2
    print(f"  p={pp}: Σ|λ_k|² (k≠0) = {sum_sq:.2f}, expected = {expected:.2f}, "
          f"per term = {sum_sq/(pp-1):.4f}, √(per term) = {(sum_sq/(pp-1))**0.5:.4f}")

print()
print("So Σ_{k≠0} |λ_k|² = (p-1)(p+1)/4, giving average |λ|² = (p+1)/4.")
print("For Paley, EVERY |λ_k|² = (p+1)/4 (Gauss sum: |g|²=p → |λ|²=(p+1)/4).")
print()
print("Wait, let me recheck for p=7:")

eigs7 = eigenvalues_circulant(7, quadratic_residues(7))
for k in range(7):
    print(f"  |λ_{k}|² = {abs(eigs7[k])**2:.4f}")

print()
print("Hmm, |λ_k|² = 2 for k≠0 at p=7, and (p+1)/4 = 2. Correct!")
print()

print("""
The key is that Σ_{k≠0} Re(λ_k^5) is maximized when all phases align constructively
AND all magnitudes are equal (concentrating mass at the "right" phases).

By the power-mean inequality:
  (Σ|λ_k|^5)/(p-1) ≥ ((Σ|λ_k|²)/(p-1))^{5/2}   — WRONG direction!

Actually, the power-mean inequality says:
  (Σ|λ_k|^5/(p-1))^{1/5} ≥ (Σ|λ_k|²/(p-1))^{1/2}

which gives: Σ|λ_k|^5 ≥ (p-1)·((p+1)/4)^{5/2}

with EQUALITY when all |λ_k| are equal, i.e., for Paley.

Wait — this says Paley MINIMIZES Σ|λ_k|^5, not maximizes!
That's because power-mean goes the other direction for p>2.

Let me think again...
""")

# Let's actually check: does Paley maximize or minimize |Σλ^5|?
print("Direct check: which circulant maximizes Σ|λ_k|^5?")
for pp in [5, 7, 11]:
    tours = all_circulant_tournaments(pp)
    qr = quadratic_residues(pp)
    results = []
    for S in tours:
        eigs = eigenvalues_circulant(pp, S)
        sum_mag5 = sum(abs(e)**5 for e in eigs[1:])
        sum_re5 = sum(e**5 for e in eigs).real
        c5 = round(sum_re5 / 5)
        is_paley = (set(S) == set(qr))
        results.append((sum_mag5, sum_re5, c5, S, is_paley))

    results.sort(key=lambda x: x[0])
    print(f"\n  p={pp}: Σ|λ_k|^5 range: {results[0][0]:.2f} to {results[-1][0]:.2f}")

    # Find Paley
    for r in results:
        if r[4]:
            rank = results.index(r) + 1
            print(f"    Paley: Σ|λ|^5 = {r[0]:.2f} (rank {rank}/{len(results)}), "
                  f"Re(Σλ^5) = {r[1]:.2f}, c_5 = {r[2]}")
            break

    # Show by c_5
    results.sort(key=lambda x: x[2], reverse=True)
    print(f"    Max c_5 = {results[0][2]} (Paley? {results[0][4]})")
    print(f"    Min c_5 = {results[-1][2]}")

print("""
KEY FINDING: For Paley (QR sets), equal magnitudes |λ_k| = √((p+1)/4)
MINIMIZE Σ|λ_k|^5 among all circulants (by power-mean inequality).

But Paley MAXIMIZES Re(Σλ^5) = 5·c_5!

This means the PHASE ALIGNMENT is the crucial factor, not the magnitudes.
The Gauss sum phases for QR sets are arithmetically special — they
constructively interfere when raised to odd powers.
""")

# ============================================================
#  SECTION 9: Phase analysis
# ============================================================
print_separator("SECTION 9: Phase Analysis — Why QR Phases Align")

for pp in [7, 11]:
    qr = quadratic_residues(pp)
    tours = all_circulant_tournaments(pp)

    print(f"\np = {pp}:")
    # Get unique classes
    seen = set()
    for S in tours:
        S_tuple = tuple(S)
        found = False
        for S2 in seen:
            if is_isomorphic_circulant(pp, S, list(S2)):
                found = True
                break
        if found:
            continue
        seen.add(S_tuple)

        eigs = eigenvalues_circulant(pp, S)
        is_paley = (set(S) == set(qr))
        label = "PALEY" if is_paley else "      "

        print(f"  {label} S={S}:")
        for k in range(1, pp):
            lam = eigs[k]
            lam5 = lam**5
            phase = cmath.phase(lam)
            phase5 = cmath.phase(lam5)
            print(f"    k={k}: λ phase={phase:+.4f} ({phase/np.pi:.4f}π), "
                  f"λ^5 phase={phase5:+.4f} ({phase5/np.pi:.4f}π), "
                  f"|λ|={abs(lam):.4f}, Re(λ^5)={lam5.real:+.4f}")

        total_re5 = sum(eigs[k]**5 for k in range(1, pp)).real
        print(f"    Σ_{'{k≠0}'} Re(λ^5) = {total_re5:.4f}")

# ============================================================
#  SECTION 10: Jacobi sum connection
# ============================================================
print_separator("SECTION 10: Jacobi Sum Connection")

print("""
For Paley tournaments, λ_k = (χ(k)g - 1)/2 where g = quadratic Gauss sum.

The 5th power sum involves g^5. Since g² = χ(-1)·p:
  g^4 = p², g^5 = p²·g

So (χg-1)^5 = -1 + 5χg - 10χ(-1)p + 10χ(-1)·χ·pg - 5p² + χp²g

Summing over k≠0 (using Σχ(k) = 0):
  Σ(χg-1)^5 = -(p-1) - 10(p-1)χ(-1)p - 5(p-1)p²
            = -(p-1)(1 + 10χ(-1)p + 5p²)

So: c_5(Paley_p) = (1/5)[((p-1)/2)^5 - (p-1)(1 + 10χ(-1)p + 5p²)/32]

This is purely in terms of p and χ(-1) = (-1)^{(p-1)/2}.
No Jacobi sums needed — the fifth power already closes over {1, g, g², g³, g⁴, g⁵}
because g² = χ(-1)·p reduces everything to ℤ[g].
""")

# Verify the closed form
print("Verification of closed form:")
for pp in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    chi_neg1 = (-1)**((pp-1)//2)
    lam0_5 = ((pp-1)//2)**5
    bracket = 1 + 10*chi_neg1*pp + 5*pp**2
    formula = (lam0_5 - (pp-1)*bracket//32)
    if formula % 5 != 0:
        # Need to be more careful with integer arithmetic
        num = 32*lam0_5 - (pp-1)*bracket
        c5 = num // (32*5)
        check = num % (32*5)
    else:
        c5 = formula // 5
        check = 0

    # Direct computation
    qr = quadratic_residues(pp)
    eigs = eigenvalues_circulant(pp, qr)
    c5_dir = round(sum(e**5 for e in eigs).real / 5)

    ok = "✓" if c5 == c5_dir else "✗"
    print(f"  p={pp:2d}: c_5(Paley) = {c5_dir:10d}  formula = {c5:10d}  {ok}")

# ============================================================
#  SECTION 11: The 3/2 ratio question
# ============================================================
print_separator("SECTION 11: Is c_5(Paley)/c_5(other) = 3/2 Special to p=7?")

for pp in [5, 7, 11, 13]:
    tours = all_circulant_tournaments(pp)
    qr = quadratic_residues(pp)

    eigs_paley = eigenvalues_circulant(pp, qr)
    c5_paley = round(cycle_count_from_eigenvalues(eigs_paley, 5).real)

    c5_set = set()
    for S in tours:
        eigs = eigenvalues_circulant(pp, S)
        c5 = round(cycle_count_from_eigenvalues(eigs, 5).real)
        c5_set.add(c5)

    c5_list = sorted(c5_set)
    print(f"p={pp}: c_5 values = {c5_list}")
    print(f"  Paley c_5 = {c5_paley}")
    for v in c5_list:
        if v != c5_paley and v > 0:
            ratio = Fraction(c5_paley, v)
            print(f"  c_5(Paley)/c_5({v}) = {c5_paley}/{v} = {ratio} = {float(ratio):.6f}")

# ============================================================
#  SECTION 12: Higher primes — does Paley always maximize c_5?
# ============================================================
print_separator("SECTION 12: Does Paley Always Maximize c_5 Among Circulants?")

print("For larger primes, enumerating all circulants is expensive.")
print("Instead, sample random circulants and compare.\n")

import random
random.seed(42)

for pp in [11, 13, 17, 19, 23, 29]:
    half = (pp - 1) // 2
    pairs = [(s, pp - s) for s in range(1, half + 1)]
    qr = quadratic_residues(pp)

    eigs_paley = eigenvalues_circulant(pp, qr)
    c5_paley = round(cycle_count_from_eigenvalues(eigs_paley, 5).real)

    n_samples = min(2**half, 500)
    max_c5_other = -float('inf')
    max_S_other = None

    if 2**half <= 500:
        # Enumerate all
        for bits in range(2**half):
            S = []
            for i in range(half):
                if bits & (1 << i):
                    S.append(pairs[i][1])
                else:
                    S.append(pairs[i][0])
            S.sort()
            if set(S) == set(qr):
                continue
            eigs = eigenvalues_circulant(pp, S)
            c5 = round(cycle_count_from_eigenvalues(eigs, 5).real)
            if c5 > max_c5_other:
                max_c5_other = c5
                max_S_other = S
    else:
        for _ in range(n_samples):
            S = []
            for i in range(half):
                if random.random() < 0.5:
                    S.append(pairs[i][1])
                else:
                    S.append(pairs[i][0])
            S.sort()
            eigs = eigenvalues_circulant(pp, S)
            c5 = round(cycle_count_from_eigenvalues(eigs, 5).real)
            if c5 > max_c5_other:
                max_c5_other = c5
                max_S_other = S

    is_max = "PALEY MAX" if c5_paley >= max_c5_other else f"PALEY NOT MAX (other: {max_c5_other})"
    ratio = c5_paley / max_c5_other if max_c5_other > 0 else float('inf')
    method = "exhaustive" if 2**half <= 500 else f"sampled {n_samples}"
    print(f"  p={pp:2d}: c_5(Paley)={c5_paley:8d}, max other={max_c5_other:8d}, "
          f"ratio={ratio:.4f}  [{method}] — {is_max}")

# ============================================================
#  SECTION 13: Summary
# ============================================================
print_separator("SUMMARY")

print("""
WHY PALEY TOURNAMENTS MAXIMIZE 5-CYCLES AMONG Z_p CIRCULANTS
=============================================================

1. EIGENVALUE STRUCTURE
   For any circulant tournament on Z_p with connection set S:
   - λ_k = Σ_{s∈S} ω^{ks} (Gauss-type sum)
   - c_m = (1/m) Σ_{k=0}^{p-1} λ_k^m  (cycle count)

2. WHY c_3 IS UNIVERSAL
   Σ λ_k^3 involves Σ_{k=0}^{p-1} (Σ_{s∈S} ω^{ks})^3
   = Σ_{a,b,c ∈ S} Σ_k ω^{k(a+b+c)}
   = p · |{(a,b,c) ∈ S³ : a+b+c ≡ 0 mod p}|
   This count is the same for ALL tournament connection sets S with |S|=(p-1)/2,
   because it depends only on the additive structure of half-sets of Z_p^*.

3. WHY c_5 IS NOT UNIVERSAL
   The analogous 5-tuple count |{(a,b,c,d,e) ∈ S⁵ : a+b+c+d+e ≡ 0}|
   DOES depend on the specific choice of S.

4. WHY PALEY MAXIMIZES c_5
   For S = QR_p: λ_k = (χ(k)g - 1)/2 where g is the quadratic Gauss sum.

   Key properties:
   (a) |λ_k| = √((p+1)/4) for ALL k≠0 (constant magnitude)
   (b) The phases of λ_k are distributed according to the Legendre symbol
   (c) g² = χ(-1)·p gives a clean algebraic closure under powers

   The PHASE ALIGNMENT (not magnitude equality) is the critical factor.
   Equal magnitudes actually MINIMIZE Σ|λ|^5 by power-mean.
   But the QR phases are arithmetically special: when raised to the 5th power,
   the real parts constructively interfere.

5. CLOSED FORM
   c_5(Paley_p) = [32·((p-1)/2)^5 - (p-1)(1 + 10χ(-1)p + 5p²)] / 160

   where χ(-1) = (-1)^{(p-1)/2}.

6. THE 3/2 RATIO
   c_5(Paley_7)/c_5(other_7) = 42/28 = 3/2 is SPECIFIC to p=7.
   The ratio varies with p and is NOT always 3/2.
   However, Paley appears to ALWAYS give the maximum c_5 among circulants.
""")
