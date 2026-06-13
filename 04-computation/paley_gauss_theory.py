#!/usr/bin/env python3
"""
THEORETICAL ANALYSIS: GAUSS SUMS AND PALEY PATH HOMOLOGY

For P_p (p ≡ 3 mod 4 prime), the Fourier eigenspaces decompose
path homology into p pieces. All non-trivial eigenspaces (k=1,...,p-1)
give IDENTICAL Ω dimensions (verified for p=7).

This means the topology is controlled by the COMMON structure of
all non-trivial eigenspaces, which is related to Gauss sums.

KEY QUESTION: What dimension d has β_d^(k) = 1 for k ≠ 0?

For P_7: d = p - 3 = 4.
Conjecture: For P_p, d = p - 3?

The Ω dims for P_7 are [1, 3, 6, 9, 9, 6, 3] for each eigenspace.
These are |A_m| for the step-sequence space. Since all eigenspaces
have the same Ω dims, the junk kernel dimension is the same at each λ.

Why? Because for Paley, the Legendre symbol χ(k) determines the
character sum Σ_{s∈QR} λ^s, and the junk matrix J_m(λ) has entries
that are sums of powers of λ.

For QR connection sets, the multiplicative structure of F_p* implies
that the junk matrices for different eigenvalues are related by the
Legendre symbol. Specifically, if χ(k) = χ(k'), then the junk
matrices are identical up to conjugation.

Since there are only TWO values of χ (±1), and p ≡ 3 mod 4 means
both QR and QNR have the same size (p-1)/2, the eigenspaces split
into two groups of equal size — but the Ω dims are the SAME for both!

This means the Legendre symbol DOESN'T affect Ω dims.
Why? Because QR is closed under multiplication, so the map
s → ks (for k ∈ QR) permutes S = QR, giving an isomorphism between eigenspaces.

And for k ∈ QNR? The map s → ks sends S to QNR. But Q + QNR structure...

Let's investigate computationally.
"""
import numpy as np
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

def qr(p):
    return set((a*a) % p for a in range(1, p))

def legendre(a, p):
    """Legendre symbol (a/p)."""
    if a % p == 0: return 0
    return 1 if (a % p) in qr(p) else -1

# ===== Gauss sum analysis =====
print("=" * 70)
print("GAUSS SUM ANALYSIS FOR PALEY PRIMES")
print("=" * 70)

for p in [3, 7, 11, 19, 23, 31, 43, 47]:
    # Check p ≡ 3 mod 4 and prime
    if p % 4 != 3: continue
    is_prime = all(p % d != 0 for d in range(2, int(p**0.5)+1))
    if not is_prime: continue

    S = qr(p)
    omega = np.exp(2j * np.pi / p)

    # Gauss sum
    g = sum(legendre(a, p) * omega**a for a in range(1, p))

    # Character sums for each k
    char_sums = []
    for k in range(p):
        cs = sum(omega**(k*s) for s in S)
        char_sums.append(cs)

    # Group by Legendre symbol value
    qr_sums = [char_sums[k] for k in range(1, p) if legendre(k, p) == 1]
    qnr_sums = [char_sums[k] for k in range(1, p) if legendre(k, p) == -1]

    print(f"\nP_{p}: |QR|={len(S)}")
    print(f"  Gauss sum g = {g.real:.4f} + {g.imag:.4f}i")
    print(f"  |g| = {abs(g):.4f}, √p = {p**0.5:.4f}")
    print(f"  g² = {(g*g).real:.1f} (expected {(-1)**((p-1)//2) * p} = -p since p≡3 mod 4)")

    # Are all QR character sums the same? All QNR the same?
    qr_vals = set(round(cs.real, 6) + 1j*round(cs.imag, 6) for cs in qr_sums)
    qnr_vals = set(round(cs.real, 6) + 1j*round(cs.imag, 6) for cs in qnr_sums)

    print(f"  QR char sums: {len(qr_vals)} distinct values (out of {len(qr_sums)})")
    print(f"  QNR char sums: {len(qnr_vals)} distinct values (out of {len(qnr_sums)})")

    if len(qr_sums) > 0:
        print(f"  QR[0] = {qr_sums[0]:.4f}")
        print(f"  QNR[0] = {qnr_sums[0]:.4f}")

    # Verify: Σ_{s∈QR} ω^{ks} = (χ(k)·g - 1)/2
    print(f"  Formula verification: Σ_QR ω^(ks) = (χ(k)·g - 1)/2:")
    for k in range(1, min(p, 4)):
        actual = char_sums[k]
        expected = (legendre(k, p) * g - 1) / 2
        match = abs(actual - expected) < 1e-8
        print(f"    k={k}: χ(k)={legendre(k,p)}, actual={actual:.4f}, formula={expected:.4f}, match={match}")

# ===== Why do all eigenspaces give the same Ω dims? =====
print("\n\n" + "=" * 70)
print("WHY ALL EIGENSPACES GIVE SAME Ω DIMS")
print("=" * 70)

print("""
For C_p^S with S = QR(p), the junk matrix J_m(λ_k) depends on λ_k = ω^k.

The entries of J_m involve sums of the form λ^shift for various shifts.
The shifts come from merged steps: (s_i + s_{i+1}) mod p for s_i, s_{i+1} ∈ S.

KEY OBSERVATION: For k ∈ QR, the map φ_k: s → ks is a bijection on QR.
So φ_k permutes the step sequences, and the junk matrix J_m(ω^k) is
obtained from J_m(ω) by permuting rows and columns.
Hence ker(J_m(ω^k)) has the same dimension as ker(J_m(ω)).

For k ∈ QNR: φ_k sends QR → QNR. But the step sequences use steps from QR,
so this doesn't directly permute them. However, the junk face types
involve merged steps in (QR+QR)∖QR, and for Paley, (QR+QR)∖QR = QNR.

The map φ_k for k ∈ QNR sends:
  QR → QNR (steps change)
  QNR → QR (junk types change)

This swaps the roles of steps and junk types, but the STRUCTURE of J_m
is preserved because the additive structure of F_p respects this swap.

More precisely: multiplication by k ∈ F_p* is an automorphism of (Z/pZ, +),
and it maps QR-step-sequences bijectively to QNR-step-sequences.
But QNR-step-sequences would be step sequences for the COMPLEMENT tournament!

Since P_p is self-complementary (P_p ≅ P_p^op), the complement has the
same path homology. So the QNR eigenspaces should give the same Ω dims.

This is why ALL eigenspaces give the same Ω dims: QR eigenspaces are
related by QR-multiplication (permutation), and QNR eigenspaces are
related to QR eigenspaces by QNR-multiplication (complement swap +
self-complementarity).
""")

# ===== Predict the topology dimension for general p =====
print("=" * 70)
print("TOPOLOGY DIMENSION PREDICTION")
print("=" * 70)

print("""
For P_p, the Ω dims are the same at all eigenspaces.
The total Betti numbers are β_d^total = p · β_d^(k=1) for d > 0
(since k=0 contributes only to β_0), or more precisely:
  β_0 = 1 (from k=0)
  β_d = (p-1) · β_d^(k=1) for d > 0

Wait, this assumes β_0^(k) = 0 for k ≠ 0. Let me verify.

For P_7:
  k=0: β = [1,0,0,0,0,0,0]
  k≠0: β = [0,0,0,0,1,0,0]
  Total: β = [1,0,0,0,6,0,0]

So β_0^(k=0) = 1, β_0^(k≠0) = 0.
And β_4^(k≠0) = 1 for each of 6 eigenspaces, giving β_4 = 6 = p-1.

Prediction for P_11:
  If the pattern holds, β_d = 10 for some d, and all other β = 0 (except β_0=1).
  What d? For P_7, d = p-3 = 4. So d = 11-3 = 8?

Alternative: d = (p-1)/2 - 1 = (p-3)/2.
  P_7: d = 2 (no, it's d=4)
  P_3: d = 0 (no, β_1=1 for k≠0, but wait P_3 has β=(1,1,0))

Let me reconsider. For P_3: β = (1,1,0). That's β_0=1, β_1=1.
  k=0: β=[1,0,0], k=1: β=[0,1,0], k=2: β=[0,1,0]? No...
  P_3 has n=3, so only k=0,1,2. If β_1 = 1, then either one eigenspace
  contributes β_1=1 or... let me check.

Actually P_3: β = (1,1,0), so β_1 = 1, not p-1=2.
So the P_3 → P_7 jump is: β_1=1 → β_4=6.

P_3 is special: it's just the 3-cycle, which is topologically S^1.
""")

# Let me also check: for P_3, what is the per-eigenspace breakdown?
p = 3
S = qr(p)
print(f"P_3: QR = {sorted(S)}")

# P_3 is C_3^{1}, the directed 3-cycle
# A_0 = 1, A_1 = 1 step sequence (just (1,)), A_2 = 0
# Eigenspaces: ω = e^{2πi/3}

omega = np.exp(2j * np.pi / 3)
for k in range(3):
    lam = omega ** k
    # M_1(λ) = [λ - 1] (since S = {1})
    m1 = lam - 1
    print(f"  k={k}: λ={lam:.4f}, M_1 = [{m1:.4f}], |M_1| = {abs(m1):.4f}")
    # Ω_0 = 1 (trivially)
    # Ω_1: only step (1,), boundary ∂(1) = λ·() - () = (λ-1)·()
    # Junk: none (all faces are valid)
    # So Ω_1 = A_1 = 1
    # ker(∂_1) = 1 if λ=1, else 0
    # β_0 = Ω_0 - rank(∂_1) = 1 - (1 if λ≠1 else 0) = (1 if λ=1 else 0)
    # Wait, that's wrong. Let me think again.
    # ∂_1: Ω_1 → Ω_0 given by ∂_1(s) = λ^s · () - ()
    # For S={1}: ∂_1 = λ - 1
    # rank = 1 if λ≠1, else 0
    # β_0^(k) = dim(ker ∂_1 on Ω_0) / im(∂_2)
    #         = Ω_0 - rank(∂_1) = 1 - (1 if λ≠1 else 0)
    # Wait, we need to think in terms of the chain complex Ω_1 → Ω_0
    # β_0 = ker(∂_0)/im(∂_1). ∂_0: Ω_0 → 0, so ker(∂_0) = Ω_0 = 1.
    # im(∂_1) = rank of ∂_1: Ω_1 → Ω_0 = rank(λ-1) = (1 if λ≠1 else 0)
    # So β_0^(k=0) = 1 - 0 = 1, β_0^(k≠0) = 1 - 1 = 0.
    # β_1 = ker(∂_1)/im(∂_2). Since A_2 = 0, Ω_2 = 0, im(∂_2) = 0.
    # ker(∂_1): if λ=1, ∂_1=0, ker=1. If λ≠1, ∂_1≠0, ker=0.
    # So β_1^(k=0) = 1, β_1^(k≠0) = 0. Total β_1 = 1. ✓

    # Wait, that says β_1 = 1 from k=0, not from k≠0!
    # And β_0 = 1 from k=0 too. So χ = 1-1 = 0 from k=0, and 0 from k≠0.
    # Total χ = 0. ✓ (1-1=0)

    # Hmm, so P_3 is different from P_7: the β_1 comes from k=0.

print("""
P_3 analysis: β_1 = 1 from k=0 eigenspace (NOT from k≠0).
This is because S = {1} has only 1 generator, so A_1 has only 1 step sequence.
The Ω_1 = A_1 = {(1)} for each eigenspace, but:
  k=0: λ=1, ∂_1 = 0, so ker(∂_1)=1, β_1=1
  k≠0: λ≠1, ∂_1 = λ-1 ≠ 0, so ker(∂_1)=0, β_1=0

So P_3 is fundamentally different from P_7.
At P_7, the β comes from k≠0 eigenspaces.
At P_3, the β comes from k=0.

The transition happens because |S| grows: P_3 has |S|=1, P_7 has |S|=3.
With more generators, higher-dimensional step sequences proliferate,
and the non-trivial eigenspaces develop nontrivial homology.
""")

# ===== Dimension formula: where does β appear? =====
print("=" * 70)
print("DIMENSION FORMULA INVESTIGATION")
print("=" * 70)

# For P_7: Ω dims = [1,3,6,9,9,6,3]. These are:
# dim 0: 1 = C(3,0)? No, 1.
# dim 1: 3 = |S| = 3
# dim 2: 6 = A_2 - rank(J_2)
# dim 3: 9 = ?
# dim 4: 9 = ?
# dim 5: 6 = ?
# dim 6: 3 = A_6 = 27, Ω_6 = 3

# The pattern [1,3,6,9,9,6,3] sums to 1+3+6+9+9+6+3 = 37
# χ(Ω^(k)) = 1-3+6-9+9-6+3 = 1 for each k.
# Total χ = 7·1 = 7 = p ✓

# The Ω dims look like they might be related to Pascal's triangle or
# some combinatorial object.

# 1, 3, 6, 9, 9, 6, 3
# C(3,0), C(3,1)?, not C(3,2)=3≠6

# Actually: 1, 3, 6, 9... looks like partial sums:
# 1, 1+2=3, 1+2+3=6, 1+2+3+3=9? No.
# Or: C(3+k-1, k) = C(2+k, k): 1, 3, 6, 10... no, 10≠9.

# Let me think about this differently.
p = 7
dims = [1, 3, 6, 9, 9, 6, 3]
print(f"P_7 Ω dims: {dims}")
print(f"  Sum = {sum(dims)}")
print(f"  Alt sum = {sum((-1)**i * dims[i] for i in range(len(dims)))}")
print(f"  Ratios: {[dims[i+1]/dims[i] if dims[i]>0 else 'inf' for i in range(len(dims)-1)]}")

# The boundary ranks for P_7 (one eigenspace):
# β = [0,0,0,0,1,0,0] at k≠0
# So: Ω dims = [1,3,6,9,9,6,3]
# ∂ ranks: r_1, r_2, r_3, r_4, r_5, r_6
# β_0 = 1 - r_1 = 0 → r_1 = 1
# β_1 = 3 - r_1 - r_2 = 0 → r_2 = 3 - 1 = 2
# β_2 = 6 - r_2 - r_3 = 0 → r_3 = 6 - 2 = 4
# β_3 = 9 - r_3 - r_4 = 0 → r_4 = 9 - 4 = 5
# β_4 = 9 - r_4 - r_5 = 1 → r_5 = 9 - 5 - 1 = 3
# β_5 = 6 - r_5 - r_6 = 0 → r_6 = 6 - 3 = 3
# β_6 = 3 - r_6 = 0 → r_6 = 3 ✓

print(f"\n  Boundary ranks (per eigenspace k≠0):")
ranks = [1, 2, 4, 5, 3, 3]
print(f"  r = {ranks}")
print(f"  Verify: dims = [{dims[0]}]", end="")
for i in range(len(ranks)):
    bd = dims[i+1] - ranks[i]
    print(f", [{dims[i+1]}-{ranks[i]}={bd}]", end="")
print()

# These ranks are: 1, 2, 4, 5, 3, 3
# Alternately: the map ∂_m has rank r_m, and the "excess" dim(ker)-r_{m+1} = β_m.

# For β_4 = 1: ker(∂_4) = 9 - 5 = 4, im(∂_5) = 3, so β_4 = 4 - 3 = 1.

# Fascinating: the ONE non-killed dimension at β_4 is the "mismatch"
# between kernel and image at that level.

print(f"""
KEY PATTERN:
  Ω_m dims are symmetric around the midpoint (palindromic-ish):
  [1, 3, 6, 9 | 9, 6, 3]

  The boundary ranks are NOT symmetric:
  [1, 2, 4, 5 | 3, 3]

  The asymmetry at the "equator" (dim 3-4) creates β_4 = 1.

  If Ω dims were perfectly symmetric and boundary ranks
  were determined by exact sequences, we'd need Poincaré duality.
  But we DON'T have a manifold — this is a directed graph!

  The asymmetry of boundary ranks at the middle dimension is what
  creates the non-trivial homology.
""")

# ===== What are the Ω dims for P_11? =====
# We know |A_m| for P_11:
# [1, 5, 25, 110, 430, 1430, 3970, 8735, 14395, 15745, 8645]
# The Ω dims depend on junk kernel dimensions.

# Predict: if the pattern matches P_7, the Ω dims should peak in the middle
# and be roughly palindromic.

# For P_7: |A_m|/Ω_m ratios:
a7 = [1, 3, 9, 21, 39, 45, 27]
o7 = [1, 3, 6, 9, 9, 6, 3]
print(f"P_7: |A_m|/Ω_m = {[f'{a/o:.1f}' for a,o in zip(a7, o7)]}")
print(f"P_7: Ω_m/|A_m| = {[f'{o/a:.3f}' for a,o in zip(a7, o7)]}")

print("\nDone.")
