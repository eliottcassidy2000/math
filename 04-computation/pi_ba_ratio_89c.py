#!/usr/bin/env python3
"""
pi_ba_ratio_89c.py — The b/a ratio and Gauss sum connection
opus-2026-03-14-S89c

b/a = (NQR endpoint count) / (QR endpoint count) per vertex

Data:
  p=3: a=0, b=1 (b/a = ∞)
  p=7: a=1, b=8 (b/a = 8)
  p=11: a=628, b=1101 (b/a ≈ 1.753)
  p=19: a=3040229428, b=3817640437 (b/a ≈ 1.256)

The Gauss sum g = Σ (a/p) e^{2πia/p} satisfies g² = -p (for p ≡ 3 mod 4).
So g = i√p.

HYPOTHESIS: b/a = 1 + f(p)/g for some function f, giving b/a → 1 as p → ∞.

Or more precisely: (b-a)/(a+b) involves 1/√p.
"""

from fractions import Fraction
from math import sqrt, pi, log

# Known data
data = {
    7: (1, 8),
    11: (628, 1101),
    19: (3040229428, 3817640437),
}

# We need a and b for p=23. Let me compute it.
def paley_tournament(p):
    qr = set()
    for a in range(1, p):
        qr.add((a*a) % p)
    adj = [[] for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                adj[i].append(j)
    return adj, qr

# p=23 computation
p = 23
adj, qr = paley_tournament(p)
nqr = set(range(1, p)) - qr

dp = [dict() for _ in range(1 << p)]
dp[1][0] = 1
import time
t0 = time.time()
for mask in range(1, 1 << p):
    for v in dp[mask]:
        if dp[mask][v] == 0:
            continue
        for u in adj[v]:
            if mask & (1 << u) == 0:
                new_mask = mask | (1 << u)
                dp[new_mask][u] = dp[new_mask].get(u, 0) + dp[mask][v]

full = (1 << p) - 1
a23 = dp[full].get(min(qr), 0)
b23 = dp[full].get(min(nqr), 0)
h0_23 = sum(dp[full].values())
t1 = time.time()
print(f"P_23: a={a23}, b={b23}, a+b={a23+b23}, time={t1-t0:.1f}s")
print(f"  h₀ = {h0_23}, d×(a+b) = {11*(a23+b23)}")
print(f"  Match: {'✓' if h0_23 == 11*(a23+b23) else '✗'}")

# Verify all within-class values are equal
qr_vals = set(dp[full].get(v, 0) for v in qr)
nqr_vals = set(dp[full].get(v, 0) for v in nqr)
print(f"  QR all equal: {'✓' if len(qr_vals) == 1 else '✗'}")
print(f"  NQR all equal: {'✓' if len(nqr_vals) == 1 else '✗'}")

data[23] = (a23, b23)

print()
print("=" * 70)
print("PART 1: The b/a ratio analysis")
print("=" * 70)

for p in sorted(data.keys()):
    a, b = data[p]
    print(f"\n  p={p}:")
    print(f"    a = {a}, b = {b}")
    print(f"    b/a = {Fraction(b, a) if a > 0 else 'inf'} = {b/a if a > 0 else float('inf'):.8f}")
    print(f"    (b-a)/(b+a) = {Fraction(b-a, b+a)} = {(b-a)/(b+a):.8f}")
    print(f"    (b-a)/(b+a) × √p = {(b-a)/(b+a) * sqrt(p):.8f}")

# Check if (b-a)/(b+a) ≈ C/√p
print(f"\n  Fitting (b-a)/(b+a) ≈ C/√p:")
for p in sorted(data.keys()):
    a, b = data[p]
    delta = (b-a)/(b+a)
    C = delta * sqrt(p)
    print(f"    p={p}: C = {C:.6f}")

# C values: try to see if they converge
# The fact that (b-a)/(b+a) × √p varies suggests the correction is not purely 1/√p.

print()
print("=" * 70)
print("PART 2: Character sum interpretation")
print("=" * 70)

# b - a = (NQR end count per vertex) - (QR end count per vertex)
# b + a = total per vertex pair = (a+b)
#
# (b-a) counts paths from 0 weighted by +1 if ending at NQR, -1 if ending at QR.
# Since NQR vertices v have (v/p) = -1 and QR vertices have (v/p) = +1,
# this is related to: Σ_{v≠0} -(v/p) × h(0→v)
#
# where h(0→v) = number of Hamiltonian paths from 0 to v.
#
# So b - a = -Σ_{v≠0} (v/p) × h(0→v) / d  (dividing by class size d = (p-1)/2)
# Wait: total QR endpoints = d × a, total NQR endpoints = d × b.
# Σ_{v ∈ QR} h(0→v) = d × a
# Σ_{v ∈ NQR} h(0→v) = d × b
# Σ_{v≠0} h(0→v) = d × (a + b) = h₀
# Σ_{v≠0} (v/p) × h(0→v) = d × a - d × b = d(a - b)
#
# So the "character-weighted" path count is:
# χ_path := Σ_{v≠0} (v/p) × h(0→v) = d(a-b)
#
# And (a-b) = χ_path / d
# (b-a)/(b+a) = -(a-b)/(a+b) = χ_path / (d × (a+b)) = χ_path / h₀

# This χ_path is a natural "spectral" quantity!

for p in sorted(data.keys()):
    a, b = data[p]
    d = (p-1)//2
    chi_path = d * (a - b)  # = Σ (v/p) h(0→v)
    print(f"\n  p={p}:")
    print(f"    χ_path = d(a-b) = {chi_path}")
    print(f"    h₀ = {d * (a + b)}")
    print(f"    χ_path / h₀ = {Fraction(chi_path, d*(a+b))} = {chi_path / (d*(a+b)):.8f}")
    # This should be (a-b)/(a+b)
    print(f"    (a-b)/(a+b) = {Fraction(a-b, a+b)}")

print()
print("=" * 70)
print("PART 3: χ_path and the Gauss sum")
print("=" * 70)

# The adjacency matrix A of P_p in the DFT basis:
# A = F^* diag(λ_0, λ_1, ..., λ_{p-1}) F
# where F_{jk} = ω^{jk}/√p, ω = e^{2πi/p}
#
# The "Legendre character" χ(v) = (v/p) is a DFT eigenvector!
# Specifically: Σ_{v} (v/p) ω^{vk} = (k/p) × g
# where g = Σ_v (v/p) ω^v is the Gauss sum.
#
# So the character acts like a "projection" onto the Gauss sum direction.
#
# χ_path = Σ_v (v/p) h(0→v) relates to the component of h(0,·)
# in the Legendre character direction.

# For a circulant tournament: h(0, v) depends only on v (mod symmetry).
# But h(0, v) is NOT circulant-invariant! (It depends on the full path structure.)
# However, we KNOW that h(0, v) = a for v ∈ QR and b for v ∈ NQR.
# So h(0, ·) is a function of the Legendre symbol:
# h(0, v) = (a+b)/2 + (a-b)/2 × (v/p) for v ≠ 0
# (since (v/p) = +1 for QR and -1 for NQR)

# Wait, that's not quite right. We need to also handle v=0.
# For v=0: h(0,0) = 0 (can't start and end at same vertex in a Ham path
# of length p-1... actually you CAN if you return to 0, but in a
# Ham PATH from 0, the last vertex is different from 0.)

# So h(0, v) = α + β × (v/p) for v ≠ 0
# where α = (a+b)/2 and β = (a-b)/2

# The DFT of h(0, ·):
# ĥ(0, k) = Σ_{v=0}^{p-1} h(0, v) ω^{-vk}
# = α × Σ_{v≠0} ω^{-vk} + β × Σ_{v≠0} (v/p) ω^{-vk}
# = α × (-1 + Σ_v ω^{-vk}) + β × (k/p) × g*   [if k ≠ 0]
# = α × (-1 + p δ_{k,0}) + β × (k/p) × ḡ

# For k = 0: ĥ(0, 0) = Σ h(0,v) = h₀ = d(a+b)
# For k ≠ 0: ĥ(0, k) = -α + β × (k/p) × ḡ
# where ḡ = conjugate Gauss sum = (−1/p) g = -g (for p ≡ 3 mod 4)
# So ḡ = -g

# ĥ(0, k) = -α - β × (k/p) × g   for k ≠ 0
# = -(a+b)/2 - (a-b)/2 × (k/p) × g

# Since g = i√p (for p ≡ 3 mod 4):
# ĥ(0, k) = -(a+b)/2 - (a-b)/2 × (k/p) × i√p

print("\n  The DFT of h(0, ·):")
print("  ĥ(0, 0) = h₀ = d(a+b)")
print("  ĥ(0, k) = -(a+b)/2 - (a-b)/2 × (k/p) × i√p    for k ≠ 0")

for p in sorted(data.keys()):
    a, b = data[p]
    d = (p-1)//2
    alpha = Fraction(a+b, 2)
    beta = Fraction(a-b, 2)

    print(f"\n  p={p}: α = (a+b)/2 = {alpha}, β = (a-b)/2 = {beta}")
    print(f"    ĥ(0, k) = {-alpha} + {-beta}×(k/p)×i√{p}")
    print(f"    |ĥ(0, k)|² = {float(alpha**2 + beta**2 * p):.4f}")
    print(f"    = α² + β²p = {alpha**2} + {beta**2 * p} = {alpha**2 + beta**2 * p}")

print()
print("=" * 70)
print("PART 4: Parseval check")
print("=" * 70)

# Parseval: Σ_k |ĥ(0,k)|² = p × Σ_v |h(0,v)|²
# LHS: |h₀|² + (p-1) × (α² + β²p)
# = d²(a+b)² + (p-1)(α² + β²p)
# RHS: p × [d × a² + d × b²]
# = pd(a² + b²)

for p in sorted(data.keys()):
    a, b = data[p]
    d = (p-1)//2
    alpha = (a+b)/2
    beta = (a-b)/2

    lhs = d**2 * (a+b)**2 + (p-1) * (alpha**2 + beta**2 * p)
    rhs = p * d * (a**2 + b**2)
    print(f"  p={p}: LHS = {lhs:.0f}, RHS = {rhs:.0f}, match: {'✓' if abs(lhs-rhs) < 1 else '✗'}")

print()
print("=" * 70)
print("PART 5: What constrains a and b?")
print("=" * 70)

# We have:
# 1. a + b = H / (pd) (known)
# 2. Parseval gives a constraint on a² + b²
# 3. Are there more constraints?
#
# From the Parseval relation:
# d²(a+b)² + (p-1)((a+b)²/4 + (a-b)²p/4) = pd(a²+b²)
# d²(a+b)² + (p-1)(a+b)²/4 + (p-1)(a-b)²p/4 = pd(a²+b²)
#
# Let s = a+b, t = a-b. Then a = (s+t)/2, b = (s-t)/2.
# a² + b² = (s² + t²)/2
#
# d²s² + (p-1)s²/4 + (p-1)pt²/4 = pd(s²+t²)/2
#
# d²s² + (p-1)s²/4 + p(p-1)t²/4 = pds²/2 + pdt²/2
#
# With d = (p-1)/2:
# (p-1)²s²/4 + (p-1)s²/4 + p(p-1)t²/4 = p(p-1)s²/4 + p(p-1)t²/4
#
# LHS: [(p-1)² + (p-1)]s²/4 + p(p-1)t²/4 = (p-1)p s²/4 + p(p-1)t²/4
# RHS: p(p-1)s²/4 + p(p-1)t²/4
#
# LHS = RHS identically! So Parseval is automatically satisfied.
# No extra constraint from Parseval.

print("  Parseval is automatically satisfied — gives NO extra constraint on a, b")
print("  (This is because h(0,·) is a 2-level function of the Legendre symbol)")

print()
print("=" * 70)
print("PART 6: Ratio b/a as a function of p — limit analysis")
print("=" * 70)

# b/a values:
# p=7: 8/1 = 8
# p=11: 1101/628 ≈ 1.7531
# p=19: 3817640437/3040229428 ≈ 1.2558
# p=23: b23/a23

print("\n  b/a ratios:")
for p in sorted(data.keys()):
    a, b = data[p]
    if a > 0:
        print(f"    p={p}: b/a = {b/a:.8f}")
        print(f"      (b/a - 1) = {b/a - 1:.8f}")
        print(f"      (b/a - 1) × √p = {(b/a - 1) * sqrt(p):.8f}")
        print(f"      (b/a - 1) × p = {(b/a - 1) * p:.8f}")

# (b/a - 1) × √p:
# p=7: 7 × 2.646 = 18.52
# p=11: 0.753 × 3.317 = 2.497
# p=19: 0.256 × 4.359 = 1.115
# p=23: check

# Not clean. Try (b/a - 1) × p:
# p=7: 7 × 7 = 49
# p=11: 0.753 × 11 = 8.28
# p=19: 0.256 × 19 = 4.86
# p=23: check

# Not clean either. Try log(b/a) / log(p):
print(f"\n  ln(b/a) and scaling:")
for p in sorted(data.keys()):
    a, b = data[p]
    if a > 0 and b > a:
        print(f"    p={p}: ln(b/a) = {log(b/a):.6f}, ln(b/a)/ln(p) = {log(b/a)/log(p):.6f}")

# Hmm, the b/a ratio is converging to 1 but the rate is unclear.
# Maybe it's 1 + C/(p-1) for some constant C?
print(f"\n  Fit b/a = 1 + C/(p-1):")
for p in sorted(data.keys()):
    a, b = data[p]
    if a > 0:
        C = (b/a - 1) * (p - 1)
        print(f"    p={p}: C = {C:.6f}")

# C values: 42, 7.53, 4.60, 3.62
# Still decreasing. Try b/a = 1 + C/p^α and fit α.

print()
print("=" * 70)
print("PART 7: The (a-b) sequence")
print("=" * 70)

# b - a gives the "character excess"
for p in sorted(data.keys()):
    a, b = data[p]
    print(f"  p={p}: b-a = {b-a}")
    if b - a != 0:
        import sympy
        print(f"    factorization: {sympy.factorint(abs(b-a))}")
    print(f"    (b-a)/(a+b) = {Fraction(b-a, a+b)}")

# b-a values: 7, 473, 777411009, ?
# 7 = 7
# 473 = 11 × 43
# 777411009 = ?

print()
print("=" * 70)
print("PART 8: Connection to directed Hamiltonian cycles")
print("=" * 70)

# From THM-214: directed Ham cycles ≡ (p-1)/2 mod p
# And hc = number of directed Ham cycles
# The ending NQR bias in paths should relate to the cycle structure.

# For Ham CYCLES: by symmetry, each vertex is equally likely to be
# "adjacent to vertex 0" in the cycle. But the cycle has two neighbors
# of 0: one predecessor and one successor. The predecessor has an
# arc TO 0 (so pred→0 means 0-pred ∈ QR, i.e., -pred ∈ QR iff pred ∈ NQR
# since -1 is NQR for p ≡ 3 mod 4).
# The successor has 0→succ, meaning succ ∈ QR.
#
# So in a Ham cycle through 0:
# - The successor of 0 is always a QR vertex
# - The predecessor of 0 is always a NQR vertex
#
# For Ham PATHS from 0: the first vertex after 0 is always QR (since 0→v requires v ∈ QR).
# The last vertex is preferentially NQR for small p but approaches balance for large p.

print("  In a Ham cycle through 0:")
print("    successor of 0 is always QR (0→succ ⟹ succ ∈ QR)")
print("    predecessor of 0 is always NQR (pred→0 ⟹ -pred ∈ QR ⟹ pred ∈ NQR)")
print()
print("  In a Ham path from 0:")
print("    first step goes to QR (0→v₁ ⟹ v₁ ∈ QR)")
print("    last step: NQR bias decreases as p → ∞")
print()
print("  The NQR excess in path endpoints comes from the")
print("  ASYMMETRIC first step: we start by entering QR territory,")
print("  which means we're more likely to END in NQR territory")
print("  (since we've 'used up' more QR vertices early on).")

print()
print("=" * 70)
print("SUMMARY — The QR/NQR Decomposition")
print("=" * 70)

print("""
  H(P_p) = p × (p-1)/2 × (a + b)

  where a = paths from 0 to each QR vertex
        b = paths from 0 to each NQR vertex

  ╔═══╦═══════════╦══════════════╦═════════╦══════════╗
  ║ p ║     a     ║      b       ║   a+b   ║   b/a    ║
  ╠═══╬═══════════╬══════════════╬═════════╬══════════╣
  ║ 7 ║         1 ║            8 ║       9 ║ 8.000000 ║
  ║11 ║       628 ║        1,101 ║   1,729 ║ 1.753185 ║
  ║19 ║3040229428 ║ 3,817,640,437║         ║ 1.255807 ║
  ╚═══╩═══════════╩══════════════╩═════════╩══════════╝

  The b/a ratio → 1 as p → ∞ (QR/NQR become balanced)

  CHARACTER INTERPRETATION:
  h(0, v) = α + β × (v/p)  where α = (a+b)/2, β = (a-b)/2
  The DFT: ĥ(0, k) = -(a+b)/2 - (a-b)/2 × (k/p) × i√p
""")

print("Done!")
