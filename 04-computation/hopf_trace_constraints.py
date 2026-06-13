"""
Hopf trace formula constraints on Paley tournament homology.

Z_m = QR(p)* acts on diff-seqs by multiplication: g_q(s_1,...,s_d) = (qs_1,...,qs_d).
This action is FREE on d-diff-seqs for d >= 1 (since s_i != 0).

The Lefschetz number L(g) = sum (-1)^d tr(g | H_d) for g in Z_m, g != id.
By the Lefschetz fixed-point theorem for chain complexes with free action:
  L(g) = chi(fixed-point set of g)

For a free action: fixed-point set is empty for d >= 1, so:
  L(g) = chi(fixed) = 1 (just the vertex 0 is fixed... wait)

Actually for diff-seqs: the action g_q acts on the CHAIN COMPLEX, not on a space.
The Lefschetz number is:
  L(g_q) = sum_{d=0}^{2m} (-1)^d tr(g_q | Omega_d)

Since g_q acts freely on A_d for d >= 1 (all diff-seq entries are nonzero QR):
  Every orbit has size m = |QR*|/(... wait, need to think about this more carefully)

Actually QR* is the multiplicative group of QR under multiplication mod p.
|QR*| = m since QR has m elements and forms a group under multiplication.
Wait: QR is closed under multiplication (product of QR = QR).

The action of q in QR on diff-seqs sends (s_1,...,s_d) -> (qs_1,...,qs_d).
This is free because if qs_i = s_i for all i, then q = 1.

The trace of g_q on the chain space C_d = span(A_d) is:
  tr(g_q | C_d) = #{sigma in A_d : g_q(sigma) = sigma} = 0 for d >= 1
  (since the action is free).

But we need tr(g_q | Omega_d), not tr on all of C_d.
Since g_q preserves the constraint equations, it maps Omega_d to Omega_d.

The trace on Omega_d:
  Omega_d = ker(C_d : A_d -> junk_space)
  g_q acts on both A_d and junk_space, commuting with C_d.

So: tr(g_q | Omega_d) = tr(g_q | A_d) - tr(g_q | im(C_d^T))
                       = 0 - tr(g_q | im(C_d^T))  for d >= 1

Hmm, this is getting complicated. Let me use the character theory approach.

The Z_m action on Omega_d gives a representation. Since the action is free,
each irreducible appears equally often (regular representation, up to correction).

For a free action of Z_m on a vector space V:
  dim V = m * (number of orbits)
  tr(g | V) = 0 for all g != id

So tr(g_q | Omega_d) = 0 for d >= 1 and g_q != id.

Therefore: L(g_q) = sum (-1)^d tr(g_q | Omega_d)
                   = tr(g_q | Omega_0) + sum_{d>=1} 0
                   = 1  (since Omega_0 = span{()} has g_q acting trivially)

And: L(g_q) = sum (-1)^d tr(g_q | H_d)
This means: sum (-1)^d tr(g_q | H_d) = 1 for all g_q != id.

Since H_0 = QQ (one copy, trivial action), tr(g_q | H_0) = 1.
For d >= 1: sum (-1)^d tr(g_q | H_d) = 0.

So: the alternating trace sum of the Z_m action on higher homology is ZERO.
This constrains the Z_m representation on H_m and H_{m+1}.

For d = m: H_m has dimension beta_m = m(m-3)/2.
  Z_m acts on H_m. Since the action is free on Omega_m (hence on H_m),
  m | beta_m.  beta_m = m(m-3)/2.
  For odd m: (m-3)/2 is integer iff m is odd. ✓ (m always odd for p ≡ 3 mod 4)
  m | m(m-3)/2 iff (m-3)/2 is integer iff m is odd. ✓
  So H_m = direct sum of (m-3)/2 copies of the regular representation of Z_m.

For d = m+1: H_{m+1} has dimension beta_{m+1} = m(m+1)/2.
  m | beta_{m+1} iff (m+1)/2 is integer iff m is odd. ✓
  So H_{m+1} = direct sum of (m+1)/2 copies of the regular representation.

TRACE CHECK:
  tr(g_q | H_m) = 0 (regular rep of Z_m has trace 0 for g != id)
  tr(g_q | H_{m+1}) = 0

  (-1)^m tr(g_q | H_m) + (-1)^{m+1} tr(g_q | H_{m+1}) = (-1)^m * 0 + (-1)^{m+1} * 0 = 0 ✓

The Hopf trace is automatically satisfied.

Now: what constrains beta_m = m(m-3)/2 specifically?
"""
from math import comb

print("Boundary rank analysis - what determines beta_m?")
print()

# For the k=0 eigenspace, the chain complex has Omega dimensions
# that we know. The question is: what is rank(∂_m | Omega_m)?

# Known: Omega = [1, m, m(m-1), m(m-1)(2m-3)/2, ...]
# The first few boundary ranks are determined by acyclicity (beta=0):
#   R_0 = 0, R_1 = 0, R_2 = m-1 (wait, this was 5 for m=5)

# Let me recheck. For P_7 k=0:
# Omega = [1, 3, 6, 9, 9, 6, 3]
# R = [0, 0, 3, 3, 6, 3, 3, 0]  (from omega_recursive.py)
# Verification: R_2 = 3 = m, not m-1.

# For P_11 k=0:
# Omega = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]
# R = [0, 0, 5, 15, 55, 150, 305, 390, 300, 150, 30, 0]
# R_2 = 5 = m.

# Let me find the pattern in R for both primes.
print("P_7 (m=3):")
R7 = [0, 0, 3, 3, 6, 3, 3, 0]
O7 = [1, 3, 6, 9, 9, 6, 3]
for d in range(7):
    print(f"  R_{d}={R7[d]:5d}  R_{d+1}={R7[d+1]:5d}  O_{d}={O7[d]:5d}  beta={O7[d]-R7[d]-R7[d+1]:5d}")

print()
print("P_11 (m=5):")
R11 = [0, 0, 5, 15, 55, 150, 305, 390, 300, 150, 30, 0]
O11 = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]
for d in range(11):
    print(f"  R_{d}={R11[d]:5d}  R_{d+1}={R11[d+1]:5d}  O_{d}={O11[d]:5d}  beta={O11[d]-R11[d]-R11[d+1]:5d}")

# Pattern: R_d for the ACYCLIC portion (d < m and d > m+1) is determined.
# For d <= m: R_{d+1} = Omega_d - R_d (since beta_d = 0)
# R_0 = 0
# R_1 = O_0 - 1 = 0 (since beta_0 = 1)
# R_2 = O_1 - R_1 = m - 0 = m
# R_3 = O_2 - R_2 = m(m-1) - m = m(m-2)
# R_4 = O_3 - R_3 = m(m-1)(2m-3)/2 - m(m-2)

print()
print("Boundary rank formula attempt:")
for m_val in [3, 5, 9]:
    R = [0, 0]  # R_0, R_1
    # Omega formulas: O_0=1, O_1=m, O_2=m(m-1), O_3=m(m-1)(2m-3)/2
    O_formulas = [1, m_val, m_val*(m_val-1), m_val*(m_val-1)*(2*m_val-3)//2]
    print(f"  m={m_val}:")
    print(f"    O = {O_formulas}")
    R.append(O_formulas[1])  # R_2 = m
    R.append(O_formulas[2] - R[-1])  # R_3 = O_2 - R_2
    if len(O_formulas) > 3:
        R.append(O_formulas[3] - R[-1])  # R_4 = O_3 - R_3

    print(f"    R = {R}")
    # Predicted beta_m:
    # From acyclicity + the boundary ranks
    # R_d + R_{d+1} = O_d for d != m, m+1
    # R_m + R_{m+1} + beta_m = O_m
    # R_{m+1} + R_{m+2} + beta_{m+1} = O_{m+1}
    # With beta_m = beta_{m+1} and all d > m+1 acyclic
    print(f"    R_2 = {m_val}")
    print(f"    R_3 = {m_val*(m_val-1) - m_val} = {m_val*(m_val-2)}")

# Key: if we know all Omega_d and use acyclicity everywhere except d=m,m+1,
# and use beta_m = beta_{m+1}, we can determine everything.
# The unknown is beta_m itself, which equals m(m-3)/2.
# This must come from a specific property of the chain complex.

print()
print("CONJECTURE: beta_m^{(0)} = m(m-3)/2 for all Paley P_p, m >= 3.")
print("  m=3: beta_3 = 0 (contractible)")
print("  m=5: beta_5 = 5")
print("  m=9: beta_9 = 27 (predicted)")
print("  m=11: beta_11 = 44 (predicted)")
print()
print("What is m(m-3)/2?")
for m_val in [3, 5, 7, 9, 11]:
    b = m_val*(m_val-3)//2
    print(f"  m={m_val}: {b} = C({m_val-1},2) - {comb(m_val-1,2) - b}")
    print(f"    Also: {b} = sum_{{k=1}}^{{(m-3)/2}} (m-2k) = ?")
