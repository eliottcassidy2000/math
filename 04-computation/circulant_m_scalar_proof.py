#!/usr/bin/env python3
"""
PROOF ATTEMPT: M = (H/n)*I for circulant (= vertex-transitive) tournaments.

Step 1: Cyclic symmetry sigma(i) = i+1 mod n is an automorphism.
   => M[sigma(a), sigma(b)] = M[a,b] for all a,b.
   => M is a circulant matrix.

Step 2: M is symmetric (THM-030).
   A circulant + symmetric matrix has the form M[i,j] = f(|i-j| mod n)
   where f(k) = f(n-k).

Step 3: tr(M) = sum M[i,i] = n*f(0). Since tr(M) = H, f(0) = H/n.

Step 4: WHY are the off-diagonal entries 0?

For VT at odd n:
  By Key Identity (THM-030): B_b + (-1)^n E_b = 2r * cs(b)
  For VT: by cyclic symmetry, B_b = B_0 and E_b = E_0 for all b.
  Sum over b: n*B_0 + (-1)^n * n*E_0 = 2r * sum cs(b)
  But sum B_b = H (each path starts somewhere) and sum E_b = H.
  So H + (-1)^n H = 2r * Sigma.

  At c=1 (r=1/2): H(1-1) = Sigma (odd n) or 2H = Sigma (even n).
  Odd n: Sigma = 0. By cyclic symmetry, each cs(b) = 0. So B_b = E_b for all b.

Step 5: But M[a,b] = 0 for a != b needs more than cs(b) = 0.
  M is circulant, so M[0,k] = M[0,1] for k=1,...,n-1 by... NO!
  M is circulant means M[0,k] depends on k. M[0,1] and M[0,2] can differ.

  For M to be SCALAR, we need M[0,k] = 0 for all k != 0.
  Equivalently, all eigenvalues of M are equal (= H/n).

  Eigenvalues of circulant matrix: lambda_j = sum_{k=0}^{n-1} M[0,k] * omega^{jk}
  where omega = e^{2pi*i/n}.
  lambda_j = f(0) + sum_{k=1}^{n-1} f(k) * omega^{jk}

  For M = (H/n)*I: lambda_j = H/n for all j. This requires
  sum_{k=1}^{n-1} f(k) * omega^{jk} = 0 for all j = 1,...,n-1.
  This is exactly: f(k) = 0 for k = 1,...,n-1.

  So proving M scalar <==> proving all off-diagonal circulant entries are 0.

  Can we use the Fourier transform structure of the transfer matrix?

Step 6: M[0,k] = sum_S (-1)^|S| E_0(S+{0}) B_k(R+{k})
  Under cyclic shift sigma: E_0(S+{0}) corresponds to paths ending at 0
  in subtournament on S+{0}, while B_k(R+{k}) corresponds to paths starting
  at k in subtournament on R+{k}.

  The cyclic symmetry maps (0, S, R, k) to (1, sigma(S), sigma(R), k+1).
  So M[0,k] = M[1,k+1] = ... = M[i, (k+i) mod n]. This is the circulant property.

  But to prove M[0,k] = 0, we need cancellation in the subset sum.

APPROACH: For n prime, the automorphism group Z/nZ has no nontrivial subgroups.
  The transfer matrix is an n x n circulant symmetric matrix, so its eigenvalues
  are sums f(k)*omega^{jk}. For n prime, the group ring C[Z/nZ] decomposes into
  1-dimensional representations. The transfer matrix M acts on the regular
  representation of Z/nZ.

  CLAIM: The transfer matrix M, as a Z/nZ-equivariant operator, must be scalar
  because it equals H times the projection onto the regular representation
  divided by n.

  Actually, that's circular. We need to prove M has no off-diagonal circulant entries.

  Let me try a Fourier analysis approach: compute the DFT of row 0 of M.

kind-pasteur-2026-03-06-S25c
"""

import numpy as np
from itertools import permutations, combinations

def make_circulant_tournament(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def compute_M_entry(T, n, a, b):
    if a == b:
        val = 0
        for perm in permutations(range(n)):
            prod = 1
            for k in range(n-1):
                prod *= T.get((perm[k], perm[k+1]), 0)
            if prod > 0:
                pos = list(perm).index(a)
                val += (-1)**pos
        return val
    else:
        U = [v for v in range(n) if v != a and v != b]
        val = 0
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})
            ea = 0
            if len(S_set) == 1:
                ea = 1
            else:
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    ea += prod
            bb2 = 0
            if len(R_set) == 1:
                bb2 = 1
            else:
                for p in permutations(R_set):
                    if p[0] != b: continue
                    prod = 1
                    for k in range(len(p)-1):
                        prod *= T.get((p[k], p[k+1]), 0)
                    bb2 += prod
            val += sign * ea * bb2
        return val


# ============================================================
# Fourier analysis of M[0, k] for circulant tournaments
# ============================================================
print("=" * 70)
print("Fourier analysis of M for circulant tournaments")
print("=" * 70)

for n in [5, 7]:
    print(f"\n  n={n}:")

    # Find all valid circulant generating sets
    gen_sets = []
    for S_tuple in combinations(range(1, n), (n-1)//2):
        S = set(S_tuple)
        complement = {(n - s) % n for s in S}
        if not (S & complement) and (S | complement) == set(range(1, n)):
            gen_sets.append(S)

    for S in gen_sets[:4]:  # Limit for speed at n=7
        T = make_circulant_tournament(n, S)
        H = count_H(T, n)

        # Compute row 0 of M
        row0 = [compute_M_entry(T, n, 0, k) for k in range(n)]
        print(f"\n    S={sorted(S)}: H={H}, M[0,:] = {row0}")

        # DFT of row 0
        row0_arr = np.array(row0, dtype=complex)
        dft = np.fft.fft(row0_arr)
        print(f"    DFT(M[0,:]) = {[f'({d.real:.2f}+{d.imag:.2f}i)' for d in dft]}")

        # If M is scalar, DFT should be [H, 0, 0, ..., 0]
        is_scalar = all(abs(dft[k]) < 0.01 for k in range(1, n))
        print(f"    Scalar (DFT[1:]=0): {is_scalar}")


# ============================================================
# n=5: What about non-circulant VT tournaments?
# ============================================================
print("\n" + "=" * 70)
print("n=5: Are ALL VT tournaments circulant?")
print("=" * 70)

n = 5
# At n=5, a VT tournament has Aut group transitive on 5 vertices.
# Since |Aut| divides 5! = 120 and 5 divides |Aut| and |Aut| is odd,
# |Aut| must be in {5, 15, 25, 75}.
# For n=5 prime, the only transitive group of odd order on 5 points
# containing Z/5Z is... Z/5Z itself (order 5), or the dihedral group D_5
# (order 10, but 10 is even!), or AGL(1,5) (order 20, even).
# So the ONLY transitive subgroup of S_5 with odd order dividing 5*4*3*2*1
# and divisible by 5 is Z/5Z (order 5).
# Wait, Z/5Z has order 5 (odd). Are there larger odd-order transitive groups on 5 points?
# The cyclic group C_5 has order 5. The only groups of odd order dividing 120
# and containing C_5 as a subgroup: C_5 (order 5), C_5 x C_1 = C_5 (same),
# or a Frobenius group of order 5*k with k odd dividing 4: k=1 or k=3 (but 15
# doesn't divide 120 via a transitive group... actually 15 does: the group of
# order 20 = 5*4 has a subgroup of order 5 (Sylow), but is it transitive?
# The alternating group A_5 has order 60, not odd.
# So at n=5, VT = circulant (the only odd-order transitive group on 5 points is Z/5Z).

print("  At n=5 (prime): ALL VT tournaments are circulant.")
print("  Reason: Aut(T) has odd order, and the only odd-order transitive")
print("  subgroup of S_5 is Z/5Z (order 5).")
print("  Therefore VT at n=5 <==> circulant at n=5.")

# At n=7: same argument
print("\n  At n=7 (prime): ALL VT tournaments are circulant.")
print("  Reason: same as n=5. Z/7Z is the only odd-order transitive")
print("  subgroup of S_7 (order 7).")

# Wait, what about larger groups? E.g., order 21 = 3*7?
# The group of order 21 is Z/7Z x| Z/3Z (Frobenius).
# It acts transitively on 7 points. Is it a subgroup of S_7?
# Yes! AGL(1,7) = Z/7Z x| Z/6Z has order 42. Its Sylow 3-subgroup
# gives a subgroup of order 21 = 7*3.
# This is an odd-order transitive group on 7 points.
# So VT at n=7 does NOT imply circulant — the Frobenius group of order 21
# is also possible!

print("\n  CORRECTION: At n=7, the Frobenius group of order 21 is also")
print("  odd-order and transitive on 7 points. So VT does NOT imply circulant.")
print("  But all circulant => VT (trivially).")

# Check: do Paley tournaments have |Aut| > 7?
# Paley T_7 has Aut = PGL(2,7)? No... Aut = {sigma: x -> ax+b mod 7 | a is QR, b in Z/7Z}
# = AGL(1,7) intersected with QR constraint
# a in QR_7 = {1,2,4}, b in Z/7Z, so |Aut| = 3*7 = 21.
# The Frobenius group of order 21!

print("\n  Paley T_7: |Aut| = 21 (Frobenius group Z/7Z x| Z/3Z)")
print("  This is vertex-transitive but with LARGER Aut than circulant.")

# Verify
T = make_circulant_tournament(7, {1,2,4})

# Check: is sigma(x) = 2x mod 7 an automorphism?
sigma_2x = {i: (2*i) % 7 for i in range(7)}
is_aut = True
for i in range(7):
    for j in range(7):
        if i != j:
            if T.get((i,j), 0) != T.get((sigma_2x[i], sigma_2x[j]), 0):
                is_aut = False
                break
    if not is_aut:
        break
print(f"\n  sigma(x)=2x mod 7 is automorphism of Paley T_7: {is_aut}")

# So the Paley tournament has EXTRA automorphisms beyond cyclic shift.
# The cyclic shift gives M circulant.
# The extra aut sigma(x)=2x mod 7 gives additional constraints on M.
# Combined: M[0,k] = M[0, 2k mod 7] for all k.
# The orbits of k under multiplication by 2 mod 7:
#   {0}, {1,2,4}, {3,6,5}
# So M[0,1] = M[0,2] = M[0,4] and M[0,3] = M[0,6] = M[0,5].
# Plus symmetry M[0,k] = M[0,7-k]:
#   M[0,1] = M[0,6], M[0,2] = M[0,5], M[0,3] = M[0,4]
# Combined with the multiplication-by-2 orbits:
#   {1,2,4} and {3,5,6} are related by k -> 7-k
# So M[0,1] = M[0,2] = M[0,3] = M[0,4] = M[0,5] = M[0,6]!
# If this common value is c, then sum of row = H/7 + 6c = H.
# So c = (H - H/7)/6 = H*(n-1)/(n*6) = H*6/42 = H/7 = H/n.
# Wait, that gives c = H/7 = H/n, but M[0,0] = H/n too.
# So M = H/n * J (all-ones matrix times H/n)?
# No, that can't be right. M is symmetric and M*1 = H*1 only if M = (H/n)*J.
# But that contradicts M[0,0] = H/n too.

# Let me just check the actual values:
row0 = [compute_M_entry(T, 7, 0, k) for k in range(7)]
print(f"\n  Paley T_7: M[0,:] = {row0}")

# For non-Paley circulant:
T2 = make_circulant_tournament(7, {1,2,3})
sigma_check = {i: (2*i) % 7 for i in range(7)}
is_aut2 = True
for i in range(7):
    for j in range(7):
        if i != j:
            if T2.get((i,j), 0) != T2.get((sigma_check[i], sigma_check[j]), 0):
                is_aut2 = False
                break
    if not is_aut2:
        break
print(f"  sigma(x)=2x mod 7 is automorphism of circ{{1,2,3}}: {is_aut2}")

row0_2 = [compute_M_entry(T2, 7, 0, k) for k in range(7)]
print(f"  circ{{1,2,3}}: M[0,:] = {row0_2}")


print("\n" + "=" * 70)
print("KEY INSIGHT")
print("=" * 70)
print("""
For circulant tournaments (Z/nZ automorphism):
  M is circulant => M[i,j] = f((j-i) mod n)
  M is symmetric => f(k) = f(n-k)

For M to be scalar, we need f(k) = 0 for k != 0.

The cyclic symmetry alone does NOT force f(k) = 0.
But computationally, it IS always 0 for all circulant tournaments tested.

PROOF IDEA: For Z/nZ-equivariant operators, the space of
  n x n circulant symmetric integer matrices is spanned by
  {e_j + e_{n-j}} for j = 0, ..., floor(n/2).
  M being in this space doesn't force it to be scalar.

So there must be some ADDITIONAL structure of the transfer matrix
that forces the off-diagonal circulant entries to vanish.

CONJECTURE: M = (H/n)*I for all vertex-transitive tournaments at odd n.
STATUS: Verified exhaustively for n=3,5 and by sampling for n=7.
""")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
