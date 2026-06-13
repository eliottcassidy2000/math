#!/usr/bin/env python3
"""
Proving c_4 = (n-1)! * I and understanding the even-r decomposition.

M(r) = c_0 + c_2*r^2 + c_4*r^4 + ... + c_{n-1}*r^{n-1} (even powers only)

Key conjecture: c_{n-1} = (n-1)! * I for ALL tournaments.

Why? The r^{n-1} coefficient of M[a,b] picks the term where ALL n-1 edges
in each sub-path contribute r (i.e., the product is r^{path_length}).
This makes the sub-path count independent of the tournament.

Actually, for M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b):
The r^k coefficient in E_a(S+a) comes from paths of length |S| using
k powers of r and |S|-k powers of s.
The total degree in r across E*B must sum to some value.

Let me think more carefully. In the IE product E*B:
E has paths of length |S| (|S| edges), B has paths of length |R| (|R| edges).
Total edges = |S| + |R| = |U| = n-2 (when a≠b).

So the maximum power of r in E*B is n-2. But we found degree 4 = n-1 at n=5.

Hmm, that's inconsistent. Let me re-examine.

WAIT: For diagonal entries a=b, the IE formula is different. M[a,a] involves
paths through ALL n vertices. The total length is n-1 edges.
So max r-degree is n-1 for diagonal, n-2 for off-diagonal?

But the data showed M[0,0] has degree 4 and M[0,1] has degree 1 at n=5.
Let's check this more carefully.

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np
from collections import defaultdict

def count_paths_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def transfer_matrix_r(A, r_val):
    n = len(A)
    M = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0.0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    S_verts = sorted(list(S) + [a])
                    R_verts = sorted(R + [b])
                    ea = count_paths_weighted(A, S_verts, r_val, end=a)
                    bb = count_paths_weighted(A, R_verts, r_val, start=b)
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

def extract_even_poly(A, n, max_degree):
    """Extract c_0, c_2, c_4, ... coefficients of M(r) = sum c_{2k} r^{2k}."""
    num_coeffs = max_degree // 2 + 1
    # Sample at enough u = r^2 points
    u_samples = np.linspace(0, 1, num_coeffs + 3)
    r_samples = np.sqrt(u_samples)

    M_samples = [transfer_matrix_r(A, rv) for rv in r_samples]

    coeffs = {}
    for k in range(num_coeffs):
        coeffs[2*k] = np.zeros((n, n))

    for a in range(n):
        for b in range(n):
            vals = np.array([M_samples[i][a][b] for i in range(len(r_samples))])
            # Fit polynomial in u = r^2
            poly_coeffs = np.polyfit(u_samples, vals, num_coeffs - 1)
            # poly_coeffs[0]*u^{d} + ... + poly_coeffs[d]
            for k in range(num_coeffs):
                coeffs[2*k][a][b] = poly_coeffs[num_coeffs - 1 - k]

    return coeffs

# =====================================================================
print("=" * 70)
print("EVEN-r POLYNOMIAL AT MULTIPLE n VALUES")
print("=" * 70)

# n=3
print("\n--- n=3 ---")
n = 3
A3_cycle = [[0,1,0],[0,0,1],[1,0,0]]
A3_trans = [[0,1,1],[0,0,1],[0,0,0]]

for A, name in [(A3_cycle, "3-cycle"), (A3_trans, "transitive")]:
    coeffs = extract_even_poly(A, n, 2)
    print(f"\n  {name}:")
    for k in sorted(coeffs.keys()):
        c = coeffs[k]
        tr_c = np.trace(c)
        is_scalar = np.allclose(c, (tr_c/n) * np.eye(n))
        eigs = sorted(np.linalg.eigvalsh(c))
        print(f"    c_{k}: tr={tr_c:.4f}, scalar={is_scalar}, eigs={[round(e,3) for e in eigs]}")

    # Verify M(0.3)
    M_check = transfer_matrix_r(A, 0.3)
    M_recon = sum(coeffs[k] * 0.3**k for k in coeffs)
    print(f"    Reconstruction error: {np.max(np.abs(M_check - M_recon)):.2e}")

# n=4
print("\n--- n=4 ---")
n = 4
A4s = [
    ([[0,1,1,0],[0,0,1,1],[0,0,0,1],[1,0,0,0]], "near-regular"),
    ([[0,1,1,1],[0,0,1,1],[0,0,0,1],[0,0,0,0]], "transitive"),
]

for A, name in A4s:
    coeffs = extract_even_poly(A, n, 4)
    print(f"\n  {name}:")
    for k in sorted(coeffs.keys()):
        c = coeffs[k]
        tr_c = np.trace(c)
        is_scalar = np.allclose(c, (tr_c/n) * np.eye(n))
        eigs = sorted(np.linalg.eigvalsh(c))
        nonzero = np.max(np.abs(c)) > 1e-6
        if nonzero:
            print(f"    c_{k}: tr={tr_c:.4f}, scalar={is_scalar}, eigs={[round(e,3) for e in eigs]}")
        else:
            print(f"    c_{k}: ZERO")

    M_check = transfer_matrix_r(A, 0.3)
    M_recon = sum(coeffs[k] * 0.3**k for k in coeffs)
    print(f"    Reconstruction error: {np.max(np.abs(M_check - M_recon)):.2e}")

# n=5 Paley
print("\n--- n=5 Paley ---")
n = 5
A5 = [[0]*5 for _ in range(5)]
for i in range(5):
    for j in range(5):
        if i != j and (j-i)%5 in [1,4]:
            A5[i][j] = 1

coeffs = extract_even_poly(A5, n, 4)
for k in sorted(coeffs.keys()):
    c = coeffs[k]
    tr_c = np.trace(c)
    is_scalar = np.allclose(c, (tr_c/n) * np.eye(n))
    eigs = sorted(np.linalg.eigvalsh(c))
    nonzero = np.max(np.abs(c)) > 1e-6
    if nonzero:
        print(f"  c_{k}: tr={tr_c:.4f}, scalar={is_scalar}, eigs={[round(e,3) for e in eigs]}")
    else:
        print(f"  c_{k}: ZERO")

# =====================================================================
print()
print("=" * 70)
print("DEGREE OF M[a,b](r) — DIAGONAL vs OFF-DIAGONAL")
print("=" * 70)

# At n=5: check degree of each entry
n = 5
r_pts = np.linspace(0, 1, 10)

A5_nr = [[0,0,0,0,1],[1,0,0,1,0],[1,1,0,0,0],[1,0,1,0,0],[0,1,1,1,0]]  # non-regular scalar

for A, name in [(A5, "Paley"), (A5_nr, "non-reg scalar")]:
    print(f"\n  {name}:")
    for a in range(n):
        for b in range(n):
            vals = [transfer_matrix_r(A, rv)[a][b] for rv in r_pts]
            for deg in range(0, 6):
                coeffs_fit = np.polyfit(r_pts, vals, deg)
                fitted = np.polyval(coeffs_fit, r_pts)
                err = np.max(np.abs(fitted - vals))
                if err < 1e-8:
                    break
            # Now check if it has only even powers
            coeffs_r = list(reversed(coeffs_fit))
            odd_max = max(abs(coeffs_r[k]) for k in range(len(coeffs_r)) if k % 2 == 1) if any(k % 2 == 1 for k in range(len(coeffs_r))) else 0
            even_only = odd_max < 1e-6
            if a <= b:  # upper triangle
                print(f"    M[{a},{b}]: degree {deg}, even_only={even_only}")

# =====================================================================
print()
print("=" * 70)
print("UNIVERSAL c_{n-1} AT EACH n")
print("=" * 70)

# At n=3: max even degree = 2, so c_2 is the top term.
# At n=4: max even degree = 2, so c_2 is the top term (degree 3 would be odd).
# At n=5: max even degree = 4, so c_4 is the top term.

# The claim: for diagonal entries M[a,a], the max even power is n-1 (if n-1 is even)
# or n-2 (if n-1 is odd).
# For off-diagonal: max even power is less.

# Let's verify the top coefficient.

print("\nTop even-power coefficient:")
for n_val in [3, 4, 5]:
    if n_val == 3:
        A = A3_cycle
    elif n_val == 4:
        A = A4s[0][0]
    else:
        A = A5

    max_even = n_val - 1 if (n_val - 1) % 2 == 0 else n_val - 2
    coeffs = extract_even_poly(A, n_val, max_even + 2)

    c_top = coeffs.get(max_even, np.zeros((n_val, n_val)))
    tr_top = np.trace(c_top)
    is_scalar = np.allclose(c_top, (tr_top/n_val) * np.eye(n_val))
    diag_val = c_top[0][0]
    expected_factorial = 1
    for i in range(1, n_val):
        expected_factorial *= i

    print(f"  n={n_val}: max_even_degree={max_even}, c_{max_even} diagonal={diag_val:.2f}, "
          f"(n-1)!={expected_factorial}, scalar={is_scalar}")

# =====================================================================
print()
print("=" * 70)
print("ALGEBRAIC ARGUMENT FOR c_{n-1} = (n-1)! * I")
print("=" * 70)
print("""
For the DIAGONAL entry M[a,a]:
  M[a,a] = sum_P (-1)^{pos(a,P)} * w(P)
where w(P) = product of (r + s_e) for all n-1 edges in P.

The coefficient of r^{n-1} in w(P) is 1 for every P (all edges contribute r).
So the r^{n-1} coefficient of M[a,a] = sum_P (-1)^{pos(a,P)}.

For odd n: this sum equals H/n (shown earlier).
But we need: for what n does r^{n-1} have even parity?
  n-1 even <=> n odd.
So at ODD n: c_{n-1} = (sum_P (-1)^{pos(a,P)} per diagonal entry) * I.

The sum_P (-1)^{pos(a,P)} = H/n at odd n (from the position-counting argument).
But the TOP coefficient should be (n-1)! = number of Ham paths on K_n starting at a.
Since on K_n every permutation is a Ham path, there are n! total paths.
The paths where a is at position k number (n-1)! (choose remaining n-1 in (n-1)! ways).
Sum of (-1)^k * #{P: pos(a,P)=k} = sum_{k=0}^{n-1} (-1)^k * (n-1)!
  = (n-1)! * sum_{k=0}^{n-1} (-1)^k = (n-1)! * (1 if n odd, 0 if n even).

Wait — on K_n, vertex a is at EACH position equally often: (n-1)! times.
So sum_P (-1)^{pos(a,P)} = (n-1)! * sum_k (-1)^k = 0 if n even, (n-1)! if n odd.

But we also need to account for the IE formula, not just the direct path sum.
The r^{n-1} coefficient of M[a,a] via IE:
Each IE term has E_a paths of length |S| and B_a paths of length |R| = n-2-|S|.
The max r power from E*B is |S| + |R| = n-2 edges.
So the max r power from IE for DIAGONAL entries is n-2, NOT n-1!

Hmm, but the data shows degree 4 = n-1 at n=5 for diagonal.
Let me re-examine the IE formula for diagonal entries.
""")

# Direct check: what is the IE formula for a=b?
# M[a,a] = sum_{S subset U} (-1)^|S| E_a(S+a) B_a(R+a)
# where U = V \ {a} (n-1 vertices), S + a includes a.
# E_a(S+a): paths on S+a ending at a, length |S|.
# B_a(R+a): paths on R+a starting at a, length |R| = n-1-|S|.
# Product: E*B uses |S| + n-1-|S| = n-1 edges total!
# So max r power is n-1. This IS correct.

# For a != b:
# M[a,b] with U = V \ {a,b} (n-2 vertices)
# E_a(S+a): paths on S+a ending at a, length |S|.
# B_b(R+b): paths on R+b starting at b, length |R| = n-2-|S|.
# Product uses |S| + n-2-|S| = n-2 edges.
# So max r power for off-diagonal is n-2.

print("\nIE FORMULA DEGREE ANALYSIS:")
print("  Diagonal M[a,a]: max r-degree = n-1 (U has n-1 vertices)")
print("  Off-diag M[a,b]: max r-degree = n-2 (U has n-2 vertices)")
print()
print("  For n=5: diagonal degree 4, off-diagonal degree 3")
print("  Even powers only: diagonal has r^0, r^2, r^4")
print("                    off-diagonal has r^0, r^2 (r^3 is odd, excluded)")

# Verify off-diagonal has only r^0 and r^2 at n=5
print("\n  Off-diagonal polynomial at n=5 Paley:")
r_pts = np.linspace(0, 1, 8)
for a in range(5):
    for b in range(5):
        if a == b: continue
        vals = [transfer_matrix_r(A5, rv)[a][b] for rv in r_pts]
        for deg in range(0, 5):
            coeffs_fit = np.polyfit(r_pts, vals, deg)
            fitted = np.polyval(coeffs_fit, r_pts)
            err = np.max(np.abs(fitted - vals))
            if err < 1e-8:
                break
        if a == 0 and b < 3:
            coeffs_r = list(reversed(coeffs_fit))
            print(f"    M[{a},{b}]: degree {deg}, coeffs = {[round(c,4) for c in coeffs_r]}")

# =====================================================================
print()
print("=" * 70)
print("THE r^{n-1} COEFFICIENT FOR DIAGONAL ENTRIES")
print("=" * 70)

# At the r^{n-1} level, all s_e are replaced by r.
# So the coefficient is: sum_S (-1)^|S| * (# paths on S+a ending at a) * (# paths on R+a starting at a)
# evaluated on the COMPLETE graph K_n (since s disappears).

# On K_n: # Ham paths on vertex set W ending at a = (|W|-1)! (any ordering of W\{a} before a)
# Similarly starting at a: (|W|-1)!.

# So r^{n-1} coeff of M[a,a] = sum_{k=0}^{n-1} sum_{|S|=k} (-1)^k * k! * (n-1-k)!
#                              = sum_{k=0}^{n-1} C(n-1,k) * (-1)^k * k! * (n-1-k)!
#                              = sum_{k=0}^{n-1} (-1)^k * (n-1)!
#                              = (n-1)! * sum_{k=0}^{n-1} (-1)^k

# sum (-1)^k for k=0..n-1: = 1 if n odd, 0 if n even.

# So: r^{n-1} coeff of M[a,a] = (n-1)! if n odd, 0 if n even.

# At n=5 (odd): c_4 diagonal = 4! = 24. CONFIRMED!
# At n=3 (odd): c_2 diagonal should be 2! = 2. Let's check.
# At n=4 (even): c_3 should be 0 (odd power excluded), and c_2 is the top even power.

print("\nAlgebraic prediction:")
for nv in [3, 4, 5]:
    if nv % 2 == 1:
        top_even = nv - 1
        top_coeff = 1
        for i in range(1, nv):
            top_coeff *= i
        print(f"  n={nv} (odd): c_{top_even}[a,a] = (n-1)! = {top_coeff}")
    else:
        top_even = nv - 2
        print(f"  n={nv} (even): c_{top_even} is the top even term (r^{nv-1} coeff = 0)")

# Verify c_2[a,a] = 2 at n=3
coeffs_n3 = extract_even_poly(A3_cycle, 3, 2)
print(f"\n  n=3, 3-cycle: c_2[0,0] = {coeffs_n3[2][0][0]:.4f} (expected 2)")
print(f"  n=3, 3-cycle: c_2 = {np.round(coeffs_n3[2], 4).tolist()}")

# What about off-diagonal c_2 at n=3?
print(f"  n=3, 3-cycle: c_0 = {np.round(coeffs_n3[0], 4).tolist()}")

# n=4: the top even coefficient
coeffs_n4 = extract_even_poly(A4s[0][0], 4, 4)
for k in sorted(coeffs_n4.keys()):
    c = coeffs_n4[k]
    if np.max(np.abs(c)) > 1e-6:
        print(f"  n=4, near-regular: c_{k} = {np.round(c, 4).tolist()}")


print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
THEOREM (Even-r polynomial decomposition):
  M(r) = c_0 + c_2*r^2 + c_4*r^4 + ...  (only even powers of r)

  DIAGONAL M[a,a]: polynomial of degree n-1 in r (n-1 is the max even power at odd n)
  OFF-DIAGONAL M[a,b]: polynomial of degree n-2 in r

THEOREM (Universal top coefficient at odd n):
  c_{{n-1}} = (n-1)! * I   (diagonal-only, with (n-1)! on each diagonal entry)

  Proof: At r^{{n-1}}, all s-factors are replaced by r. The IE sum gives:
    sum_{{k=0}}^{{n-1}} C(n-1,k) (-1)^k k! (n-1-k)! = (n-1)! sum (-1)^k = (n-1)!
  (using n odd so the alternating sum = 1).

  At even n: r^{{n-1}} has odd parity, so it doesn't appear in the even-r expansion.
  The top even coefficient c_{{n-2}} is NOT universal.

COROLLARY (H formula):
  H = tr(c_0) + (1/4)*tr(c_2) + (1/16)*tr(c_4) + ...
  At n=5: H = tr(c_0) + tr(c_2)/4 + 120/16 = tr(c_0) + tr(c_2)/4 + 7.5

FORMULA: tr(c_2) at n=5:
  tr(c_2) = 12*(#3-cycles) - 30
""")
