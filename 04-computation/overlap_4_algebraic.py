#!/usr/bin/env python3
"""
ALGEBRAIC DERIVATION OF THE (4,) CONTRIBUTION TO tr(c_{n-5})

THEOREM: At n=7, the (4,) contribution to tr(c_2) = sum_P e_4(s_P) is:
  126 - 36*t_3 + 12*t_5

PROOF:
=====

1. A (4,) position pattern is {i, i+1, i+2, i+3} (4 consecutive).
   At n=7: 3 such patterns ({0,1,2,3}, {1,2,3,4}, {2,3,4,5}).

2. For each pattern, the product s_i*s_{i+1}*s_{i+2}*s_{i+3} involves
   5 consecutive vertices: p_i, p_{i+1}, p_{i+2}, p_{i+3}, p_{i+4}.

3. sum_P s_i*s_{i+1}*s_{i+2}*s_{i+3}
   = (n-5)! * sum_{5-vertex subsets S} G(S)

   where G(S) = sum_{orderings of S} prod_{k=0}^3 (A[p_k,p_{k+1}]-1/2)
   = tr(c_0) of the 5-vertex subtournament T[S].

4. At n=5 (the subtournament): tr(c_0) = sum_P prod s_i.
   By the coefficient relation: c_0 + c_2/4 + c_4/16 = H(S).
   c_4 = (5-1)! = 24 (universal).
   c_2 = 2*(5-2)!*(t_3(S) - C(5,3)/4) = 12*t_3(S) - 30.

   So: c_0(S) = H(S) - (12*t_3(S)-30)/4 - 24/16
              = H(S) - 3*t_3(S) + 30/4 - 24/16
              = H(S) - 3*t_3(S) + 7.5 - 1.5
              = H(S) - 3*t_3(S) + 6.

   Wait, that gives c_0(S) = H(S) - 3*t_3(S) + 6, but earlier we had
   G(S) = H(S) - 3*t_3(S). Let me recheck.

   Actually: c_0 + c_2/4 + c_4/16 = H.
   c_0 = H - c_2/4 - c_4/16.
   c_2/4 = (12*t_3(S) - 30)/4 = 3*t_3(S) - 30/4.
   c_4/16 = 24/16 = 3/2.
   c_0 = H - 3*t_3(S) + 30/4 - 3/2 = H - 3*t_3(S) + 15/2 - 3/2 = H - 3*t_3(S) + 6.

   Hmm, but the numerical verification showed tr(c_0)(S) = H(S) - 3*t_3(S).
   Let me recheck the c_2 formula at n=5.

   THM-054: tr(c_{n-3}) = 2*(n-2)!*(t_3 - C(n,3)/4).
   At n=5: tr(c_2) = 2*3!*(t_3 - 10/4) = 12*(t_3 - 5/2) = 12*t_3 - 30.
   tr(c_4) = (n-1)! = 24.

   H = c_0 + c_2/4 + c_4/16.

   Wait, the powers of 4: H = sum_{k} c_{2k} * (1/2)^{2k} where 2k ranges...
   No. At n=5: H = c_0 + c_2*(1/4) + c_4*(1/16)? Or H = c_0 + c_2*(1/2)^2 + c_4*(1/2)^4?

   Actually: M(r) = c_0 + c_2*r^2 + c_4*r^4.
   H = M(1/2) = c_0 + c_2/4 + c_4/16.

   c_4 = tr(c_4) = 24.
   c_2 = tr(c_2) = 12*t_3 - 30.
   c_0 = H - c_2/4 - c_4/16 = H - (12*t_3-30)/4 - 24/16
       = H - 3*t_3 + 15/2 - 3/2 = H - 3*t_3 + 6.

   But at n=5 with OCF: H = 1 + 2*(t_3+t_5).
   c_0 = 1 + 2*t_3 + 2*t_5 - 3*t_3 + 6 = 7 - t_3 + 2*t_5.

   For the "transitive" tournament at n=5 (t_3=0, t_5=0, H=1):
   c_0 = 1 - 0 + 6 = 7.
   But tr(c_0) should equal sum_P prod s_i.

   Let me compute directly for the transitive tournament on 5 vertices:
   At n=5 transitive: all A[i][j]=1 for i<j.
   H = 1 (only the natural ordering works).

   sum_P prod_{k=0}^3 s_k = sum_P prod_{k=0}^3 (A[p_k,p_{k+1}]-1/2)

   For the natural ordering (0,1,2,3,4): all s_k = 1/2, prod = 1/16.
   For the reverse (4,3,2,1,0): all s_k = -1/2, prod = 1/16 (even number of factors).

   Hmm, this doesn't look like it sums to 7. Let me actually compute.

   Actually I realize the formula G(S) in the (4,) derivation should be
   tr(c_0)(S) = sum_P prod s_i. And the numerical verification showed
   sum_S tr(c_0)(S) = C(n,5) - C(n-3,2)*t_3 + 2*t_5 at n=7,
   matching 21 - 6*t_3 + 2*t_5.

   If tr(c_0)(S) = H(S) - 3*t_3(S) + 6 = 1 - 0 + 6 = 7 for transitive at n=5,
   then sum_S over C(7,5)=21 subsets would be much larger than 21.

   Something's wrong. Let me recheck the formula. The issue might be the
   convention for the r-polynomial.

   Hmm: M(r) at n=5 is a matrix whose trace gives a polynomial in r^2.
   Let me verify: at r=0, tr(M(0)) = c_0 = sum_P prod s_i.
   At r=1/2, tr(M(1/2)) = H.

   For the 5-vertex transitive tournament:
   sum_P prod_{k=0}^3 s_k = (1/2)^4 * sum_P (-1)^{f_P} = (1/16)*sum_P(-1)^f

   sum_P (-1)^f: for each permutation of 5 vertices, count forward arcs.
   Transitive: A[i][j]=1 iff i<j.

   For perm (0,1,2,3,4): f=4, (-1)^4=1.
   For perm (4,3,2,1,0): f=0, (-1)^0=1.
   In general: (-1)^f always equals... let me think.

   f(P) = #{i : p_i < p_{i+1}} (since A[a][b]=1 iff a<b for transitive).
   This is the number of ascents in the permutation.

   sum_P (-1)^{ascents(P)}: this is related to the Euler numbers.
   At n=5 (permutations of 5 elements with 4 positions):
   #{ascents} = number of i where p_i < p_{i+1}.

   The signed sum over all permutations by number of ascents:
   sum_P (-1)^{asc(P)} = alternating Euler polynomial evaluation.

   E_n(x) = sum_{P ∈ S_n} x^{asc(P)}.
   At x=-1: E_n(-1) = sum_P (-1)^{asc(P)}.

   For n>=2: E_n(-1) = 0 if n is even, (-1)^{(n-1)/2} * E_{n-1} if n is odd?
   Actually, E_n(-1) is related to tangent/secant numbers.

   Let me just compute numerically for the transitive tournament.

   Sorry, this script has gotten too theoretical. Let me just VERIFY the claim.
"""

from itertools import permutations, combinations

# Transitive tournament at n=5
n = 5
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        A[i][j] = 1

# tr(c_0) = sum_P prod s_i
total = 0
for p in permutations(range(n)):
    prod = 1
    for k in range(n-1):
        prod *= (A[p[k]][p[k+1]] - 0.5)
    total += prod

print(f"Transitive at n=5: tr(c_0) = {total}")
# (-1)^f signed sum
sign_sum = sum((-1)**sum(A[p[k]][p[k+1]] for k in range(n-1)) for p in permutations(range(n)))
print(f"  sign_sum = {sign_sum}, (1/16)*sign_sum = {sign_sum/16}")

# H and cycle counts
H = sum(1 for p in permutations(range(n)) if all(A[p[k]][p[k+1]] for k in range(n-1)))
t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
         if (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i]) > 0)
print(f"  H={H}, t3={t3}")
print(f"  H - 3*t3 = {H - 3*t3}")

# Also verify c_2 and c_4
c2 = 12*t3 - 30
c4 = 24
c0_from_H = H - c2/4 - c4/16
print(f"  c_2={c2}, c_4={c4}, c_0 (from H)={c0_from_H}")
print(f"  Direct c_0={total}")

# The issue: c_0 from the formula = H - c_2/4 - c_4/16 = 1 + 30/4 - 24/16 = 1+7.5-1.5 = 7.
# But direct computation gives total = ?

# For a RANDOM tournament at n=5
import random
random.seed(42)
B = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            B[i][j] = 1
        else:
            B[j][i] = 1

total_B = sum(
    eval("*".join(f"({B[p[k]][p[k+1]]}-0.5)" for k in range(n-1)))
    for p in permutations(range(n))
)
H_B = sum(1 for p in permutations(range(n)) if all(B[p[k]][p[k+1]] for k in range(n-1)))
t3_B = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
           if (B[i][j]*B[j][k]*B[k][i] + B[i][k]*B[k][j]*B[j][i]) > 0)
c0_B = H_B - (12*t3_B-30)/4 - 24/16
print(f"\nRandom at n=5: tr(c_0)={total_B}, formula c_0={c0_B}, H={H_B}, t3={t3_B}")

# Hmm, something is definitely wrong. Let me double-check the convention.
# tr(M(r)) = sum_P prod(r + s_i)
# = sum_P sum_k e_k(s_P) * r^{n-1-k}
# So the coefficient of r^{n-1-2k} is sum_P e_{2k}(s_P).
# c_{n-1-2k} = sum_P e_{2k}(s_P).
# At n=5: c_4 = sum_P e_0 = sum_P 1 = n! = 120. NOT 24!

print(f"\nOh! c_4 = n! = {120} at n=5, not (n-1)! = 24.")
print("The convention is c_{n-1-2k} = sum_P e_{2k}(s_P).")
print("At n=5: c_4 corresponds to 2k=0, so c_4 = e_0 sum = n!.")
print("c_2 corresponds to 2k=2, so c_2 = sum_P e_2(s_P).")
print("c_0 corresponds to 2k=4, so c_0 = sum_P e_4(s_P) = sum_P prod s_i.")

# So the CORRECT relation is:
# tr(M(r)) = c_0 + c_2*r^2 + c_4*r^4
# H = tr(M(1/2)) = c_0 + c_2/4 + c_4/16
# c_4 = 120 (at n=5)
# c_2 = sum_P e_2(s_P) = ?

# THM-054: tr(c_{n-3}) = 2*(n-2)!*(t_3 - C(n,3)/4)
# At n=5: n-3=2, so tr(c_2) = 2*3!*(t_3 - 10/4) = 12*(t_3 - 5/2) = 12*t_3 - 30.

# Wait, but c_4 at n=5 = tr(c_4) = sum_P e_0 = n! = 120.
# And THM-054 says tr(c_{n-3}) at n=5 = tr(c_2) = 12*t_3 - 30.
# But that seems too small. For transitive: 12*0 - 30 = -30.
# c_0 = H - c_2/4 - c_4/16 = 1 + 30/4 - 120/16 = 1 + 7.5 - 7.5 = 1. That makes sense!

# So c_4 = n! = 120, NOT (n-1)! = 24.
print(f"\nCORRECTED: c_4 = n! = 120 at n=5.")
c4_correct = 120
c2_val = 12*t3 - 30
c0_correct = H - c2_val/4 - c4_correct/16
print(f"c_0 = {H} - {c2_val}/4 - {c4_correct}/16 = {c0_correct}")
print(f"Direct c_0 = {total}")
print(f"Match: {abs(c0_correct - total) < 0.001}")

# Also for the random tournament
c0_B_correct = H_B - (12*t3_B-30)/4 - 120/16
print(f"\nRandom: c_0 = {H_B} - {12*t3_B-30}/4 - 120/16 = {c0_B_correct}")
print(f"Direct c_0 = {total_B}")
print(f"Match: {abs(c0_B_correct - total_B) < 0.001}")

# So the issue was: c_{n-1} = n! (not (n-1)!).
# This is because c_{n-1} = sum_P e_0(s_P) = n! (e_0 = 1 for all P).
# But in the THM-055 table: tr(c_6) = 720 at n=7. 720 = 6! = (n-1)!.
# That contradicts c_{n-1} = n!.

# Wait: at n=7, the r-polynomial has c_0, c_2, c_4, c_6.
# c_6 = sum_P e_0 = n! = 5040 = 7!. But the verified formula said c_6 = 720 = 6!.

# There must be a normalization issue. Let me check the convention in THM-055.

# From THM-055: tr(c_{n-1-2k}) = sum_P e_{2k}(s_P).
# At n=7: c_6 = tr(c_6) = sum_P e_0(s_P) = sum_P 1 = 7! = 5040.
# But the verified table says c_6 = 720.

# This means the table uses a DIFFERENT convention.
# Maybe: M(r) = (1/n!) * sum_P prod(r + s_i)?
# Or: the "c_k" in the table are the TRACE of the transfer matrix
# as defined in the project (which may have a different normalization).

# Actually, looking at the verify script: the c_k are extracted from
# tr(M(r)) where M(r) is the actual transfer matrix (not normalized by n!).
# The transfer matrix is defined via inclusion-exclusion, not as sum over permutations.

# The relation: tr(M(r)) = (sum_P prod(r+s_i)) / ???

# Hmm, I think the issue is that tr(M(r)) ≠ sum_P prod(r+s_i) directly.
# The transfer matrix involves an alternating sign IE formula, not a plain sum.

# From THM-053: tr(M(r)) = sum_P [sum_v (-1)^{pos(v,P)}] * prod(r+s_i).
# At odd n: sum_v (-1)^{pos(v,P)} = 1 for all P.
# So tr(M(r)) = sum_P prod(r+s_i).

# But that gives c_{n-1} = sum_P 1 = n!. And the verify gave c_6 = 720 = 6! at n=7.
# 720 ≠ 5040 = 7!.

# Unless the actual formula involves a normalization factor of 1/n?
# Or the diagonal formula involves 1/n somehow...

# Let me just recheck by direct computation.

print(f"\n\n=== CONVENTION CHECK ===")
print(f"Computing tr(M(r)) at n=5 for the transitive tournament")

# Build the actual transfer matrix at n=5 with r=0 and r=1
import numpy as np

n = 5
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        A[i][j] = 1

# Transfer matrix via IE
def compute_transfer_matrix(A, n, r):
    """Compute M[a,b](r) via the inclusion-exclusion formula."""
    full = (1 << n) - 1
    M = [[0.0]*n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            # IE sum over subsets not containing a or b
            # M[a][b] = sum_{S: a,b not in S} (-1)^|S| * (path from a to b through V\S)
            # This is complex. Let me use the DP approach instead.
            pass

    # Actually, use the direct DP method from the project
    # Forward DP
    dp_fwd = [[0.0]*n for _ in range(1 << n)]
    for v in range(n):
        dp_fwd[1 << v][v] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp_fwd[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r + (A[v][u] - 0.5)
                dp_fwd[mask | (1 << u)][u] += dp_fwd[mask][v] * wt
    # Backward DP
    dp_bwd = [[0.0]*n for _ in range(1 << n)]
    for v in range(n):
        dp_bwd[1 << v][v] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp_bwd[mask][v] == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r + (A[u][v] - 0.5)
                dp_bwd[mask | (1 << u)][u] += dp_bwd[mask][v] * wt
    # Trace
    tr_val = 0.0
    for v in range(n):
        for mask_before in range(1 << n):
            if mask_before & (1 << v): continue
            mask_with_v = mask_before | (1 << v)
            if dp_fwd[mask_with_v][v] == 0: continue
            mask_after = full ^ mask_before
            if not (mask_after & (1 << v)): continue
            if dp_bwd[mask_after][v] == 0: continue
            k = bin(mask_before).count('1')
            tr_val += ((-1)**k) * dp_fwd[mask_with_v][v] * dp_bwd[mask_after][v]
    return tr_val

# Compute at several r values
r_values = [0.0, 0.3, 0.5]
for r in r_values:
    tr = compute_transfer_matrix(A, n, r)
    perm_sum = sum(
        eval("*".join(f"({r}+{A[p[k]][p[k+1]]}-0.5)" for k in range(n-1)))
        for p in permutations(range(n))
    )
    print(f"  r={r}: tr(M)={tr:.4f}, sum_P prod(r+s)={perm_sum:.4f}, ratio={perm_sum/tr if tr != 0 else 'N/A'}")

print("\nDONE")
