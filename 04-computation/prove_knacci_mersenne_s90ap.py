#!/usr/bin/env python3
"""
prove_knacci_mersenne_s90ap.py — Prove Tr(M_k^k) = 2^k - 1
opus-2026-03-16-S90ap

THEOREM: For the k-nacci companion matrix M_k (k x k), Tr(M_k^k) = 2^k - 1.

This connects:
  - Powers of 2 on the Cayley line (tanh at half-bit spacing)
  - Mersenne numbers as Cayley numerators
  - k-nacci transfer matrix traces
  - Forbidden H values of k-nacci theories
"""

import numpy as np
from math import log2

# ========================================================================
# Part 1: Define the k-nacci companion matrix
# ========================================================================

def knacci_matrix(k):
    """
    The k-nacci companion matrix M_k is k x k:
      Row 0: [1, 1, 1, ..., 1]  (all ones)
      Row i (i>0): e_{i-1}^T (shift matrix)

    Equivalently:
      M[0][j] = 1 for all j
      M[i][i-1] = 1 for i = 1,...,k-1
      All other entries 0.

    Characteristic polynomial: lambda^k - lambda^{k-1} - ... - lambda - 1
    """
    M = np.zeros((k, k), dtype=np.int64)
    for j in range(k):
        M[0][j] = 1
    for i in range(1, k):
        M[i][i-1] = 1
    return M

print("""


        PROOF: Tr(M_k^k) = 2^k - 1


Part 1: The k-nacci companion matrix
""")

for k in [2, 3, 4, 5]:
    M = knacci_matrix(k)
    print(f"M_{k} =")
    for row in M:
        print("  ", list(row))
    Mk = np.linalg.matrix_power(M, k)
    tr = int(round(np.trace(Mk)))
    print(f"  Tr(M_{k}^{k}) = {tr}, 2^{k}-1 = {2**k - 1}, match: {tr == 2**k - 1}")
    print()

# ========================================================================
# Part 2: Build the trace sequence and verify the recurrence
# ========================================================================

print("=" * 70)
print("Part 2: Full trace sequences Tr(M_k^n) for n = 0, 1, ..., 2k")
print("=" * 70)

for k in [2, 3, 4, 5, 6, 7]:
    M = knacci_matrix(k)
    traces = []
    for n in range(2*k + 1):
        Mn = np.linalg.matrix_power(M, n)
        traces.append(int(round(np.trace(Mn))))
    print(f"\nk={k}: Tr(M^n) for n=0..{2*k}:")
    print(f"  {traces}")
    print(f"  Tr(M^0) = {traces[0]} = k = {k}")
    print(f"  Tr(M^1) = {traces[1]} = 1")
    print(f"  Tr(M^k) = {traces[k]} = 2^k-1 = {2**k-1}")

    # Verify k-nacci recurrence: a(n) = a(n-1) + a(n-2) + ... + a(n-k)
    ok = True
    for n in range(k, 2*k + 1):
        s = sum(traces[n-j] for j in range(1, k+1))
        if s != traces[n]:
            ok = False
            print(f"  RECURRENCE FAILS at n={n}: sum={s}, actual={traces[n]}")
    if ok:
        print(f"  k-nacci recurrence verified for n=k..{2*k}")

# ========================================================================
# Part 3: Algebraic proof via Newton's identity
# ========================================================================

print(f"""

{"=" * 70}
Part 3: ALGEBRAIC PROOF via Newton's identity
{"=" * 70}

The k-nacci companion matrix M_k has characteristic polynomial:

  p(x) = x^k - x^(k-1) - x^(k-2) - ... - x - 1.

Let e_j = (-1)^j * (elementary symmetric polynomial of degree j in eigenvalues).
From the char poly: e_1 = -1, e_2 = -1, ..., e_(k-1) = -1, e_k = -1.

Wait -- let's be more careful. If p(x) = x^k + c_(k-1)*x^(k-1) + ... + c_0,
then e_j = (-1)^j * c_(k-j).

For k-nacci: p(x) = x^k - x^(k-1) - x^(k-2) - ... - x - 1.
So c_(k-1) = -1, c_(k-2) = -1, ..., c_1 = -1, c_0 = -1.

Elementary symmetric polynomials (with signs from Vieta):
  e_1 = sum of eigenvalues = -c_(k-1) = 1    (= Tr(M) = 1)
  e_2 = sum of products of pairs = c_(k-2) = -1
  e_3 = -c_(k-3) = 1
  ...
  e_j = (-1)^(j+1) * c_(k-j) = (-1)^(j+1) * (-1) = (-1)^j

Actually, let me use Newton's identities directly.

Newton's identities relate power sums p_n = Tr(M^n) to the
coefficients of the characteristic polynomial.
""")

# Verify Newton's identities for each k
from sympy import symbols, Poly, Matrix as SMatrix

print("Verifying via exact symbolic computation (sympy):")
print()

for k in [2, 3, 4, 5, 7, 11, 13, 17, 19]:
    # Build M_k as a sympy matrix
    M_data = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == 0:
                row.append(1)
            elif j == i - 1:
                row.append(1)
            else:
                row.append(0)
        M_data.append(row)

    M_sym = SMatrix(M_data)
    Mk_sym = M_sym ** k
    tr = Mk_sym.trace()
    mersenne = 2**k - 1
    match = (tr == mersenne)
    print(f"  k={k:2d}: Tr(M_{k}^{k}) = {tr}, 2^{k}-1 = {mersenne}, match: {match}")

# ========================================================================
# Part 4: The proof
# ========================================================================

print(f"""

{"=" * 70}
Part 4: THE PROOF
{"=" * 70}

THEOREM: Tr(M_k^k) = 2^k - 1 for all k >= 1.

PROOF (by direct computation using the recurrence):

The trace sequence t(n) = Tr(M_k^n) satisfies:
  (i)   t(0) = k  (trace of identity = dimension)
  (ii)  t(n) = t(n-1) + t(n-2) + ... + t(n-k)  for n >= 1
        (Newton's identity + the fact that ALL char poly coefficients are -1)

The second initial condition and early values come from Newton's identity:
  t(1) = e_1 = 1
  t(2) = e_1*t(1) - 2*e_2 = 1*1 - 2*(-1) = 3

In general, Newton's identities for this char poly give:
  t(n) = t(n-1) + t(n-2) + ... + t(n-k)  for n >= k+1  (standard)

But ALSO for 1 <= n <= k, using the initial condition t(0) = k and
the recursion with the convention that t(j) = 0 for j < 0 EXCEPT t(0) = k.

CLAIM: Define s(n) = 2^n - 1 for n >= 1, s(0) = k.
Then s(n) satisfies the same recurrence for 1 <= n <= k.

Check s(0) = k. YES (by definition).
Check s(1) = s(0) = k? No, s(1) = 1, and the recurrence gives
  s(1) = s(0) = k, which only works if k = 1.

So the proof needs more care. Let me compute the ACTUAL initial values.
""")

# Compute initial values carefully
print("Actual trace sequences (exact):")
for k in [2, 3, 4, 5, 6, 7, 8]:
    M_data = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == 0:
                row.append(1)
            elif j == i - 1:
                row.append(1)
            else:
                row.append(0)
        M_data.append(row)
    M_sym = SMatrix(M_data)

    traces = []
    for n in range(k+1):
        traces.append(int(M_sym.pow(n).trace()))
    print(f"  k={k}: t(0..{k}) = {traces}")

print(f"""

KEY OBSERVATION: For the k-nacci companion matrix,

  t(0) = k
  t(1) = 1
  t(n) = t(n-1) + t(n-2) + ... + t(max(n-k, 0))  for n >= 1,
         but with Newton's correction for n <= k.

Let me instead use a DIRECT approach: compute M^k explicitly.
""")

# ========================================================================
# Part 5: Direct matrix power computation
# ========================================================================

print("=" * 70)
print("Part 5: DIRECT PROOF via (M+I) and binomial-like expansion")
print("=" * 70)

# Key insight: M_k = J + N where J has first row all 1s, N is the nilpotent shift
# Actually, let's try: sum of eigenvalues of M^k

# Better: note that M_k has char poly x^k - x^{k-1} - ... - x - 1
# Multiply by (x-1): (x-1)(x^k - x^{k-1} - ... - 1) = x^{k+1} - 2x^k + 1
# So the eigenvalues of M_k satisfy x^{k+1} - 2x^k + 1 = 0, OR x = 1.

# Cayley-Hamilton: M_k^k = M_k^{k-1} + M_k^{k-2} + ... + M_k + I
# So Tr(M_k^k) = Tr(M_k^{k-1}) + ... + Tr(M_k) + Tr(I) = Tr(M_k^{k-1}) + ... + 1 + k

print("""
APPROACH: Use the augmented characteristic polynomial.

The k-nacci char poly is p(x) = x^k - x^{k-1} - ... - x - 1.

Multiply by (x-1):
  (x-1) * p(x) = x^{k+1} - 2x^k + 1.

So the eigenvalues L_1, ..., L_k satisfy L^{k+1} - 2*L^k + 1 = 0 (plus L=1).

Therefore: L_i^{k+1} = 2*L_i^k - 1 for each eigenvalue.

Summing: Tr(M^{k+1}) = 2*Tr(M^k) - k  (NO! sum of 1's is number of eigenvalues = k)

Wait: sum_i L_i^{k+1} = 2 * sum_i L_i^k - sum_i 1 = 2*Tr(M^k) - k.

So: Tr(M^{k+1}) = 2*Tr(M^k) - k.

This is a RECURRENCE IN k (not in n)!

If Tr(M_k^k) = 2^k - 1, then Tr(M_k^{k+1}) = 2(2^k - 1) - k = 2^{k+1} - 2 - k.

Hmm, that's Tr(M_k^{k+1}), still for the SAME matrix M_k (dimension k).
We want Tr(M_{k+1}^{k+1}) for the DIFFERENT matrix M_{k+1} (dimension k+1).

These are different matrices! So the augmented poly trick doesn't directly
give a recurrence across different k's.
""")

# Let me try a completely different approach: explicit formula for Tr(M^n)
# via the generating function.

print("=" * 70)
print("Part 6: PROOF BY GENERATING FUNCTION")
print("=" * 70)

print("""
The trace Tr(M_k^n) has generating function:

  sum_{n>=0} Tr(M_k^n) * x^n = k - d/dx log(det(I - x*M_k)) ...

Actually, sum_{n>=0} Tr(M^n) x^n = Tr((I - xM)^{-1}) = sum of 1/(1 - lambda_i * x)

For the k-nacci matrix, the char poly is p(L) = L^k - L^{k-1} - ... - L - 1.

The trace generating function is:
  f(x) = sum_i 1/(1 - L_i * x) = -x * p'(1/x) / p(1/x)  ... (up to signs)

Actually, let's use: det(I - xM) = x^k * p(1/x) (the reversed char poly).

p(L) = L^k - L^{k-1} - ... - 1.
det(I - xM) = x^k * p(1/x) = x^k * (x^{-k} - x^{-k+1} - ... - 1)
            = 1 - x - x^2 - ... - x^k.

So det(I - xM) = 1 - x - x^2 - ... - x^k = 1 - x(1 + x + ... + x^{k-1})
              = 1 - x * (x^k - 1)/(x - 1).

And -d/dx log det(I - xM) = (1 + 2x + 3x^2 + ... + k*x^{k-1}) / (1 - x - ... - x^k).

The generating function for Tr(M^n) is:
  T(x) = k / (1 - x - x^2 - ... - x^k)  ...

No. The standard formula is:
  sum_{n>=0} Tr(M^n) x^n = d/dx [-log det(I - xM)] ... (incorrect)

Let me just verify: Tr((I-xM)^{-1}) = sum Tr(M^n) x^n.
And Tr((I-xM)^{-1}) = sum_i 1/(1-L_i x).

Also: sum_i 1/(1-L_i x) = -sum_i L_i / (L_i - 1/x) ... not helpful.

Better: use sum_i 1/(1-L_i x) = (d/dx)[sum_i (-1/L_i) log(1-L_i x)] ... also awkward.

Simpler: sum_i 1/(1-L_i x) = [derivative of log prod(1-L_i x)] / (-1)
Nope. Let me just use the recurrence.
""")

# ========================================================================
# Part 7: CLEAN PROOF
# ========================================================================

print("=" * 70)
print("Part 7: CLEAN PROOF by induction on n within fixed k")
print("=" * 70)

print("""
THEOREM: For the k-nacci companion matrix M_k (k x k, k >= 2),
         Tr(M_k^k) = 2^k - 1.

PROOF:

Step 1: The trace sequence t(n) = Tr(M_k^n) satisfies the k-nacci recurrence
  t(n) = t(n-1) + t(n-2) + ... + t(n-k)  for n >= k
with initial conditions determined by Newton's identities from the char poly.

Step 2: The char poly of M_k is p(L) = L^k - L^{k-1} - ... - L - 1.
The coefficients are c_{k-1} = -1, c_{k-2} = -1, ..., c_0 = -1.

Step 3: Newton's identities give t(n) for n = 0, ..., k-1:
  t(0) = k
  t(1) = 1 (= -c_{k-1} = -(-1) = 1)
  For n = 2: t(2) = t(1)*(-c_{k-1}) + (-2*c_{k-2}) = 1*1 + 2 = 3
  For n = 3: t(3) = t(2)*1 + t(1)*1 + (-3*c_{k-3}) = 3 + 1 + 3 = 7  (if k >= 4)
  ...
  For general n <= k-1:
    t(n) = sum_{j=1}^{n-1} t(n-j)*(-c_{k-j}) + n*(-c_{k-n})
         = sum_{j=1}^{n-1} t(n-j)*1 + n*1  (since -c_{k-j} = 1 for all j)
         = sum_{j=1}^{n-1} t(n-j) + n.
""")

# Let me verify this formula and check that it gives 2^n - 1
print("Verification that t(n) = 2^n - 1 for n = 1, ..., k:")
print()

for k in [3, 5, 7, 10, 15]:
    # Build trace sequence using Newton's identities
    t = [0] * (k + 1)
    t[0] = k

    for n in range(1, k + 1):
        if n <= k - 1:
            # Newton's identity for n < k:
            # t(n) = sum_{j=1}^{n-1} t(n-j) + n
            # (using -c_{k-j} = 1 for j = 1, ..., n, and n*(-c_{k-n}) = n)
            s = sum(t[n - j] for j in range(1, n))
            t[n] = s + n
        else:
            # n = k: use the recurrence
            # t(k) = sum_{j=1}^{k} t(k-j)
            # But also Newton: t(k) = sum_{j=1}^{k-1} t(k-j) + k*(-c_0)
            # c_0 = -1, so -c_0 = 1, so t(k) = sum_{j=1}^{k-1} t(k-j) + k
            s = sum(t[k - j] for j in range(1, k))
            t[n] = s + k

    expected = [k] + [2**n - 1 for n in range(1, k + 1)]
    match = (t == expected)
    print(f"  k={k:2d}: t(0..k) = {t}")
    print(f"         expected = {expected}")
    print(f"         match: {match}")
    print()

print("""
Now I need to prove: if t(0) = k and t(n) = sum_{j=1}^{n-1} t(n-j) + n for 1 <= n <= k,
then t(n) = 2^n - 1 for 1 <= n <= k.

PROOF by induction on n:

Base case: t(1) = (empty sum) + 1 = 1 = 2^1 - 1. CHECK.

Inductive step: Assume t(m) = 2^m - 1 for 1 <= m <= n-1 (and n <= k).
Then:
  t(n) = sum_{j=1}^{n-1} t(n-j) + n
       = sum_{m=1}^{n-1} t(m) + n           (substituting m = n-j)
       = sum_{m=1}^{n-1} (2^m - 1) + n
       = (2 + 4 + ... + 2^{n-1}) - (n-1) + n
       = (2^n - 2) - (n-1) + n
       = 2^n - 2 - n + 1 + n
       = 2^n - 1.

QED.

Note: We never needed the value t(0) = k in the inductive argument!
The formula t(n) = 2^n - 1 is INDEPENDENT of k, for all 1 <= n <= k.

This means:
  Tr(M_k^n) = 2^n - 1 for ALL 1 <= n <= k, not just n = k.

In particular:
  Tr(M_3^3) = 7 = 2^3 - 1.   (tribonacci)
  Tr(M_5^5) = 31 = 2^5 - 1.  (pentanacci)
  Tr(M_7^7) = 127 = 2^7 - 1. (heptanacci)
  etc.

But also:
  Tr(M_7^3) = 7 (same as tribonacci at n=3!)
  Tr(M_7^5) = 31 (same as pentanacci at n=5!)

The trace at power n depends ONLY on n, not on k, as long as k >= n.

COROLLARY: For any k >= n >= 1, Tr(M_k^n) = 2^n - 1 = the n-th Mersenne number.
""")

# Verify the corollary
print("=" * 70)
print("Verification of COROLLARY: Tr(M_k^n) = 2^n - 1 for all 1 <= n <= k")
print("=" * 70)
print()
header = 'k\\n'
print(f"{header:>4}", end="")
for n in range(1, 11):
    print(f"{n:>6}", end="")
print()
print("-" * 64)

for k in range(2, 11):
    print(f"{k:>4}", end="")
    M_data = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == 0:
                row.append(1)
            elif j == i - 1:
                row.append(1)
            else:
                row.append(0)
        M_data.append(row)
    M_sym = SMatrix(M_data)
    for n in range(1, 11):
        if n > k:
            # Beyond k, the traces diverge
            tr = int(M_sym.pow(n).trace())
            mersenne = 2**n - 1
            if tr == mersenne:
                print(f"  {tr:>4}", end="")
            else:
                print(f" *{tr:>3}", end="")  # mark divergence
        else:
            tr = int(M_sym.pow(n).trace())
            print(f"  {tr:>4}", end="")
    print()

print()
print("(Stars mark n > k where trace may diverge from 2^n - 1)")

# What happens beyond n = k?
print()
print("=" * 70)
print("What happens at n = k+1?")
print("=" * 70)

for k in range(2, 9):
    M_data = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == 0:
                row.append(1)
            elif j == i - 1:
                row.append(1)
            else:
                row.append(0)
        M_data.append(row)
    M_sym = SMatrix(M_data)
    tr_k = int(M_sym.pow(k).trace())
    tr_k1 = int(M_sym.pow(k+1).trace())
    mersenne_k1 = 2**(k+1) - 1
    deficit = mersenne_k1 - tr_k1
    print(f"  k={k}: Tr(M^k) = {tr_k} = 2^k-1, Tr(M^{{k+1}}) = {tr_k1}, 2^{{k+1}}-1 = {mersenne_k1}, deficit = {deficit}")

print("""

At n = k+1, the trace is no longer 2^{k+1}-1 because the k-nacci
recurrence kicks in instead of Newton's identity.

For n = k+1: t(k+1) = t(k) + t(k-1) + ... + t(1)  (k-nacci recurrence, NO +n term)
           = (2^k - 1) + (2^{k-1} - 1) + ... + (2^1 - 1)
           = (2^{k+1} - 2) - k
           = 2^{k+1} - k - 2.

Deficit from 2^{k+1}-1 = k + 1.

So Tr(M_k^{k+1}) = 2^{k+1} - 1 - (k+1) = 2^{k+1} - k - 2.

The Mersenne formula BREAKS by exactly (k+1) at the first step beyond k.
""")

# Verify
for k in range(2, 9):
    M_data = []
    for i in range(k):
        row = []
        for j in range(k):
            if i == 0:
                row.append(1)
            elif j == i - 1:
                row.append(1)
            else:
                row.append(0)
        M_data.append(row)
    M_sym = SMatrix(M_data)
    tr_k1 = int(M_sym.pow(k+1).trace())
    predicted = 2**(k+1) - k - 2
    print(f"  k={k}: Tr(M^{{k+1}}) = {tr_k1}, predicted 2^{{k+1}}-k-2 = {predicted}, match: {tr_k1 == predicted}")

print(f"""

{"=" * 70}
SUMMARY
{"=" * 70}

THEOREM (k-nacci Mersenne Identity):
  For the k-nacci companion matrix M_k and any 1 <= n <= k:
    Tr(M_k^n) = 2^n - 1.

  In particular, Tr(M_k^k) = 2^k - 1 (the k-th Mersenne number).

PROOF:
  Newton's identities, applied to the k-nacci characteristic polynomial
  p(x) = x^k - x^(k-1) - ... - x - 1 (all coefficients = -1),
  give the recurrence:
    t(n) = t(n-1) + t(n-2) + ... + t(1) + n    for 1 <= n <= k.

  By induction: t(n) = sum_{{m=1}}^{{n-1}} (2^m - 1) + n = 2^n - 1. QED.

COROLLARY: The trace is INDEPENDENT of k for 1 <= n <= k.
  All k-nacci matrices agree on Tr(M^n) = 2^n-1 within the "Newton range" n <= k.

COROLLARY: At n = k+1, the trace drops by exactly (k+1):
  Tr(M_k^{{k+1}}) = 2^{{k+1}} - k - 2.

CONNECTION TO TOURNAMENTS:
  The tribonacci matrix (k=3) is the transfer matrix for tournament path counting.
  Tr(M_3^3) = 7 = the first FORBIDDEN Hamiltonian path count.
  Tr(M_k^k) = 2^k-1 gives the forbidden value for the k-nacci generalization.
  When 2^k-1 is PRIME (Mersenne prime), the forbidden value is a prime number.

CONNECTION TO THE CAYLEY LINE:
  tanh(k*ln(2)/2) = (2^k-1)/(2^k+1) places 2^k at rapidity k/2 bits.
  The NUMERATOR of this Cayley address is 2^k-1 = Tr(M_k^k).
  So the forbidden k-nacci value lives at the k/2-bit mark on the Cayley line.

THE THREAD:
  tanh --> powers of 2 --> Mersenne numbers --> k-nacci traces --> forbidden values.
  One identity: Tr(M_k^k) = 2^k - 1.
  One proof: induction on Newton's identities with all-(-1) coefficients.
  The entire connection follows from the fact that the k-nacci char poly
  has ALL coefficients equal to -1.
""")

print("opus-2026-03-16-S90ap: PROVE Tr(M_k^k) = 2^k - 1")
