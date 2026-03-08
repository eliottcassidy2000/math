#!/usr/bin/env python3
"""
F(T,x) AND EHRHART THEORY / ORDER POLYTOPES

For a tournament T, the order polytope O(T) consists of points
(x_0, ..., x_{n-1}) in [0,1]^n such that x_i < x_j whenever i -> j in T.

The h*-vector of O(T) (Ehrhart h*-polynomial) should be related to
descents of linear extensions — and HPs of T are precisely the
linear extensions of the partial order induced by T (though tournaments
may have cycles, making this only valid for acyclic tournaments =
transitive tournaments).

For GENERAL tournaments (with cycles), T has no linear extensions.
But F(T,x) still makes sense: it counts HPs by ascent number.

QUESTION 1: Is F(T,x) at positive integer m related to counting
something geometric (lattice points, P-partitions, etc.)?

QUESTION 2: The Worpitzky coefficients w_k of F(T,x) are integers.
What is their combinatorial meaning?

QUESTION 3: F(T,-1) = ? For odd degree d=n-1, palindromic polynomials
satisfy p(-1) = 0. But F is NOT palindromic, so F(-1) can be nonzero.
"""
from itertools import permutations
from math import comb, factorial
import numpy as np
from collections import defaultdict, Counter
from fractions import Fraction

def all_tournaments(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

def forward_edge_poly(A, n):
    """F_k = #{HPs with k ascents}"""
    F = [0] * n
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        asc = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[asc] += 1
    return F

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))

def worpitzky_expand(F_coeffs):
    """Expand polynomial with given coefficients in Worpitzky basis.
    p(x) = sum_k w_k C(x+k, d) where d = len(F_coeffs)-1.
    Uses the standard inversion: evaluate at x=0,1,...,d and solve."""
    d = len(F_coeffs) - 1
    n = d + 1
    # Evaluate p at x = 0, 1, ..., d
    vals = [sum(F_coeffs[k] * x**k for k in range(n)) for x in range(n)]
    # Build matrix M[i,j] = C(i+j, d)
    M = [[comb(i+j, d) for j in range(n)] for i in range(n)]
    # Solve M * w = vals
    M = np.array(M, dtype=float)
    w = np.linalg.solve(M, vals)
    return [int(round(x)) for x in w]

# ===== SECTION 1: F(T, m) at integer m =====
print("=" * 70)
print("SECTION 1: F(T, m) AT POSITIVE INTEGERS")
print("=" * 70)
print("""
F(T, m) = sum_k F_k * m^k.
At m=0: F(T,0) = F_0 = #{HPs with 0 ascents} = #{decreasing HPs}
At m=1: F(T,1) = sum F_k = H = total HP count
At m=2: F(T,2) = sum F_k * 2^k
At m=-1: F(T,-1) = sum F_k * (-1)^k

For the EULERIAN polynomial A_n(x):
  A_n(0) = 1 (identity permutation has 0 ascents... wait, n=0 case)
  Actually A(n,0) = 1 always (one perm with 0 ascents: n,n-1,...,1)
  A_n(1) = n!
  A_n(-1) = sum A(n,k)(-1)^k = 0 for n >= 2 (by inclusion-exclusion)

For F(T,x), F(T,-1) = sum F_k (-1)^k. This is the "alternating sum"
of the ascent distribution. It need not be 0.
""")

for n in [3, 4, 5]:
    print(f"\nn={n}:")
    E = [eulerian_number(n, k) for k in range(n)]
    print(f"  Eulerian: A({n},k) = {E}")
    print(f"  A({n},-1) = {sum(E[k]*(-1)**k for k in range(n))}")

    F_neg1_vals = []
    F_neg1_counts = Counter()

    count = 0
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        Fm1 = sum(F[k] * (-1)**k for k in range(n))
        F_neg1_vals.append(Fm1)
        F_neg1_counts[Fm1] += 1

        if count < 5:
            vals_at = {m: sum(F[k]*m**k for k in range(n)) for m in range(-2, 5)}
            print(f"  T#{count}: F={F}, H={H}")
            print(f"    F(-1)={Fm1}, F(-2)={vals_at[-2]}, F(2)={vals_at[2]}, F(3)={vals_at[3]}")
        count += 1

    print(f"\n  F(T,-1) distribution: {dict(sorted(F_neg1_counts.items()))}")
    print(f"  Sum F(T,-1) over all T = {sum(F_neg1_vals)}")
    # Should be A_n(-1) * 2^{extra} = 0

    # Check if F(-1) is always even, or has a pattern
    print(f"  All F(-1) even? {all(v % 2 == 0 for v in F_neg1_vals)}")

# ===== SECTION 2: Worpitzky coefficients =====
print("\n\n" + "=" * 70)
print("SECTION 2: WORPITZKY COEFFICIENTS OF F(T,x)")
print("=" * 70)
print("""
Worpitzky identity: x^d = sum_{k=0}^{d} A(d,k) * C(x+k, d)

So F(T,x) = sum_k w_k C(x+k, d) where d = n-1.

Properties to check:
1. w_k are always integers (known: yes)
2. w_k >= 0? (known: no, can be negative)
3. sum w_k = ? (at x=1: sum w_k C(1+k,d))
4. w_0 = F(T,0) = F_0 (value at x=0, since C(k,d)=0 for k<d, C(d,d)=1 only for k=d)
   Wait: C(0+k, d) = C(k,d) which is 0 for k < d and 1 for k = d.
   So p(0) = w_d. Hence w_{d} = w_{n-1} = F_0.
5. w_{n-1} = F_0 = #{decreasing HPs}
""")

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}, d={d}:")

    all_w = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = worpitzky_expand(F)
        H = sum(F)
        all_w.append((H, F, w))

    # Print a few examples
    for i, (H, F, w) in enumerate(all_w[:8]):
        print(f"  T#{i}: F={F}, w={w}, H={H}")
        # Verify w_{n-1} = F_0
        assert w[n-1] == F[0], f"w_{{n-1}}={w[n-1]} != F_0={F[0]}"

    # Statistics on w
    w_arr = np.array([w for _, _, w in all_w])
    print(f"\n  w_k ranges:")
    for k in range(n):
        vals = w_arr[:, k]
        print(f"    w_{k}: min={int(vals.min())}, max={int(vals.max())}, mean={vals.mean():.4f}")

    # Check: w_0 = ?
    # At x=1: F(T,1) = H = sum_k w_k C(1+k, d)
    # C(1+k, d) = C(1+k, n-1)
    # For k=0: C(1, n-1) = 0 if n>2
    # For k=1: C(2, n-1) = 0 if n>3
    # Only k=n-2: C(n-1, n-1) = 1 and k=n-1: C(n, n-1) = n
    # So H = w_{n-2} + n * w_{n-1} = w_{n-2} + n * F_0
    # Hence w_{n-2} = H - n * F_0
    print(f"\n  Checking w_{{n-2}} = H - n*F_0:")
    for H, F, w in all_w[:5]:
        pred = H - n * F[0]
        print(f"    H={H}, F_0={F[0]}, w_{{n-2}}={w[n-2]}, H-n*F_0={pred}, match={w[n-2]==pred}")

# ===== SECTION 3: Can we get w_k from F directly? =====
print("\n\n" + "=" * 70)
print("SECTION 3: WORPITZKY INVERSION FORMULA")
print("=" * 70)
print("""
The Worpitzky coefficients w_k can be obtained by:
  w_k = sum_{j=0}^{d-k} (-1)^j C(d+1, j) * p(k-j)

where p(x) = F(T,x) and d = n-1. This is the "Worpitzky inversion".

Equivalently, w_k = Delta^d [p(x)] evaluated at x=k, where Delta is
the forward difference operator. But the standard formula uses:

  w_k = sum_{j=0}^{d} (-1)^{d-j} C(d, d-j) F(T, k+j-d+1)

Wait, let me use the proper inversion. If p(x) = sum w_k C(x+k,d), then
  w_k = sum_{j=0}^{d} (-1)^{d-j} C(d,j) * p(k-d+j)

Let me verify this computationally.
""")

for n in [4]:
    d = n - 1
    for A in [list(all_tournaments(n))[-1]]:
        F = forward_edge_poly(A, n)
        w = worpitzky_expand(F)

        print(f"n={n}: F={F}, w={w}")

        # Verify by substitution
        for x in range(n+2):
            val_F = sum(F[k] * x**k for k in range(n))
            val_w = sum(w[k] * comb(x+k, d) for k in range(n))
            print(f"  x={x}: F(x)={val_F}, sum w_k C(x+k,{d})={val_w}")

# ===== SECTION 4: F(T,x) palindromicity test =====
print("\n\n" + "=" * 70)
print("SECTION 4: IS F(T,x) PALINDROMIC?")
print("=" * 70)
print("""
CRITICAL CHECK: Earlier analysis claimed F(T,x) is palindromic.
F_k(T) = F_{n-1-k}(T) would mean #{HPs with k ascents} = #{HPs with k descents}.

This IS true for the Eulerian polynomial (Worpitzky identity gives this).
For a general tournament, this requires a bijection on HPs that swaps
ascents and descents. The map P -> reverse(P) does this, but reverse(P)
is a HP of T^op, not T. So F(T,x) = x^{n-1} F(T^op, 1/x).

F(T,x) is palindromic iff T and T^op have the same HP structure, which
happens iff T is "self-complementary" in some sense.
""")

for n in [4, 5]:
    d = n - 1
    pal_count = 0
    total = 0
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        is_pal = all(F[k] == F[d-k] for k in range(n))
        if is_pal:
            pal_count += 1
        total += 1

        if total <= 3:
            print(f"  n={n}, T#{total}: F={F}, palindromic={is_pal}")

    print(f"  n={n}: palindromic F: {pal_count}/{total}")

    # How many pairs (T, T^op) with same F?
    # Build T^op
    by_F = defaultdict(int)
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        by_F[tuple(F)] += 1

    # For each T, F(T^op, x) has F^op_k = F_{d-k}
    # So F(T^op) = reversed F(T)
    rev_match = 0
    for F_tuple, count in by_F.items():
        rev_F = tuple(reversed(F_tuple))
        if F_tuple == rev_F:
            rev_match += count
    print(f"  n={n}: F=reversed(F) count: {rev_match}/{total}")

# ===== SECTION 5: Deeper Worpitzky structure =====
print("\n\n" + "=" * 70)
print("SECTION 5: WORPITZKY COEFFICIENTS — DEEPER STRUCTURE")
print("=" * 70)

for n in [4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    all_data = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = worpitzky_expand(F)
        H = sum(F)

        # Count 3-cycles
        t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

        scores = tuple(sorted([sum(A[i]) for i in range(n)]))
        all_data.append((H, t3, scores, F, w))

    # Group by w (what determines w?)
    by_w = defaultdict(list)
    for H, t3, scores, F, w in all_data:
        by_w[tuple(w)].append((H, t3, scores, F))

    print(f"  Distinct w vectors: {len(by_w)}")
    print(f"  Distinct F vectors: {len(set(tuple(F) for _,_,_,F,_ in all_data))}")

    # w and F are related by a fixed linear transform — so #distinct should be same
    # unless there are collisions

    # What determines w beyond H and F_0?
    # We know w_{n-1} = F_0, w_{n-2} = H - n*F_0
    # What about w_{n-3}?
    if n >= 4:
        print(f"\n  w_{{n-3}} analysis (w_{n-3}):")
        # At x=2: F(T,2) = sum w_k C(2+k, d)
        # C(2+k, d) for each k:
        for k in range(n):
            print(f"    C(2+{k}, {d}) = {comb(2+k, d)}")

        # So F(2) = sum_k w_k C(2+k, d)
        # = w_0 C(2,d) + w_1 C(3,d) + ... + w_{d} C(d+2, d)
        # For n=4, d=3: C(2,3)=0, C(3,3)=1, C(4,3)=4, C(5,3)=10
        # F(2) = w_1 + 4*w_2 + 10*w_3
        # We know w_3 = F_0, w_2 = H - 4*F_0
        # F(2) = w_1 + 4*(H-4*F_0) + 10*F_0 = w_1 + 4H - 16*F_0 + 10*F_0 = w_1 + 4H - 6*F_0
        # w_1 = F(2) - 4H + 6*F_0 = F_0 + 2*F_1 + 4*F_2 + 8*F_3 - 4H + 6*F_0
        # = 7*F_0 + 2*F_1 + 4*F_2 + 8*F_3 - 4*(F_0+F_1+F_2+F_3)
        # = 3*F_0 - 2*F_1 + 0*F_2 + 4*F_3
        if n == 4:
            print(f"\n  For n=4: w_1 = F(2) - 4H + 6F_0 = 3F_0 - 2F_1 + 0F_2 + 4F_3")
            for H, t3, scores, F, w in all_data[:5]:
                pred = 3*F[0] - 2*F[1] + 0*F[2] + 4*F[3]
                print(f"    F={F}, w_1={w[1]}, pred={pred}, match={w[1]==pred}")

    # The relationship between w and (H, F_0, F(-1), etc)
    print(f"\n  w as function of F values at specific points:")
    for H, t3, scores, F, w in all_data[:8]:
        Fm1 = sum(F[k]*(-1)**k for k in range(n))
        F2 = sum(F[k]*2**k for k in range(n))
        Fm2 = sum(F[k]*(-2)**k for k in range(n))
        print(f"    F={F}, w={w}, F(-1)={Fm1}, F(2)={F2}")

# ===== SECTION 6: Connection to OCF =====
print("\n\n" + "=" * 70)
print("SECTION 6: WORPITZKY COEFFICIENTS AND OCF INVARIANTS")
print("=" * 70)
print("""
The OCF states H(T) = I(Omega(T), 2). The Worpitzky expansion gives
a finer decomposition of F(T,x) than just H = F(T,1).

QUESTION: Are the Worpitzky coefficients w_k of F(T,x) related to
the transfer matrix M or the OCF invariants?

Since w_{n-1} = F_0 and w_{n-2} = H - n*F_0, the first two "new"
coefficients are F_0 and H. The next ones (w_{n-3}, etc.) carry
additional information about the F-polynomial.
""")

for n in [5]:
    d = n - 1

    all_data = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = worpitzky_expand(F)
        H = sum(F)

        t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]))

        all_data.append((H, t3, F, w))

    # Are the w's determined by (H, t3)?
    by_Ht3 = defaultdict(set)
    for H, t3, F, w in all_data:
        by_Ht3[(H, t3)].add(tuple(w))

    print(f"\nn={n}: Is w determined by (H, t3)?")
    ambig = sum(1 for v in by_Ht3.values() if len(v) > 1)
    print(f"  Ambiguous: {ambig}/{len(by_Ht3)}")

    # Is w determined by F? (Of course — it's a linear transform)
    # But is F determined by (H, F_0, F(-1))?
    by_HF0Fm1 = defaultdict(set)
    for H, t3, F, w in all_data:
        Fm1 = sum(F[k]*(-1)**k for k in range(n))
        by_HF0Fm1[(H, F[0], Fm1)].add(tuple(F))

    print(f"\n  Is F determined by (H, F_0, F(-1))?")
    ambig2 = sum(1 for v in by_HF0Fm1.values() if len(v) > 1)
    print(f"  Ambiguous: {ambig2}/{len(by_HF0Fm1)}")

    # Correlation between w_k and H, t3
    w_arr = np.array([w for _, _, _, w in all_data])
    H_arr = np.array([H for H, _, _, _ in all_data])
    t3_arr = np.array([t3 for _, t3, _, _ in all_data])

    print(f"\n  Correlations:")
    for k in range(n):
        corr_H = np.corrcoef(w_arr[:,k], H_arr)[0,1]
        corr_t3 = np.corrcoef(w_arr[:,k], t3_arr)[0,1]
        print(f"    Corr(w_{k}, H) = {corr_H:.4f}, Corr(w_{k}, t3) = {corr_t3:.4f}")

# ===== SECTION 7: Order polytope connection =====
print("\n\n" + "=" * 70)
print("SECTION 7: ORDER POLYTOPE / P-PARTITION CONNECTION")
print("=" * 70)
print("""
For a POSET P, the order polytope O(P) = {x in [0,1]^n : x_i <= x_j if i < j in P}.
Its Ehrhart polynomial is related to P-partitions.
The h*-vector of O(P) equals the descent polynomial of linear extensions.

A tournament T defines a COMPLETE BINARY RELATION (total, antisymmetric).
If T is acyclic = transitive, it IS a total order (poset), and:
  O(T) = {x : x_{sigma(1)} <= ... <= x_{sigma(n)}} (a simplex)
  The only linear extension is sigma, so h* = x^{n-1} or x^0 depending.

For a GENERAL tournament with cycles, T is not a poset.
But we can still ask: what is F(T,m) for positive integers m?

Idea: F(T,m) = sum_k F_k m^k counts HPs weighted by m^{ascents}.
By the transfer matrix method (separate computation):
  F(T,m) = #{(P, f) : P is HP of T, f:[n]->{1,...,m} weakly increasing along P
            with strict increase at ascents}

This is exactly Stanley's P-partition theory! Even though T isn't a poset,
a HP defines a LABELED total order, and f is a (P,w)-partition.

So: F(T,m) = sum over HPs P of (multiset coefficient related to m and ascents of P).
""")

# Verify the P-partition interpretation
for n in [4]:
    d = n - 1
    print(f"\nn={n}:")

    for idx, A in enumerate(all_tournaments(n)):
        if idx >= 3: break
        F = forward_edge_poly(A, n)
        H = sum(F)

        print(f"\n  Tournament #{idx}: F={F}, H={H}")

        for m in [1, 2, 3, 4]:
            # F(T,m) by polynomial evaluation
            Fm_poly = sum(F[k] * m**k for k in range(n))

            # F(T,m) by P-partition counting:
            # For each HP P, count f:[n]->{1,...,m} such that
            #   f(P_i) <= f(P_{i+1}) always
            #   f(P_i) < f(P_{i+1}) when P_i < P_{i+1} (ascent)
            # This is equivalent to: for each HP, the number of ways
            # to assign values 1..m to positions such that:
            #   at a descent: f(P_i) >= f(P_{i+1}) (i.e., f weakly decreasing)
            #   at an ascent: f(P_i) < f(P_{i+1}) (strictly increasing)
            # Wait, I need to be more careful about the convention.

            # Actually, the standard result is:
            # sum_{sigma} x^{des(sigma)} = h*-vector =>
            # Ehrhart at m = sum_{sigma} C(m + n - 1 - des(sigma), n)
            #              = sum_k h_k C(m+n-1-k, n)

            # For our F(T,x) = sum_k F_k x^k:
            # F(T,m) should equal sum over HPs of C(m + d - asc, d)?
            # No... Let me think.

            # By Worpitzky: x^d = sum A(d,k) C(x+k,d)
            # So F(T,x) = sum_k F_k x^k = sum_k F_k sum_j A(k... no this mixes things up.

            # The direct computation: for each HP P with a ascents,
            # the contribution at m should be... let me just check numerically.

            # For a HP with a ascents and d-a descents,
            # #{weakly increasing f with strict at ascents} = C(m-1+d-a, d) ... maybe?
            # No, let me think about it as: place n values from {1,...,m} in n positions
            # such that at each step, the value weakly increases, and strictly at ascents.
            # This is equivalent to choosing n values from {1,...,m} and a non-decreasing
            # sequence with strict increases at specified positions.
            # The number of such sequences = C(m + n - 1 - a, n - 1)?

            # Stars and bars: n-1 "steps", a of which are strict, (n-1-a) weakly.
            # For weakly increasing with a strict steps from {1,...,m}:
            # Substitute y_i = f(P_i) - (# strict steps before i): get weakly increasing from reduced range
            # This gives C(m - a + n - 1, n)? Not quite...

            # Actually for a weakly increasing sequence of length n from {1,...,m} with
            # exactly a specified positions being strict: it's C(m - 1 + n - a, n - 1)
            # = C(m + d - a, d) where d = n - 1.

            # Wait no. Weakly increasing from {1,...,m}: C(m+n-1, n) total.
            # With k strict positions: subtract 1 at those -> C(m-1+n-k, n-1)?
            # Hmm, let me just verify directly for small cases.

            Fm_ppart = 0
            for P in permutations(range(n)):
                if not all(A[P[i]][P[i+1]] for i in range(n-1)):
                    continue
                asc = sum(1 for i in range(n-1) if P[i] < P[i+1])
                # Count weakly increasing f with strict at ascents
                Fm_ppart += comb(m - 1 + d - asc, d)

            match = "OK" if Fm_poly == Fm_ppart else "MISMATCH"
            print(f"    m={m}: F(T,m)={Fm_poly}, P-part={Fm_ppart} {match}")

# ===== SECTION 8: The correct P-partition interpretation =====
print("\n\n" + "=" * 70)
print("SECTION 8: FINDING THE RIGHT P-PARTITION FORMULA")
print("=" * 70)

for n in [4]:
    d = n - 1
    A = list(all_tournaments(n))[7]  # pick one
    F = forward_edge_poly(A, n)
    print(f"n={n}: F={F}, H={sum(F)}")

    # For each HP, try various formulas
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        asc = sum(1 for i in range(n-1) if P[i] < P[i+1])
        desc = d - asc
        print(f"  HP {P}, asc={asc}, desc={desc}")

        # Try each formula for various m
        for m in [1, 2, 3]:
            # Direct count by brute force
            count = 0
            for f in range(m**n):
                vals = [(f // m**i) % m + 1 for i in range(n)]
                valid = True
                for i in range(n-1):
                    if P[i] < P[i+1]:  # ascent: strict increase
                        if vals[i] >= vals[i+1]: valid = False; break
                    else:  # descent: weak increase
                        if vals[i] > vals[i+1]: valid = False; break
                if valid:
                    count += 1

            # Formula attempts
            f1 = comb(m + d - asc, d)
            f2 = comb(m - 1 + d - asc, d)
            f3 = comb(m + desc, d)
            f4 = comb(m - 1 + desc, d)
            f5 = comb(m + asc - 1, d)

            matches = []
            for name, val in [("C(m+d-a,d)", f1), ("C(m-1+d-a,d)", f2),
                              ("C(m+desc,d)", f3), ("C(m-1+desc,d)", f4),
                              ("C(m+asc-1,d)", f5)]:
                if val == count:
                    matches.append(name)

            print(f"    m={m}: brute={count}, matches: {matches}")
        break  # just first HP

print("\n\nDone.")
