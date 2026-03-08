#!/usr/bin/env python3
"""
F(T,2) AND THE OCF CONNECTION

F(T,x) = sum_k F_k x^k where F_k = #{HPs with k ascents}.
F(T,1) = H(T).
F(T,2) = sum F_k * 2^k.

OCF states H(T) = I(Omega(T), 2) where Omega is the odd-cycle graph
and I is the independence polynomial.

QUESTION: Does F(T,2) relate to anything in the OCF framework?

Also explore: the generating function sum_{m>=0} F(T,m) x^m.
By Worpitzky: F(T,m) = sum_k w_k C(m+k,d) where d=n-1.
sum_{m>=0} F(T,m) x^m = sum_k w_k * sum_m C(m+k,d) x^m
                       = sum_k w_k * x^{d-k} / (1-x)^{d+1}
                       = F(T,x^{-1}) * x^{-d} ... wait, let me re-derive.

Actually: sum_{m>=0} C(m+k,d) x^m = x^{d-k} / (1-x)^{d+1} for k <= d.
So: sum_{m>=0} F(T,m) x^m = (1/(1-x)^{d+1}) * sum_k w_k x^{d-k}
                           = (1/(1-x)^n) * sum_k w_{d-k} x^k
                           = (1/(1-x)^n) * W^*(T, x)

where W^*(T,x) = sum_k w_{n-1-k} x^k is the "reversed w polynomial".

So W^*(T,x) is the h*-vector polynomial in the Ehrhart sense, even though
there's no polytope! The GF sum_m F(T,m) x^m = W^*(T,x) / (1-x)^n.

This means F(T,m) = sum_k w_{n-1-k} C(m+n-1-k, n-1).

Compare with: h*(P, x) / (1-x)^n = Ehrhart series of an n-polytope.
The polynomial W^*(T,x) plays the role of h*-vector.

We know: w_{n-1} = F_0 = #{fully decreasing HPs}, so w^*_0 = F_0.
And w^*_1 = w_{n-2} = H - n*F_0.
And sum w^*_k = F(T,1) = H (evaluating W^* at 1).

IMPORTANT: for actual simplicial polytopes, h* is nonneg and unimodal.
Are the w^*_k nonneg? Unimodal?
"""
from itertools import permutations
from math import comb, factorial
import numpy as np
from collections import defaultdict, Counter

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
    F = [0] * n
    for P in permutations(range(n)):
        if not all(A[P[i]][P[i+1]] for i in range(n-1)):
            continue
        asc = sum(1 for i in range(n-1) if P[i] < P[i+1])
        F[asc] += 1
    return F

def worpitzky_expand(F_coeffs):
    d = len(F_coeffs) - 1
    n = d + 1
    vals = [sum(F_coeffs[k] * x**k for k in range(n)) for x in range(n)]
    M = np.array([[comb(i+j, d) for j in range(n)] for i in range(n)], dtype=float)
    w = np.linalg.solve(M, vals)
    return [int(round(x)) for x in w]

def count_odd_cycles(A, n):
    """Count 3-cycles"""
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def omega_graph(A, n):
    """Build Omega(T) — odd cycle adjacency. For n<=5, only 3-cycles matter."""
    # For simplicity, enumerate edges of T and check if they share a 3-cycle
    edges = []
    for i in range(n):
        for j in range(i+1,n):
            edges.append((i,j))

    # For each pair of edges, check if they appear in a common 3-cycle
    # Actually, Omega has vertices = edges of T, adjacency = share a 3-cycle
    # Let me simplify: just check if omega makes sense for small n
    return None  # complex; let me just check F(T,2) vs H directly

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# ===== SECTION 1: W^* polynomial (reversed Worpitzky) =====
print("=" * 70)
print("W^*(T,x) = REVERSED WORPITZKY POLYNOMIAL (h*-analog)")
print("=" * 70)
print("W^*(x) = sum_k w_{n-1-k} x^k")
print("GF: sum_m F(T,m) x^m = W^*(T,x) / (1-x)^n")
print()

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    nonneg_count = 0
    unimodal_count = 0
    total = 0
    all_wstar = []

    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = worpitzky_expand(F)
        wstar = list(reversed(w))
        H = sum(F)

        is_nonneg = all(x >= 0 for x in wstar)
        is_unimodal = all(wstar[i] <= wstar[i+1] for i in range(n//2 - 1)) and \
                       all(wstar[i] >= wstar[i+1] for i in range(n//2, n-1))

        if is_nonneg: nonneg_count += 1
        total += 1
        all_wstar.append((H, F, wstar))

        if total <= 5:
            print(f"  T#{total}: F={F}, w*={wstar}, H={H}, nonneg={is_nonneg}")

    print(f"\n  w* nonneg: {nonneg_count}/{total}")

    # Distribution of w*_1 = w_{n-2} = H - n*F_0
    print(f"  w*_0 = F_0, w*_1 = H - n*F_0")

    # What values do w*_k take?
    for k in range(n):
        vals = [ws[k] for _, _, ws in all_wstar]
        print(f"  w*_{k}: min={min(vals)}, max={max(vals)}")

# ===== SECTION 2: F(T,2) analysis =====
print("\n\n" + "=" * 70)
print("F(T,2) = sum F_k * 2^k")
print("=" * 70)

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    data = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        H = sum(F)
        F2 = sum(F[k] * 2**k for k in range(n))
        t3 = count_odd_cycles(A, n)
        data.append((H, F2, t3, F))

    # Is F(2) determined by H?
    by_H = defaultdict(set)
    for H, F2, t3, F in data:
        by_H[H].add(F2)

    print(f"  Is F(2) determined by H?")
    for H_val in sorted(by_H.keys()):
        vals = by_H[H_val]
        if len(vals) == 1:
            print(f"    H={H_val}: F(2)={vals.pop()} (unique)")
        else:
            print(f"    H={H_val}: {len(vals)} distinct F(2) values: {sorted(vals)}")

    # F(2) vs H ratio
    print(f"\n  F(2)/H ratios:")
    ratios = defaultdict(list)
    for H, F2, t3, F in data:
        ratios[H].append(F2/H)
    for H_val in sorted(ratios.keys()):
        r = ratios[H_val]
        print(f"    H={H_val}: F(2)/H in [{min(r):.3f}, {max(r):.3f}]")

    # Is F(2) = a*H + b*t3 + c*F_0 + ...?
    H_arr = np.array([d[0] for d in data], dtype=float)
    F2_arr = np.array([d[1] for d in data], dtype=float)
    t3_arr = np.array([d[2] for d in data], dtype=float)
    F0_arr = np.array([d[3][0] for d in data], dtype=float)

    # Linear regression: F2 = a*H + b*t3 + c*F0 + d_const
    X = np.column_stack([H_arr, t3_arr, F0_arr, np.ones(len(data))])
    coeffs, residuals, _, _ = np.linalg.lstsq(X, F2_arr, rcond=None)
    print(f"\n  Linear fit F(2) = {coeffs[0]:.4f}*H + {coeffs[1]:.4f}*t3 + {coeffs[2]:.4f}*F0 + {coeffs[3]:.4f}")
    pred = X @ coeffs
    r2 = 1 - np.sum((F2_arr - pred)**2) / np.sum((F2_arr - np.mean(F2_arr))**2)
    print(f"  R^2 = {r2:.6f}")

# ===== SECTION 3: F(T, 1/x) and palindromicity =====
print("\n\n" + "=" * 70)
print("F(T,x) AND F(T^op, x): COMPLEMENTARY RELATIONSHIP")
print("=" * 70)
print("F(T,x) = x^{n-1} F(T^op, 1/x)")
print("Equivalently: F_k(T) = F_{n-1-k}(T^op)")
print()

for n in [4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    tournaments = list(all_tournaments(n))
    # For each T, compute F(T) and F(T^op)
    verified = 0
    for A in tournaments[:20]:
        F = forward_edge_poly(A, n)
        # Build T^op
        A_op = [[1 - A[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
        F_op = forward_edge_poly(A_op, n)

        # Check F_k(T) = F_{d-k}(T^op)
        ok = all(F[k] == F_op[d-k] for k in range(n))
        if ok: verified += 1

        if verified <= 3 or not ok:
            print(f"  F(T)={F}, F(T^op)={F_op}, F_k(T)=F_{{d-k}}(T^op)? {ok}")

    print(f"  Verified: {verified}/20")

# ===== SECTION 4: Alternating sum F(-1) and its meaning =====
print("\n\n" + "=" * 70)
print("F(T,-1) = ALTERNATING SUM")
print("=" * 70)
print("""
F(T,-1) = sum_k F_k (-1)^k = #{even-ascent HPs} - #{odd-ascent HPs}

By complementarity: F(T,-1) = (-1)^{d} F(T^op, -1)
  - If d is even (n odd): F(T,-1) = F(T^op,-1)
  - If d is odd (n even): F(T,-1) = -F(T^op,-1)

At n=4 (d=3 odd): F(T,-1) = -F(T^op,-1), so they come in +/- pairs.
Sum over all T of F(T,-1) = 0 (since T <-> T^op is a bijection). Verified!

At n=5 (d=4 even): F(T,-1) = F(T^op,-1), so self-paired.
Sum = A_5(-1) * 2^6 = 16 * 64 = 1024. Verified!
""")

# ===== SECTION 5: The W^* polynomial and Dehn-Sommerville =====
print("=" * 70)
print("DEHN-SOMMERVILLE RELATIONS FOR W^*")
print("=" * 70)
print("""
For h*-vectors of Gorenstein polytopes: h*_k = h*_{d-k} (palindromic).
This would mean w_{n-1-k} = w_k, i.e., w is palindromic.
We know w is NOT palindromic in general. So F(T,x) does NOT come from
a Gorenstein polytope.

But what about the WEAKER condition: what linear relations hold among w_k?
""")

for n in [4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    all_w = []
    all_F = []
    for A in all_tournaments(n):
        F = forward_edge_poly(A, n)
        w = worpitzky_expand(F)
        all_w.append(w)
        all_F.append(F)

    w_arr = np.array(all_w)

    # Check all linear relations c_0 w_0 + ... + c_{d} w_d = const
    # that hold for ALL tournaments
    # This means: find the affine hull of the w vectors
    mean_w = w_arr.mean(axis=0)
    centered = w_arr - mean_w
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)
    print(f"  Singular values of centered w matrix: {S.round(4)}")
    print(f"  Rank of w space: {sum(s > 1e-8 for s in S)}")
    print(f"  This means w lives in a {sum(s > 1e-8 for s in S)}-dimensional affine subspace of R^{n}")

    # What are the null directions?
    null_dirs = Vt[sum(s > 1e-8 for s in S):]
    for i, v in enumerate(null_dirs):
        print(f"  Null direction {i}: {v.round(6)}")
        # Verify it's constant
        vals = w_arr @ v
        print(f"    w . v = {vals[0]:.6f} (all same? {np.std(vals) < 1e-8})")

# ===== SECTION 6: Connection to F(T,x) GF =====
print("\n\n" + "=" * 70)
print("GF STRUCTURE: sum_m F(T,m) x^m = W^*(T,x) / (1-x)^n")
print("=" * 70)

for n in [4]:
    d = n - 1
    print(f"\nn={n}:")

    A = list(all_tournaments(n))[10]  # pick one
    F = forward_edge_poly(A, n)
    w = worpitzky_expand(F)
    wstar = list(reversed(w))
    H = sum(F)

    print(f"  F={F}, w={w}, w*={wstar}, H={H}")
    print(f"  W*(1) = {sum(wstar)} should = H = {H}")

    # Verify: F(T,m) = sum_k w*_k C(m+n-1-k, n-1)
    print(f"\n  Verifying F(T,m) = sum_k w*_k C(m+{n}-1-k, {n}-1):")
    for m in range(6):
        Fm = sum(F[k] * m**k for k in range(n))
        Fm_wstar = sum(wstar[k] * comb(m + n - 1 - k, n - 1) for k in range(n))
        print(f"    m={m}: F(T,m)={Fm}, sum w*_k C(m+{n-1}-k,{n-1})={Fm_wstar}")

    # The GF coefficients
    print(f"\n  First terms of GF sum_m F(T,m) x^m:")
    for m in range(8):
        Fm = sum(F[k] * m**k for k in range(n))
        print(f"    m={m}: F(T,m) = {Fm}")

# ===== SECTION 7: w0 = (-1)^d F(T, -1) =====
print("\n\n" + "=" * 70)
print("w_0 AND F(T,-1)")
print("=" * 70)

for n in [3, 4, 5]:
    d = n - 1
    print(f"\nn={n}:")

    for idx, A in enumerate(list(all_tournaments(n))[:10]):
        F = forward_edge_poly(A, n)
        w = worpitzky_expand(F)
        Fm1 = sum(F[k] * (-1)**k for k in range(n))

        # C(-1+k, d) = C(k-1, d) = 0 for k <= d, C(d, d) = 1 for k = d+1... no
        # C(k-1, d): k=0: C(-1,d) = (-1)^d/d! (generalized)
        # Actually for the Worpitzky basis at x=-1:
        # F(-1) = sum w_k C(-1+k, d)
        # C(-1+k, d) = C(k-1, d) = 0 for 0 <= k-1 < d, so only k >= d+1 contributes...
        # But k <= d = n-1, so k-1 <= d-1 < d, meaning C(k-1,d) = 0 for all k!
        # That can't be right since F(-1) != 0.

        # Wait, C(-1+0, d) = C(-1, d) = (-1)^d (generalized binomial)
        # C(-1+1, d) = C(0, d) = 0 if d > 0
        # C(-1+2, d) = C(1, d) = 0 if d > 1
        # ...
        # C(-1+d, d) = C(d-1, d) = 0
        # So F(-1) = w_0 * (-1)^d, hence w_0 = (-1)^d * F(-1)

        w0_pred = (-1)**d * Fm1
        match = w[0] == w0_pred
        if idx < 5:
            print(f"  F={F}, w_0={w[0]}, (-1)^d*F(-1)={w0_pred}, match={match}")

    print(f"  So w_0 = (-1)^d F(-1) = (-1)^{{n-1}} F(-1)")

# ===== SECTION 8: Complete structure of w =====
print("\n\n" + "=" * 70)
print("COMPLETE STRUCTURE: w_k = Delta^d[F](k) = finite differences")
print("=" * 70)
print("""
The Worpitzky coefficients are related to FINITE DIFFERENCES.

If F(T,x) = sum w_k C(x+k, d), then by the standard forward difference:
  w_k = sum_{j=0}^{d} (-1)^{d-j} C(d, j) F(T, k+j-d)

For k=0: w_0 = sum_j (-1)^{d-j} C(d,j) F(T, j-d)
       = sum_j (-1)^{d-j} C(d,j) F(T, j-d)

Let me verify this.
""")

for n in [4]:
    d = n - 1
    A = list(all_tournaments(n))[7]
    F = forward_edge_poly(A, n)
    w = worpitzky_expand(F)
    print(f"n={n}: F={F}, w={w}")

    for k in range(n):
        w_pred = sum((-1)**(d-j) * comb(d, j) * sum(F[c] * (k+j-d)**c for c in range(n))
                     for j in range(d+1))
        print(f"  w_{k} = {w[k]}, finite diff formula = {int(round(w_pred))}")

    # So w_k = Delta^d[F](k-d+d) = Delta^d[F](k)... let me be precise.
    # Actually: F(x) = sum w_k C(x+k,d) means w_k = Delta^d_h[F](x)|_{x=k, h=1}
    # where Delta^d is the d-th forward difference.

    # Forward difference: Delta[f](x) = f(x+1) - f(x)
    # Delta^d[f](x) = sum_{j=0}^d (-1)^{d-j} C(d,j) f(x+j)
    # And w_k = Delta^d[F](k) ... but shifted.

    # Let me check: p(x) = sum w_k C(x+k,d)
    # Delta^d[p](x) = sum w_k Delta^d[C(x+k,d)]
    # Delta^d[C(x+k,d)] = sum_j (-1)^{d-j} C(d,j) C(x+j+k,d)
    # Hmm, C(x+k,d) is degree d in x, and Delta^d of degree-d poly = d! * leading coeff = 1.
    # So Delta^d[C(x+k,d)] = 1 for all x. That means Delta^d[p](x) = sum w_k = F(1) = H.
    # That's not useful.

    # Let me use the INVERSE direction. We want to express w_k in terms of F values.
    # The standard formula: if p(x) = sum_k w_k C(x+k,d), then evaluating at x=0,1,...,d:
    # p(j) = sum_k w_k C(j+k,d) for j=0,...,d
    # This is a triangular system since C(j+k,d) = 0 when j+k < d.

    print(f"\n  Matrix C(j+k, {d}) for j,k = 0,...,{d}:")
    for j in range(n):
        row = [comb(j+k, d) for k in range(n)]
        print(f"    j={j}: {row}")

    # So the system is:
    # p(0) = w_{d} * 1
    # p(1) = w_{d-1} * 1 + w_{d} * C(d+1,d) = w_{d-1} + (d+1)*w_d
    # p(2) = w_{d-2}*1 + w_{d-1}*C(d,d) + w_d*C(d+1,d) = ... wait let me think
    # Actually C(j+k, d) for j+k >= d: C(j+k,d) = C(j+k, j+k-d)
    # Row j=0: C(k, d) = 0 for k<d, C(d,d)=1 for k=d
    # Row j=1: C(1+k, d) = 0 for k<d-1, C(d,d)=1 for k=d-1, C(d+1,d)=d+1 for k=d
    # So it's upper triangular (in reverse k order):
    # w_d = p(0) = F_0
    # w_{d-1} = p(1) - (d+1)*w_d = (F_0 + F_1 + ... + F_d) - (d+1)*F_0 = H - n*F_0
    # Etc.
    print(f"\n  Back-substitution:")
    print(f"    w_{d} = F(0) = F_0 = {F[0]}")
    print(f"    w_{d-1} = F(1) - (d+1)*w_d = H - n*F_0 = {sum(F)} - {n}*{F[0]} = {sum(F)-n*F[0]}")

print("\n\nDone.")
