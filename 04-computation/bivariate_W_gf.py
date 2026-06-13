#!/usr/bin/env python3
"""
Bivariate generating function for tournament W(r).

We know: W(r) = sum_k a_k(T) * p^k * q^{n-1-k} where p=r+1/2, q=r-1/2.

This means: W(r) = E_T(p/q) * q^{n-1} where E_T(x) = sum_k a_k x^k.

Question: what is the bivariate GF?
G_T(x, y) = sum_k a_k(T) x^k y^{n-1-k} = W at (r -> (x-y)/2 + (x+y)/2) = ...

Actually more natural: define the *tournament Eulerian generating function*:
  Phi_T(x, y) = sum_{P in S_n} x^{fwd(P)} y^{bwd(P)}

where fwd(P) = #forward edges, bwd(P) = #backward edges, fwd + bwd = n-1.

So Phi_T(x, y) = sum_k a_k x^k y^{n-1-k}.

For the transitive tournament: Phi = classical Eulerian polynomial in (x, y).

Key identity: W(r) = Phi_T(r+1/2, r-1/2).

Now, from the OCF decomposition:
  Phi_T(x, y) = A_n(x, y) + sum_I 2^{parts} * c_I(x, y) * I(T)

where A_n(x, y) = sum_k A(n,k) x^k y^{n-1-k} (classical Eulerian)
and c_I(x, y) = inflated F_{f_I} in the (x, y) basis.

The inflation formula: F_f(x, y) = sum_j A(f+1, j) x^j y^{f-j}
inflates to degree d=n-1 via multiplication by (x-y)^{d-f} = ?

Wait, we used (p-q)^{d-f} = 1 because p-q = (r+1/2)-(r-1/2) = 1.
But for general (x, y): x - y is NOT 1! We need to be more careful.

The correct statement: W(r) = Phi_T(r+1/2, r-1/2). Setting x=r+1/2, y=r-1/2:
  x - y = 1 (always!)

So Phi_T restricted to the line x - y = 1 gives W(r).
On this line, the inflation formula works because (x-y)^{d-f} = 1.

But Phi_T is defined on all (x, y). Off the line x-y=1, the inflation
formula needs the FULL (x-y)^{d-f} factor:

  Phi_T(x, y) = A_n(x, y) + sum_I 2^{parts} * F_{f_I}(x, y) * (x-y)^{n-1-f_I} * I(T)

where F_f(x, y) = sum_j A(f+1, j) x^j y^{f-j}.

This is the COMPLETE bivariate formula!

Let me verify this.

opus-2026-03-07-S32
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb
import random

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def random_tournament(n, seed=42):
    rng = random.Random(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_t3(A, n):
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    return t3

def count_directed_cycles(A, n, cl):
    if n < cl: return 0
    total = 0
    for verts in combinations(range(n), cl):
        sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
        dp = [[0]*cl for _ in range(1 << cl)]
        dp[1][0] = 1
        for m in range(1, 1 << cl):
            for v in range(cl):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(cl):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << cl) - 1
        total += sum(dp[full][v] for v in range(1, cl) if sub[v][0])
    return total

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def phi_actual(A, n, x, y):
    """Compute Phi_T(x, y) = sum_P x^{fwd(P)} y^{bwd(P)}."""
    total = 0.0
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        bwd = (n-1) - fwd
        total += x**fwd * y**bwd
    return total

def phi_formula(n, invariants, x, y):
    """
    Phi_T(x, y) = A_n(x,y) + sum_I 2^parts * F_{f_I}(x,y) * (x-y)^{n-1-f_I} * I(T)

    where A_n(x,y) = sum_k A(n,k) x^k y^{n-1-k}
    and F_f(x,y) = sum_j A(f+1,j) x^j y^{f-j}
    """
    d = n - 1

    # Baseline
    result = sum(eulerian_number(n, k) * x**k * y**(d-k) for k in range(n))

    # Invariant corrections
    for name, val, f, parts in invariants:
        if val == 0:
            continue
        # F_f(x, y)
        F_f = sum(eulerian_number(f+1, j) * x**j * y**(f-j) for j in range(f+1))
        # Inflation factor
        inflation = (x - y) ** (d - f)
        result += 2**parts * F_f * inflation * val

    return result

# ====================================================================
# Verify at n=5
# ====================================================================
print("BIVARIATE VERIFICATION at n=5")
print("=" * 70)

n = 5
all_match = True
for trial in range(15):
    A = random_tournament(n, n * 888 + trial)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)

    inv = [('t3', t3, 2, 1), ('t5', t5, 0, 1)]

    # Test at several (x, y) points
    for x, y in [(2.0, 1.0), (1.5, 0.5), (3.0, 2.0), (1.0, -1.0), (0.7, 0.3)]:
        actual = phi_actual(A, n, x, y)
        predicted = phi_formula(n, inv, x, y)
        if abs(actual - predicted) > 0.001:
            print(f"  FAIL: trial={trial}, (x,y)=({x},{y}), actual={actual:.2f}, pred={predicted:.2f}")
            all_match = False

print(f"  n=5: {'ALL PASS' if all_match else 'FAIL'} (15 tournaments x 5 points = 75 checks)")

# Also verify W(r) = Phi(r+1/2, r-1/2)
print(f"\n  W(r) = Phi(r+1/2, r-1/2) check:")
for trial in range(5):
    A = random_tournament(n, n * 888 + trial)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)

    # W(r) at r = 0.7
    r = 0.7
    phi_val = phi_actual(A, n, r + 0.5, r - 0.5)

    # Also compute W(r) directly via product formula
    W_val = 0.0
    for perm in permutations(range(n)):
        prod = 1.0
        for i in range(n-1):
            s = A[perm[i]][perm[i+1]] - 0.5
            prod *= (r + s)
        W_val += prod

    print(f"    trial {trial}: Phi({r+0.5:.1f}, {r-0.5:.1f}) = {phi_val:.4f}, W({r}) = {W_val:.4f}, match={abs(phi_val-W_val)<0.001}")

# ====================================================================
# Verify at n=7
# ====================================================================
print(f"\n{'=' * 70}")
print("BIVARIATE VERIFICATION at n=7")
print("=" * 70)

n = 7
all_match = True
for trial in range(10):
    A = random_tournament(n, n * 999 + trial)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)

    inv = [('t3', t3, 4, 1), ('t5', t5, 2, 1), ('t7', t7, 0, 1), ('bc', bc, 2, 2)]

    for x, y in [(2.0, 1.0), (1.5, 0.5), (1.0, -1.0)]:
        actual = phi_actual(A, n, x, y)
        predicted = phi_formula(n, inv, x, y)
        if abs(actual - predicted) > 0.1:  # larger tolerance for n=7 (bigger numbers)
            print(f"  FAIL: trial={trial}, (x,y)=({x},{y}), actual={actual:.2f}, pred={predicted:.2f}, diff={abs(actual-predicted):.4f}")
            all_match = False

print(f"  n=7: {'ALL PASS' if all_match else 'FAIL'} (10 tournaments x 3 points = 30 checks)")

# ====================================================================
# Key observation
# ====================================================================
print(f"\n{'=' * 70}")
print("THEOREM: BIVARIATE DEFORMED EULERIAN FORMULA")
print("=" * 70)
print("""
Phi_T(x, y) = sum_P x^{fwd(P)} y^{bwd(P)}
            = A_n(x,y) + sum_I 2^{parts(I)} * A_{f_I+1}(x,y) * (x-y)^{n-1-f_I} * I(T)

where:
  A_n(x,y) = sum_k A(n,k) x^k y^{n-1-k}   [bivariate Eulerian polynomial]
  A_{f+1}(x,y) = sum_j A(f+1,j) x^j y^{f-j}

Special evaluations:
  Phi_T(1, 0) = n!                           [total permutations; only k=n-1 survives]
  Phi_T(1, 1) = n!                           [sum of all a_k = n!]
  Phi_T(x, 0) = x^{n-1} * H(T)              [only Hamiltonian paths contribute]
  Phi_T(0, x) = x^{n-1} * H(T)              [palindromy: a_0 = a_{n-1} = H(T)]

The formula on the line x - y = 1 reduces to W(r) (with x = r+1/2, y = r-1/2).
Off this line, the (x-y)^{n-1-f_I} factor is the key new ingredient.
""")
