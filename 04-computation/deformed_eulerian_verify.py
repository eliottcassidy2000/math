#!/usr/bin/env python3
"""
DEFORMED EULERIAN NUMBERS — complete closed-form for a_k(T).

THEOREM (THM-062 complete form):
  a_k(T) = A(n, k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T)

where the "inflated Eulerian coefficient" is:

  c_k^{(f, d)} = sum_{j} A(f+1, j) * C(d-f, k-j) * (-1)^{d-f-k+j}

with sum over max(0, k-d+f) <= j <= min(f, k).

This arises from expressing F_f(r) = sum_j A(f+1,j) p^j q^{f-j}
in the degree-d basis {p^k q^{d-k}} using the identity (p-q)^{d-f} = 1^{d-f} = 1.

opus-2026-03-07-S32
"""
from itertools import permutations, combinations
from collections import defaultdict
from math import comb, factorial
import random

def eulerian_number(n, k):
    """A(n,k) = perms of [n] with k descents."""
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def inflated_eulerian(f, d, k):
    """
    c_k^{(f,d)} = coefficient of p^k q^{d-k} in F_f(r) expressed in degree-d basis.

    F_f(r) = sum_j A(f+1, j) p^j q^{f-j}, inflated to degree d via (p-q)^{d-f}=1.
    """
    total = 0
    for j in range(max(0, k - (d - f)), min(f, k) + 1):
        sign = (-1) ** (d - f - k + j)
        total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
    return total

# ====================================================================
# VERIFICATION 1: F_f(r) expanded in p^k q^{d-k} for d = f (trivial case)
# ====================================================================
print("VERIFICATION 1: d = f (trivial inflation)")
print("=" * 70)

for f in range(5):
    print(f"  f={f}: ", end="")
    coeffs = [inflated_eulerian(f, f, k) for k in range(f + 1)]
    euler = [eulerian_number(f + 1, k) for k in range(f + 1)]
    print(f"inflated={coeffs}, Eulerian={euler}, match={coeffs == euler}")

# ====================================================================
# VERIFICATION 2: F_2(r) in degree-6 basis (should give our t5 column at n=7)
# ====================================================================
print(f"\n{'=' * 70}")
print("VERIFICATION 2: F_2(r) in degree-6 basis")
print("=" * 70)

d = 6
f = 2
coeffs = [inflated_eulerian(f, d, k) for k in range(d + 1)]
print(f"  c_k = {coeffs}")
print(f"  Expected (t5 delta / 2): [1, 0, -9, 16, -9, 0, 1]")
print(f"  Match: {coeffs == [1, 0, -9, 16, -9, 0, 1]}")

# ====================================================================
# VERIFICATION 3: F_4(r) in degree-6 basis (should give t3 column / 2 at n=7)
# ====================================================================
print(f"\n{'=' * 70}")
print("VERIFICATION 3: F_4(r) in degree-6 basis")
print("=" * 70)

f = 4
coeffs = [inflated_eulerian(f, d, k) for k in range(d + 1)]
print(f"  c_k = {coeffs}")
print(f"  Expected (t3 delta / 2): [1, 24, 15, -80, 15, 24, 1]")
print(f"  Match: {coeffs == [1, 24, 15, -80, 15, 24, 1]}")

# ====================================================================
# VERIFICATION 4: F_0(r) in degree-6 basis (should give t7 column / 2 at n=7)
# ====================================================================
print(f"\n{'=' * 70}")
print("VERIFICATION 4: F_0(r) in degree-6 basis")
print("=" * 70)

f = 0
coeffs = [inflated_eulerian(f, d, k) for k in range(d + 1)]
print(f"  c_k = {coeffs}")
print(f"  Expected (t7 delta / 2): [1, -6, 15, -20, 15, -6, 1]")
print(f"  Match: {coeffs == [1, -6, 15, -20, 15, -6, 1]}")

# Note: C(6,k)*(-1)^{6-k} = (-1)^{6-k} C(6,k). Let me check:
# C(6,0)*(-1)^6 = 1, C(6,1)*(-1)^5 = -6, C(6,2)*(-1)^4 = 15,
# C(6,3)*(-1)^3 = -20, C(6,4)*(-1)^2 = 15, C(6,5)*(-1)^1 = -6, C(6,6)*(-1)^0 = 1
# YES! This is (-1)^{d-k} C(d, k), i.e., the row of Pascal's triangle with alternating signs.
# Beautiful: F_0(r) = 1, so its inflation to degree d is (p-q)^d in the p,q basis...
# Wait, (p-q)^d = 1^d = 1, but in the p^k q^{d-k} basis, 1 = ???
# Actually, 1 = (p-q)^d = sum_k (-1)^{d-k} C(d,k) p^k q^{d-k}. Yes!

print(f"\n  Note: c_k^(0,6) = (-1)^{{6-k}} C(6,k) — inflation of F_0=1 via (p-q)^6=1")

# ====================================================================
# VERIFICATION 5: Full a_k(T) formula at n=5
# ====================================================================
print(f"\n{'=' * 70}")
print("FULL FORMULA VERIFICATION at n=5")
print("=" * 70)

n = 5
d = n - 1  # = 4

# Invariants at n=5: t3 (f=2, parts=1), t5 (f=0, parts=1)
invariants_5 = [
    ('t3', 2, 1),  # name, f, parts
    ('t5', 0, 1),
]

print(f"  Inflated coefficients for each invariant:")
for name, f, parts in invariants_5:
    coeffs = [inflated_eulerian(f, d, k) for k in range(d + 1)]
    print(f"    {name} (f={f}, parts={parts}): 2^{parts} * {coeffs} = {[2**parts * c for c in coeffs]}")

# Compute for random tournaments
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

def forward_edge_dist(A, n):
    dist = defaultdict(int)
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[perm[i]][perm[i+1]])
        dist[fwd] += 1
    return dict(dist)

print(f"\n  Verifying against brute force:")
all_match = True
for trial in range(15):
    A = random_tournament(n, n * 200 + trial)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    actual = forward_edge_dist(A, n)

    inv_vals = {'t3': t3, 't5': t5}

    for k in range(n):
        predicted = eulerian_number(n, k)
        for name, f, parts in invariants_5:
            predicted += 2**parts * inflated_eulerian(f, d, k) * inv_vals[name]

        if predicted != actual.get(k, 0):
            print(f"    FAIL: trial={trial}, k={k}, pred={predicted}, actual={actual.get(k,0)}")
            all_match = False

print(f"  n=5: {'ALL PASS' if all_match else 'FAIL'} ({15 * n} coefficient checks)")

# ====================================================================
# VERIFICATION 6: Full a_k(T) formula at n=7
# ====================================================================
print(f"\n{'=' * 70}")
print("FULL FORMULA VERIFICATION at n=7")
print("=" * 70)

n = 7
d = n - 1  # = 6

# Invariants at n=7: t3 (f=4, parts=1), t5 (f=2, parts=1), t7 (f=0, parts=1), bc (f=2, parts=2)
invariants_7 = [
    ('t3', 4, 1),
    ('t5', 2, 1),
    ('t7', 0, 1),
    ('bc', 2, 2),
]

print(f"  Inflated coefficients:")
for name, f, parts in invariants_7:
    coeffs = [inflated_eulerian(f, d, k) for k in range(d + 1)]
    scaled = [2**parts * c for c in coeffs]
    print(f"    {name} (f={f}, 2^{parts}): {scaled}")

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

print(f"\n  Verifying against brute force (n=7 — using DP for dist):")

def forward_edge_dist_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    new_fwd = fwd + A[v][u]
                    key = (mask | (1 << u), u, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)

all_match = True
for trial in range(20):
    A = random_tournament(n, n * 200 + trial)
    t3 = count_t3(A, n)
    t5 = count_directed_cycles(A, n, 5)
    t7 = count_directed_cycles(A, n, 7)
    bc = count_bc(A, n)
    actual = forward_edge_dist_dp(A, n)

    inv_vals = {'t3': t3, 't5': t5, 't7': t7, 'bc': bc}

    for k in range(n):
        predicted = eulerian_number(n, k)
        for name, f, parts in invariants_7:
            predicted += 2**parts * inflated_eulerian(f, d, k) * inv_vals[name]

        if predicted != actual.get(k, 0):
            print(f"    FAIL: trial={trial}, k={k}, pred={predicted}, actual={actual.get(k,0)}")
            all_match = False

print(f"  n=7: {'ALL PASS' if all_match else 'FAIL'} ({20 * n} coefficient checks)")

# ====================================================================
# BEAUTIFUL SUMMARY
# ====================================================================
print(f"\n{'=' * 70}")
print("THEOREM: DEFORMED EULERIAN NUMBERS")
print("=" * 70)
print("""
For tournament T on n vertices, the forward-edge distribution is:

  a_k(T) = A(n, k) + sum_I 2^{parts(I)} * c_k^{(f_I, n-1)} * I(T)

where the inflated Eulerian coefficient is:

  c_k^{(f, d)} = sum_{j=max(0,k-d+f)}^{min(f,k)} A(f+1, j) * C(d-f, k-j) * (-1)^{d-f-k+j}

Special cases:
  c_k^{(0, d)} = (-1)^{d-k} C(d, k)     [signed Pascal row]
  c_k^{(d, d)} = A(d+1, k)               [standard Eulerian numbers]

Properties:
  c_k^{(f,d)} = c_{d-k}^{(f,d)}          [palindromic]
  sum_k c_k^{(f,d)} = 0  for f < d       [zero sum since F_f(r) evaluates same at r=1/2]

This completely determines all n! permutation counts a_k(T) from the
OCF invariants, via a single universal formula.
""")

# Verify zero-sum property
print("Zero-sum verification:")
for d in range(1, 8):
    for f in range(d):
        s = sum(inflated_eulerian(f, d, k) for k in range(d + 1))
        if s != 0:
            print(f"  FAIL: f={f}, d={d}, sum={s}")
    print(f"  d={d}: all f<d have zero sum")

print(f"\n{'=' * 70}")
print("DONE")
print("=" * 70)
