#!/usr/bin/env python3
"""
Coefficient matrix analysis at n=9.

At n=9, d=8, the OCF invariants are the coefficients alpha_k of the
independence polynomial I(Omega(T), x) = sum alpha_k x^k.

The conflict graph Omega(T) has one vertex per directed odd cycle.
Two cycles are adjacent iff they share a vertex.

alpha_0 = 1
alpha_1 = number of directed odd cycles = t3 + t5 + t7 + t9
alpha_2 = number of independent pairs of directed odd cycles
alpha_3 = number of independent triples of directed odd cycles

At n=9, independent pairs can be:
  (3,3): two VD 3-cycles [bc33]
  (3,5): VD 3-cycle + 5-cycle [bc35]
  (3,7): needs 10 vertices -- IMPOSSIBLE
  (5,5): needs 10 vertices -- IMPOSSIBLE
  Others need even more vertices -- IMPOSSIBLE

So alpha_2 = bc33 + bc35.

Independent triples can only be three VD 3-cycles (9 vertices total): a3.
No other triple fits in 9 vertices.

The formula: d-f = sum over cycles of (L_i - 1).
  t3:   f=6, t5: f=4, t7: f=2, t9: f=0
  bc33: f=4, bc35: f=2, a3: f=2

IMPORTANT: cycle counts here count DIRECTED cycles, and VD pair/triple counts
must account for multiplicity (a 5-vertex set may have multiple directed 5-cycles).
"""

from itertools import combinations
from collections import defaultdict
from math import comb, factorial
from fractions import Fraction
import random
import numpy as np


def eulerian_number(n, k):
    if n == 0:
        return 1 if k == 0 else 0
    if k < 0 or k >= n:
        return 0
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+2))


def get_ck(f, d):
    """Inflated Eulerian coefficient c_k^{(f,d)} for k = 0, ..., d."""
    result = []
    for k in range(d + 1):
        total = 0
        for j in range(f + 1):
            if k - j < 0 or k - j > d - f:
                continue
            sign = (-1) ** (d - f - k + j)
            total += eulerian_number(f + 1, j) * comb(d - f, k - j) * sign
        result.append(total)
    return result


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


def count_directed_hamcycles(A, verts):
    """Count directed Hamiltonian cycles on vertex subset verts.
    Returns the number of distinct directed cycles (each cyclic ordering counted once).
    We fix the first vertex and count Hamiltonian paths from it that return to it."""
    cl = len(verts)
    if cl < 3:
        return 0
    sub = [[A[verts[i]][verts[j]] for j in range(cl)] for i in range(cl)]
    dp = [[0]*cl for _ in range(1 << cl)]
    dp[1][0] = 1  # start from vertex 0 in subset
    for m in range(1, 1 << cl):
        for v in range(cl):
            if not (m & (1 << v)) or dp[m][v] == 0:
                continue
            for u in range(cl):
                if m & (1 << u):
                    continue
                if sub[v][u]:
                    dp[m | (1 << u)][u] += dp[m][v]
    full = (1 << cl) - 1
    return sum(dp[full][v] for v in range(1, cl) if sub[v][0])


def count_directed_cycles_total(A, n, cl):
    """Count total number of directed cl-cycles in tournament A on n vertices."""
    if n < cl:
        return 0
    total = 0
    for verts in combinations(range(n), cl):
        total += count_directed_hamcycles(A, list(verts))
    return total


def get_cycle_data(A, n, cl):
    """For each cl-element subset, return (frozenset(verts), num_directed_cycles).
    Only include subsets with at least one directed cycle."""
    data = []
    for verts in combinations(range(n), cl):
        count = count_directed_hamcycles(A, list(verts))
        if count > 0:
            data.append((frozenset(verts), count))
    return data


def forward_edge_dist_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v, 0)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            for fwd in range(n):
                c = dp.get((mask, v, fwd), 0)
                if c == 0:
                    continue
                for u_node in range(n):
                    if mask & (1 << u_node):
                        continue
                    new_fwd = fwd + A[v][u_node]
                    key = (mask | (1 << u_node), u_node, new_fwd)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    dist = defaultdict(int)
    for v in range(n):
        for fwd in range(n):
            dist[fwd] += dp.get((full, v, fwd), 0)
    return dict(dist)


def count_all_invariants_n9(A):
    """Count all OCF invariants for n=9, properly handling multiplicity."""
    n = 9

    # Get cycle data: for each vertex set, how many directed cycles it supports
    cyc3_data = get_cycle_data(A, n, 3)  # each has count 1 (tournament on 3 vertices has 0 or 1 directed 3-cycle)
    cyc5_data = get_cycle_data(A, n, 5)
    cyc7_data = get_cycle_data(A, n, 7)
    cyc9_data = get_cycle_data(A, n, 9)

    # Total directed cycle counts (= vertices of Omega by cycle length)
    t3 = sum(c for _, c in cyc3_data)
    t5 = sum(c for _, c in cyc5_data)
    t7 = sum(c for _, c in cyc7_data)
    t9 = sum(c for _, c in cyc9_data)

    # VD pairs: bc33 - pairs of VD directed 3-cycles
    # Since each 3-vertex cyclic triple has exactly 1 directed 3-cycle,
    # bc33 = number of pairs of disjoint cyclic triples
    bc33 = 0
    for i in range(len(cyc3_data)):
        for j in range(i+1, len(cyc3_data)):
            if cyc3_data[i][0].isdisjoint(cyc3_data[j][0]):
                # Each contributes 1 directed cycle, so 1*1 = 1 pair in Omega
                bc33 += cyc3_data[i][1] * cyc3_data[j][1]

    # VD pairs: bc35 - pairs (directed 3-cycle, directed 5-cycle) vertex-disjoint
    bc35 = 0
    for s3, c3 in cyc3_data:
        for s5, c5 in cyc5_data:
            if s3.isdisjoint(s5):
                bc35 += c3 * c5

    # VD triples: a3 - triples of mutually VD directed 3-cycles
    a3_count = 0
    for i in range(len(cyc3_data)):
        for j in range(i+1, len(cyc3_data)):
            if not cyc3_data[i][0].isdisjoint(cyc3_data[j][0]):
                continue
            for k in range(j+1, len(cyc3_data)):
                if (cyc3_data[i][0].isdisjoint(cyc3_data[k][0]) and
                    cyc3_data[j][0].isdisjoint(cyc3_data[k][0])):
                    a3_count += cyc3_data[i][1] * cyc3_data[j][1] * cyc3_data[k][1]

    return {
        't3': t3, 't5': t5, 't7': t7, 't9': t9,
        'bc33': bc33, 'bc35': bc35, 'a3': a3_count
    }


# ========================================================================
# First: verify the formula at n=7 to make sure our counting is right
# ========================================================================

print("=" * 70)
print("VERIFICATION AT n=7 (sanity check)")
print("=" * 70)

n = 7
An7 = [eulerian_number(7, k) for k in range(7)]
print(f"A(7,k) = {An7}")

# n=7 invariants: t3, t5, t7, bc33
# f values: t3->f=4, t5->f=2, t7->f=0, bc33->f=2
inv_n7 = [
    ('t3',  1, 4),  # (name, parts, f)
    ('t5',  1, 2),
    ('t7',  1, 0),
    ('bc33', 2, 2),
]

for seed in range(5):
    A = random_tournament(7, 7000 + seed)
    dist = forward_edge_dist_dp(A, 7)
    a_actual = [dist.get(k, 0) for k in range(7)]

    cyc3 = get_cycle_data(A, 7, 3)
    cyc5 = get_cycle_data(A, 7, 5)
    cyc7 = get_cycle_data(A, 7, 7)

    t3 = sum(c for _, c in cyc3)
    t5 = sum(c for _, c in cyc5)
    t7 = sum(c for _, c in cyc7)
    bc33 = 0
    for i in range(len(cyc3)):
        for j in range(i+1, len(cyc3)):
            if cyc3[i][0].isdisjoint(cyc3[j][0]):
                bc33 += cyc3[i][1] * cyc3[j][1]

    inv_vals = {'t3': t3, 't5': t5, 't7': t7, 'bc33': bc33}

    # Predict a_k
    a_pred = []
    for k in range(7):
        val = An7[k]
        for name, parts, f in inv_n7:
            ck = get_ck(f, 6)
            val += (2**parts) * ck[k] * inv_vals[name]
        a_pred.append(val)

    match = (a_actual == a_pred)
    print(f"  seed={seed}: t3={t3}, t5={t5}, t7={t7}, bc33={bc33}")
    print(f"    actual:  {a_actual}")
    print(f"    predict: {a_pred}")
    print(f"    {'OK' if match else 'FAIL'}")

# ========================================================================
# Now n=9
# ========================================================================

print(f"\n{'=' * 70}")
print("n=9 COEFFICIENT MATRIX ANALYSIS")
print("=" * 70)

n = 9
d = 8
An9 = [eulerian_number(9, k) for k in range(9)]
print(f"A(9,k) = {An9}")
print(f"Sum = {sum(An9)} = 9! = {factorial(9)}: {'OK' if sum(An9) == factorial(9) else 'FAIL'}")

invariants_n9 = [
    ('t3',   1, 6),  # (name, parts, f)
    ('t5',   1, 4),
    ('t7',   1, 2),
    ('t9',   1, 0),
    ('bc33', 2, 4),
    ('bc35', 2, 2),
    ('a3',   3, 2),
]

inv_names = [name for name, _, _ in invariants_n9]
n_inv = len(inv_names)

print(f"\nInflated Eulerian coefficients:")
for name, parts, f in invariants_n9:
    ck = get_ck(f, 8)
    print(f"  {name:6s}: parts={parts}, f={f}, 2^p*c_k = {[2**parts * c for c in ck]}")

# Build coefficient matrix
M_full = np.zeros((9, n_inv), dtype=int)
for j, (name, parts, f) in enumerate(invariants_n9):
    ck = get_ck(f, 8)
    for k in range(9):
        M_full[k][j] = (2**parts) * ck[k]

print(f"\nCoefficient matrix (rows a_0..a_8, cols {inv_names}):")
for k in range(9):
    print(f"  a_{k}: {list(M_full[k])}")

# Palindromy check
print(f"\nPalindromy: {all(all(M_full[k][j] == M_full[8-k][j] for j in range(n_inv)) for k in range(4))}")

# Independent rows: a_0, a_1, a_2, a_3, a_4
M_indep = M_full[:5].astype(int)

# ========================================================================
# Rank and null space
# ========================================================================

print(f"\n{'=' * 70}")
print("RANK AND NULL SPACE")
print("=" * 70)

# Exact computation using fractions
M_frac = [[Fraction(int(M_indep[i][j])) for j in range(n_inv)] for i in range(5)]

def rref_fraction(mat):
    m = len(mat)
    n_cols = len(mat[0])
    A = [row[:] for row in mat]
    pivot_cols = []
    row_idx = 0
    for col in range(n_cols):
        pivot = None
        for r in range(row_idx, m):
            if A[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        pivot_cols.append(col)
        A[row_idx], A[pivot] = A[pivot], A[row_idx]
        scale = A[row_idx][col]
        A[row_idx] = [x / scale for x in A[row_idx]]
        for r in range(m):
            if r != row_idx and A[r][col] != 0:
                factor = A[r][col]
                A[r] = [A[r][j] - factor * A[row_idx][j] for j in range(n_cols)]
        row_idx += 1
    return A, pivot_cols

rref, pivot_cols = rref_fraction(M_frac)
free_cols = [j for j in range(n_inv) if j not in pivot_cols]
rank = len(pivot_cols)
null_dim = len(free_cols)

print(f"Matrix dimensions: {5} x {n_inv}")
print(f"Rank: {rank}")
print(f"Null space dimension: {null_dim}")
print(f"Pivot columns: {pivot_cols} = {[inv_names[j] for j in pivot_cols]}")
print(f"Free columns: {free_cols} = {[inv_names[j] for j in free_cols]}")

print(f"\nRREF:")
for i in range(rank):
    entries = [f"{rref[i][j]}" for j in range(n_inv)]
    print(f"  [{', '.join(entries)}]")

# Extract null vectors
from math import gcd

null_vectors_exact = []
for fc in free_cols:
    vec = [Fraction(0)] * n_inv
    vec[fc] = Fraction(1)
    for i, pc in enumerate(pivot_cols):
        vec[pc] = -rref[i][fc]
    null_vectors_exact.append(vec)

print(f"\n{'=' * 70}")
print("NULL VECTORS (integer form)")
print("=" * 70)

for idx, nv in enumerate(null_vectors_exact):
    # Convert to integers
    def lcm(a, b):
        return a * b // gcd(a, b)
    dl = 1
    for x in nv:
        if x != 0:
            dl = lcm(dl, abs(x.denominator))
    int_vec = [int(x * dl) for x in nv]
    g = 0
    for x in int_vec:
        g = gcd(g, abs(x))
    if g > 0:
        int_vec = [x // g for x in int_vec]

    print(f"\n  Null vector {idx+1}: {int_vec}")
    terms_pos = []
    terms_neg = []
    for j in range(n_inv):
        if int_vec[j] > 0:
            terms_pos.append(f"{int_vec[j]}*{inv_names[j]}")
        elif int_vec[j] < 0:
            terms_neg.append(f"{-int_vec[j]}*{inv_names[j]}")
    if terms_pos and terms_neg:
        print(f"    Meaning: {' + '.join(terms_pos)} = {' + '.join(terms_neg)}")
    else:
        terms = []
        for j in range(n_inv):
            if int_vec[j] != 0:
                terms.append(f"{int_vec[j]}*{inv_names[j]}")
        print(f"    Meaning: {' + '.join(terms)} = 0")

    # Verify
    for k in range(5):
        val = sum(M_frac[k][j] * nv[j] for j in range(n_inv))
        assert val == 0, f"Null vector check failed at row {k}: {val}"
    print(f"    Verification: M * v = 0  [PASSED]")

# ========================================================================
# Numerical verification at n=9
# ========================================================================

print(f"\n{'=' * 70}")
print("NUMERICAL VERIFICATION AT n=9")
print("=" * 70)

print("\nTesting OCF formula (may take a few minutes for n=9)...")

n_test = 8
all_ok = True
for seed in range(n_test):
    A = random_tournament(9, 9000 + seed)

    # Actual forward-edge distribution
    dist = forward_edge_dist_dp(A, 9)
    a_actual = [dist.get(k, 0) for k in range(9)]

    # Count invariants
    inv = count_all_invariants_n9(A)

    # Predict a_k
    a_pred = []
    for k in range(9):
        val = An9[k]
        for j, (name, parts, f) in enumerate(invariants_n9):
            ck = get_ck(f, 8)
            val += (2**parts) * ck[k] * inv[name]
        a_pred.append(val)

    match = (a_actual == a_pred)
    if not match:
        all_ok = False
    print(f"\n  seed={seed}: {'OK' if match else 'FAIL'}")
    print(f"    t3={inv['t3']}, t5={inv['t5']}, t7={inv['t7']}, t9={inv['t9']}")
    print(f"    bc33={inv['bc33']}, bc35={inv['bc35']}, a3={inv['a3']}")
    if not match:
        diff = [a_actual[k] - a_pred[k] for k in range(9)]
        print(f"    actual:  {a_actual}")
        print(f"    predict: {a_pred}")
        print(f"    diff:    {diff}")
    else:
        print(f"    a_k = {a_actual}")

print(f"\n{'=' * 70}")
print(f"Overall verification: {'ALL PASS' if all_ok else 'SOME FAILURES'}")
print(f"{'=' * 70}")

# ========================================================================
# SUMMARY
# ========================================================================

print(f"\n{'=' * 70}")
print("FINAL SUMMARY")
print("=" * 70)

print(f"""
At n=9, the deformed Eulerian numbers a_k(T) are determined by
{n_inv} cycle-structure invariants: {', '.join(inv_names)}.

Due to palindromy (a_k = a_{{8-k}}), there are 5 independent coefficients.

Coefficient matrix M (5 x {n_inv}):
  Rank = {rank}
  Null space dimension = {null_dim}

This means: from the 5 deformed Eulerian numbers, we can recover only
{rank} of the {n_inv} invariants. {null_dim} linear combinations of
cycle structures are invisible to the forward-edge distribution.
""")

print("Null vectors (invisible combinations):")
for idx, nv in enumerate(null_vectors_exact):
    def _lcm(a, b):
        return a * b // gcd(a, b)
    dl = 1
    for x in nv:
        if x != 0:
            dl = _lcm(dl, abs(x.denominator))
    iv = [int(x * dl) for x in nv]
    g = 0
    for x in iv:
        g = gcd(g, abs(x))
    if g > 0:
        iv = [x // g for x in iv]
    terms_pos, terms_neg = [], []
    for j in range(n_inv):
        if iv[j] > 0:
            terms_pos.append(f"{iv[j]}*{inv_names[j]}")
        elif iv[j] < 0:
            terms_neg.append(f"{-iv[j]}*{inv_names[j]}")
    if terms_pos and terms_neg:
        print(f"  {idx+1}. {' + '.join(terms_pos)} = {' + '.join(terms_neg)}")

print(f"""
Compare with n=7:
  4 invariants (t3, t5, t7, bc33), 4 independent a_k values
  Rank = 3, null dim = 1
  Null vector: 1*bc33 = 2*t5 (bc33 and 2*t5 are indistinguishable)

At n=9, the same phenomenon persists and expands: the "multi-cycle"
invariants (bc33, bc35, a3) have their inflated Eulerian coefficients
proportional to those of single cycles (t5, t7), so they cannot be
separately resolved from the forward-edge distribution alone.
""")
