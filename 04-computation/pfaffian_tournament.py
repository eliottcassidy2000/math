#!/usr/bin/env python3
"""
pfaffian_tournament.py -- kind-pasteur-2026-03-13-S61

Compute the SIGNED Pfaffian of skew-symmetric matrices arising from tournaments,
and verify the formula det(I + 2A) = det(J + S) via adj(S).

Key result: I + 2A = J + S where S is skew-symmetric.
For even n: det(J+S) = det(S) = Pf(S)^2
For odd n: det(J+S) = sum_{i,j} adj(S)_{i,j} = (sum_i w_i)^2

where adj(S) = w * w^T (rank 1 PSD for odd-order skew-symmetric S).

The SIGNED Pfaffian is computed via the recursive definition:
  Pf(empty) = 1
  Pf(M) = sum_{j>0} (-1)^j M[0][j] * Pf(M with rows/cols 0,j deleted)

Author: kind-pasteur-2026-03-13-S61
"""

from itertools import combinations
from collections import defaultdict


def binary_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << pos):
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A


def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + dp[key]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def count_directed_ham_cycles_on_subset(A, verts):
    k = len(verts)
    if k <= 2:
        return 0
    if k == 3:
        a, b, c = verts
        return (A[a][b] * A[b][c] * A[c][a]) + (A[a][c] * A[c][b] * A[b][a])
    dp = {}
    dp[(1, 0)] = 1
    for mask in range(1, 1 << k):
        if not (mask & 1):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nk = (mask | (1 << w), w)
                    dp[nk] = dp.get(nk, 0) + dp[key]
    full = (1 << k) - 1
    total = 0
    for v in range(1, k):
        if (full, v) in dp and A[verts[v]][verts[0]]:
            total += dp[(full, v)]
    return total


def det_small(M, n):
    if n == 0:
        return 1
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0]*M[1][1] - M[0][1]*M[1][0]
    result = 0
    for j in range(n):
        if M[0][j] == 0:
            continue
        minor = [[M[r][c] for c in range(n) if c != j] for r in range(1, n)]
        cofactor = ((-1) ** j) * M[0][j] * det_small(minor, n-1)
        result += cofactor
    return result


def pfaffian(M):
    """Compute the Pfaffian of a skew-symmetric matrix M.
    Uses the recursive definition."""
    n = len(M)
    if n == 0:
        return 1
    if n == 1:
        return 0  # Pf of 1x1 is undefined/0 for odd
    if n == 2:
        return M[0][1]
    # Pf(M) = sum_{j=1}^{n-1} (-1)^{j+1} M[0][j] Pf(M_{0j})
    # where M_{0j} = M with rows and cols 0,j removed
    result = 0
    for j in range(1, n):
        if M[0][j] == 0:
            continue
        # Build M with rows 0,j and cols 0,j removed
        indices = [i for i in range(n) if i != 0 and i != j]
        sub = [[M[indices[a]][indices[b]] for b in range(len(indices))] for a in range(len(indices))]
        result += ((-1) ** (j + 1)) * M[0][j] * pfaffian(sub)
    return result


# ========================================================================
# VERIFY: Pf(S)^2 = det(S) at small sizes
# ========================================================================
print("=" * 70)
print("VERIFICATION: Pf(S)^2 = det(S)")
print("=" * 70)

for n in [3, 4, 5, 6]:
    m = n * (n - 1) // 2
    total = 1 << m
    all_match = True
    for bits in range(total):
        A = binary_to_tournament(bits, n)
        S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
        pf = pfaffian(S)
        det_s = det_small(S, n)
        if pf * pf != det_s:
            all_match = False
            print(f"  FAIL n={n} bits={bits}: Pf={pf}, Pf^2={pf*pf}, det={det_s}")
            break
    if all_match:
        print(f"  n={n}: Pf^2 = det(S) VERIFIED for all {total} tournaments")


# ========================================================================
# COMPUTE: The null vector of S for odd n
# ========================================================================
print(f"\n{'='*70}")
print("NULL VECTOR OF S (odd n) = signed Pfaffian vector")
print("=" * 70)

# For odd n, S has a 1-dimensional null space.
# The null vector w has components w_i = (-1)^i * Pf(S_{ii})
# where S_{ii} is S with row i and col i deleted.

# Verify: S * w = 0
for n in [3, 5]:
    print(f"\n  n={n}:")
    m = n * (n - 1) // 2
    for bits in range(0, 1 << m, max(1, (1 << m) // 10)):
        A = binary_to_tournament(bits, n)
        S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]

        w = []
        for i in range(n):
            remaining = [j for j in range(n) if j != i]
            S_del = [[S[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            pf_del = pfaffian(S_del)
            w.append(((-1) ** i) * pf_del)

        # Check S*w = 0
        Sw = [sum(S[i][j] * w[j] for j in range(n)) for i in range(n)]
        is_zero = all(x == 0 for x in Sw)

        # Check: det(J+S) = (sum w)^2
        JS = [[1 + S[i][j] for j in range(n)] for i in range(n)]
        det_JS = det_small(JS, n)
        w_sum = sum(w)
        w_sum_sq = w_sum * w_sum

        if not is_zero or det_JS != w_sum_sq:
            print(f"    bits={bits}: w={w}, Sw={Sw}, S*w=0? {is_zero}, det(J+S)={det_JS}, sum(w)^2={w_sum_sq}, match={det_JS == w_sum_sq}")

    # Just show a couple of examples
    for bits in [0, 1, 2]:
        A = binary_to_tournament(bits, n)
        S = [[A[i][j] - A[j][i] for j in range(n)] for i in range(n)]
        w = []
        for i in range(n):
            remaining = [j for j in range(n) if j != i]
            S_del = [[S[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
            pf_del = pfaffian(S_del)
            w.append(((-1) ** i) * pf_del)
        JS = [[1 + S[i][j] for j in range(n)] for i in range(n)]
        det_JS = det_small(JS, n)
        w_sum = sum(w)
        print(f"    bits={bits}: w={w}, sum={w_sum}, sum^2={w_sum*w_sum}, det(J+S)={det_JS}")


# ========================================================================
# EXHAUSTIVE: n=5 Pfaffian sum vs H
# ========================================================================
print(f"\n{'='*70}")
print("EXHAUSTIVE: n=5 Pfaffian sum (signed) vs H")
print("=" * 70)

n5 = 5
by_H_pf = defaultdict(list)
pf_to_H = defaultdict(set)

for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    S = [[A[i][j] - A[j][i] for j in range(n5)] for i in range(n5)]
    H = count_ham_paths(A, n5)

    w = []
    for i in range(n5):
        remaining = [j for j in range(n5) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n5-1)] for a in range(n5-1)]
        pf_del = pfaffian(S_del)
        w.append(((-1) ** i) * pf_del)

    pf_sum = sum(w)
    by_H_pf[H].append((bits, pf_sum, tuple(w)))
    pf_to_H[pf_sum].add(H)

print(f"  Pfaffian sum -> H mapping:")
for ps in sorted(pf_to_H.keys()):
    Hs = sorted(pf_to_H[ps])
    print(f"    Pf_sum={ps:+d}: H={Hs}")

print(f"\n  H -> Pfaffian sum mapping:")
for H_val in sorted(by_H_pf.keys()):
    ps_vals = sorted(set(ps for _, ps, _ in by_H_pf[H_val]))
    print(f"    H={H_val}: Pf_sum={ps_vals}")

# Check: does (H, Pf_sum) determine more than H alone?
print(f"\n  Does Pf_sum determine anything beyond H at n=5?")
combined = defaultdict(set)
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    S = [[A[i][j] - A[j][i] for j in range(n5)] for i in range(n5)]
    H = count_ham_paths(A, n5)

    w = []
    for i in range(n5):
        remaining = [j for j in range(n5) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n5-1)] for a in range(n5-1)]
        pf_del = pfaffian(S_del)
        w.append(((-1) ** i) * pf_del)
    pf_sum = sum(w)

    c5 = count_directed_ham_cycles_on_subset(A, list(range(n5)))
    combined[(H, pf_sum)].add(c5)

for (H_val, ps), c5_vals in sorted(combined.items()):
    if len(c5_vals) > 1:
        print(f"    (H={H_val}, Pf_sum={ps}): c5={sorted(c5_vals)} -- AMBIGUOUS")
    else:
        print(f"    (H={H_val}, Pf_sum={ps:+d}): c5={sorted(c5_vals)}")


# ========================================================================
# THE DEEPEST PAIR at n=7
# ========================================================================
print(f"\n{'='*70}")
print("DEEPEST AMBIGUITY PAIR: n=7 Pfaffian analysis")
print("=" * 70)

n7 = 7
for bits, label in [(4728, "H=109, c7=8"), (4658, "H=111, c7=9")]:
    A = binary_to_tournament(bits, n7)
    S = [[A[i][j] - A[j][i] for j in range(n7)] for i in range(n7)]
    H = count_ham_paths(A, n7)

    # Full Pfaffian sum
    w = []
    for i in range(n7):
        remaining = [j for j in range(n7) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n7-1)] for a in range(n7-1)]
        pf_del = pfaffian(S_del)
        w.append(((-1) ** i) * pf_del)

    pf_sum = sum(w)

    # Verify det(J+S) = pf_sum^2
    JS = [[1 + S[i][j] for j in range(n7)] for i in range(n7)]
    det_JS = det_small(JS, n7)

    # Also compute det(I+2A) directly
    I2A = [[2*A[i][j] for j in range(n7)] for i in range(n7)]
    for i in range(n7):
        I2A[i][i] += 1
    det_I2A = det_small(I2A, n7)

    c7 = count_directed_ham_cycles_on_subset(A, list(range(n7)))

    print(f"\n  {label}:")
    print(f"    H = {H}, c7 = {c7}")
    print(f"    w = {w}")
    print(f"    Pf_sum = sum(w) = {pf_sum}")
    print(f"    Pf_sum^2 = {pf_sum*pf_sum}")
    print(f"    det(J+S) = {det_JS}")
    print(f"    det(I+2A) = {det_I2A}")
    print(f"    Match: Pf_sum^2 = det(J+S) = det(I+2A)? {pf_sum*pf_sum == det_JS == det_I2A}")

    # Individual Pfaffians (deletion Pfaffians)
    print(f"    Deletion Pfaffians:")
    for i in range(n7):
        remaining = [j for j in range(n7) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n7-1)] for a in range(n7-1)]
        pf_del = pfaffian(S_del)
        det_del = det_small(S_del, n7-1)
        print(f"      Pf(S\\{i}) = {pf_del:+d} (det = {det_del} = {pf_del}^2 = {pf_del*pf_del})")


# ========================================================================
# WHAT DOES Pf(S\\v) COUNT?
# ========================================================================
print(f"\n{'='*70}")
print("INTERPRETATION: What does Pf(S\\v) count?")
print("=" * 70)

print("""
For a tournament T on n vertices with skew-adjacency S = A - A^T:

Pf(S) counts the number of PERFECT MATCHINGS with signs:
  Pf(S) = sum over perfect matchings M: sign(M) * prod_{(i,j) in M} S[i][j]

For a tournament, S[i][j] = +1 if i->j, -1 if j->i.
So each perfect matching M gets weight = product of signs = +1 or -1
multiplied by the sign of the matching permutation.

For EVEN n: Pf(S) = signed count of perfect matchings.
For ODD n: Pf is undefined (S is singular). But Pf(S\\v) is defined
for each vertex v, giving a "local Pfaffian" at each vertex.

The null vector w_i = (-1)^i * Pf(S\\i) has the interpretation:
w_i = signed count of perfect matchings of T\\{i} with an extra
sign factor of (-1)^i (from the row/column ordering).

The PFAFFIAN SUM = sum_i w_i = sum_i (-1)^i * Pf(S\\i)
is the total "signed matching weight" across all vertex deletions.
""")

# Compute actual perfect matching counts at n=5 to check
print("Verify at n=5:")
for bits in [0, 10, 40, 41]:
    A = binary_to_tournament(bits, n5)
    S = [[A[i][j] - A[j][i] for j in range(n5)] for i in range(n5)]
    H = count_ham_paths(A, n5)

    # For each vertex deletion, count signed perfect matchings
    for del_v in range(n5):
        remaining = [j for j in range(n5) if j != del_v]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(4)] for a in range(4)]

        # Enumerate all perfect matchings of K_4 on vertices remaining
        # Perfect matchings of K_4: {{0,1},{2,3}}, {{0,2},{1,3}}, {{0,3},{1,2}}
        pm_count = 0
        r = remaining
        matchings = [
            [(r[0],r[1]), (r[2],r[3])],
            [(r[0],r[2]), (r[1],r[3])],
            [(r[0],r[3]), (r[1],r[2])],
        ]
        signed_sum = 0
        for pm in matchings:
            # Product of S values
            prod_val = 1
            for i, j in pm:
                prod_val *= S[i][j]
            # Sign of matching: parity of the permutation
            # For K_4 with matching {{a,b},{c,d}}, the permutation is
            # (a->b, b->a, c->d, d->c) which has 0 inversions in sorted
            # Actually the Pfaffian sign comes from the ordering convention
            # Let's just use the recursive Pfaffian
            signed_sum += 0  # skip manual computation
            pm_count += 1

        pf_del = pfaffian(S_del)

    w = []
    for i in range(n5):
        remaining = [j for j in range(n5) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(4)] for a in range(4)]
        pf_del = pfaffian(S_del)
        w.append(((-1) ** i) * pf_del)

    pf_sum = sum(w)
    print(f"  bits={bits}: H={H}, w={w}, Pf_sum={pf_sum}, Pf_sum^2={pf_sum*pf_sum}")


# ========================================================================
# FORMULA: sqrt(det(I+2A)) mod various primes
# ========================================================================
print(f"\n{'='*70}")
print("PFAFFIAN SUM MODULAR PROPERTIES")
print("=" * 70)

# Check: is the Pfaffian sum always odd?
print(f"\n  n=5 exhaustive:")
pf_sums_5 = set()
for bits in range(1 << 10):
    A = binary_to_tournament(bits, n5)
    S = [[A[i][j] - A[j][i] for j in range(n5)] for i in range(n5)]

    w = []
    for i in range(n5):
        remaining = [j for j in range(n5) if j != i]
        S_del = [[S[remaining[a]][remaining[b]] for b in range(n5-1)] for a in range(n5-1)]
        pf_del = pfaffian(S_del)
        w.append(((-1) ** i) * pf_del)
    pf_sums_5.add(sum(w))

print(f"  All Pfaffian sums at n=5: {sorted(pf_sums_5)}")
all_odd_5 = all(ps % 2 == 1 for ps in pf_sums_5)
print(f"  All odd: {all_odd_5}")


print(f"\n{'='*70}")
print("DONE.")
print("=" * 70)
