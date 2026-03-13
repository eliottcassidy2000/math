#!/usr/bin/env python3
"""
constant_complement.py -- kind-pasteur-2026-03-13-S60

Why is C_odd(comp) constant for 3-cycle complements at p=11 Paley?
Uses inclusion-exclusion to decompose the answer.
"""

from itertools import combinations
from collections import defaultdict

def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A

def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total

def total_odd_cycles(A, verts):
    n = len(verts)
    total = 0
    for k in range(3, n+1, 2):
        for sub in combinations(verts, k):
            total += count_ham_cycles(A, list(sub))
    return total


# Test constancy across ALL orientations at p=11
p = 11
m = (p - 1) // 2
N = 1 << m

print(f"p={p}, m={m}, testing all {N} orientations")

# First: just test c3(8-vertex comp) constancy
print(f"\n=== c3(8-vertex complement) constancy per orientation ===")

constant_c3 = 0
for bits in range(N):
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A = build_adj(p, S)

    c3_comp_vals = set()
    seen_orbits = set()
    for subset in combinations(range(p), 3):
        canon = min(frozenset((v+t)%p for v in subset) for t in range(p))
        if canon in seen_orbits:
            continue
        seen_orbits.add(canon)

        comp = sorted(frozenset(range(p)) - frozenset(subset))
        c3_comp = sum(count_ham_cycles(A, list(s)) for s in combinations(comp, 3))
        c3_comp_vals.add(c3_comp)

    is_constant = len(c3_comp_vals) == 1
    if is_constant:
        constant_c3 += 1

    c3 = sum(count_ham_cycles(A, list(s)) for s in combinations(range(p), 3))
    if bits < 4 or not is_constant:
        print(f"  bits={bits}: S={S}, c3={c3}, c3(comp) values={sorted(c3_comp_vals)}, "
              f"{'CONST' if is_constant else 'VARIES'}")

print(f"\n  c3(comp) CONSTANT: {constant_c3}/{N}")

# Now check with inclusion-exclusion for Paley
print(f"\n\n=== INCLUSION-EXCLUSION for Paley p={p} ===")
S_qr = sorted(j for j in range(1, p) if pow(j, (p - 1) // 2, p) == 1)
A = build_adj(p, S_qr)
qr_set = set(S_qr)

c3 = sum(count_ham_cycles(A, list(s)) for s in combinations(range(p), 3))
print(f"  c3 = {c3}")

# 3-cycles per vertex (by vertex transitivity)
cyc_per_vertex = c3 * 3 // p
print(f"  3-cycles per vertex = {cyc_per_vertex}")

# 3-cycles through a QR pair vs NQR pair
pair_qr = sum(count_ham_cycles(A, [0, 1, k]) for k in range(2, p))
# For NQR: diff 2 is NQR at p=11
pair_nqr = sum(count_ham_cycles(A, [0, 2, k]) for k in range(p) if k not in (0, 2))
print(f"  3-cycles through QR pair: {pair_qr}")
print(f"  3-cycles through NQR pair: {pair_nqr}")

# Inc-exc: c3(comp) = c3 - 3*cyc_per_vertex + [sum of pair counts] - n(V)
# For 3-set V with q QR-pairs (undirected):
# c3(comp) = c3 - 3*cpv + q*pair_qr + (3-q)*pair_nqr - n(V)
# = 55 - 45 + q*pair_qr + (3-q)*pair_nqr - n(V)

print(f"\n  Formula: c3(comp) = {c3} - 3*{cyc_per_vertex} + q*{pair_qr} + (3-q)*{pair_nqr} - n(V)")
print(f"         = 10 + q*{pair_qr - pair_nqr} + 3*{pair_nqr} - n(V)")
delta = pair_qr - pair_nqr
base = 10 + 3 * pair_nqr
print(f"         = {base} + {delta}*q - n(V)")

# For Paley: relationship between q (QR pairs) and n(V):
# A 3-set {a,b,c} has 3 undirected pairs. Each pair's difference is QR or NQR.
# n(V) = number of directed 3-cycles on {a,b,c} = 0 or 2.
# For Paley, n(V) = 2 iff the edges don't form a "transitive triple".
# In a tournament, a 3-set has n=2 (cyclic) or n=0 (transitive).
# Cyclic <=> all 3 edges form a cycle.

# For Paley: the tournament is self-complementary. The number of QR differences
# in a 3-set is always... let's check.

q_dist = defaultdict(lambda: defaultdict(int))
for subset in combinations(range(p), 3):
    verts = sorted(subset)
    q = 0
    for i in range(3):
        for j in range(i+1, 3):
            d = min((verts[j] - verts[i]) % p, (verts[i] - verts[j]) % p)
            if d in qr_set:
                q += 1
    nv = count_ham_cycles(A, verts)
    q_dist[q][nv] += 1

print(f"\n  QR-pair distribution:")
for q in sorted(q_dist):
    for nv in sorted(q_dist[q]):
        count = q_dist[q][nv]
        pred_c3comp = base + delta * q - nv
        print(f"    q={q}, n(V)={nv}: {count} subsets, predicted c3(comp) = {pred_c3comp}")

# Verify predictions
print(f"\n  Verification against actual c3(comp):")
for q_val in sorted(q_dist):
    for nv_val in sorted(q_dist[q_val]):
        # Find one such subset
        for subset in combinations(range(p), 3):
            verts = sorted(subset)
            q = 0
            for i in range(3):
                for j in range(i+1, 3):
                    d = min((verts[j] - verts[i]) % p, (verts[i] - verts[j]) % p)
                    if d in qr_set:
                        q += 1
            nv = count_ham_cycles(A, verts)
            if q == q_val and nv == nv_val:
                comp = sorted(frozenset(range(p)) - frozenset(subset))
                c3_actual = sum(count_ham_cycles(A, list(s)) for s in combinations(comp, 3))
                pred = base + delta * q - nv
                print(f"    q={q}, n(V)={nv}: actual={c3_actual}, pred={pred}, "
                      f"match={c3_actual == pred}")
                break

# So c3(comp) is constant iff base + delta*q - n(V) is constant,
# i.e., delta*q = n(V) + const.
# If delta = 0 (pair_qr = pair_nqr), then n(V) must be constant.
# Otherwise, q and n(V) must satisfy a linear relation.

if delta == 0:
    print(f"\n  delta = 0: pair_qr = pair_nqr => c3(comp) constant iff n(V) constant")
else:
    print(f"\n  delta = {delta}: constancy requires {delta}*q - n(V) = constant")
    # Check what the actual values are
    vals = set()
    for q_val in q_dist:
        for nv_val in q_dist[q_val]:
            vals.add(delta * q_val - nv_val)
    print(f"  {delta}*q - n(V) takes values: {sorted(vals)}")

print("\nDONE.")
