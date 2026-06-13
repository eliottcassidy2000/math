#!/usr/bin/env python3
"""
marginal_H_verification.py -- kind-pasteur-2026-03-13-S60

Verify the effective marginal H coefficients by computing
the ACTUAL average disjoint pair count per cycle by length.

For a k-cycle C on vertex set V:
- It contributes 1 to N
- It forms disjoint pairs with C_odd(T[comp(V)]) other cycles
- It forms disjoint triples in a more complex way

So the "marginal alpha_2" from one k-cycle on V is exactly
C_odd(T[comp(V)]).

The effective b_k should satisfy:
  b_k = 2 + 4 * avg_over_V[C_odd(comp(V))] / dc_k_relationship

Wait, that's not quite right. Let me think more carefully.

dH = 2*dN + 4*d(alpha_2) + 8*d(alpha_3)

When we change orientation, MANY cycles change simultaneously.
The effective b_k captures the NET effect including all interactions.

But we can compute the AVERAGE C_odd(comp) per cycle, which gives
the "per-cycle disjoint pair contribution."

At Paley p=11:
- k=3: D = C_odd(8-vert comp) = 172 for ALL 3-cycles (constant!)
  Each 3-cycle creates 172 disjoint pairs.
  Per-cycle: dH_per = 2 + 4*172 = 690 (but this counts EACH pair twice)

Actually, let me think about this differently.
alpha_2 = (1/2) * sum_C w_disj(C) where w_disj(C) counts disjoint neighbors.
No: alpha_2 = sum_{C1 < C2} [V(C1) disj V(C2)] * n(V1) * n(V2).

The key issue is that when orientation changes, MANY c_k values change
simultaneously. The effective b_k is an AVERAGE effect, not a per-cycle one.
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


p = 11
m = (p - 1) // 2
S_qr = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
A = build_adj(p, S_qr)

print(f"PALEY p={p}: Average C_odd(comp) per cycle by length")
print(f"{'='*60}")

# For each active vertex set, compute C_odd(comp)
active_by_k = defaultdict(list)  # k -> [(fs, nc, C_odd_comp)]

for k in range(3, p+1, 2):
    for subset in combinations(range(p), k):
        fs = frozenset(subset)
        nc = count_ham_cycles(A, list(subset))
        if nc > 0:
            # Compute C_odd(comp)
            comp = sorted(frozenset(range(p)) - fs)
            c_odd_comp = 0
            if len(comp) >= 3:
                for k2 in range(3, len(comp)+1, 2):
                    for sub2 in combinations(comp, k2):
                        c_odd_comp += count_ham_cycles(A, list(sub2))
            active_by_k[k].append((fs, nc, c_odd_comp))

for k in sorted(active_by_k):
    entries = active_by_k[k]
    n_sets = len(entries)
    total_nc = sum(nc for _, nc, _ in entries)
    total_cod = sum(nc * cod for _, nc, cod in entries)  # weighted by cycle count
    avg_cod = total_cod / total_nc if total_nc > 0 else 0

    # The contribution of all k-cycles to alpha_2 is:
    # alpha_2 contribution = sum over k-cycles C, sum over C' disj C of [C' counted]
    # = sum over V active for k, nc(V) * C_odd(comp(V))
    # This double-counts: each disjoint pair {C1, C2} is counted from both sides.
    # But actually no: we're summing over ONE cycle C and counting all C' disjoint.
    # So this gives the TOTAL "disjoint degree sum" = 2 * alpha_2(involving k-cycles)
    # plus disjoint pairs where BOTH are k-cycles (counted twice).
    # Hmm, this is getting complicated. Let's just use the direct approach.

    cod_vals = sorted(set(cod for _, _, cod in entries))
    comp_size = p - k

    print(f"\n  k={k}: {n_sets} active sets, total directed cycles = {total_nc}")
    print(f"    complement size = {comp_size}")
    print(f"    C_odd(comp) values: {cod_vals}")
    print(f"    Weighted avg C_odd(comp) = {avg_cod:.2f}")
    print(f"    Per-cycle 'disjoint contribution' = {avg_cod:.2f}")
    print(f"    Per-cycle effective dH = 2 + correction from disjoint structure")

# Direct computation: alpha_2 decomposed by (k1, k2) source length
print(f"\n{'='*60}")
print(f"  ALPHA_2 DECOMPOSITION BY (k1, k2)")
print(f"{'='*60}")

all_active = []
for k in sorted(active_by_k):
    for fs, nc, cod in active_by_k[k]:
        all_active.append((fs, k, nc))

disj_kk = defaultdict(int)
for i in range(len(all_active)):
    for j in range(i+1, len(all_active)):
        V1, k1, n1 = all_active[i]
        V2, k2, n2 = all_active[j]
        if not (V1 & V2):
            key = tuple(sorted([k1, k2]))
            disj_kk[key] += n1 * n2

total_a2 = sum(disj_kk.values())
print(f"\n  Paley p={p} alpha_2 = {total_a2}")
for kk in sorted(disj_kk):
    val = disj_kk[kk]
    print(f"    disj({kk[0]},{kk[1]}) = {val:>6} ({100*val/total_a2:>5.1f}%)")

# Now the crucial insight: when orientation changes from Paley,
# which (k1,k2) terms change most?
# From the earlier data:
# Paley:  disj(3,3)=495, disj(3,5)=3630, disj(3,7)=4840, disj(5,5)=1914
# ClassB: disj(3,3)=517, disj(3,5)=3740, disj(3,7)=4862, disj(5,5)=2101
# ClassC: disj(3,3)=550, disj(3,5)=3476, disj(3,7)=5434, disj(5,5)=1650
# ClassD: disj(3,3)=506, disj(3,5)=3652, disj(3,7)=4840, disj(5,5)=1914

disj_data = {
    'A': {(3,3): 495, (3,5): 3630, (3,7): 4840, (5,5): 1914},
    'B': {(3,3): 517, (3,5): 3740, (3,7): 4862, (5,5): 2101},
    'C': {(3,3): 550, (3,5): 3476, (3,7): 5434, (5,5): 1650},
    'D': {(3,3): 506, (3,5): 3652, (3,7): 4840, (5,5): 1914},
}

print(f"\n{'='*60}")
print(f"  DISJ PAIR CHANGES FROM PALEY (per type)")
print(f"{'='*60}")
for cl in ['B', 'C', 'D']:
    print(f"\n  Class {cl} vs Paley:")
    total_diff = 0
    for kk in sorted(disj_data['A']):
        diff = disj_data[cl][kk] - disj_data['A'][kk]
        total_diff += diff
        print(f"    d_disj({kk[0]},{kk[1]}) = {diff:>+5}")
    print(f"    Total d(alpha_2) = {total_diff:>+5}")

# The disj(3,3) increase is interesting: c3 is constant (55) for all orientations!
# How can disj(3,3) change if c3 is constant?
# Because disj(3,3) = sum_{V1,V2 disj 3-sets} n(V1)*n(V2)
# The vertex SETS are all C(11,3)=165, but the n(V) varies (0 or 2).
# When n(V) changes (3-set becomes transitive vs cyclic), the disjoint
# pair weight changes.
# But c3 = 55 means exactly 55 out of 165 three-sets are cyclic, always.
# The NUMBER of cyclic 3-sets is constant, but WHICH sets are cyclic varies!
# And this changes the disjoint pair structure.

print(f"\n  c3 is always 55: same NUMBER of active 3-sets.")
print(f"  But WHICH 3-sets are active varies by orientation.")
print(f"  This changes disjoint pair count even for (3,3) pairs.")

# Actually: for CIRCULANT tournaments on Z_p, all 3-sets form orbits of
# size p under translation. There are C(p,3)/p orbits.
# At p=11: 165/11 = 15 orbits.
# Each orbit has CONSTANT n(V) within the orbit (by circulant symmetry).
# So the structure is: which of the 15 orbits are active?
# 55/11 = 5 orbits per position... wait: 55 active 3-sets / 11 per orbit = 5 orbits.
# There are 15 total orbits, 5 are active (cyclic), 10 are transitive.
# This should be the same for ALL circulant orientations (since c3=55 is constant).

# But is the SAME set of 5 orbits active for all orientations?
# Let me check!

print(f"\n{'='*60}")
print(f"  3-SET ORBIT ANALYSIS: which orbits are active per orientation?")
print(f"{'='*60}")

# Enumerate the 15 orbits
orbits_3 = []
seen = set()
for subset in combinations(range(p), 3):
    canon = min(frozenset((v+t) % p for v in subset) for t in range(p))
    if canon not in seen:
        seen.add(canon)
        orbits_3.append(canon)

print(f"  Total 3-set orbits: {len(orbits_3)}")

# For a few orientations, check which orbits are active
for bits in [29, 4, 0, 1]:  # A, B, C, D representatives
    S = []
    for j in range(m):
        if bits & (1 << j):
            S.append(j + 1)
        else:
            S.append(p - (j + 1))

    A_mat = build_adj(p, S)

    active_orbits = set()
    for orb in orbits_3:
        rep = sorted(orb)
        nc = count_ham_cycles(A_mat, rep)
        if nc > 0:
            active_orbits.add(orb)

    print(f"\n  bits={bits}: {len(active_orbits)} active orbits out of {len(orbits_3)}")
    for orb in sorted(active_orbits, key=lambda x: tuple(sorted(x))):
        rep = sorted(orb)
        # Check if ALL orbits have the same orbit representatives
        diffs = [(rep[1]-rep[0]) % p, (rep[2]-rep[0]) % p]
        print(f"    orbit {rep}: diffs = {diffs}")

print("\nDONE.")
