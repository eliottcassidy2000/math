#!/usr/bin/env python3
"""
Attempt algebraic proof of NONHAM=0 for position-uniform tournaments.

Strategy: At n=3, the only position-uniform tournament is the 3-cycle.
Let's prove M[a,b]=0 for T[a,b]=0 directly, then see if the argument
generalizes.

For the 3-cycle T = {0->1, 1->2, 2->0}:
  T[0,1]=1, T[1,2]=1, T[2,0]=1
  T[1,0]=0, T[0,2]=0, T[2,1]=0

Consider (a,b)=(1,0) where T[1,0]=0.
U = {2} (single vertex).

M[1,0] = sum_{S subset of {2}} (-1)^|S| E_1(S+{1}) * B_0((U\S)+{0})

S={}: E_1({1})=1, B_0({0,2})=#{paths in {0,2} starting at 0 with T-edges}
  = #{(0,2)} if T[0,2]=1 else 0. T[0,2]=0, so B_0({0,2})=0.
  Contribution: (+1) * 1 * 0 = 0

S={2}: E_1({1,2})=#{paths in {1,2} ending at 1 with T-edges}
  = #{(2,1)} if T[2,1]=1 else 0. T[2,1]=0, so E_1({1,2})=0.
  Contribution: (-1) * 0 * 1 = 0

M[1,0] = 0 + 0 = 0. Trivially zero because both subsets give 0!

Now (a,b)=(0,2) where T[0,2]=0.
U = {1}.

S={}: E_0({0})=1, B_2({1,2})=#{paths in {1,2} starting at 2}
  = #{(2,1)} if T[2,1]=1 else 0. T[2,1]=0, so B_2({1,2})=0.
  Contribution: 0

S={1}: E_0({0,1})=#{paths in {0,1} ending at 0}
  = #{(1,0)} if T[1,0]=1 else 0. T[1,0]=0, so E_0({0,1})=0.
  Contribution: 0

M[0,2] = 0. Again trivially zero.

Interesting: at n=3, M[a,b]=0 for T[a,b]=0 because the subtournament
paths themselves vanish. No cancellation needed!

Let's check n=5 (Paley). Here the cancellation is nontrivial.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations
import numpy as np

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def E_paths(T, verts, a):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == a else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != a: continue
        valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
        if valid: count += 1
    return count

def B_paths(T, verts, b):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == b else 0
    count = 0
    for p in permutations(verts):
        if p[0] != b: continue
        valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
        if valid: count += 1
    return count

# ============================================================
# n=5 Paley: Analyze WHY M[0,2]=0 (T[0,2]=0)
# ============================================================
print("=" * 70)
print("n=5 Paley T_5: Algebraic analysis of M[0,2]=0")
print("=" * 70)

n = 5
T = circulant_tournament(n, {1, 2})  # 3-cycle: 0->1->2->0
print(f"n=3 cycle: trivially 0 (subtournament paths vanish)")

# Paley T_5: QR mod 5 = {1,4}
T5 = circulant_tournament(5, {1, 4})
print(f"\nn=5 Paley T_5: gen_set={{1,4}}")
print(f"T[0,1]={T5[(0,1)]}, T[0,2]={T5[(0,2)]}, T[0,3]={T5[(0,3)]}, T[0,4]={T5[(0,4)]}")

# Check M[0,2] where T[0,2]=0
a, b = 0, 2
U = [v for v in range(5) if v != a and v != b]
print(f"\nM[{a},{b}] analysis (T[{a},{b}]={T5[(a,b)]}):")
print(f"U = {U}")

total = 0
for mask in range(1 << len(U)):
    S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
    sign = (-1)**len(S_list)
    S_set = sorted(set(S_list) | {a})
    R_set = sorted(set(R) | {b})

    ea = E_paths(T5, S_set, a)
    bb = B_paths(T5, R_set, b)

    contrib = sign * ea * bb
    total += contrib

    # Show the paths
    e_paths = []
    for p in permutations(S_set):
        if p[-1] != a: continue
        if all(T5.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
            e_paths.append(p)

    b_paths = []
    for p in permutations(R_set):
        if p[0] != b: continue
        if all(T5.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1)):
            b_paths.append(p)

    print(f"  S={str(sorted(S_list)):>10}, |S|={len(S_list)}, sign={sign:+d}:")
    print(f"    S+a={S_set}, E-paths ending at {a}: {e_paths} (count={ea})")
    print(f"    R+b={R_set}, B-paths starting at {b}: {b_paths} (count={bb})")
    print(f"    Contribution: {sign:+d} * {ea} * {bb} = {contrib:+d}")

print(f"\n  M[{a},{b}] = {total}")


# ============================================================
# Now analyze the pairing: which subsets cancel each other?
# ============================================================
print("\n" + "=" * 70)
print("Pairing analysis: which subsets cancel?")
print("=" * 70)

# Analyze for all non-edge pairs
for a in range(5):
    for b in range(5):
        if a == b: continue
        if T5[(a,b)] == 1: continue  # Only non-edges

        U = [v for v in range(5) if v != a and v != b]
        print(f"\n  ({a},{b}) T[{a},{b}]=0:")

        nonzero_subsets = []
        for mask in range(1 << len(U)):
            S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
            R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
            sign = (-1)**len(S_list)
            S_set = sorted(set(S_list) | {a})
            R_set = sorted(set(R) | {b})

            ea = E_paths(T5, S_set, a)
            bb = B_paths(T5, R_set, b)

            if ea > 0 and bb > 0:
                nonzero_subsets.append((sorted(S_list), ea, bb, sign * ea * bb))

        for s, ea, bb, c in nonzero_subsets:
            print(f"    S={s}: E={ea}, B={bb}, contrib={c:+d}")

        # Show the pairing
        positive = [(s, ea, bb, c) for s, ea, bb, c in nonzero_subsets if c > 0]
        negative = [(s, ea, bb, c) for s, ea, bb, c in nonzero_subsets if c < 0]
        print(f"    Positive: {sum(c for _,_,_,c in positive)}, Negative: {sum(c for _,_,_,c in negative)}")


# ============================================================
# KEY QUESTION: For T[a,b]=0, does T[a,v]*T[v,b] pattern matter?
# ============================================================
print("\n" + "=" * 70)
print("Edge pattern analysis for non-edge pairs")
print("=" * 70)

for a in range(5):
    for b in range(5):
        if a == b: continue
        if T5[(a,b)] == 1: continue

        U = [v for v in range(5) if v != a and v != b]
        patterns = {}
        for v in U:
            p = (T5.get((a,v),0), T5.get((v,b),0))
            patterns[v] = p

        print(f"  ({a},{b}): " + ", ".join(f"T[{a},{v}]={T5[(a,v)]},T[{v},{b}]={T5[(v,b)]}" for v in U))

        # For T[a,b]=0, the "mediators" (v where a->v->b) should balance
        mediators = [v for v in U if T5[(a,v)] == 1 and T5[(v,b)] == 1]
        blockers = [v for v in U if T5[(a,v)] == 0 or T5[(v,b)] == 0]
        print(f"    Mediators (a->v->b): {mediators}")
        print(f"    Non-mediators: {blockers}")


# ============================================================
# What about reversal symmetry? E_a(W,T) = B_a(W,T^op)
# ============================================================
print("\n" + "=" * 70)
print("Reversal identity: E_a(W,T) vs B_a(W,T^op)")
print("=" * 70)

T5_op = {}
for i in range(5):
    for j in range(5):
        if i == j: continue
        T5_op[(i,j)] = T5[(j,i)]

print("T5 gen_set = {1,4}")
print("T5^op gen_set = {5-1, 5-4} = {4, 1} = {1,4}")
print("T5 = T5^op! (Paley T_5 is self-converse)")

for a in range(3):
    for W_mask in range(1, 1 << 4):
        W = [v for v in range(5) if W_mask & (1 << v)]
        if a not in W: continue
        ea = E_paths(T5, W, a)
        ba_op = B_paths(T5_op, W, a)
        ba = B_paths(T5, W, a)
        if ea != ba_op and ea > 0:
            print(f"  a={a}, W={W}: E_a={ea}, B_a(T^op)={ba_op} -- MISMATCH!")
        if ea > 0:
            print(f"  a={a}, W={W}: E_a(T)={ea}, B_a(T)={ba}, B_a(T^op)={ba_op}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
