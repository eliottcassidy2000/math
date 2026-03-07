#!/usr/bin/env python3
"""
CONJECTURE: For position-uniform tournaments, ALL valid split pairs
(E_a * B_b > 0) on complementary vertex sets form Hamiltonian paths
when concatenated.

If true: M[a,b] = sum_j (-1)^j consec(a,b,j) for position-uniform T.

Then: M[a,b] = 0 iff sum_j (-1)^j consec(a,b,j) = 0.

For T[a,b]=0: no path has a->b consecutive, so consec=0 and M[a,b]=0.
For T[a,b]=1: need alternating sum of consecutive counts = 0.

NOTE: For symmetric M, M[a,b] = M[b,a]. If T[a,b]=1 then T[b,a]=0.
  M[a,b] = sum_j (-1)^j consec(a,b,j)  (nontrivial, T[a,b]=1)
  M[b,a] = sum_j (-1)^j consec(b,a,j) = 0  (trivial, T[b,a]=0)
  Symmetry M[a,b] = M[b,a] = 0 would then be automatic!

WAIT — but M[a,b] = M[b,a] is proved (THM-030), and if the formula holds,
then M[a,b] = some nontrivial sum and M[b,a] = 0. But M[a,b] = M[b,a]
means the nontrivial sum MUST also be 0.

So the formula M[a,b] = sum_j (-1)^j consec(a,b,j) + SYMMETRY gives:
  For T[a,b]=1: M[a,b] = sum(-1)^j consec(a,b,j) = M[b,a] = 0.

THIS WOULD BE A PROOF THAT POSITION-UNIFORM => M SCALAR!

Let me verify the "all split pairs are Ham paths" claim exhaustively.

kind-pasteur-2026-03-06-S25c
"""

from itertools import permutations, combinations
import numpy as np

def E_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[-1] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

def B_v(T, verts, v):
    verts = list(verts)
    if len(verts) == 1:
        return 1 if verts[0] == v else 0
    count = 0
    for p in permutations(verts):
        if p[0] != v: continue
        valid = True
        for k in range(len(p)-1):
            if T.get((p[k], p[k+1]), 0) != 1:
                valid = False; break
        if valid: count += 1
    return count

def count_H(T, n):
    count = 0
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        count += prod
    return count

def position_matrix(T, n):
    P = np.zeros((n, n), dtype=int)
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            for k in range(n):
                P[perm[k], k] += 1
    return P

def tournament_from_bits(n, bits):
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    T = {}
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1:
            T[(i,j)] = 1; T[(j,i)] = 0
        else:
            T[(i,j)] = 0; T[(j,i)] = 1
    return T


# ============================================================
# n=5: Verify "all split pairs are Ham paths" for position-uniform
# ============================================================
print("=" * 70)
print("n=5: Do all valid split pairs form Ham paths?")
print("=" * 70)

n = 5
pairs_list = [(i,j) for i in range(n) for j in range(i+1, n)]

total_uniform = 0
split_always_ham = 0
split_not_always_ham = 0

for bits in range(1 << len(pairs_list)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)

    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))
    if not is_uniform:
        continue

    total_uniform += 1
    all_ham = True

    # Check ALL pairs (a,b)
    for a in range(n):
        for b in range(n):
            if a == b:
                continue
            U = [v for v in range(n) if v != a and v != b]

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                # Find all valid (P,Q) pairs
                for p in permutations(S_set):
                    if p[-1] != a: continue
                    p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
                    if not p_valid: continue

                    for q in permutations(R_set):
                        if q[0] != b: continue
                        q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                        if not q_valid: continue

                        # Check concatenation
                        concat = list(p) + list(q)
                        is_ham = all(T.get((concat[k], concat[k+1]), 0) == 1 for k in range(len(concat)-1))
                        if not is_ham:
                            all_ham = False
                            if split_not_always_ham < 3:
                                print(f"  bits={bits}: P={p}, Q={q} NOT Ham path!")
                                print(f"    Junction T[{a},{b}]={T.get((a,b),0)}")

            if not all_ham:
                break
        if not all_ham:
            break

    if all_ham:
        split_always_ham += 1
    else:
        split_not_always_ham += 1

print(f"\n  Position-uniform tournaments: {total_uniform}")
print(f"  All split pairs are Ham paths: {split_always_ham}")
print(f"  Some split pairs are NOT Ham paths: {split_not_always_ham}")

if split_not_always_ham == 0:
    print("""
  CONFIRMED: For ALL position-uniform tournaments at n=5,
  every valid split pair (P ending at a, Q starting at b)
  on complementary vertex sets forms a Hamiltonian path
  when concatenated!

  This means: M[a,b] = sum_j (-1)^j consec(a,b,j)
  for ALL position-uniform tournaments.

  Combined with THM-030 (symmetry M[a,b] = M[b,a]):
    For T[a,b]=1: M[a,b] = sum(-1)^j consec(a,b,j)
    For T[b,a]=0: M[b,a] = 0 (trivially)
    Symmetry: M[a,b] = M[b,a] = 0.

  THEREFORE: Position-uniform => M[a,b] = 0 for all a != b.
  Combined with diagonal: M = (H/n)*I.

  Q.E.D. (modulo verifying that split pairs = Ham paths
  holds at all n, not just n=5.)
""")
else:
    print("  Some counterexamples found — need investigation.")


# ============================================================
# Also check n=3
# ============================================================
print("=" * 70)
print("n=3: Verify split pairs = Ham paths")
print("=" * 70)

n = 3
pairs_3 = [(i,j) for i in range(n) for j in range(i+1, n)]

for bits in range(1 << len(pairs_3)):
    T = tournament_from_bits(n, bits)
    H = count_H(T, n)
    P = position_matrix(T, n)
    is_uniform = H % n == 0 and all(P[v,k] == H // n for v in range(n) for k in range(n))

    if not is_uniform:
        continue

    print(f"\n  bits={bits}: H={H}, uniform_pos=True")

    all_ham = True
    for a in range(n):
        for b in range(n):
            if a == b: continue
            U = [v for v in range(n) if v != a and v != b]

            for mask in range(1 << len(U)):
                S_list = [U[k] for k in range(len(U)) if mask & (1 << k)]
                R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                S_set = sorted(set(S_list) | {a})
                R_set = sorted(set(R) | {b})

                for p in permutations(S_set):
                    if p[-1] != a: continue
                    p_valid = all(T.get((p[k], p[k+1]), 0) == 1 for k in range(len(p)-1))
                    if not p_valid: continue

                    for q in permutations(R_set):
                        if q[0] != b: continue
                        q_valid = all(T.get((q[k], q[k+1]), 0) == 1 for k in range(len(q)-1))
                        if not q_valid: continue

                        concat = list(p) + list(q)
                        is_ham = all(T.get((concat[k], concat[k+1]), 0) == 1 for k in range(len(concat)-1))
                        if not is_ham:
                            all_ham = False
                            print(f"    P={p}, Q={q} NOT Ham path! T[{a},{b}]={T.get((a,b),0)}")

    print(f"    All split pairs form Ham paths: {all_ham}")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
