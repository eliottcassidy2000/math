#!/usr/bin/env python3
"""
Verify whether tr(A^5) = 5*c5_dir at n=7, and check c5_dir vs lambda more carefully.

Key question: is c5_dir (total directed 5-cycle count) the same as tr(A^5)/5 at n>=6?

At n=5,6: CONFIRMED tr(A^5) = 5*c5_dir (no degenerate walks).
At n>=6: there MIGHT be degenerate 5-walks on 5 vertices that traverse some vertices
from outside the 5-subset. Wait, no — tr(A^5) counts ALL closed 5-walks on the full
vertex set, not just on 5-vertex subsets.

IMPORTANT DISTINCTION:
- c5_dir = number of directed 5-cycles = #{5-element subsets V ⊂ [n] with Hamiltonian directed cycle on V}
- tr(A^5) = number of closed walks of length 5 in the full graph

At n>5, a closed walk of length 5 can visit vertices from DIFFERENT 5-subsets.
Moreover, degenerate walks (with repeated vertices) could exist when n>5
because there are more routing options.

Let me check this carefully.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
from itertools import combinations, permutations
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

# At n=7, check tr(A^5) vs 5*c5_dir
print("Checking tr(A^5) vs 5*c5_dir at n=7...")
n = 7
tb = n*(n-1)//2
np.random.seed(42)

fail_count = 0
for trial in range(1000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    A5 = np.linalg.matrix_power(A, 5)
    tr5 = int(np.trace(A5))

    c5_total = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = all(A[path[i]][path[(i+1) % 5]] for i in range(5))
            if valid:
                c5_total += 1

    if tr5 != 5 * c5_total:
        fail_count += 1
        if fail_count <= 3:
            # Decompose
            w_count = defaultdict(int)
            for i0 in range(n):
                for i1 in range(n):
                    if i1==i0 or A[i0][i1]==0: continue
                    for i2 in range(n):
                        if i2==i1 or A[i1][i2]==0: continue
                        for i3 in range(n):
                            if i3==i2 or A[i2][i3]==0: continue
                            for i4 in range(n):
                                if i4==i3 or A[i3][i4]==0: continue
                                if i4==i0 or A[i4][i0]==0: continue
                                nd = len({i0,i1,i2,i3,i4})
                                w_count[nd] += 1
            print(f"  tr5={tr5}, 5*c5={5*c5_total}, diff={tr5-5*c5_total}")
            print(f"    walks by distinct vertices: {dict(sorted(w_count.items()))}")

print(f"\nResult: {fail_count}/{1000} failures")
if fail_count > 0:
    print("tr(A^5) ≠ 5*c5_dir at n=7!")
    print("The difference comes from degenerate walks — walks on 5 vertices")
    print("that reuse vertices but are NOT on a 5-vertex induced subtournament.")
else:
    print("tr(A^5) = 5*c5_dir at n=7 too!")

# Now: more carefully check c5_dir vs lambda at n=7
# Kind-pasteur used 50k samples and found 0 ambiguities.
# But my tr(A^5) test found 35 ambiguous lambda groups.
# If tr(A^5) ≠ 5*c5_dir, then the two could differ:
# c5_dir might be lambda-determined even if tr(A^5) is not.

print(f"\n{'='*60}")
print("DIRECT CHECK: c5_dir vs lambda at n=7")
print(f"{'='*60}")

lam_groups_c5 = defaultdict(set)
for trial in range(20000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lam_key = tuple(sorted(L[i][j] for i in range(n) for j in range(i+1, n)))

    c5_total = 0
    for combo in combinations(range(n), 5):
        A_sub = A[np.ix_(list(combo), list(combo))]
        tr5_sub = int(np.trace(np.linalg.matrix_power(A_sub, 5)))
        c5_total += tr5_sub // 5  # At n=5, tr(A^5) = 5*c5_dir

    lam_groups_c5[lam_key].add(c5_total)

ambig_c5 = sum(1 for v in lam_groups_c5.values() if len(v) > 1)
print(f"  Lambda multiset groups: {len(lam_groups_c5)}")
print(f"  Ambiguous: {ambig_c5}")
if ambig_c5 > 0:
    count_shown = 0
    for key, vals in sorted(lam_groups_c5.items()):
        if len(vals) > 1:
            print(f"    {key}: c5 = {sorted(vals)}")
            count_shown += 1
            if count_shown >= 5:
                break

print("\nDone.")
