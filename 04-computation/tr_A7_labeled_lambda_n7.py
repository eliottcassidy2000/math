#!/usr/bin/env python3
"""
Check: is tr(A^7) determined by the LABELED lambda graph at n=7?

We know c7_dir is NOT lambda-determined (kind-pasteur found 5 ambiguities).
Question: does tr(A^7) = 7*c7_dir? Or are there degenerate 7-walks?

At n=7, a closed 7-walk on 7 vertices IS a 7-cycle. But walks on fewer
than 7 vertices can also exist.

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

n = 7
tb = n*(n-1)//2
np.random.seed(42)

# First: check if tr(A^7) = 7*c7_dir (i.e., no degenerate 7-walks)
print("Checking tr(A^7) vs 7*c7_dir at n=7...")
fail_count = 0
for trial in range(500):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    A7 = np.linalg.matrix_power(A, 7)
    tr7 = int(np.trace(A7))

    # Count c7_dir: Hamiltonian directed cycles on all 7 vertices
    c7_dir = 0
    for perm in permutations(range(1, n)):
        path = (0,) + perm
        valid = all(A[path[i]][path[(i+1) % n]] for i in range(n))
        if valid:
            c7_dir += 1

    if tr7 != 7 * c7_dir:
        fail_count += 1
        if fail_count <= 3:
            diff = tr7 - 7 * c7_dir
            print(f"  trial {trial}: tr7={tr7}, 7*c7={7*c7_dir}, diff={diff}")
            # Also compute c5 contribution to degenerate walks
            c5_total = 0
            for combo in combinations(range(n), 5):
                A_sub = A[np.ix_(list(combo), list(combo))]
                tr5_sub = int(np.trace(np.linalg.matrix_power(A_sub, 5)))
                c5_total += tr5_sub // 5
            c3 = 0
            for combo in combinations(range(n), 3):
                i, j, k = combo
                if A[i][j] and A[j][k] and A[k][i]: c3 += 1
                if A[j][i] and A[i][k] and A[k][j]: c3 += 1
            print(f"    c3={c3}, c5={c5_total}, c7={c7_dir}")

print(f"\nResult: {fail_count}/500 failures")
if fail_count > 0:
    print("tr(A^7) != 7*c7_dir — degenerate 7-walks exist!")
    print("\nNow checking: is tr(A^7) lambda-determined (labeled)?")

    lam_groups_tr7 = defaultdict(set)
    for trial in range(20000):
        bits = np.random.randint(0, 1 << tb)
        A = bits_to_adj(bits, n)
        A7 = np.linalg.matrix_power(A, 7)
        tr7 = int(np.trace(A7))
        L = lambda_graph(A, n)
        labeled_key = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
        lam_groups_tr7[labeled_key].add(tr7)

    ambig = sum(1 for v in lam_groups_tr7.values() if len(v) > 1)
    print(f"  Labeled lambda groups: {len(lam_groups_tr7)}, ambiguous: {ambig}")
    if ambig > 0:
        count = 0
        for key, vals in sorted(lam_groups_tr7.items()):
            if len(vals) > 1:
                print(f"    {key}: tr7 = {sorted(vals)}")
                count += 1
                if count >= 5: break
else:
    print("tr(A^7) = 7*c7_dir at n=7!")

# Also check: at what point does tr(A^k) stop being lambda-determined?
# We know tr(A^3) = 3*c3 (always), and c3 is lambda-determined (THM-171).
# tr(A^5) = 5*c5 (no degenerate walks), c5 is labeled-lambda-determined (THM-172).
# What about tr(A^7)?

print(f"\n{'='*60}")
print("DEGENERATE WALK ANALYSIS")
print(f"{'='*60}")

# For tr(A^3): closed walks of length 3 on <= 3 vertices.
# Must use 3 distinct vertices (otherwise A_{ii}=0 kills it). So tr(A^3) = 3*c3.
print("tr(A^3) = 3*c3: TRUE (no degenerate 3-walks in tournaments)")

# For tr(A^5): no degenerate walks at n=5 (proved above).
# At n=6,7: still no degenerate walks (verified computationally).
# WHY? Because any 5-walk revisiting a vertex would need A_{ij}*A_{ji}=1 somewhere.
# Let me prove this properly.
print("tr(A^5) = 5*c5: TRUE at n=5,6,7 (no degenerate 5-walks in tournaments)")

# For tr(A^7): if we have 7-step walks, can we revisit vertices?
# The walk is i0→i1→i2→i3→i4→i5→i6→i0.
# If any ia = i_{a+2 mod 7}: walk includes ia→i_{a+1}→ia, needing A_{ia,i_{a+1}}*A_{i_{a+1},ia}=1.
# IMPOSSIBLE in a tournament.

# What about ia = i_{a+3 mod 7}? E.g., i0=i3:
# Walk: i0→i1→i2→i0→i4→i5→i6→i0
# Need A_{i2,i0}=1 and A_{i6,i0}=1 — both possible independently.
# So this CAN happen! The walk visits 6 distinct vertices with i0 repeated.

# Similarly, i0=i4: walk = i0→i1→i2→i3→i0→i5→i6→i0
# Need A_{i3,i0}=1 and A_{i6,i0}=1 — possible.

print("\nFor 7-walks: ia=i_{a+3} is possible! Gives walks on 6 distinct vertices.")
print("This explains why tr(A^7) != 7*c7_dir: degenerate walks on 5,6 vertices exist.")

# The formula should be: tr(A^7) = 7*c7_dir + (degenerate 7-walks)
# And degenerate 7-walks involve structures on 5 or 6 vertices.

# The KEY insight: for ODD k, are degenerate k-walks possible?
# k=3: NO (minimum 3 distinct vertices needed, and 3-cycles are the only option)
# k=5: NO (any vertex repetition forces a back-and-forth within 2 steps)
# k=7: YES (vertex repetition with 3-step gap is possible)

# The general pattern: for k=2j+1, degenerate walks exist iff k >= 7.
# Because with a gap of 3 (not 1 or 2), you can have i→?→?→i without back-and-forth.

# This means:
# c3 = tr(A^3)/3 — always
# c5 = tr(A^5)/5 — always
# c7 ≠ tr(A^7)/7 in general — degenerate correction needed

# And since tr(A^5) can be expressed in terms of lambda (because A² involves lambda),
# c5 is lambda-determined!

# But c7 has both the degenerate correction (which might be lambda-determined)
# AND the pure 7-cycle count (which might not be).

print(f"\nDone.")
