"""
hierarchy_c9_test.py -- kind-pasteur-2026-03-13-S61

Does (lambda, sigma) determine c9?
If yes: Level 1.5 captures all odd cycles up to 9.
If no: need Level 1.75 or Level 2.

Also: what is the minimal additional invariant beyond lambda
that determines c9?

At n=9: this is the first n where c9 (Hamiltonian directed cycles) exists.
tr(A^9)/9 = c9_dir.

Approach: at n=9, sample tournaments, group by (lambda, sigma), check c9.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict

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
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

# First: does lambda determine c9 at n=9? (Almost certainly not — c7 already fails at n=7)
print("=" * 60)
print("DOES LAMBDA DETERMINE c9 AT n=9?")
print("=" * 60)

n = 9
total_bits = n * (n-1) // 2

np.random.seed(42)
lambda_c9 = defaultdict(set)
lambda_c7 = defaultdict(set)
lambda_c5 = defaultdict(set)

for trial in range(10000):
    # Generate random bits (n=9 has C(9,2)=36 bits)
    bits = int(np.random.randint(0, 2**31)) | (int(np.random.randint(0, 2**5)) << 31)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    lam_key = tuple(L[i][j] for i in range(n) for j in range(n))

    c5 = count_directed_k_cycles(A, n, 5)
    c7 = count_directed_k_cycles(A, n, 7)
    c9 = count_directed_k_cycles(A, n, 9)

    lambda_c5[lam_key].add(c5)
    lambda_c7[lam_key].add(c7)
    lambda_c9[lam_key].add(c9)

    if trial % 2000 == 0 and trial > 0:
        amb5 = sum(1 for v in lambda_c5.values() if len(v) > 1)
        amb7 = sum(1 for v in lambda_c7.values() if len(v) > 1)
        amb9 = sum(1 for v in lambda_c9.values() if len(v) > 1)
        print(f"  {trial}/10000: c5 amb={amb5}, c7 amb={amb7}, c9 amb={amb9}")

amb5 = sum(1 for v in lambda_c5.values() if len(v) > 1)
amb7 = sum(1 for v in lambda_c7.values() if len(v) > 1)
amb9 = sum(1 for v in lambda_c9.values() if len(v) > 1)
print(f"\nFinal: lambda groups={len(lambda_c5)}")
print(f"  c5 ambiguous: {amb5}")
print(f"  c7 ambiguous: {amb7}")
print(f"  c9 ambiguous: {amb9}")
print(f"  Lambda determines c5? {amb5 == 0}")
print(f"  Lambda determines c7? {amb7 == 0}")
print(f"  Lambda determines c9? {amb9 == 0}")

# Now test: does (lambda, sigma) determine c9?
print(f"\n{'='*60}")
print("DOES (LAMBDA, SIGMA) DETERMINE c9 AT n=9?")
print(f"{'='*60}")

np.random.seed(42)
ls_c7 = defaultdict(set)
ls_c9 = defaultdict(set)

for trial in range(10000):
    bits = int(np.random.randint(0, 2**31)) | (int(np.random.randint(0, 2**5)) << 31)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)
    A2 = A @ A

    lam_key = tuple(L[i][j] for i in range(n) for j in range(n))
    sig_vals = []
    for u in range(n):
        for v in range(u+1, n):
            sig_vals.append(n - 2 - int(A2[u][v]) - int(A2[v][u]))
    sig_key = tuple(sig_vals)
    combined = lam_key + sig_key

    c7 = count_directed_k_cycles(A, n, 7)
    c9 = count_directed_k_cycles(A, n, 9)

    ls_c7[combined].add(c7)
    ls_c9[combined].add(c9)

    if trial % 2000 == 0 and trial > 0:
        amb7 = sum(1 for v in ls_c7.values() if len(v) > 1)
        amb9 = sum(1 for v in ls_c9.values() if len(v) > 1)
        print(f"  {trial}/10000: c7 amb={amb7}, c9 amb={amb9}")

amb7_ls = sum(1 for v in ls_c7.values() if len(v) > 1)
amb9_ls = sum(1 for v in ls_c9.values() if len(v) > 1)
print(f"\nFinal: (lambda,sigma) groups={len(ls_c7)}")
print(f"  c7 ambiguous: {amb7_ls}")
print(f"  c9 ambiguous: {amb9_ls}")
print(f"  (Lambda,sigma) determines c7? {amb7_ls == 0}")
print(f"  (Lambda,sigma) determines c9? {amb9_ls == 0}")

# If (lambda, sigma) doesn't determine c9, try A^2
if amb9_ls > 0:
    print(f"\n{'='*60}")
    print("DOES A^2 DETERMINE c9 AT n=9?")
    print(f"{'='*60}")

    np.random.seed(42)
    a2_c9 = defaultdict(set)

    for trial in range(10000):
        bits = int(np.random.randint(0, 2**31)) | (int(np.random.randint(0, 2**5)) << 31)
        A = bits_to_adj(bits, n)
        A2 = A @ A
        a2_key = tuple(int(A2[i][j]) for i in range(n) for j in range(n))
        c9 = count_directed_k_cycles(A, n, 9)
        a2_c9[a2_key].add(c9)

    amb9_a2 = sum(1 for v in a2_c9.values() if len(v) > 1)
    print(f"  A^2 groups: {len(a2_c9)}")
    print(f"  c9 ambiguous: {amb9_a2}")
    print(f"  A^2 determines c9? {amb9_a2 == 0}")

# Also: does A^2 UNIQUELY determine A at n=7?
print(f"\n{'='*60}")
print("DOES A^2 DETERMINE A UNIQUELY?")
print(f"{'='*60}")

for n_test in [6, 7]:
    total_bits_test = n_test * (n_test - 1) // 2

    if n_test <= 7:
        a2_tours = defaultdict(list)
        # For n=7, 2^21 = 2M — we can sample
        np.random.seed(42)
        max_trials = min(1 << total_bits_test, 200000)
        for trial in range(max_trials):
            if total_bits_test <= 15:
                bits = trial
            else:
                bits = np.random.randint(0, 1 << total_bits_test)
            A = bits_to_adj(bits, n_test)
            A2 = A @ A
            key = tuple(int(A2[i][j]) for i in range(n_test) for j in range(n_test))
            a2_tours[key].append(bits)

        multi = sum(1 for v in a2_tours.values() if len(v) > 1)
        max_fiber = max(len(v) for v in a2_tours.values())
        print(f"\n  n={n_test}: {len(a2_tours)} A^2 classes from {max_trials} tournaments")
        print(f"  Multiple-tournament classes: {multi}")
        print(f"  Max fiber size: {max_fiber}")

        if multi > 0:
            # Show example
            for key, tours in a2_tours.items():
                if len(tours) > 1:
                    print(f"  Example: {len(tours)} tournaments with same A^2")
                    for b in tours[:3]:
                        A = bits_to_adj(b, n_test)
                        c7 = count_directed_k_cycles(A, n_test, min(7, n_test)) if n_test >= 7 else 0
                        c5 = count_directed_k_cycles(A, n_test, 5)
                        print(f"    bits={b}: scores={sorted(int(sum(A[i])) for i in range(n_test))}, c5={c5}, c7={c7}")
                    break

print("\nDone.")
