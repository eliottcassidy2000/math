#!/usr/bin/env python3
"""
At n=7, labeled lambda does NOT uniquely determine the tournament.
Question: do all tournaments with the same labeled lambda have the same c7?

We know tr7 has 6 ambiguous labeled lambda groups (from earlier test).
But the sampled non-unique pairs all had matching tr7.
Let's check specifically for c7 differences.

opus-2026-03-13-S71c
"""
import sys, time
import numpy as np
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

def lambda_key(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
    return tuple(L[i][j] for i in range(n) for j in range(i+1, n))

n = 7
tb = n*(n-1)//2
np.random.seed(42)

# Store full data per labeled lambda
lam_data = defaultdict(list)

t0 = time.time()
for trial in range(200000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    key = lambda_key(A, n)
    tr5 = int(np.trace(np.linalg.matrix_power(A, 5)))
    tr7 = int(np.trace(np.linalg.matrix_power(A, 7)))
    lam_data[key].append((tr5, tr7, bits))

dt = time.time() - t0
print(f"200k samples, {len(lam_data)} labeled lambda groups, {dt:.1f}s")

# Check for c5 ambiguities
c5_ambig = 0
c7_ambig = 0
tr7_ambig = 0
for key, entries in lam_data.items():
    tr5_set = set(e[0] for e in entries)
    tr7_set = set(e[1] for e in entries)
    if len(tr5_set) > 1:
        c5_ambig += 1
    if len(tr7_set) > 1:
        tr7_ambig += 1
        c7_ambig += 1  # Since tr7 != 7*c7 at n=7, tr7 ambiguity implies structural difference

print(f"\nLabeled lambda ambiguities:")
print(f"  c5 (tr5): {c5_ambig} ambiguous groups")
print(f"  tr7: {tr7_ambig} ambiguous groups")

if tr7_ambig > 0:
    print(f"\ntr7 ambiguous examples:")
    count = 0
    for key, entries in sorted(lam_data.items()):
        tr7_set = set(e[1] for e in entries)
        if len(tr7_set) > 1:
            print(f"  lambda={key[:10]}...: tr7 = {sorted(tr7_set)}")
            count += 1
            if count >= 10: break

# Now let's check: is c7_dir (actual 7-cycle count) determined by labeled lambda?
# We need to compute c7_dir directly, not from tr7 (since tr7 ≠ 7*c7).
from itertools import permutations

print(f"\nComputing actual c7_dir for ambiguous groups...")

for key, entries in lam_data.items():
    tr7_set = set(e[1] for e in entries)
    if len(tr7_set) > 1:
        # Compute c7_dir for each entry
        c7_vals = set()
        for tr5, tr7, bits in entries:
            A = bits_to_adj(bits, n)
            c7 = 0
            for perm in permutations(range(1, n)):
                path = (0,) + perm
                valid = all(A[path[i]][path[(i+1) % n]] for i in range(n))
                if valid: c7 += 1
            c7_vals.add(c7)
        if len(c7_vals) > 1:
            print(f"  lambda={key[:10]}...: c7 = {sorted(c7_vals)} (AMBIGUOUS!)")
        else:
            print(f"  lambda={key[:10]}...: c7 = {c7_vals} (tr7 differs but c7 same)")
        break  # Just check one example

# Actually, we already know c7_dir is NOT labeled-lambda-determined from kind-pasteur's
# finding of 5 ambiguous groups. Let me instead focus on understanding
# WHAT information beyond lambda is needed for c7.

# The "sigma" statistic: sigma(u,v) = #{common successors} + #{common predecessors}
# This was mentioned in THM-173 as part of the sigma-algebra hierarchy.

print(f"\nChecking sigma(u,v) = common_succ + common_pred:")
print(f"Does (lambda, sigma) determine c7?")

def lambda_sigma_key(A, n):
    L = np.zeros((n, n), dtype=int)
    Sig = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
                # Common successors: A_{uw}=1 AND A_{vw}=1
                if A[u][w] and A[v][w]:
                    Sig[u][v] += 1; Sig[v][u] += 1
                # Common predecessors: A_{wu}=1 AND A_{wv}=1
                if A[w][u] and A[w][v]:
                    Sig[u][v] += 1; Sig[v][u] += 1
    lam = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    sig = tuple(Sig[i][j] for i in range(n) for j in range(i+1, n))
    return lam + sig

lam_sig_data = defaultdict(set)
for trial in range(50000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    key = lambda_sigma_key(A, n)
    tr7 = int(np.trace(np.linalg.matrix_power(A, 7)))
    lam_sig_data[key].add(tr7)

ambig_ls = sum(1 for v in lam_sig_data.values() if len(v) > 1)
print(f"  (lambda, sigma) groups: {len(lam_sig_data)}, tr7 ambiguous: {ambig_ls}")

# Also check: does sigma add info beyond lambda?
# Relationship: for a tournament, let d_u = out-degree of u.
# Common successors of (u,v) = #{w: u→w, v→w} = co(u,v)
# Common predecessors = #{w: w→u, w→v} = ci(u,v)
# sigma(u,v) = co(u,v) + ci(u,v) = (n-2) - P_{uv} - P_{vu}
# where P_{uv} = (A²)_{uv}.

# And: lambda(u,v) = A_{uv}·P_{vu} + A_{vu}·P_{uv}
# So: P_{uv} + P_{vu} = (n-2) - sigma(u,v)
# And: lambda(u,v) = A_{uv}·P_{vu} + (1-A_{uv})·P_{uv}
#                   = A_{uv}·(P_{vu} - P_{uv}) + P_{uv}
#                   = A_{uv}·(n-2-sigma-2P_{uv}) + P_{uv}  [using P_{vu}=n-2-sigma-P_{uv}]
# Hmm, this gives:
# lambda = A_{uv}·(n-2-sigma) - A_{uv}·2P_{uv} + P_{uv} + A_{uv}·P_{uv}
# Ugh, getting messy. Let me just check computationally.

# P_{uv} = (A²)_{uv} = #{k: u→k→v}
# For a tournament with i→j: P_{ij} = d_i - 1 - co(i,j)
# Wait: P_{ij} = #{k≠i,j: i→k, k→j}. Among k≠i,j:
# - k is out-neighbor of i: there are d_i - A_{ij} such k
# - Among these, k→j: P_{ij}, k←j (i.e., j→k): co(i,j) if i→k too (which it is)
# Actually: #{k≠i,j: i→k} = d_i - A_{ij}. Among these:
# #{k: i→k, k→j} = P_{ij}
# #{k: i→k, j→k} = co(i,j) (common out, not counting mutual direction effects)
# Wait I need to be more careful. co(i,j) = #{k: i→k AND j→k}. But for k such that i→k,
# either k→j or j→k. So: d_i - A_{ij} = P_{ij} + co(i,j).
# Therefore: P_{ij} = d_i - A_{ij} - co(i,j).
# Similarly: P_{ji} = d_j - A_{ji} - co(j,i) = d_j - (1-A_{ij}) - co(i,j)  [co is symmetric]

# So: P_{ij} + P_{ji} = d_i + d_j - 1 - 2·co(i,j)
# Also: P_{ij} + P_{ji} + co(i,j) + ci(i,j) = n - 2
# So: d_i + d_j - 1 - 2·co(i,j) + co(i,j) + ci(i,j) = n-2
# d_i + d_j - 1 - co(i,j) + ci(i,j) = n - 2
# co(i,j) - ci(i,j) = d_i + d_j - n + 1
# sigma(i,j) = co + ci = (n-2) - P_{ij} - P_{ji} = (n-2) - (d_i + d_j - 1 - 2co) = n-1 - d_i - d_j + 2co
# From co - ci = d_i + d_j - n + 1:
# co = (sigma + d_i + d_j - n + 1) / 2
# ci = (sigma - d_i - d_j + n - 1) / 2

# So sigma = co + ci, and co - ci = d_i + d_j - n + 1.
# Lambda(i,j) = A_{ij}·P_{ji} + A_{ji}·P_{ij}
# If i→j: lambda = P_{ji} = d_j - co(i,j) [since A_{ji}=0]
#                = d_j - (sigma + d_i + d_j - n + 1)/2

# So: lambda(i,j) and sigma(i,j) together with d_i, d_j determine P_{ij} and P_{ji}!
# And (lambda, sigma, scores) determines A² completely.
# And tr(A^5) = tr(A·A^4) = tr(A·(A²)²).

# But does lambda alone determine sigma? Let's check.
print(f"\n  Does lambda determine sigma at n=7?")
lam_only_data = defaultdict(set)
for trial in range(50000):
    bits = np.random.randint(0, 1 << tb)
    A = bits_to_adj(bits, n)
    L = np.zeros((n, n), dtype=int)
    Sig = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v: continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1; L[v][u] += 1
                if A[u][w] and A[v][w]: Sig[u][v] += 1; Sig[v][u] += 1
                if A[w][u] and A[w][v]: Sig[u][v] += 1; Sig[v][u] += 1
    lam = tuple(L[i][j] for i in range(n) for j in range(i+1, n))
    sig = tuple(Sig[i][j] for i in range(n) for j in range(i+1, n))
    lam_only_data[lam].add(sig)

ambig_sig = sum(1 for v in lam_only_data.values() if len(v) > 1)
print(f"  Lambda groups: {len(lam_only_data)}, sigma ambiguous: {ambig_sig}")

print(f"\nDone.")
