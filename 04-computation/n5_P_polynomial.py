"""
Exhaustive computation of P(u, x) at n=5.

G_T(t) = sum_{sigma in S_n} t^{fwd(sigma)}
where fwd(sigma) = # forward edges in Hamiltonian path sigma(0)->sigma(1)->...->sigma(n-1).

This is a polynomial of degree n-1=4 in t, palindromic.

Relation to P(u): G_T(t) = t^m * P(u) where u = t + 1/t, m = (n-1)/2 = 2.

Let's first just compute everything and examine the structure, without
assuming p_2 = I(Omega, 2). We'll check that relationship empirically.
"""

from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def tournament_from_bits(bits, n=5):
    """Create adjacency matrix from bit encoding."""
    adj = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            k += 1
    return adj

def count_3cycles(adj, n=5):
    """Count directed 3-cycles."""
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count

def count_5cycles(adj, n=5):
    """Count directed 5-cycles (Hamiltonian cycles)."""
    count = 0
    # Fix vertex 0 to avoid counting rotations
    for perm in permutations(range(1, 5)):
        cycle = [0] + list(perm)
        is_cycle = True
        for idx in range(5):
            if not adj[cycle[idx]][cycle[(idx+1)%5]]:
                is_cycle = False
                break
        if is_cycle:
            count += 1
    # Each undirected cycle counted twice (two directions), each directed once
    # With vertex 0 fixed, we get each directed cycle exactly once
    return count

def independence_poly(adj, n=5):
    """Compute I(Omega(T), x) coefficients for conflict graph."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles.append(frozenset({i,j,k}))
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles.append(frozenset({i,j,k}))
    c3 = len(cycles)
    # At n=5, no two 3-cycles can be vertex-disjoint
    # But they CAN be non-conflicting if they share no vertex? No - 3+3=6>5
    # So alpha_2 = 0 always at n=5
    return [1, c3]

def forward_edge_distribution(adj, n=5):
    """Compute a_k = number of permutations with exactly k forward edges."""
    dist = [0] * n
    for perm in permutations(range(n)):
        fwd = 0
        for idx in range(n-1):
            if adj[perm[idx]][perm[idx+1]]:
                fwd += 1
        dist[fwd] += 1
    return dist

# Main computation
n = 5
num_edges = n*(n-1)//2  # = 10
num_tournaments = 2**num_edges  # = 1024

print(f"Enumerating all {num_tournaments} tournaments on n={n} vertices")
print()

# Group by (c3, c5) -> forward edge distribution and counts
results = defaultdict(lambda: {'count': 0, 'p_values': None, 'a_dist': None, 'ip': None, 'bits_example': None})

for bits in range(num_tournaments):
    adj = tournament_from_bits(bits, n)
    c3 = count_3cycles(adj, n)
    c5 = count_5cycles(adj, n)
    a = forward_edge_distribution(adj, n)

    # Verify palindromic
    assert a[0] == a[4], f"Not palindromic: a={a}"
    assert a[1] == a[3], f"Not palindromic: a={a}"

    # Extract P coefficients
    p2 = a[0]
    p1 = a[1]
    p0 = a[2] - 2*a[0]

    # Verify total = 120
    total = sum(a)
    assert total == 120

    # Verify P(2) = 120
    assert p0 + 2*p1 + 4*p2 == 120

    ip = independence_poly(adj, n)
    H_T = ip[0] + 2*ip[1]  # = 1 + 2*c3

    key = (c3, c5)
    if results[key]['p_values'] is None:
        results[key]['p_values'] = (p0, p1, p2)
        results[key]['a_dist'] = tuple(a)
        results[key]['ip'] = tuple(ip)
        results[key]['bits_example'] = bits
    else:
        if results[key]['p_values'] != (p0, p1, p2):
            print(f"WARNING: (c3,c5)={key} has DIFFERENT P values!")
            print(f"  Previous: {results[key]['p_values']}")
            print(f"  Current:  {(p0, p1, p2)}, bits={bits}")
    results[key]['count'] += 1

print("=" * 90)
print("FULL TABLE: (c3, c5) -> forward-edge distribution -> P-coefficients")
print("=" * 90)
header = f"{'c3':>3} {'c5':>3} | {'a0':>4} {'a1':>4} {'a2':>4} {'a3':>4} {'a4':>4} | {'p0':>5} {'p1':>5} {'p2':>5} | {'H=I(O,2)':>9} | {'#tours':>6}"
print(header)
print("-" * 90)

for key in sorted(results.keys()):
    c3, c5 = key
    r = results[key]
    p0, p1, p2 = r['p_values']
    a = r['a_dist']
    H_T = 1 + 2*c3
    print(f"{c3:>3} {c5:>3} | {a[0]:>4} {a[1]:>4} {a[2]:>4} {a[3]:>4} {a[4]:>4} | {p0:>5} {p1:>5} {p2:>5} | {H_T:>9} | {r['count']:>6}")

print()
print("=" * 90)
print("ANALYSIS: p_j as functions of c3 and c5")
print("=" * 90)

data_points = [(key[0], key[1], results[key]['p_values']) for key in sorted(results.keys())]

print("\nRaw data:")
for c3, c5, (p0, p1, p2) in data_points:
    print(f"  c3={c3}, c5={c5}: p0={p0}, p1={p1}, p2={p2}")

c3_arr = np.array([d[0] for d in data_points], dtype=float)
c5_arr = np.array([d[1] for d in data_points], dtype=float)
p0_arr = np.array([d[2][0] for d in data_points], dtype=float)
p1_arr = np.array([d[2][1] for d in data_points], dtype=float)
p2_arr = np.array([d[2][2] for d in data_points], dtype=float)

# Linear fit: p_j = alpha + beta*c3 + gamma*c5
A_lin = np.column_stack([np.ones_like(c3_arr), c3_arr, c5_arr])

print("\nLinear fits (p_j = alpha + beta*c3 + gamma*c5):")
for name, vals in [("p0", p0_arr), ("p1", p1_arr), ("p2", p2_arr)]:
    coeffs, residuals, rank, sv = np.linalg.lstsq(A_lin, vals, rcond=None)
    resid = np.sum((A_lin @ coeffs - vals)**2)
    print(f"  {name} = {coeffs[0]:.4f} + {coeffs[1]:.4f}*c3 + {coeffs[2]:.4f}*c5  (residual={resid:.6f})")

# Since linear in c3, c5 might not be exact, try with more terms
A_quad = np.column_stack([np.ones_like(c3_arr), c3_arr, c5_arr, c3_arr**2, c5_arr**2, c3_arr*c5_arr])
print("\nQuadratic fits:")
for name, vals in [("p0", p0_arr), ("p1", p1_arr), ("p2", p2_arr)]:
    coeffs, residuals, rank, sv = np.linalg.lstsq(A_quad, vals, rcond=None)
    resid = np.sum((A_quad @ coeffs - vals)**2)
    terms = ["1", "c3", "c5", "c3^2", "c5^2", "c3*c5"]
    expr = " + ".join(f"{c:.4f}*{t}" for c, t in zip(coeffs, terms) if abs(c) > 0.0001)
    print(f"  {name} = {expr}  (residual={resid:.6f})")

print()
print("=" * 90)
print("RELATIONSHIP BETWEEN p2 AND I(Omega(T), 2)")
print("=" * 90)

for key in sorted(results.keys()):
    c3, c5 = key
    r = results[key]
    p0, p1, p2 = r['p_values']
    H_T = 1 + 2*c3
    print(f"  c3={c3}, c5={c5}: p2={p2}, I(Omega,2)=1+2*{c3}={H_T}, ratio p2/H_T={p2/H_T:.4f}, diff={p2-H_T}")

print()
print("=" * 90)
print("ROOTS OF P(u, 2) — QUADRATIC IN u")
print("=" * 90)

print(f"\nP(u) = p2*u^2 + p1*u + p0")
print(f"Roots: u = (-p1 +/- sqrt(p1^2 - 4*p0*p2)) / (2*p2)")
print()

for key in sorted(results.keys()):
    c3, c5 = key
    p0, p1, p2 = results[key]['p_values']

    disc = p1**2 - 4*p0*p2
    if disc >= 0:
        sqrt_d = np.sqrt(disc)
        r1 = (-p1 + sqrt_d) / (2*p2)
        r2 = (-p1 - sqrt_d) / (2*p2)
        root_str = f"u = {r1:.6f}, {r2:.6f}"
        real_str = "REAL"
    else:
        re = -p1 / (2*p2)
        im = np.sqrt(-disc) / (2*p2)
        root_str = f"u = {re:.6f} +/- {im:.6f}i"
        real_str = "COMPLEX"

    product = p0/p2
    sum_r = -p1/p2
    print(f"c3={c3}, c5={c5}: disc={disc:>8}, {real_str:>7}, {root_str}")
    print(f"  sum=-p1/p2={sum_r:.4f}, product=p0/p2={product:.4f}")

print()
print("=" * 90)
print("P(0, 2) = p_0 ANALYSIS (t = i)")
print("=" * 90)

print(f"\nWhen u=0: t+1/t=0 => t=i.")
print(f"G_T(i) = sum a_k * i^k = a0 - a2 + a4 + (a1-a3)i = 2a0 - a2 + 0i")
print(f"G_T(i) = i^2 * P(0) = -P(0) = -p0")
print(f"So p0 = a2 - 2a0 = -(2a0 - a2) = -G_T(i)")
print()

for key in sorted(results.keys()):
    c3, c5 = key
    p0, p1, p2 = results[key]['p_values']
    a = results[key]['a_dist']
    G_i = 2*a[0] - a[2]
    print(f"  c3={c3}, c5={c5}: p0={p0:>5}, G_T(i)={G_i:>5}, -G_T(i)={-G_i:>5}, match={p0==-G_i}")

print()
print("=" * 90)
print("CHECKING: Is c5 determined by c3 at n=5?")
print("=" * 90)

c3_to_c5 = defaultdict(set)
for key in results:
    c3_to_c5[key[0]].add(key[1])
for c3 in sorted(c3_to_c5):
    print(f"  c3={c3}: c5 in {sorted(c3_to_c5[c3])}")

print()
print("=" * 90)
print("SCORE SEQUENCES AND CYCLE COUNTS")
print("=" * 90)

score_to_cycles = defaultdict(lambda: defaultdict(int))
for bits in range(num_tournaments):
    adj = tournament_from_bits(bits, n)
    scores = tuple(sorted([sum(adj[i]) for i in range(n)]))
    c3 = count_3cycles(adj, n)
    c5 = count_5cycles(adj, n)
    score_to_cycles[scores][(c3, c5)] += 1

for scores in sorted(score_to_cycles):
    print(f"  Score {scores}: {dict(score_to_cycles[scores])}")

print()
print("=" * 90)
print("DEEPER: What does p2 actually count?")
print("=" * 90)

print("""
p2 = a_0 = number of permutations with 0 forward edges.
A path with 0 forward edges means EVERY consecutive pair goes backward:
  sigma(1)->sigma(0), sigma(2)->sigma(1), ..., sigma(4)->sigma(3).
This means sigma(0)<-sigma(1)<-sigma(2)<-sigma(3)<-sigma(4) is a Hamiltonian path
in the REVERSE tournament T^op.
So a_0(T) = a_4(T) = H(T^op) restricted to fully-forward paths in T^op.
Wait - a_0(T) = # permutations where ALL edges go backward = # Hamiltonian paths in T^op
that go fully forward. Actually:

A permutation (v0,v1,v2,v3,v4) with 0 forward edges means v1->v0, v2->v1, v3->v2, v4->v3
in T, i.e., v4->v3->v2->v1->v0 is a Hamiltonian PATH in T (all forward).
So a_0 = # Hamiltonian paths in T. But H(T) is the total number of Hamiltonian paths...
No wait, H(T) = sum of all a_k? No, H(T) = I(Omega(T), 2).

Actually, H(T) as the number of Hamiltonian PATHS = total # permutations where
the sequence forms a directed path = sum over all orderings of those that are
Hamiltonian paths. But ANY permutation of vertices defines a sequence, and
we're counting forward edges in that sequence.

a_k = # permutations (v0,...,v4) with exactly k forward edges.
sum a_k = 5! = 120.

The number of Hamiltonian PATHS in T = a_4 (all forward) + those reversed = a_0.
Wait no. a_4 = # permutations where all 4 edges go forward = # Hamiltonian paths in T.
And a_0 = # permutations where all edges go backward = # Hamiltonian paths in T
(just read the permutation backwards).

Since reversing a permutation with 0 forward edges gives one with 4 forward edges,
a_0 = a_4 = number of Hamiltonian paths in T.

So: H(T) = a_0 = a_4 = p_2.

Wait but H(T) = I(Omega(T), 2) and we saw p2 != I(Omega,2). Let me recheck.
""")

# Let's directly count Hamiltonian paths for a specific tournament
bits = 8  # The one that failed
adj = tournament_from_bits(bits, n)
a = forward_edge_distribution(adj, n)
c3 = count_3cycles(adj, n)
ip = independence_poly(adj, n)

print(f"Tournament bits=8:")
print(f"  Adjacency matrix:")
for i in range(n):
    print(f"    {adj[i]}")
print(f"  c3={c3}, I(Omega,2)=1+2*{c3}={1+2*c3}")
print(f"  a = {a}")
print(f"  a[0] = a[4] = {a[0]} (# Hamiltonian paths in T)")

# Count Hamiltonian paths directly
hp_count = 0
for perm in permutations(range(n)):
    is_path = True
    for idx in range(n-1):
        if not adj[perm[idx]][perm[idx+1]]:
            is_path = False
            break
    if is_path:
        hp_count += 1

print(f"  Direct Hamiltonian path count: {hp_count}")
print(f"  So a[4] = hp_count = {hp_count}, and a[0] should also = {a[0]}")
print(f"  H(T) from OCF = I(Omega,2) = {1+2*c3}")
print()
print(f"  => p_2 = a_0 = H(T) (Hamiltonian path count)")
print(f"  => H(T) via direct counting = {hp_count}")
print(f"  => H(T) via OCF (I(Omega,2)) = {1+2*c3}")
print(f"  => These {'AGREE' if hp_count == 1+2*c3 else 'DISAGREE!'}")

if hp_count != 1 + 2*c3:
    print(f"\n  DISCREPANCY! Let me check c3 counting...")
    # Recount 3-cycles carefully
    cycles_list = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    cycles_list.append((i,j,k, "i->j->k->i"))
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    cycles_list.append((i,j,k, "i->k->j->i"))
    print(f"  3-cycles found: {len(cycles_list)}")
    for c in cycles_list:
        print(f"    {c}")

    # Let me also check: what are the Hamiltonian paths?
    print(f"\n  Hamiltonian paths (all-forward):")
    count = 0
    for perm in permutations(range(n)):
        is_path = True
        for idx in range(n-1):
            if not adj[perm[idx]][perm[idx+1]]:
                is_path = False
                break
        if is_path:
            count += 1
            if count <= 20:
                print(f"    {perm}")

print()
print("=" * 90)
print("RE-EXAMINING: H(T) = I(Omega(T), 2) vs direct path count")
print("=" * 90)

# Check for ALL distinct (c3,c5) classes
mismatches = 0
for key in sorted(results.keys()):
    c3, c5 = key
    p2 = results[key]['p_values'][2]
    H_ocf = 1 + 2*c3
    bits_ex = results[key]['bits_example']

    # Directly count HP for the example tournament
    adj = tournament_from_bits(bits_ex, n)
    hp = 0
    for perm in permutations(range(n)):
        ok = all(adj[perm[idx]][perm[idx+1]] for idx in range(n-1))
        if ok:
            hp += 1

    match_p2 = (hp == p2)
    match_ocf = (hp == H_ocf)

    if not match_p2 or not match_ocf:
        mismatches += 1
    print(f"  c3={c3}, c5={c5}: H(direct)={hp}, p2=a0={p2}, I(Omega,2)={H_ocf}, "
          f"p2=H:{match_p2}, OCF=H:{match_ocf}")

print(f"\nMismatches: {mismatches}")
