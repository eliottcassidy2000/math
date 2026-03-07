"""
H21_gap_search.py -- Investigate whether H=21 is achievable at n=8,9,10,11,12.

H(T) = number of Hamiltonian paths in tournament T.
Uses OCF: H(T) = I(Omega(T), 2) where Omega is the conflict graph on directed
odd cycles and I is the independence polynomial.

Strategy:
1. Random sampling at each n to collect achieved H values
2. Report which odd values <= 25 are missing
3. Targeted construction attempts for H=21
4. Analysis of why H=21 might be impossible
"""

import random
import sys
from itertools import combinations
from collections import defaultdict

def count_H(A, n):
    """Count Hamiltonian paths in tournament with adjacency matrix A."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            c = dp.get((mask, v), 0)
            if c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def random_tournament(n):
    """Generate a random tournament on n vertices."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def transitive_tournament(n):
    """Generate the transitive tournament: i beats j iff i < j."""
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    return A

def flip_edge(A, n, i, j):
    """Flip the edge between i and j. Returns a new matrix."""
    B = [row[:] for row in A]
    B[i][j], B[j][i] = B[j][i], B[i][j]
    return B

def find_directed_3cycles(A, n):
    """Find all directed 3-cycles in tournament A."""
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check all 2 orientations of the 3-cycle on {i,j,k}
                if A[i][j] and A[j][k] and A[k][i]:
                    cycles.append((i,j,k))
                elif A[i][k] and A[k][j] and A[j][i]:
                    cycles.append((i,k,j))
    return cycles

def find_odd_cycles_up_to(A, n, max_len=None):
    """Find all directed odd cycles (as vertex sets) in tournament A.
    Returns list of frozensets of vertices."""
    if max_len is None:
        max_len = n
    cycle_sets = set()

    # For each subset of odd size 3,5,7,...
    for size in range(3, max_len+1, 2):
        for subset in combinations(range(n), size):
            # Check if there's a directed Hamiltonian cycle on this subset
            sub_nodes = list(subset)
            # Build sub-adjacency
            idx = {v: i for i, v in enumerate(sub_nodes)}
            sub_A = [[0]*size for _ in range(size)]
            for a in sub_nodes:
                for b in sub_nodes:
                    if a != b and A[a][b]:
                        sub_A[idx[a]][idx[b]] = 1
            # Count Hamiltonian cycles using DP
            # dp[(mask, v)] = number of Hamiltonian paths from node 0 to v using mask
            s = size
            dp = {}
            dp[(1, 0)] = 1  # start at node 0
            for mask in range(1, 1 << s):
                if not (mask & 1):
                    continue  # must include node 0
                for v in range(s):
                    if not (mask & (1 << v)):
                        continue
                    c = dp.get((mask, v), 0)
                    if c == 0:
                        continue
                    for u in range(s):
                        if mask & (1 << u):
                            continue
                        if sub_A[v][u]:
                            dp[(mask | (1 << u), u)] = dp.get((mask | (1 << u), u), 0) + c
            full = (1 << s) - 1
            # Check if any path from 0 to v closes back to 0
            has_cycle = False
            for v in range(1, s):
                if dp.get((full, v), 0) > 0 and sub_A[v][0]:
                    has_cycle = True
                    break
            if has_cycle:
                cycle_sets.add(frozenset(subset))
    return cycle_sets

def build_conflict_graph(cycle_sets):
    """Build conflict graph: vertices are cycle sets, edges between non-disjoint ones."""
    cycles = list(cycle_sets)
    m = len(cycles)
    adj = [[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1, m):
            if cycles[i] & cycles[j]:  # share a vertex
                adj[i][j] = True
                adj[j][i] = True
    return cycles, adj

def independence_polynomial_at_2(cycles_list, adj, max_size=None):
    """Compute I(Omega, 2) = sum over independent sets S of 2^|S|."""
    m = len(cycles_list)
    if max_size is None:
        max_size = m

    # For small m, enumerate all independent sets
    if m > 25:
        # Too large, use approximation or skip
        return None

    total = 0
    for mask in range(1 << m):
        # Check if mask is an independent set
        is_indep = True
        nodes_in = []
        for i in range(m):
            if mask & (1 << i):
                nodes_in.append(i)
        for idx_a in range(len(nodes_in)):
            if not is_indep:
                break
            for idx_b in range(idx_a+1, len(nodes_in)):
                if adj[nodes_in[idx_a]][nodes_in[idx_b]]:
                    is_indep = False
                    break
        if is_indep:
            total += (1 << len(nodes_in))  # 2^|S|
    return total

def analyze_tournament_ocf(A, n):
    """Compute H via OCF decomposition: find cycle structure."""
    cycle_sets = find_odd_cycles_up_to(A, n)
    cycles_list, adj = build_conflict_graph(cycle_sets)

    # Count 3-cycles
    c3 = sum(1 for cs in cycle_sets if len(cs) == 3)
    # Count 5-cycles
    c5 = sum(1 for cs in cycle_sets if len(cs) == 5)

    # alpha_1 = total odd cycle vertex-sets
    alpha_1 = len(cycle_sets)

    # alpha_2 = number of disjoint pairs (independent edges in conflict graph complement)
    # = pairs that are disjoint (no shared vertex)
    alpha_2 = 0
    for i in range(len(cycles_list)):
        for j in range(i+1, len(cycles_list)):
            if not adj[i][j]:  # disjoint
                alpha_2 += 1

    ip = independence_polynomial_at_2(cycles_list, adj)
    return {
        'c3': c3, 'c5': c5, 'alpha_1': alpha_1, 'alpha_2': alpha_2,
        'I_Omega_2': ip, 'num_cycles': len(cycle_sets)
    }

# ============================================================
# PART 1: Random sampling
# ============================================================
print("=" * 70)
print("PART 1: Random sampling to find achieved H values")
print("=" * 70)

target_values = set(range(1, 28, 2))  # odd values 1,3,...,27

for n in range(8, 13):
    num_samples = 10000
    if n >= 11:
        num_samples = 5000  # slower for larger n
    if n >= 12:
        num_samples = 2000

    achieved = set()
    h_counts = defaultdict(int)
    h21_found = False

    print(f"\nn={n}: Sampling {num_samples} random tournaments...")
    sys.stdout.flush()

    for trial in range(num_samples):
        A = random_tournament(n)
        h = count_H(A, n)
        achieved.add(h)
        if h <= 30:
            h_counts[h] += 1
        if h == 21:
            h21_found = True
            # Save the tournament
            print(f"  *** H=21 FOUND at n={n}, trial {trial}! ***")
            print(f"  Adjacency: {[A[i] for i in range(n)]}")
            sys.stdout.flush()

    # Also add transitive tournament (H=1)
    A_trans = transitive_tournament(n)
    h_trans = count_H(A_trans, n)
    achieved.add(h_trans)

    small_achieved = sorted([h for h in achieved if h <= 30])
    missing_odd = sorted([h for h in target_values if h not in achieved])

    print(f"  Achieved H values <= 30: {small_achieved}")
    print(f"  Missing odd values <= 27: {missing_odd}")
    print(f"  H counts for small values: {dict(sorted(h_counts.items()))}")
    if not h21_found:
        print(f"  H=21 NOT found in {num_samples} samples")
    sys.stdout.flush()

# ============================================================
# PART 2: Targeted construction for H=21
# ============================================================
print("\n" + "=" * 70)
print("PART 2: Targeted edge-flip construction for H=21")
print("=" * 70)

def try_targeted_construction(n, target_h=21, max_attempts=5000):
    """Try to construct a tournament with H=target_h by flipping edges."""
    best_diff = float('inf')
    best_h = None
    found = False

    for attempt in range(max_attempts):
        # Start with transitive tournament
        A = transitive_tournament(n)
        h = 1  # transitive has H=1

        # Randomly flip edges and track H
        edges = [(i,j) for i in range(n) for j in range(i+1, n)]
        random.shuffle(edges)

        for i, j in edges:
            B = flip_edge(A, n, i, j)
            new_h = count_H(B, n)

            if new_h == target_h:
                print(f"  FOUND H={target_h} at n={n}!")
                print(f"  Adjacency: {[B[r] for r in range(n)]}")
                return B, True

            if abs(new_h - target_h) < abs(h - target_h):
                A = B
                h = new_h
            elif random.random() < 0.1:  # sometimes accept worse
                A = B
                h = new_h

            if abs(h - target_h) < best_diff:
                best_diff = abs(h - target_h)
                best_h = h

    print(f"  n={n}: Best H found near 21: {best_h} (diff={best_diff})")
    return None, False

for n in [8, 9, 10]:
    print(f"\nTargeted search at n={n}...")
    sys.stdout.flush()
    try_targeted_construction(n, target_h=21, max_attempts=2000)
    sys.stdout.flush()

# ============================================================
# PART 3: Exhaustive search at n=8 for H=21
# ============================================================
print("\n" + "=" * 70)
print("PART 3: Near-exhaustive search at n=8 for small H values")
print("=" * 70)

# At n=8, there are 2^28 tournaments. We can't enumerate all,
# but we can do a large sample and also systematically explore
# tournaments with few 3-cycles (which tend to have small H).

# Strategy: start from transitive, flip 1,2,3,... edges systematically
print("\nn=8: Systematic edge-flip from transitive tournament")
sys.stdout.flush()

n = 8
A_trans = transitive_tournament(n)
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
num_edges = len(edges)  # 28

achieved_h = set()
achieved_h.add(1)

# Flip 1 edge
print("  Flipping 1 edge...")
for e1 in range(num_edges):
    i, j = edges[e1]
    A = flip_edge(A_trans, n, i, j)
    h = count_H(A, n)
    achieved_h.add(h)

print(f"  After 1-flip: achieved H <= 30: {sorted([h for h in achieved_h if h <= 30])}")
sys.stdout.flush()

# Flip 2 edges
print("  Flipping 2 edges...")
for e1 in range(num_edges):
    A1 = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
    for e2 in range(e1+1, num_edges):
        A2 = flip_edge(A1, n, edges[e2][0], edges[e2][1])
        h = count_H(A2, n)
        achieved_h.add(h)

print(f"  After 2-flip: achieved H <= 30: {sorted([h for h in achieved_h if h <= 30])}")
sys.stdout.flush()

# Flip 3 edges (this is C(28,3)=3276, feasible)
print("  Flipping 3 edges...")
for e1 in range(num_edges):
    A1 = flip_edge(A_trans, n, edges[e1][0], edges[e1][1])
    for e2 in range(e1+1, num_edges):
        A2 = flip_edge(A1, n, edges[e2][0], edges[e2][1])
        for e3 in range(e2+1, num_edges):
            A3 = flip_edge(A2, n, edges[e3][0], edges[e3][1])
            h = count_H(A3, n)
            achieved_h.add(h)

small_h = sorted([h for h in achieved_h if h <= 30])
print(f"  After 3-flip: achieved H <= 30: {small_h}")
missing = sorted([h for h in range(1, 28, 2) if h not in achieved_h])
print(f"  Missing odd values <= 27: {missing}")
sys.stdout.flush()

# ============================================================
# PART 4: OCF analysis of near-21 tournaments
# ============================================================
print("\n" + "=" * 70)
print("PART 4: OCF analysis -- what cycle structures give H near 21?")
print("=" * 70)

# For small n where we can compute the full OCF
n = 7
print(f"\nn={n}: Sampling tournaments with H near 21 and analyzing cycle structure...")
sys.stdout.flush()

near_21 = []
for trial in range(20000):
    A = random_tournament(n)
    h = count_H(A, n)
    if 15 <= h <= 27:
        near_21.append((h, A))

# Group by H value
by_h = defaultdict(list)
for h, A in near_21:
    by_h[h].append(A)

print(f"  H value distribution near 21 at n={n}:")
for h in sorted(by_h.keys()):
    count = len(by_h[h])
    print(f"    H={h}: {count} tournaments found")

# Analyze cycle structure for H values around 21
for target in [19, 21, 23]:
    if target in by_h and len(by_h[target]) > 0:
        A = by_h[target][0]
        info = analyze_tournament_ocf(A, n)
        print(f"\n  OCF analysis for H={target} tournament at n={n}:")
        print(f"    c3={info['c3']}, c5={info['c5']}, alpha_1={info['alpha_1']}")
        print(f"    alpha_2={info['alpha_2']}, I(Omega,2)={info['I_Omega_2']}")
    else:
        print(f"\n  H={target}: NOT FOUND at n={n}")
sys.stdout.flush()

# ============================================================
# PART 5: Algebraic constraint analysis
# ============================================================
print("\n" + "=" * 70)
print("PART 5: Algebraic constraint analysis for H=21")
print("=" * 70)

print("""
OCF: H(T) = I(Omega(T), 2) = sum_{k>=0} i_k * 2^k

where i_k = number of independent sets of size k in Omega(T).

H = 21 = 10101 in binary = 1 + 4 + 16 = 2^0 + 2^2 + 2^4

So we need: i_0=1 (always), and the sum i_1*2 + i_2*4 + i_3*8 + ... = 20

Possible decompositions of 20 = sum i_k * 2^k (k>=1):
  i_1=10, rest=0: 10*2 = 20  => H = 1 + 20 = 21
  i_1=8, i_2=1: 16+4 = 20    => H = 1 + 16 + 4 = 21
  i_1=6, i_2=2: 12+8 = 20    => H = 1 + 12 + 8 = 21
  i_1=4, i_2=3: 8+12 = 20    => H = 1 + 8 + 12 = 21
  i_1=2, i_2=4: 4+16 = 20    => H = 1 + 4 + 16 = 21
  i_1=0, i_2=5: 0+20 = 20    => H = 1 + 0 + 20 = 21
  i_1=8, i_3=0.5: impossible (i_3 must be integer)
  i_1=6, i_2=1, i_3=0.25: impossible
  i_1=4, i_2=1, i_3=1: 8+4+8=20  => H = 1+8+4+8 = 21
  i_1=2, i_2=1, i_3=2: 4+4+16=24 != 20
  etc.

Key question: Can any tournament have alpha_1 (= i_1) = 10 with all cycles
pairwise conflicting (i_2 = 0)?  That gives H = 1 + 20 = 21.

Or: alpha_1=8, i_2=1 (exactly one disjoint pair)?  H = 1+16+4 = 21.
""")
sys.stdout.flush()

# ============================================================
# PART 6: Check if H=21 is achievable by analyzing i_k decomposition
# ============================================================
print("=" * 70)
print("PART 6: Verify OCF decomposition for achieved H values near 21")
print("=" * 70)

# At n=7, do detailed OCF for all achieved H values
n = 7
h_to_ocf = {}
random.seed(42)

print(f"\nn={n}: Computing OCF decomposition for tournaments with various H values...")
sys.stdout.flush()

count = 0
for trial in range(30000):
    A = random_tournament(n)
    h = count_H(A, n)
    if h not in h_to_ocf and h <= 50:
        info = analyze_tournament_ocf(A, n)
        h_to_ocf[h] = info
        count += 1
        if count % 5 == 0:
            sys.stdout.flush()

print(f"\n  OCF decompositions at n={n}:")
for h in sorted(h_to_ocf.keys()):
    info = h_to_ocf[h]
    print(f"    H={h}: alpha_1={info['alpha_1']}, alpha_2={info['alpha_2']}, "
          f"c3={info['c3']}, c5={info['c5']}, I(Omega,2)={info['I_Omega_2']}")

# Check which H values are achieved
all_h_n7 = set()
random.seed(123)
for trial in range(50000):
    A = random_tournament(7)
    h = count_H(A, 7)
    all_h_n7.add(h)

print(f"\n  All achieved H values at n=7 (up to 50): {sorted([h for h in all_h_n7 if h <= 50])}")
print(f"  Missing odd values <= 27 at n=7: {sorted([h for h in range(1, 28, 2) if h not in all_h_n7])}")
sys.stdout.flush()

# ============================================================
# PART 7: Summary
# ============================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
The key constraint for H=21:
  H = I(Omega(T), 2) = 1 + 2*i_1 + 4*i_2 + 8*i_3 + ...

  21 = 1 + 20, so we need sum_{k>=1} i_k * 2^k = 20.

  Since i_k >= 0, the simplest decomposition is i_1=10, i_k=0 for k>=2.
  This means 10 odd-cycle vertex sets, all pairwise sharing a vertex.

  The question is whether the combinatorial constraints on tournaments
  allow exactly these independence polynomial values.
""")
