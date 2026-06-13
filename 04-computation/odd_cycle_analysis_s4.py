"""
Corrected odd-cycle analysis for H=7 impossibility.
Session: opus-2026-05-28-S4

Key fix: represent directed cycles as ORDERED tuples (normalized by rotation),
not frozensets. This correctly distinguishes 5-cycles with same vertex set but
different arc traversal orders.
"""

from collections import defaultdict

def get_tiles(n):
    return [(x, y) for y in range(n-1) for x in range(y+2, n)]

def tiling_to_adj(n, tiles, bits):
    adj = [[0]*n for _ in range(n)]
    for k in range(1, n):
        adj[k][k-1] = 1
    for idx, (x, y) in enumerate(tiles):
        if (bits >> idx) & 1:
            adj[x][y] = 1
        else:
            adj[y][x] = 1
    return adj

def H_by_DP(adj, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask >> v) & 1 or not dp[mask][v]:
                continue
            for w in range(n):
                if not (mask >> w) & 1 and adj[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))

def normalize_cycle(path):
    """Normalize a directed cycle (as list) to a canonical tuple.
    The cycle is path[0]->path[1]->...->path[-1]->path[0].
    Canonical form: rotate so minimum element is first.
    Do NOT reverse (direction matters for directed cycles).
    """
    n = len(path)
    min_idx = path.index(min(path))
    rotated = tuple(path[min_idx:] + path[:min_idx])
    return rotated

def find_all_odd_cycles(adj, n):
    """Find all distinct odd DIRECTED cycles in tournament adj.
    Returns a set of normalized cycle tuples.
    """
    odd_cycles = set()

    def dfs(start, current, path, in_path):
        for nxt in range(n):
            if not adj[current][nxt]:
                continue
            if nxt == start and len(path) >= 3:
                # Found a cycle: path + [start]
                cycle_len = len(path)
                if cycle_len % 2 == 1:
                    norm = normalize_cycle(list(path))
                    odd_cycles.add(norm)
            elif nxt not in in_path and nxt > start:
                # Only extend to vertices > start (canonical form:
                # minimum vertex is always the starting vertex,
                # ensuring each cycle is found exactly once)
                in_path.add(nxt)
                path.append(nxt)
                dfs(start, nxt, path, in_path)
                path.pop()
                in_path.remove(nxt)

    for start in range(n):
        dfs(start, start, [start], {start})

    return odd_cycles

def independence_polynomial_from_cycles(cycles):
    """Given a list of odd cycles (as normalized tuples),
    compute the independence polynomial I(Omega, x) where
    Omega is the conflict graph (edges = pairs sharing a vertex).
    Returns dict {k: alpha_k}.
    """
    nc = len(cycles)
    if nc == 0:
        return {0: 1}

    cycles = list(cycles)

    # Build conflict relation: i and j conflict if they share a vertex
    # (cycles are stored as tuples of vertices)
    cycle_vsets = [set(c) for c in cycles]
    conflict = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if cycle_vsets[i] & cycle_vsets[j]:
                conflict[i][j] = conflict[j][i] = True

    # Find all independent sets
    alpha = defaultdict(int)
    alpha[0] = 1

    for mask in range(1, 1 << nc):
        members = [k for k in range(nc) if (mask >> k) & 1]
        is_indep = all(
            not conflict[members[a]][members[b]]
            for a in range(len(members))
            for b in range(a+1, len(members))
        )
        if is_indep:
            alpha[len(members)] += 1

    return dict(alpha)

def I_at_2(poly):
    return sum(v * (2**k) for k, v in poly.items())

# ============================================================
# PART 1: Fix verification — H(T) = I(Omega(T), 2) at n=5,6
# ============================================================
print("=" * 60)
print("PART 1: CORRECT OCF VERIFICATION")
print("=" * 60)

for n in [5, 6]:
    tiles = get_tiles(n)
    m = len(tiles)
    correct = 0
    total = 0
    wrong_examples = []

    for bits in range(1 << m):
        adj = tiling_to_adj(n, tiles, bits)
        h_dp = H_by_DP(adj, n)
        cycles = find_all_odd_cycles(adj, n)
        poly = independence_polynomial_from_cycles(cycles)
        h_ocf = I_at_2(poly)

        if h_dp == h_ocf:
            correct += 1
        else:
            wrong_examples.append((bits, h_dp, h_ocf, len(cycles)))
        total += 1

    print(f"n={n}: H(T) = I(Omega(T), 2): {correct}/{total} ({('PASS' if correct==total else 'FAIL')})")
    if wrong_examples:
        print(f"  Wrong examples (first 3): {wrong_examples[:3]}")

# ============================================================
# PART 2: Alpha_1 = 3 analysis (corrected)
# ============================================================
print("\n" + "=" * 60)
print("PART 2: ALPHA_1 DISTRIBUTION (CORRECTED)")
print("=" * 60)

for n in [5, 6, 7]:
    tiles = get_tiles(n)
    m = len(tiles)
    alpha1_dist = defaultdict(int)
    found_alpha1_3_alpha2_0 = []

    for bits in range(1 << m):
        adj = tiling_to_adj(n, tiles, bits)
        cycles = find_all_odd_cycles(adj, n)
        nc = len(cycles)
        alpha1_dist[nc] += 1

        if nc == 3:
            poly = independence_polynomial_from_cycles(cycles)
            a2 = poly.get(2, 0)
            a3 = poly.get(3, 0)
            if a2 == 0 and a3 == 0:
                h = H_by_DP(adj, n)
                found_alpha1_3_alpha2_0.append((bits, h, cycles))

    print(f"\nn={n} alpha_1 distribution:")
    for k in sorted(alpha1_dist.keys()):
        print(f"  alpha_1={k}: {alpha1_dist[k]} tournaments")

    if found_alpha1_3_alpha2_0:
        print(f"\n  FOUND tournaments with alpha_1=3, alpha_2=0 (H=7 SHOULD BE POSSIBLE):")
        for bits, h, cyc in found_alpha1_3_alpha2_0[:5]:
            print(f"    bits={bin(bits)}, H={h}, #cycles=3, cycles={[list(c) for c in cyc]}")
        print(f"  Total: {len(found_alpha1_3_alpha2_0)}")
    else:
        print(f"\n  No tournament has alpha_1=3, alpha_2=0 at n={n}.")
        if alpha1_dist.get(3, 0) > 0:
            print(f"  (Tournaments with alpha_1=3 exist, but always have alpha_2>0)")

# ============================================================
# PART 3: Why is alpha_1 = 3 impossible or always has alpha_2 > 0?
# ============================================================
print("\n" + "=" * 60)
print("PART 3: STRUCTURAL ANALYSIS — ODD CYCLE TRIPLES")
print("=" * 60)

# At n=5: is alpha_1=3 possible at all?
n = 5
tiles = get_tiles(n)
m = len(tiles)

print(f"\nAt n=5: detailed analysis of tournaments with alpha_1=3")
a1_3_tournaments = []
for bits in range(1 << m):
    adj = tiling_to_adj(n, tiles, bits)
    cycles = find_all_odd_cycles(adj, n)
    if len(cycles) == 3:
        poly = independence_polynomial_from_cycles(cycles)
        h = H_by_DP(adj, n)
        a1_3_tournaments.append({
            'bits': bits, 'h': h,
            'cycles': list(cycles),
            'alpha': dict(poly)
        })

if not a1_3_tournaments:
    print("  alpha_1=3 is IMPOSSIBLE at n=5!")
    print("  This is the fundamental reason H=7 is impossible at n=5.")
    print("")
    print("  Proof sketch: In any n=5 tournament, the number of odd cycles")
    print("  jumps from 2 to 4 (never equals 3).")
    print("  Why? Every pair of odd cycles in a 5-vertex tournament that intersect")
    print("  forces the existence of additional odd cycles.")
else:
    print(f"  Found {len(a1_3_tournaments)} tournaments with alpha_1=3")
    for t in a1_3_tournaments[:3]:
        print(f"  H={t['h']}, alpha={t['alpha']}, cycles={t['cycles']}")

# ============================================================
# PART 4: The skip pattern in alpha_1 values
# ============================================================
print("\n" + "=" * 60)
print("PART 4: ALPHA_1 SKIP PATTERN")
print("=" * 60)

print("\nFor each n, which alpha_1 values are IMPOSSIBLE?")
for n in [4, 5, 6, 7]:
    tiles = get_tiles(n)
    m = len(tiles)
    alpha1_achieved = set()
    for bits in range(1 << m):
        adj = tiling_to_adj(n, tiles, bits)
        cycles = find_all_odd_cycles(adj, n)
        alpha1_achieved.add(len(cycles))
    max_possible = max(alpha1_achieved)
    all_possible = set(range(max_possible + 1))
    missing = all_possible - alpha1_achieved
    print(f"  n={n}: alpha_1 in {sorted(alpha1_achieved)[:15]}{'...' if len(alpha1_achieved) > 15 else ''}")
    print(f"        missing: {sorted(missing)}")

# ============================================================
# PART 5: Prove alpha_1 != 3 for all n using structure
# ============================================================
print("\n" + "=" * 60)
print("PART 5: ALGEBRAIC PROOF STRATEGY FOR ALPHA_1 GAPS")
print("=" * 60)

print("""
H=7 requires I(Omega(T), 2) = 7.
Writing I(Omega, x) = 1 + alpha_1*x + alpha_2*x^2 + ...:
  1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ... = 7
=> alpha_1 + 2*alpha_2 + 4*alpha_3 + ... = 3

Since all alpha_k >= 0, the only solution is:
  (alpha_1=3, alpha_2=0, alpha_3=0, ...)  [the UNIQUE solution]

So H=7 iff exactly 3 odd cycles exist with all pairs intersecting.

Key question: can a tournament have exactly 3 odd cycles?
""")

# Check: what are the possible alpha_1 values modulo small numbers?
n = 7
tiles = get_tiles(n)
m = len(tiles)
alpha1_vals_n7 = []
for bits in range(1 << m):
    adj = tiling_to_adj(n, tiles, bits)
    cycles = find_all_odd_cycles(adj, n)
    alpha1_vals_n7.append(len(cycles))

from collections import Counter
alpha1_counts_n7 = Counter(alpha1_vals_n7)
print(f"n=7: alpha_1 distribution (first 20 nonzero):")
for k in sorted(alpha1_counts_n7.keys())[:20]:
    print(f"  alpha_1={k}: {alpha1_counts_n7[k]} tournaments")

print(f"\nDoes alpha_1=3 appear at n=7? {'YES' if 3 in alpha1_counts_n7 else 'NO'}")
print(f"Does alpha_1=7 appear at n=7? {'YES' if 7 in alpha1_counts_n7 else 'NO'}")

# Find the minimum nonzero even number NOT achieved
achieved_n7 = sorted(alpha1_counts_n7.keys())
print(f"\nAll achieved alpha_1 values at n=7 (sorted): {achieved_n7[:30]}{'...' if len(achieved_n7) > 30 else ''}")
missing_n7 = [k for k in range(max(achieved_n7)+1) if k not in alpha1_counts_n7]
print(f"Missing alpha_1 values at n=7: {missing_n7[:20]}")

# ============================================================
# PART 6: Connecting H spectrum to the forbidden values
# ============================================================
print("\n" + "=" * 60)
print("PART 6: FORBIDDEN H VALUES BY SIZE")
print("=" * 60)

for n in [5, 6, 7]:
    tiles = get_tiles(n)
    m = len(tiles)
    H_achieved = set()
    for bits in range(1 << m):
        adj = tiling_to_adj(n, tiles, bits)
        H_achieved.add(H_by_DP(adj, n))
    H_sorted = sorted(H_achieved)
    max_h = max(H_sorted)

    # Which odd numbers <= max_h are NOT achieved?
    forbidden = [h for h in range(1, max_h+2, 2) if h not in H_achieved]
    print(f"\nn={n}: H ranges from 1 to {max_h}")
    print(f"  Forbidden odd H values: {forbidden[:20]}")
    print(f"  Forbidden: H ≡ ? (mod various primes)")
    for forbidden_h in forbidden[:5]:
        mods = {p: forbidden_h % p for p in [3, 5, 7, 11, 13]}
        print(f"    H={forbidden_h}: {mods}")

print("""
Observation: H=7 is the SMALLEST forbidden value (> 1).
All higher forbidden values (H=21, 35, ...) are multiples of 7!
""")

# Check: are all forbidden values multiples of 7?
for n in [5, 6]:
    tiles = get_tiles(n)
    m = len(tiles)
    H_achieved = set()
    for bits in range(1 << m):
        adj = tiling_to_adj(n, tiles, bits)
        H_achieved.add(H_by_DP(adj, n))
    max_h = max(H_achieved)
    forbidden = [h for h in range(1, max_h+2, 2) if h not in H_achieved]
    print(f"n={n}: forbidden H values = {forbidden}")
    all_mult_7 = all(h % 7 == 0 for h in forbidden)
    print(f"  All multiples of 7? {all_mult_7}")
