#!/usr/bin/env python3
"""
recursive_tiling_s20cv.py — The recursive decomposition of the staircase
kind-pasteur-2026-03-23-S20cv

THE KEY RECURSION:
  delta_{n-2} = bottom ∪ top ∪ {apex}
  where:
    bottom = {(r,c) : r+c <= n-2} = Grid(n-1)  [arcs NOT involving vertex n-1]
    top    = {(r,c) : c >= 2}     ≅ Grid(n-1)  [arcs NOT involving vertex 0]
    apex   = (n-2, 1)                            [THE arc between vertices 0 and n-1]
    bottom ∩ top = {(r,c) : c>=2 and r+c<=n-2} ≅ Grid(n-2)

  In tournament terms:
    T on n vertices = T\(n-1) ∪ T\0 ∪ {arc(0, n-1)}
    The overlap T\(n-1) ∩ T\0 = T\{0, n-1} on n-2 vertices

  The apex arc determines whether there's a Hamiltonian CYCLE through P_0.

CELL COUNT: m(n) = 2*m(n-1) - m(n-2) + 1
  C(n-1,2) = 2*C(n-2,2) - C(n-3,2) + 1  ✓

RECURSIVE META-GRAPH STRUCTURE:
  How do the meta-graph edges of G_n/Z_2 decompose across bottom, top, overlap, and apex?
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  RECURSIVE TILING DECOMPOSITION")
print("  delta_{n-2} = bottom(n-1) U top(n-1) U {apex}")
print("  kind-pasteur-2026-03-23-S20cv")
print("=" * 80)

# ============================================================================
# HELPERS (compact)
# ============================================================================
def tadj(n, bits, pairs):
    a = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(pairs):
        if bits & (1 << k): a[i][j] = 1
        else: a[j][i] = 1
    return a

def canon(a, n):
    sc = [sum(a[i][j] for j in range(n)) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(set(sc))]
    if all(len(g)==1 for g in gs):
        p = [g[0] for g in gs]
        return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        if best is None or f < best: best = f
    return best

def Hdp(a, n):
    dp = {}
    for v in range(n): dp[(1<<v, v)] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            if not (S & (1<<v)): continue
            val = dp.get((S,v), 0)
            if val == 0: continue
            for u in range(n):
                if S & (1<<u): continue
                if a[v][u]:
                    dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
    return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))


# ============================================================================
# STAIRCASE DECOMPOSITION
# ============================================================================

def grid(n):
    """All cells (r,c) in Grid(n)."""
    return [(r,c) for r in range(1, n-1) for c in range(1, n-r)]

def cell_to_arc(r, c):
    """Cell (r,c) -> 0-indexed arc (i,j) where i=r+c, j=c-1."""
    return (r + c, c - 1)

def decompose(n):
    """Decompose Grid(n) into bottom, top, overlap, and apex."""
    cells = grid(n)
    m = len(cells)

    bottom = [(r,c) for r,c in cells if r+c <= n-2]   # Grid(n-1)
    top = [(r,c) for r,c in cells if c >= 2]           # ≅ Grid(n-1)
    overlap = [(r,c) for r,c in cells if c >= 2 and r+c <= n-2]  # ≅ Grid(n-2)
    apex = (n-2, 1)  # The single apex cell

    # Cells ONLY in bottom (not top, not apex)
    bottom_only = [(r,c) for r,c in bottom if c < 2 and (r,c) != apex]
    # Cells ONLY in top (not bottom, not apex)
    top_only = [(r,c) for r,c in top if r+c > n-2]

    return {
        'cells': cells, 'm': m,
        'bottom': bottom, 'top': top, 'overlap': overlap,
        'apex': apex,
        'bottom_only': bottom_only, 'top_only': top_only
    }


# ============================================================================
# ANALYSIS
# ============================================================================

for n in [4, 5, 6, 7]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")

    d = decompose(n)
    m = d['m']
    apex_r, apex_c = d['apex']
    apex_arc = cell_to_arc(apex_r, apex_c)

    print(f"\n  Grid(n={n}): {m} cells = C({n-1},2) = {comb(n-1,2)}")
    print(f"  Bottom (Grid({n-1})): {len(d['bottom'])} cells = C({n-2},2) = {comb(n-2,2)}")
    print(f"  Top (Grid({n-1}) shifted): {len(d['top'])} cells = C({n-2},2) = {comb(n-2,2)}")
    print(f"  Overlap (Grid({n-2})): {len(d['overlap'])} cells = C({n-3},2) = {comb(n-3,2)}")
    print(f"  Apex: ({apex_r},{apex_c}) = arc {apex_arc} (vertices 0 and {n-1})")
    print(f"  Bottom-only (column c=1, not apex): {len(d['bottom_only'])} cells")
    print(f"  Top-only (strip s=n-1, not apex): {len(d['top_only'])} cells")

    # Verify: m = |bottom| + |top| - |overlap| + 1
    check = len(d['bottom']) + len(d['top']) - len(d['overlap']) + 1
    print(f"\n  Inclusion-exclusion: |B|+|T|-|O|+1 = {len(d['bottom'])}+{len(d['top'])}-{len(d['overlap'])}+1 = {check} = m={m} {'OK' if check==m else 'FAIL'}")

    # Partition: cells = bottom_only ∪ overlap ∪ top_only ∪ {apex}
    partition_size = len(d['bottom_only']) + len(d['overlap']) + len(d['top_only']) + 1
    print(f"  Partition: |B_only|+|overlap|+|T_only|+1 = {len(d['bottom_only'])}+{len(d['overlap'])}+{len(d['top_only'])}+1 = {partition_size} = m={m} {'OK' if partition_size==m else 'FAIL'}")

    # Visualize
    print(f"\n  STAIRCASE VISUALIZATION (B=bottom-only, O=overlap, T=top-only, A=apex):")
    max_r = n - 3
    for r in range(max_r, 0, -1):
        row = f"  r={r}: "
        for c in range(1, n-r):
            cell = (r, c)
            if cell == d['apex']:
                row += " A"
            elif cell in d['overlap']:
                row += " O"
            elif cell in d['bottom_only']:
                row += " B"
            elif cell in d['top_only']:
                row += " T"
            else:
                row += " ?"
        print(row)
    cols = "       " + "".join(f" {c}" for c in range(1, n-1))
    print(cols)

    # ================================================================
    # THE RECURSIVE STRUCTURE OF H
    # ================================================================
    if n <= 6:
        print(f"\n  --- RECURSIVE H DECOMPOSITION ---")

        tiles = grid(n)
        tile_to_idx = {t: i for i, t in enumerate(tiles)}
        tile_arcs = [cell_to_arc(r, c) for r, c in tiles]

        # Build all tilings and classify
        def build_tournament(tile_bits):
            a = [[0]*n for _ in range(n)]
            for j in range(1, n): a[j][j-1] = 1  # base path
            for k in range(m):
                i, j = tile_arcs[k]
                if tile_bits & (1 << k): a[i][j] = 1
                else: a[j][i] = 1
            return a

        # For each tiling, record:
        # - The apex bit
        # - The bottom sub-tiling (restriction to bottom cells)
        # - The top sub-tiling (restriction to top cells)
        # - The overlap sub-tiling

        apex_idx = tile_to_idx[d['apex']]
        bottom_idxs = [tile_to_idx[c] for c in d['bottom']]
        top_idxs = [tile_to_idx[c] for c in d['top']]
        overlap_idxs = [tile_to_idx[c] for c in d['overlap']]
        bottom_only_idxs = [tile_to_idx[c] for c in d['bottom_only']]
        top_only_idxs = [tile_to_idx[c] for c in d['top_only']]

        # H of each tiling
        H_by_apex = {0: [], 1: []}

        for tb in range(1 << m):
            T = build_tournament(tb)
            H = Hdp(T, n)
            apex_bit = (tb >> apex_idx) & 1
            H_by_apex[apex_bit].append(H)

        avg_H_0 = sum(H_by_apex[0]) / len(H_by_apex[0])
        avg_H_1 = sum(H_by_apex[1]) / len(H_by_apex[1])
        print(f"    Apex=0 ({len(H_by_apex[0])} tilings): avg H = {avg_H_0:.2f}")
        print(f"    Apex=1 ({len(H_by_apex[1])} tilings): avg H = {avg_H_1:.2f}")
        print(f"    Overall avg H = {(avg_H_0 + avg_H_1) / 2:.2f}")

        # Does the apex bit correlate with H?
        all_H = [Hdp(build_tournament(tb), n) for tb in range(1 << m)]
        total_H = sum(all_H)
        print(f"    Total sum H = {total_H} (= n! = {factorial(n)} ? {'YES' if total_H == factorial(n) else 'NO'})")

        # THE REDEI SUM: sum over all tilings of H(T) = n!
        # Because each Hamiltonian path appears in exactly one tiling
        # (the tiling whose base path is that Hamiltonian path... wait, no)
        # Actually: sum over all 2^m tilings of H(T_tiling) counts
        # the total number of (tiling, Hamiltonian path) pairs.
        # Each labeled tournament T appears H(T)/|Aut(T)| times as a tiling.
        # So sum = sum_C (H(C)/|Aut(C)|) * H(C) = sum_C H(C)^2/|Aut(C)|

        print(f"    sum H^2/|Aut| would be: sum of H(tiling)... let me compute")
        print(f"    Actually sum_tilings H = {total_H}")

        # The TRUE identity: sum over all LABELED tournaments of H(T) = n! * avg_H
        # = n! * (avg paths per tournament)
        # With 2^m tilings, sum of H over tilings is NOT n! in general

        # Check Redei's theorem in tiling model:
        # Each tiling gives a tournament with ODD number of Hamiltonian paths.
        odd_count = sum(1 for h in all_H if h % 2 == 1)
        print(f"    All H odd (Redei)? {odd_count}/{len(all_H)} = {100*odd_count/len(all_H):.0f}%")

        # ================================================================
        # APEX CONTROLS HAMILTONIAN CYCLE
        # ================================================================
        print(f"\n  --- APEX AND HAMILTONIAN CYCLES ---")

        # The arc between 0 and n-1: if it goes n-1 -> 0, then
        # P_0 = (n-1)->...->0 can be closed to form a cycle n-1->...->0->n-1
        # Wait, P_0 goes n-1 -> n-2 -> ... -> 0. If arc 0->n-1 exists (apex=1),
        # then there's a Hamiltonian cycle: n-1 -> n-2 -> ... -> 0 -> n-1.

        has_cycle_0 = 0; has_cycle_1 = 0; no_cycle_0 = 0; no_cycle_1 = 0
        for tb in range(1 << m):
            T = build_tournament(tb)
            apex_bit = (tb >> apex_idx) & 1
            # Check if there's a Hamiltonian cycle
            # P_0 closes iff arc from vertex 0 to vertex n-1 exists
            # In our encoding: apex_bit=1 means a[n-1][0]=1? Or a[0][n-1]=1?
            # cell (n-2, 1) = arc (r+c, c-1) = arc (n-1, 0)
            # apex_bit=1 means a[n-1][0]=1, i.e., n-1 beats 0
            # P_0 = n-1 -> ... -> 0. For cycle need 0 -> n-1, i.e., a[0][n-1]=1
            # That happens when apex_bit = 0 (since bit=0 means lower beats upper: a[0][n-1]=1)
            p0_closes = (apex_bit == 0)
            if p0_closes:
                has_cycle_0 += 1
            else:
                no_cycle_1 += 1

        print(f"    Apex=0 (vertex 0 beats vertex {n-1}): P_0 closes into cycle: {has_cycle_0} tilings")
        print(f"    Apex=1 (vertex {n-1} beats vertex 0): P_0 does NOT close: {1<<(m-1)} tilings")

        # Avg H when P_0 is a cycle vs not
        print(f"    Avg H when P_0-cycle exists: {avg_H_0:.2f}")
        print(f"    Avg H when P_0-cycle absent: {avg_H_1:.2f}")

        # ================================================================
        # RECURSIVE H: H(T) from sub-tournaments
        # ================================================================
        print(f"\n  --- H RECURSION: T = T\\(n-1) + T\\0 + apex ---")
        print(f"    Checking: H(T) = H(T\\(n-1)) + H(T\\0) when apex closes P_0")
        print(f"    (This would be the tiling recursion)")

        match_count = 0; total_check = 0
        for tb in range(min(1 << m, 256)):  # sample for speed
            T = build_tournament(tb)
            H_full = Hdp(T, n)

            # T\(n-1): remove vertex n-1
            T_bot = [[T[i][j] for j in range(n-1)] for i in range(n-1)]
            H_bot = Hdp(T_bot, n-1)

            # T\0: remove vertex 0
            T_top = [[T[i][j] for j in range(1, n)] for i in range(1, n)]
            H_top = Hdp(T_top, n-1)

            # Check various formulas
            apex_bit = (tb >> apex_idx) & 1

            # Claim A says: H(T) - H(T\v) = 2 * sum cycles through v
            # For v = n-1: H(T) = H(T\(n-1)) + 2 * (sum of odd cycles through n-1)
            # For v = 0: H(T) = H(T\0) + 2 * (sum of odd cycles through 0)

            # But what about H(T) in terms of BOTH deletions?
            # H(T\(n-1)) + H(T\0) - H(T\{0,n-1}) might relate to H(T)?
            T_mid = [[T[i][j] for j in range(1, n-1)] for i in range(1, n-1)]
            H_mid = Hdp(T_mid, n-2) if n-2 >= 2 else 1

            total_check += 1
            # Test: H(T) = H_bot + H_top - H_mid + correction?
            diff = H_full - (H_bot + H_top - H_mid)
            if total_check <= 10 or diff == 0:
                print(f"    tb={tb:0{m}b}: H={H_full}, H_bot={H_bot}, H_top={H_top}, "
                      f"H_mid={H_mid}, H_bot+H_top-H_mid={H_bot+H_top-H_mid}, "
                      f"diff={diff:+d}, apex={apex_bit}")

        # Compute the correction for ALL tilings
        print(f"\n    FULL CORRECTION ANALYSIS:")
        corrections = Counter()
        for tb in range(1 << m):
            T = build_tournament(tb)
            H_full = Hdp(T, n)
            T_bot = [[T[i][j] for j in range(n-1)] for i in range(n-1)]
            T_top = [[T[i][j] for j in range(1, n)] for i in range(1, n)]
            T_mid = [[T[i][j] for j in range(1, n-1)] for i in range(1, n-1)]
            H_bot = Hdp(T_bot, n-1)
            H_top = Hdp(T_top, n-1)
            H_mid = Hdp(T_mid, n-2) if n-2 >= 2 else 1
            apex_bit = (tb >> apex_idx) & 1
            diff = H_full - (H_bot + H_top - H_mid)
            corrections[(diff, apex_bit)] += 1

        print(f"    {'diff':>6} {'apex':>5} {'count':>6}")
        for (diff, apex), count in sorted(corrections.items()):
            print(f"    {diff:+6d} {apex:5d} {count:6d}")

        # Is the correction related to the apex?
        # If correction = +1 when apex creates cycle, -1 otherwise?
        print(f"\n    Summary:")
        for apex in [0, 1]:
            diffs = [d for (d, a), c in corrections.items() if a == apex for _ in range(c)]
            if diffs:
                print(f"      Apex={apex}: corrections = {sorted(set(diffs))}, "
                      f"avg = {sum(diffs)/len(diffs):.2f}")


print(f"\n\n{'='*80}")
print("  THE RECURSIVE IDENTITIES")
print(f"{'='*80}")
print("""
  STAIRCASE DECOMPOSITION:
    delta_{n-2} = bottom(n-3) U top(n-3) U {apex}
    |bottom| = |top| = C(n-2,2)
    |overlap| = C(n-3,2)
    |apex| = 1
    m(n) = 2*m(n-1) - m(n-2) + 1

  IN TOURNAMENT TERMS:
    T = T\\(n-1) U T\\0 U {arc(0,n-1)}
    T\\(n-1) ∩ T\\0 = T\\{0,n-1}

  VISUALIZATION:
    For n=6, Grid(6) has 10 cells:
         c=1  c=2  c=3  c=4
    r=4:  A
    r=3:  B    T
    r=2:  B    O    T
    r=1:  B    O    O    T

    B = bottom-only = left column (arcs FROM vertex 0)
    T = top-only = top anti-diagonal (arcs TO vertex n-1)
    O = overlap = the n-2 staircase in the middle
    A = apex = THE arc between vertices 0 and n-1

  The overlap region O is a SMALLER COPY of the full staircase,
  encoding the tournament T\\{0, n-1} on the n-2 internal vertices.

  This gives rise to the META-GRAPH RECURSION:
  G_n/Z_2 is built from two copies of G_{n-1}/Z_2
  sharing G_{n-2}/Z_2 in common, plus one apex bit
  that controls Hamiltonian cycle existence.
""")

print(f"  DONE.")
print("=" * 80)
