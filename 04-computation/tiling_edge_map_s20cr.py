#!/usr/bin/env python3
"""
tiling_edge_map_s20cr.py — Map meta-graph edges to staircase cells
kind-pasteur-2026-03-23-S20cr

KEY IDEA: Each edge in G_n/Z_2 corresponds to flipping ONE arc in a
tournament. Each arc corresponds to ONE CELL in the staircase delta_{n-2}.

QUESTION: Which cells generate BLUE edges (SC-preserving) and which
generate BLACK edges (SC-breaking)? Does this pattern have the
diagonal/perpendicular symmetry of the staircase triangle?

PIN GRID COORDINATES:
  Cell (r,c): r = row (1-indexed from bottom), c = column (1-indexed from left)
  Constraint: r >= 1, c >= 1, r + c <= n-1
  Arc: cell (r,c) corresponds to arc from vertex (c+r) to vertex (c-1)
       [using 0-indexed vertices, the "non-base-path" arc]

  Actually using our convention: pairs P = [(i,j) for i<j], bit k = P[k].
  Cell in staircase: arc (i,j) with i<j maps to row r=j-i-1, col c=i
  (where r=0 is NOT in the staircase — r starts at 0 for adjacent pairs)

  STRIPS: cells with same r+c value (anti-diagonal lines parallel to hypotenuse)
  DIAGONAL: cells with r = c (principal diagonal of the triangle)
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TILING EDGE MAP: STAIRCASE CELLS -> META-GRAPH EDGES")
print("  kind-pasteur-2026-03-23-S20cr")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================
def tadj(n, bits):
    a = [[0]*n for _ in range(n)]; idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): a[i][j] = 1
            else: a[j][i] = 1
            idx += 1
    return a

def Hdp(a, n):
    dp = [0]*((1<<n)*n)
    for v in range(n): dp[(1<<v)*n+v] = 1
    for S in range(1, 1<<n):
        if bin(S).count('1') >= n: continue
        for v in range(n):
            if not(S&(1<<v)): continue
            val = dp[S*n+v]
            if val == 0: continue
            for u in range(n):
                if S&(1<<u): continue
                if a[v][u]: dp[(S|(1<<u))*n+u] += val
    return sum(dp[((1<<n)-1)*n+v] for v in range(n))

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

def comp(a, n):
    return [[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

# ============================================================================
# BUILD AND MAP EDGES TO CELLS
# ============================================================================

def analyze_tiling_edges(n):
    m = comb(n,2); total = 1<<m; t0 = time.time()
    PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Map arc index to staircase cell (row, col)
    # Arc k = (i,j) with i<j. In staircase: row = j-i-1, col = i
    # Row 0 = adjacent pairs (base path neighbors), not in "strict" staircase
    arc_to_cell = {}
    for k, (i, j) in enumerate(PAIRS):
        row = j - i - 1   # distance-1 = row 0, distance-2 = row 1, etc.
        col = i            # left vertex index
        arc_to_cell[k] = (row, col)

    # Build iso classes
    iso = defaultdict(list)
    for b in range(total):
        a = tadj(n, b); iso[canon(a, n)].append(b)

    cl = []; c2c = {}
    for idx, (cn, mem) in enumerate(sorted(iso.items())):
        a = tadj(n, mem[0])
        cl.append({'cid':idx,'canon':cn,'mem':mem,'adj':a,
                   'H':Hdp(a,n),'sc':canon(comp(a,n),n)==cn})
        c2c[cn] = idx
    for d in cl:
        cc = canon(comp(d['adj'],n),n); d['comp_cid'] = c2c.get(cc,-1)

    # Merge map
    mi = {}; mid = 0
    for d in cl:
        ci = d['cid']
        if ci in mi: continue
        cp = d['comp_cid']; mi[ci] = mid
        if cp != ci and cp >= 0: mi[cp] = mid
        mid += 1
    Vm = mid

    mc = {}
    for d in cl:
        mv = mi[d['cid']]
        if mv not in mc:
            mc[mv] = {'H':d['H'],'sc':d['sc']}

    # Track which ARC (cell) generates each edge, and its color
    # For each tournament in each class, flip each arc and see where it goes
    cell_blue = defaultdict(int)   # cell -> count of blue edges generated
    cell_black = defaultdict(int)  # cell -> count of black edges generated
    cell_total = defaultdict(int)

    # Also track per-class: which cells generate blue/black from SC vs NS classes
    cell_blue_from_sc = defaultdict(int)
    cell_black_from_sc = defaultdict(int)
    cell_blue_from_ns = defaultdict(int)
    cell_black_from_ns = defaultdict(int)

    # For edge deduplication: track unique merged edges with their cell
    edge_cells = defaultdict(set)  # (merged_a, merged_b) -> set of cells that generate it

    for d in cl:
        ci = d['cid']
        src_merged = mi[ci]
        src_sc = d['sc']

        for b in d['mem']:
            for k in range(m):
                fb = b ^ (1 << k)
                fa = tadj(n, fb)
                fc = canon(fa, n)
                nb = c2c.get(fc)
                if nb is None: continue
                tgt_merged = mi[nb]
                if tgt_merged == src_merged: continue  # self-loop

                cell = arc_to_cell[k]
                tgt_sc = cl[nb]['sc']

                # Edge color
                is_blue = (mc[src_merged]['sc'] == mc[tgt_merged]['sc'])

                if is_blue:
                    cell_blue[cell] += 1
                    if src_sc: cell_blue_from_sc[cell] += 1
                    else: cell_blue_from_ns[cell] += 1
                else:
                    cell_black[cell] += 1
                    if src_sc: cell_black_from_sc[cell] += 1
                    else: cell_black_from_ns[cell] += 1
                cell_total[cell] += 1

                e = (min(src_merged, tgt_merged), max(src_merged, tgt_merged))
                edge_cells[e].add(cell)

    print(f"  Built n={n}: V={Vm}, {len(edge_cells)} merged edges [{time.time()-t0:.1f}s]")
    return {
        'n': n, 'V': Vm, 'm': m, 'PAIRS': PAIRS,
        'arc_to_cell': arc_to_cell,
        'cell_blue': cell_blue, 'cell_black': cell_black, 'cell_total': cell_total,
        'cell_blue_from_sc': cell_blue_from_sc, 'cell_black_from_sc': cell_black_from_sc,
        'cell_blue_from_ns': cell_blue_from_ns, 'cell_black_from_ns': cell_black_from_ns,
        'edge_cells': edge_cells, 'mc': mc
    }


# ============================================================================
# VISUALIZE STAIRCASE HEAT MAP
# ============================================================================

def print_staircase(n, data, title, value_fn):
    """Print the staircase with values from value_fn(row, col)."""
    max_row = n - 2  # rows 0 to n-3
    print(f"\n  {title}")
    print(f"  Staircase delta_{n-2} (row = arc distance - 1, col = left vertex)")
    print(f"  Row 0 = adjacent arcs (base path), Row {max_row-1} = longest arcs")
    print()

    for r in range(max_row - 1, -1, -1):  # top to bottom
        max_c = n - 2 - r  # columns 0 to n-2-r-1? Actually r+c <= n-2-1 = n-3
        # In our convention: row r, col c with r >= 0, c >= 0, r+c <= n-3
        line = f"  r={r}: "
        for c in range(max_c):
            val = value_fn(r, c)
            line += f" {val:>6}"
        print(line)

    # Column labels
    line = "       "
    for c in range(n-1):
        line += f" {'c='+str(c):>6}"
    print(line)


# ============================================================================
# MAIN
# ============================================================================

for n in [3, 4, 5, 6]:
    print(f"\n{'#'*80}")
    print(f"  n = {n}")
    print(f"{'#'*80}")

    data = analyze_tiling_edges(n)

    # ================================================================
    # 1. CELL HEAT MAP: BLUE FRACTION
    # ================================================================
    def blue_frac(r, c):
        b = data['cell_blue'].get((r,c), 0)
        t = data['cell_total'].get((r,c), 0)
        if t == 0: return "  .   "
        return f"{b/t:.3f}"

    print_staircase(n, data, "BLUE FRACTION per cell (blue_edges / total_edges)", blue_frac)

    # ================================================================
    # 2. CELL HEAT MAP: TOTAL EDGE COUNT
    # ================================================================
    def total_count(r, c):
        t = data['cell_total'].get((r,c), 0)
        return f"{t:6d}" if t > 0 else "     ."

    print_staircase(n, data, "TOTAL EDGES per cell", total_count)

    # ================================================================
    # 3. DIAGONAL vs OFF-DIAGONAL
    # ================================================================
    print(f"\n  --- DIAGONAL vs OFF-DIAGONAL ANALYSIS ---")

    diag_blue = 0; diag_black = 0; diag_total = 0
    offdiag_blue = 0; offdiag_black = 0; offdiag_total = 0

    for (r, c), b in data['cell_blue'].items():
        if r == c:
            diag_blue += b
        else:
            offdiag_blue += b

    for (r, c), b in data['cell_black'].items():
        if r == c:
            diag_black += b
        else:
            offdiag_black += b

    for (r, c), t in data['cell_total'].items():
        if r == c:
            diag_total += t
        else:
            offdiag_total += t

    if diag_total > 0:
        print(f"    Diagonal (r=c): blue={diag_blue}, black={diag_black}, total={diag_total}, "
              f"blue_frac={diag_blue/diag_total:.4f}")
    if offdiag_total > 0:
        print(f"    Off-diagonal: blue={offdiag_blue}, black={offdiag_black}, total={offdiag_total}, "
              f"blue_frac={offdiag_blue/offdiag_total:.4f}")

    # ================================================================
    # 4. STRIP ANALYSIS (anti-diagonal: cells with same r+c)
    # ================================================================
    print(f"\n  --- STRIP ANALYSIS (r+c = const = anti-diagonal) ---")

    strips = defaultdict(lambda: {'blue':0, 'black':0, 'total':0})
    for (r, c), b in data['cell_blue'].items():
        strips[r+c]['blue'] += b
    for (r, c), b in data['cell_black'].items():
        strips[r+c]['black'] += b
    for (r, c), t in data['cell_total'].items():
        strips[r+c]['total'] += t

    for s in sorted(strips.keys()):
        d = strips[s]
        cells_in_strip = min(s+1, n-2-s)  # number of cells in this strip
        bf = d['blue']/d['total'] if d['total'] > 0 else 0
        print(f"    Strip {s} ({cells_in_strip} cells): blue={d['blue']:6d} black={d['black']:6d} "
              f"total={d['total']:6d} blue_frac={bf:.4f}")

    # ================================================================
    # 5. ROW ANALYSIS (constant arc distance)
    # ================================================================
    print(f"\n  --- ROW ANALYSIS (constant row = constant arc distance) ---")

    rows = defaultdict(lambda: {'blue':0, 'black':0, 'total':0})
    for (r, c), b in data['cell_blue'].items():
        rows[r]['blue'] += b
    for (r, c), b in data['cell_black'].items():
        rows[r]['black'] += b
    for (r, c), t in data['cell_total'].items():
        rows[r]['total'] += t

    for r in sorted(rows.keys()):
        d = rows[r]
        bf = d['blue']/d['total'] if d['total'] > 0 else 0
        cells_in_row = n - 2 - r
        print(f"    Row {r} ({cells_in_row} cells, dist={r+1}): blue={d['blue']:6d} black={d['black']:6d} "
              f"total={d['total']:6d} blue_frac={bf:.4f}")

    # ================================================================
    # 6. COLUMN ANALYSIS (constant left vertex)
    # ================================================================
    print(f"\n  --- COLUMN ANALYSIS (constant col = constant source vertex) ---")

    cols = defaultdict(lambda: {'blue':0, 'black':0, 'total':0})
    for (r, c), b in data['cell_blue'].items():
        cols[c]['blue'] += b
    for (r, c), b in data['cell_black'].items():
        cols[c]['black'] += b
    for (r, c), t in data['cell_total'].items():
        cols[c]['total'] += t

    for c in sorted(cols.keys()):
        d = cols[c]
        bf = d['blue']/d['total'] if d['total'] > 0 else 0
        cells_in_col = n - 2 - c
        print(f"    Col {c} ({cells_in_col} cells, vertex {c}): blue={d['blue']:6d} black={d['black']:6d} "
              f"total={d['total']:6d} blue_frac={bf:.4f}")

    # ================================================================
    # 7. REFLECTION SYMMETRY TEST
    # ================================================================
    print(f"\n  --- REFLECTION SYMMETRY: cell(r,c) vs cell(c,r) ---")

    max_rc = n - 3
    sym_match = 0; sym_total = 0
    for r in range(max_rc + 1):
        for c in range(r + 1, n - 2 - r):
            t1 = data['cell_total'].get((r, c), 0)
            t2 = data['cell_total'].get((c, r), 0)
            b1 = data['cell_blue'].get((r, c), 0)
            b2 = data['cell_blue'].get((c, r), 0)
            if t1 > 0 and t2 > 0:
                sym_total += 1
                bf1 = b1/t1; bf2 = b2/t2
                if abs(bf1 - bf2) < 0.001:
                    sym_match += 1
                if n <= 5 or abs(bf1 - bf2) > 0.01:
                    print(f"    ({r},{c}): bf={bf1:.4f}  vs  ({c},{r}): bf={bf2:.4f}  "
                          f"{'MATCH' if abs(bf1-bf2)<0.001 else f'DIFF={bf1-bf2:+.4f}'}")

    if sym_total > 0:
        print(f"    Symmetric pairs: {sym_match}/{sym_total} exact matches ({100*sym_match/sym_total:.1f}%)")

    # ================================================================
    # 8. SC vs NS SOURCE ANALYSIS
    # ================================================================
    print(f"\n  --- SC vs NS SOURCE ANALYSIS ---")
    print(f"    Blue edges from SC tournaments vs NS tournaments per cell:")

    def sc_ns_ratio(r, c):
        bs = data['cell_blue_from_sc'].get((r,c), 0)
        bn = data['cell_blue_from_ns'].get((r,c), 0)
        ks = data['cell_black_from_sc'].get((r,c), 0)
        kn = data['cell_black_from_ns'].get((r,c), 0)
        total_sc = bs + ks
        total_ns = bn + kn
        if total_sc + total_ns == 0: return "  .   "
        sc_bf = bs/total_sc if total_sc > 0 else 0
        ns_bf = bn/total_ns if total_ns > 0 else 0
        return f"{sc_bf:.2f}/{ns_bf:.2f}"

    print_staircase(n, data, "SC_blue_frac / NS_blue_frac per cell", sc_ns_ratio)


# ============================================================================
# CROSS-n SYNTHESIS
# ============================================================================

print(f"\n\n{'='*80}")
print("  CROSS-n PATTERNS")
print(f"{'='*80}")

print("""
  KEY QUESTIONS:
  1. Is the blue fraction UNIFORM across cells, or does it vary?
  2. Do diagonal cells (r=c) behave differently from off-diagonal?
  3. Is there (r,c) <-> (c,r) reflection symmetry?
  4. Do strips (anti-diagonals) have monotone blue fraction?
  5. Does the hypotenuse (highest strip) generate more blue than the legs?
""")

print(f"\n  DONE.")
print("=" * 80)
