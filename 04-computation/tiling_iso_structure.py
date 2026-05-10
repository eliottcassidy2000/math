"""
tiling_iso_structure.py — Full analysis of the tiling model and isomorphism class structure.

CONVENTIONS (from tournament-tiling-explorer.html):
  VERTS = [N, N-1, ..., 1]
  Base Hamiltonian path: N → N-1 → ... → 2 → 1  (vertex N beats vertex N-1 etc.)
  TILES: for y in 1..N-2: for x in N down to y+2: tile [x, y]
    bit=0: A[xi][yi]=1 (x beats y, "forward")
    bit=1: A[yi][xi]=1 (y beats x, "backward/active")
  TRANS_MAP: tile (x,y) ↔ tile (N-y+1, N-x+1)  [reflect along diagonal x+y=N+1]

BLUE LINE: flip-pair {T, T'} where T' = all tile bits flipped, BOTH are grid-symmetric
BLACK LINE: flip-pair where BOTH are NOT grid-symmetric (if T is grid-sym, T' is too)

SC/NS PAIRING: iso classes {C, C^op} where C^op = full arc complement tournament
  NS pairs: C ≠ C^op, merged node has EVEN #tilings (= 2 × #tilings(C))
  SC classes: C = C^op, has ODD #tilings (since H odd and H = |Aut|×#tilings implies both odd)

oracle-2026-05-10
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(encoding='utf-8')

# ══════════════════════════════════════════════════════════════════════════
# CORE TILING MODEL (exact explorer convention)
# ══════════════════════════════════════════════════════════════════════════

def make_tiles(n):
    """Tiles in explorer order: y from 1..n-2, x from n down to y+2."""
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    return tiles

def make_trans_map(n, tiles):
    tile_idx = {(x,y): i for i,(x,y) in enumerate(tiles)}
    return [tile_idx[(n-y+1, n-x+1)] for x,y in tiles]

def bits_to_adj(bits, n, tiles):
    verts = list(range(n, 0, -1))  # [N, N-1, ..., 1]
    A = [[0]*n for _ in range(n)]
    for k in range(n-1): A[k][k+1] = 1  # base path N→N-1→...→1
    for i, (xl, yl) in enumerate(tiles):
        xi = verts.index(xl); yi = verts.index(yl)
        if bits[i] == 0: A[xi][yi] = 1  # x beats y
        else: A[yi][xi] = 1             # y beats x
    return A

def adj_to_str(A): return ''.join(''.join(map(str,r)) for r in A)

def is_grid_sym(bits, trans_map):
    for i, ti in enumerate(trans_map):
        if ti != i and bits[i] != bits[ti]: return False
    return True

def full_complement_adj(A, n):
    """Full arc complement: reverse all arcs."""
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def canonicalize(A, n, perms):
    best = None
    for p in perms:
        s = adj_to_str([[A[p[i]][p[j]] for j in range(n)] for i in range(n)])
        if best is None or s < best: best = s
    return best

def count_ham_paths(A, n):
    dp = {}
    for v in range(n): dp[(1<<v, v)] = 1
    for ms in range(2, n+1):
        for mask in range(1<<n):
            if bin(mask).count('1') != ms: continue
            for v in range(n):
                if not (mask>>v&1): continue
                pm = mask ^ (1<<v)
                t = sum(dp.get((pm,u),0) for u in range(n) if (pm>>u&1) and A[u][v])
                if t: dp[(mask,v)] = t
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

def get_automorphism_size(A, n, perms):
    count = 0
    for p in perms:
        if all(A[p[i]][p[j]] == A[i][j] for i in range(n) for j in range(n)):
            count += 1
    return count

# ══════════════════════════════════════════════════════════════════════════
# BUILD ISO CLASS DATA
# ══════════════════════════════════════════════════════════════════════════

def build_iso_data(n):
    tiles = make_tiles(n)
    m = len(tiles)
    trans_map = make_trans_map(n, tiles)
    perms = list(permutations(range(n)))
    num_T = 1 << m

    tilings = []
    for mask in range(num_T):
        bits = [(mask>>k)&1 for k in range(m)]
        A = bits_to_adj(bits, n, tiles)
        canon = canonicalize(A, n, perms)
        flip_mask = sum((1-b)<<k for k,b in enumerate(bits))
        Acomp = full_complement_adj(A, n)
        comp_canon = canonicalize(Acomp, n, perms)
        grid_sym = is_grid_sym(bits, trans_map)
        tilings.append({
            'mask': mask, 'bits': bits, 'A': A,
            'canon': canon, 'comp_canon': comp_canon,
            'flip_mask': flip_mask, 'grid_sym': grid_sym
        })

    # Group by canonical signature → iso classes
    groups = defaultdict(list)
    for t in tilings: groups[t['canon']].append(t)
    canon_sigs = sorted(groups.keys())

    # For each iso class: compute H and |Aut|
    classes = []
    for sig in canon_sigs:
        members = groups[sig]
        A = members[0]['A']
        H = count_ham_paths(A, n)
        aut_size = get_automorphism_size(A, n, perms)
        comp_sig = members[0]['comp_canon']  # canonical sig of complement
        sc = (comp_sig == sig)  # self-complementary
        classes.append({
            'sig': sig,
            'members': members,
            'num_tilings': len(members),
            'H': H,
            'aut_size': aut_size,
            'comp_sig': comp_sig,
            'sc': sc,
        })

    return classes, tiles, m, trans_map, tilings

# ══════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS
# ══════════════════════════════════════════════════════════════════════════

for n in range(3, 7):
    print("=" * 70)
    print(f"n={n}: COMPLETE ISO CLASS ANALYSIS")
    print("=" * 70)

    classes, tiles, m, trans_map, tilings = build_iso_data(n)
    total_tilings = 1 << m

    print(f"\n  Tiles (explorer order): {tiles}")
    print(f"  TRANS_MAP: {trans_map}")
    print(f"  m={m} tiles, 2^m={total_tilings} total tilings")
    print(f"  Iso classes: {len(classes)}")

    # Count grid-symmetric tilings
    grid_sym_count = sum(1 for t in tilings if t['grid_sym'])
    print(f"  Grid-symmetric tilings: {grid_sym_count} = 2^{m - (m-grid_sym_count.bit_length()+1)}?")
    # Actually: grid-sym tilings = 2^(free bits). Free bits = self-paired tiles + 1 per paired-tile-pair
    self_paired = sum(1 for i,ti in enumerate(trans_map) if ti==i)
    non_self_pairs = (m - self_paired)//2
    free_bits = self_paired + non_self_pairs  # self-paired: free; each non-self-pair: one free bit
    print(f"  Self-paired tiles: {self_paired}, Non-self pairs: {non_self_pairs}")
    print(f"  Free bits for grid-sym: {free_bits} → 2^{free_bits} = {1<<free_bits} grid-sym tilings")

    # Blue and black flip-pairs
    blue_pairs = 0; black_pairs = 0
    seen_fp = set()
    mask_to_t = {t['mask']: t for t in tilings}
    for t in tilings:
        fm = t['flip_mask']
        pk = (min(t['mask'],fm), max(t['mask'],fm))
        if pk in seen_fp: continue
        seen_fp.add(pk)
        t2 = mask_to_t[fm]
        if t['grid_sym'] and t2['grid_sym']: blue_pairs += 1
        else: black_pairs += 1
    print(f"  Blue flip-pairs (both grid-sym): {blue_pairs}")
    print(f"  Black flip-pairs (neither grid-sym): {black_pairs}")
    print(f"  Total flip-pairs: {blue_pairs + black_pairs} = 2^m / 2 = {total_tilings//2}")

    # SC vs NS analysis
    sc_classes = [c for c in classes if c['sc']]
    ns_classes = [c for c in classes if not c['sc']]
    print(f"\n  SC (self-complementary) classes: {len(sc_classes)}")
    print(f"  NS (non-self-complementary) classes: {len(ns_classes)}")
    ns_pairs = []
    seen_ns = set()
    for c in ns_classes:
        if c['sig'] in seen_ns: continue
        # Find partner
        partner = next(x for x in ns_classes if x['sig'] == c['comp_sig'])
        seen_ns.add(c['sig']); seen_ns.add(partner['sig'])
        ns_pairs.append((c, partner))
    print(f"  NS pairs: {len(ns_pairs)}")

    # Verify tiling counts are all ODD
    all_odd = all(c['num_tilings'] % 2 == 1 for c in classes)
    print(f"  All tiling counts odd? {all_odd} (expected: True, since H=odd=|Aut|×#tilings)")

    # Build the merged sum
    merged_items = []
    for c1, c2 in ns_pairs:
        merged_items.append(c1['num_tilings'] + c2['num_tilings'])  # even
    for c in sc_classes:
        merged_items.append(c['num_tilings'])  # odd
    merged_items.sort()
    print(f"\n  MERGED tiling counts (sorted): {merged_items}")
    print(f"  Sum: {sum(merged_items)} = 2^{m} ✓" if sum(merged_items)==total_tilings else f"  Sum: {sum(merged_items)} ✗")

    # Detailed class table
    print(f"\n  Iso class table:")
    print(f"  {'ci':>3} | {'#tilings':>8} | {'H':>6} | {'|Aut|':>6} | {'H/|Aut|':>8} | {'SC':>4} | {'grid-sym tilings'}")
    print("  " + "-"*70)
    for i,c in enumerate(classes):
        gs_count = sum(1 for t in c['members'] if t['grid_sym'])
        ck = "SC" if c['sc'] else "NS"
        H_check = c['H'] // c['aut_size']
        ok = "✓" if H_check == c['num_tilings'] else "✗"
        print(f"  {i:>3} | {c['num_tilings']:>8} | {c['H']:>6} | {c['aut_size']:>6} | {H_check:>8} {ok} | {ck:>4} | {gs_count}")

    # Zeckendorf tilings in each class
    print(f"\n  Zeckendorf tilings (no two consecutive active tiles, explorer order) per class:")
    def is_zeckendorf(bits):
        return not any(bits[k]==1 and bits[k+1]==1 for k in range(len(bits)-1))
    zeck_tilings = [t for t in tilings if is_zeckendorf(t['bits'])]
    zeck_per_class = defaultdict(int)
    for t in zeck_tilings:
        zeck_per_class[t['canon']] += 1
    from math import isqrt
    from fractions import Fraction
    def fib(k):
        if k<=0: return 0
        a,b=1,1
        for _ in range(k-2): a,b=b,a+b
        return b if k>=2 else 1
    Zm = fib(m+2)
    print(f"  Total Zeckendorf tilings: {len(zeck_tilings)} = F_{{{m+2}}} = {Zm}")
    print(f"  {'ci':>3} | {'#Z-tilings':>10} | {'#tilings':>8} | {'Z/#tot':>8} | {'H':>5} | SC")
    for i,c in enumerate(classes):
        zc = zeck_per_class.get(c['sig'], 0)
        ratio = Fraction(zc, c['num_tilings']) if c['num_tilings'] else Fraction(0)
        print(f"  {i:>3} | {zc:>10} | {c['num_tilings']:>8} | {ratio!s:>8} | {c['H']:>5} | {'SC' if c['sc'] else 'NS'}")

    # Zeckendorf connection to Fibonacci
    print(f"\n  Sum of Zeckendorf counts: {sum(zeck_per_class.values())} = F_{{{m+2}}} = {Zm}")

    # H distribution over Zeckendorf tilings
    H_zeck = Counter(t['H'] for t in zeck_tilings for t in [{'H': count_ham_paths(t['A'],n)}] if True)
    # Actually compute H per zeckendorf tiling
    H_zeck_dist = Counter()
    for t in zeck_tilings:
        H_zeck_dist[count_ham_paths(t['A'],n)] += 1
    print(f"  H-distribution over Zeckendorf tilings: {dict(sorted(H_zeck_dist.items()))}")

    print()

# ══════════════════════════════════════════════════════════════════════════
# PATTERN EXTENSION: tiling counts per n
# ══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("PATTERN SUMMARY ACROSS n")
print("=" * 70)

for n in range(3, 7):
    classes, tiles, m, trans_map, tilings = build_iso_data(n)
    sc_classes = [c for c in classes if c['sc']]
    ns_classes = [c for c in classes if not c['sc']]

    # NS pairs
    seen_ns = set()
    ns_pairs = []
    for c in ns_classes:
        if c['sig'] in seen_ns: continue
        partner = next(x for x in ns_classes if x['sig'] == c['comp_sig'])
        seen_ns.add(c['sig']); seen_ns.add(partner['sig'])
        ns_pairs.append((c, partner))

    merged_items = []
    for c1,c2 in ns_pairs:
        merged_items.append(c1['num_tilings']+c2['num_tilings'])
    for c in sc_classes:
        merged_items.append(c['num_tilings'])
    merged_items.sort()

    # Tiling counts per class (unmerged)
    all_counts = sorted(c['num_tilings'] for c in classes)

    print(f"\nn={n}: m={m} tiles, 2^m={1<<m}, #classes={len(classes)}, #SC={len(sc_classes)}, #NS-pairs={len(ns_pairs)}")
    print(f"  Merged sum: {'+'.join(map(str,merged_items))} = {sum(merged_items)}")
    print(f"  Unmerged counts: {all_counts}")
    print(f"  Fibonacci |Z_n| = F_{{{m+2}}} = {fib(m+2)}")

    # Connection: do the tiling counts appear in Fibonacci/h-value sequences?
    h_vals = {1+(1<<(r-1)) for r in range(1,10)}
    fib_vals = {fib(k) for k in range(1,20)}
    for cnt in all_counts:
        tags = []
        if cnt in h_vals: tags.append(f"h({list(h_vals).index(cnt)+1})")
        if cnt in fib_vals: tags.append(f"F_{list(fib_vals).index(cnt)+1}")
        if tags: print(f"  Count {cnt} = {', '.join(tags)}")
