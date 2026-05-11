"""
n7_four_category.py — Efficient four-category analysis for n=3..7.

KEY ALGORITHM: orbit-based (vs per-tiling canonicalization).
For each unprocessed mask, compute the full S_n-orbit in O(n! × m) time.
This gives orbit_size = #tilings, |Aut| = n!/orbit_size, SC/NS, blue count — all at once.
Total: O(num_classes × n! × m) instead of O(2^m × n! × n²).

For n=7: 191 classes × 5040 × 15 ≈ 14M ops instead of 32768 × 5040 × 49 ≈ 8B ops.

oracle-2026-05-11
"""

import sys, time
from itertools import permutations
from collections import Counter
sys.stdout.reconfigure(encoding='utf-8')

def make_tiles(n):
    return [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]

def make_trans_map(n, tiles):
    idx = {t:i for i,t in enumerate(tiles)}
    return [idx[(n-y+1,n-x+1)] for x,y in tiles]

def H_count(A, n):
    dp = {}
    for v in range(n): dp[(1<<v, v)] = 1
    full = (1<<n) - 1
    masks_by_pop = [[] for _ in range(n+1)]
    for mask in range(1<<n): masks_by_pop[bin(mask).count('1')].append(mask)
    for ms in range(2, n+1):
        for mask in masks_by_pop[ms]:
            for v in range(n):
                if not (mask>>v&1): continue
                pm = mask ^ (1<<v)
                t = sum(dp.get((pm,u),0) for u in range(n) if (pm>>u&1) and A[u][v])
                if t: dp[(mask,v)] = t
    return sum(dp.get((full,v),0) for v in range(n))

def analyze_n(n):
    t_start = time.time()
    tiles = make_tiles(n)
    m = len(tiles)
    tm = make_trans_map(n, tiles)
    verts = list(range(n, 0, -1))
    xi_arr = [verts.index(x) for x,y in tiles]
    yi_arr = [verts.index(y) for x,y in tiles]
    perms = list(permutations(range(n)))
    np_count = len(perms)

    # Precompute flat index lookups: for each (perm, tile_k): sigma[xi_k]*n + sigma[yi_k]
    perm_xy = [[perms[p][xi_arr[k]]*n + perms[p][yi_arr[k]] for k in range(m)] for p in range(np_count)]
    perm_yx = [[perms[p][yi_arr[k]]*n + perms[p][xi_arr[k]] for k in range(m)] for p in range(np_count)]
    tile_masks = [1<<k for k in range(m)]

    # Grid-sym: pairs (i,j) with i<j and tm[i]=j
    gs_pairs = [(i, tm[i]) for i in range(m) if tm[i] > i]
    # Precompute all grid-sym masks
    all_gs = set()
    for mask in range(1<<m):
        if all(((mask>>i)&1) == ((mask>>j)&1) for i,j in gs_pairs):
            all_gs.add(mask)

    # Free bits = self-paired + non-self-pairs
    self_paired = sum(1 for i,ti in enumerate(tm) if ti==i)
    free_bits = self_paired + (m - self_paired)//2

    def mask_to_Aflat(mask):
        A = [0]*(n*n)
        for k in range(n-1): A[k*n+k+1] = 1
        for k in range(m):
            if (mask>>k)&1 == 0: A[xi_arr[k]*n+yi_arr[k]] = 1
            else: A[yi_arr[k]*n+xi_arr[k]] = 1
        return A

    def get_orbit(Aflat):
        orbit = set()
        for p in range(np_count):
            new_mask = 0
            idxs = perm_xy[p]
            for k in range(m):
                if Aflat[idxs[k]] == 0:
                    new_mask |= tile_masks[k]
            orbit.add(new_mask)
        return orbit

    def get_comp_orbit(Aflat):
        # Orbit of T^op (full arc complement): A^op[i][j]=A[j][i], so use perm_yx
        orbit = set()
        for p in range(np_count):
            new_mask = 0
            idxs = perm_yx[p]
            for k in range(m):
                if Aflat[idxs[k]] == 0:
                    new_mask |= tile_masks[k]
            orbit.add(new_mask)
        return orbit

    def Aflat_to_2d(Aflat):
        return [[Aflat[i*n+j] for j in range(n)] for i in range(n)]

    # Main processing
    processed = set()
    classes = []
    mask_to_ci = {}

    for start_mask in range(1<<m):
        if start_mask in processed:
            continue
        Aflat = mask_to_Aflat(start_mask)
        orbit = get_orbit(Aflat)
        comp_orbit = get_comp_orbit(Aflat)

        processed |= orbit
        orbit_size = len(orbit)
        aut_size = np_count // orbit_size
        sc = bool(orbit & comp_orbit)
        gs_count = sum(1 for mk in orbit if mk in all_gs)

        ci = len(classes)
        for mk in orbit: mask_to_ci[mk] = ci
        classes.append({
            'ci': ci, 'Aflat': Aflat, 'comp_orbit': comp_orbit,
            'nt': orbit_size, 'au': aut_size, 'sc': sc, 'gs': gs_count,
            'H': None, 'comp_partner': None
        })

    t_orbit = time.time() - t_start

    # Complement partners
    for c in classes:
        if c['sc']: c['comp_partner'] = c['ci']
        else:
            pm = next(iter(c['comp_orbit']))
            c['comp_partner'] = mask_to_ci[pm]

    # H computation
    t_h = time.time()
    for c in classes:
        c['H'] = H_count(Aflat_to_2d(c['Aflat']), n)
    t_h = time.time() - t_h

    # NS pairs
    sc_classes = [c for c in classes if c['sc']]
    ns_classes = [c for c in classes if not c['sc']]
    seen_ns = set(); ns_pairs = []
    for c in ns_classes:
        if c['ci'] in seen_ns: continue
        partner = next(x for x in ns_classes if x['ci'] == c['comp_partner'])
        seen_ns.add(c['ci']); seen_ns.add(partner['ci'])
        ns_pairs.append((c, partner))

    # SC data
    sc_data = []
    for c in sc_classes:
        blue = c['gs']; black = c['nt'] - blue
        cat = 'pure_blue' if black == 0 else 'mixed'
        sc_data.append({'ci':c['ci'],'H':c['H'],'au':c['au'],'nt':c['nt'],
                        'blue':blue,'black':black,'cat':cat})

    # NS data
    ns_data = []
    for ca, cb in ns_pairs:
        nt = ca['nt'] + cb['nt']
        ns_data.append({'a':ca['ci'],'b':cb['ci'],'nt':nt,'black':nt,
                        'H_a':ca['H'],'H_b':cb['H'],'au_a':ca['au'],'au_b':cb['au']})

    total_t = 1<<m; total_blue = 1<<free_bits; total_black = total_t - total_blue

    return {
        'n': n, 'm': m, 'free_bits': free_bits,
        'total_t': total_t, 'total_blue': total_blue, 'total_black': total_black,
        'classes': classes, 'sc_classes': sc_classes, 'ns_classes': ns_classes,
        'ns_pairs': ns_pairs, 'sc_data': sc_data, 'ns_data': ns_data,
        't_orbit': t_orbit, 't_h': t_h,
        'num_classes': len(classes), 'num_sc': len(sc_classes), 'num_ns_pairs': len(ns_pairs)
    }

# ═══════════════════════════════════════════════════════════════
# RUN FOR n=3..7 AND PRINT RESULTS
# ═══════════════════════════════════════════════════════════════

results = {}
for n in range(3, 8):
    print(f"\nAnalyzing n={n}...", flush=True)
    t0 = time.time()
    data = analyze_n(n)
    results[n] = data
    elapsed = time.time() - t0
    print(f"  Done in {elapsed:.1f}s (orbit:{data['t_orbit']:.1f}s, H:{data['t_h']:.1f}s)", flush=True)
    print(f"  Classes: {data['num_classes']}, SC: {data['num_sc']}, NS-pairs: {data['num_ns_pairs']}")

print("\n" + "="*70)
print("FULL RESULTS")
print("="*70)

for n in range(3, 8):
    d = results[n]
    sc = d['sc_data']; ns = d['ns_data']
    pb = [x for x in sc if x['cat']=='pure_blue']
    mx = [x for x in sc if x['cat']=='mixed']

    pb_vals = sorted(x['blue'] for x in pb)
    mx_blue = sorted(x['blue'] for x in mx)
    mx_black = sorted(x['black'] for x in mx)
    ns_black = sorted(x['black'] for x in ns)

    print(f"\n{'='*70}")
    print(f"n={n}  m={d['m']}  2^m={d['total_t']}  free_bits={d['free_bits']}")
    print(f"  blue={d['total_blue']} black={d['total_black']}")
    print(f"  Classes={d['num_classes']} SC={d['num_sc']} NS-pairs={d['num_ns_pairs']}")
    print(f"{'='*70}")

    print(f"\nFOUR CATEGORIES:")
    print(f"  Pure blue:  {pb_vals}  sum={sum(pb_vals)}")
    print(f"  Mixed blue: {mx_blue}  sum={sum(mx_blue)}")
    print(f"  Mixed blk:  {mx_black}  sum={sum(mx_black)}")
    print(f"  Pure blk:   {ns_black[:30]}{'...' if len(ns_black)>30 else ''}  sum={sum(ns_black)}")
    b_chk = '✓' if sum(pb_vals)+sum(mx_blue)==d['total_blue'] else '✗'
    bk_chk = '✓' if sum(mx_black)+sum(ns_black)==d['total_black'] else '✗'
    print(f"  Check: blue {sum(pb_vals)}+{sum(mx_blue)}={sum(pb_vals)+sum(mx_blue)}/{d['total_blue']}{b_chk}  "
          f"black {sum(mx_black)}+{sum(ns_black)}={sum(mx_black)+sum(ns_black)}/{d['total_black']}{bk_chk}")

    print(f"\nMIXED PAIRS (black, blue) [{len(mx)} total]:")
    mx_pairs = sorted([(x['black'],x['blue'],x['H'],x['au'],x['nt']) for x in mx])
    for bk,bl,H,au,nt in mx_pairs[:30]:
        print(f"  ({bk:>5},{bl:>3})  #t={nt:>5}  H={H:>7}  |Aut|={au}")
    if len(mx_pairs) > 30:
        print(f"  ... ({len(mx_pairs)} total)")

    print(f"\nPURE BLUE CLASSES [{len(pb)}]:")
    for x in sorted(pb, key=lambda x:x['blue']):
        print(f"  ci{x['ci']}: blue={x['blue']}, H={x['H']}, |Aut|={x['au']}, #t={x['nt']}")

    print(f"\nPURE BLACK (NS pairs) [{len(ns)}]:")
    for x in sorted(ns, key=lambda x:x['black'])[:20]:
        print(f"  [{x['a']},{x['b']}]: black={x['black']}, H_a={x['H_a']}(au={x['au_a']}), H_b={x['H_b']}(au={x['au_b']})")
    if len(ns) > 20:
        print(f"  ... ({len(ns)} total)")

    print(f"\nH-DISTRIBUTION BY CATEGORY:")
    H_pb = Counter(x['H'] for x in pb)
    H_mx = Counter(x['H'] for x in mx)
    H_ns = Counter()
    for x in ns: H_ns[x['H_a']] += 1; H_ns[x['H_b']] += 1
    all_H = sorted(set(list(H_pb)+list(H_mx)+list(H_ns)))
    print(f"  {'H':>7} | {'pure_blue':>9} | {'mixed':>6} | {'NS'}")
    for H in all_H:
        pb_c = H_pb.get(H,0); mx_c = H_mx.get(H,0); ns_c = H_ns.get(H,0)
        if pb_c+mx_c+ns_c > 0:
            print(f"  {H:>7} | {pb_c:>9} | {mx_c:>6} | {ns_c}")

# ═══════════════════════════════════════════════════════════════
# CROSS-n PATTERNS
# ═══════════════════════════════════════════════════════════════
print("\n" + "="*70)
print("CROSS-n PATTERNS")
print("="*70)

print("\n1. EVOLUTION TABLE:")
print(f"{'n':>3} | {'m':>3} | {'2^m':>7} | {'free_b':>7} | {'blue':>6} | {'#SC':>5} | {'#NS-pr':>7} | {'#mixed':>7} | {'pure_B':>7} | {'pure_bk':>8}")
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    pb = [x for x in sc if x['cat']=='pure_blue']
    mx = [x for x in sc if x['cat']=='mixed']
    print(f"{n:>3} | {d['m']:>3} | {d['total_t']:>7} | {d['free_bits']:>7} | {d['total_blue']:>6} | "
          f"{d['num_sc']:>5} | {d['num_ns_pairs']:>7} | {len(mx):>7} | {len(pb):>7} | {len(d['ns_data']):>8}")

print("\n2. PURE BLUE VALUES:")
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    pb = sorted(x['blue'] for x in sc if x['cat']=='pure_blue')
    H_pb = sorted(x['H'] for x in sc if x['cat']=='pure_blue')
    print(f"  n={n}: vals={pb} sum={sum(pb)} H={H_pb}")

print("\n3. MIXED BLUE VALUES (sorted):")
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    mx_b = sorted(x['blue'] for x in sc if x['cat']=='mixed')
    print(f"  n={n}: {mx_b} sum={sum(mx_b)}")

print("\n4. MIXED BLACK VALUES (sorted):")
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    mx_bk = sorted(x['black'] for x in sc if x['cat']=='mixed')
    print(f"  n={n}: {mx_bk} sum={sum(mx_bk)}")

print("\n5. MIXED PAIR RATIOS (bk/bl) sorted:")
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    mx = sorted([x for x in sc if x['cat']=='mixed'], key=lambda y:(y['black'],y['blue']))
    ratios = [round(x['black']/x['blue'],3) for x in mx]
    print(f"  n={n}: {ratios[:20]}{'...' if len(ratios)>20 else ''}")

print("\n6. PURE BLACK DISTRIBUTION (histogram of values):")
for n in range(3,8):
    d = results[n]
    ns_bk = sorted(x['black'] for x in d['ns_data'])
    c = Counter(ns_bk)
    print(f"  n={n}: {len(ns_bk)} NS-pairs, top counts: {sorted(c.items())[:10]}")

print("\n7. UNIVERSAL (H, |Aut|) → (bk, bl) MAP (across n):")
# Build dictionary of (H, au) → set of (bk, bl) across all n
pair_map = {}
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    for x in sc:
        key = (x['H'], x['au'])
        val = (x['black'], x['blue'])
        if key not in pair_map: pair_map[key] = set()
        pair_map[key].add(val)
print(f"  Entries with unique (bk,bl): {sum(1 for v in pair_map.values() if len(v)==1)}/{len(pair_map)}")
print(f"  Entries with MULTIPLE (bk,bl):")
for key, vals in sorted(pair_map.items()):
    if len(vals) > 1:
        print(f"    (H={key[0]}, |Aut|={key[1]}) → {sorted(vals)}")
print(f"  Universal entries (sample):")
univ = [(k,v) for k,v in sorted(pair_map.items()) if len(v)==1]
for k,v in univ[:15]:
    print(f"    (H={k[0]:>7}, |Aut|={k[1]:>4}) → (bk={list(v)[0][0]:>5}, bl={list(v)[0][1]:>3})")
if len(univ) > 15:
    print(f"    ... ({len(univ)} universal entries total)")
