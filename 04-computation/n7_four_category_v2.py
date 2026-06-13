"""
n7_four_category_v2.py — Correct four-category analysis for n=3..7.

ALGORITHM: Standard canonicalization using all n! vertex permutations,
but implemented with numpy for speed.

For each of the 2^m tilings, find the canonical form = lexicographically
minimum adjacency matrix under all n! relabelings. Group by canonical form
to get iso classes.

KEY CORRECTION from v1: permuting a tournament T^σ changes ALL arc directions
(including base-path arcs), so we cannot use an orbit shortcut — we must
canonicalize the full adjacency matrix.

NUMPY SPEEDUP:
  - A_all: (2^m, n²) uint8 matrix of all tiling adjacencies
  - For each of n! permutations p: compute perm_A = A_all[:, perm_indices[p]]
  - Update canonical = lexicographic min over all permutations
  - n=7: 5040 × 32768 × 49 = 8B operations, numpy speed ~40 seconds

oracle-2026-05-11
"""

import sys, time
import numpy as np
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(encoding='utf-8')

def make_tiles(n):
    return [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]

def make_trans_map(n, tiles):
    idx = {t:i for i,t in enumerate(tiles)}
    return [idx[(n-y+1,n-x+1)] for x,y in tiles]

def H_count(A_flat, n):
    """Count Hamiltonian paths using DP on subsets."""
    A = [[A_flat[i*n+j] for j in range(n)] for i in range(n)]
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

def analyze_n(n, verbose=True):
    t_start = time.time()
    tiles = make_tiles(n)
    m = len(tiles)
    tm = make_trans_map(n, tiles)
    verts = list(range(n, 0, -1))
    xi_arr = [verts.index(x) for x,y in tiles]
    yi_arr = [verts.index(y) for x,y in tiles]
    perms_list = list(permutations(range(n)))
    np_count = len(perms_list)

    # Free bits (grid-symmetric tilings = 2^free_bits)
    self_paired = sum(1 for i,ti in enumerate(tm) if ti==i)
    free_bits = self_paired + (m - self_paired)//2

    # ── Build A_all: (2^m, n²) uint8 ──────────────────────────────────────
    num_t = 1 << m
    A_all = np.zeros((num_t, n*n), dtype=np.uint8)
    for mask in range(num_t):
        for k in range(n-1): A_all[mask, k*n+k+1] = 1  # base path
        for k in range(m):
            if (mask>>k)&1 == 0: A_all[mask, xi_arr[k]*n+yi_arr[k]] = 1
            else: A_all[mask, yi_arr[k]*n+xi_arr[k]] = 1

    # Complement: A^op[i][j] = A[j][i], i.e., transpose off-diagonal
    trans_idx = np.array([j*n+i for i in range(n) for j in range(n)], dtype=np.int32)
    comp_A_all = A_all[:, trans_idx]

    # ── Permutation flat indices: for each perm p, perm_idx[i*n+j] = p[i]*n+p[j] ──
    perm_flat = np.zeros((np_count, n*n), dtype=np.int32)
    for p_idx, sigma in enumerate(perms_list):
        for i in range(n):
            for j in range(n):
                perm_flat[p_idx, i*n+j] = sigma[i]*n + sigma[j]

    if verbose:
        print(f"  Setup done in {time.time()-t_start:.2f}s  "
              f"(n={n}, m={m}, 2^m={num_t}, n!={np_count})", flush=True)

    # ── Compute canonical and comp_canonical arrays ───────────────────────
    def compute_canonical(src, label=""):
        """src: (num_t, n²) uint8. Returns (num_t, n²) uint8 canonical forms."""
        canon = np.full((num_t, n*n), 255, dtype=np.uint8)
        t0 = time.time()
        for p_idx in range(np_count):
            perm_A = src[:, perm_flat[p_idx]]      # (num_t, n²) fancy indexing
            # Lexicographic compare: update where perm_A < canon
            diff = perm_A != canon                  # (num_t, n²) bool
            has_diff = diff.any(axis=1)             # (num_t,) bool
            first = np.argmax(diff, axis=1)         # (num_t,) int
            ridx = np.arange(num_t)
            is_less = (perm_A[ridx, first] < canon[ridx, first]) & has_diff
            canon[is_less] = perm_A[is_less]
            if verbose and (p_idx+1) % 1000 == 0:
                elapsed = time.time()-t0
                print(f"    {label}: {p_idx+1}/{np_count} perms  {elapsed:.1f}s", flush=True)
        return canon

    t0 = time.time()
    if verbose: print(f"  Computing canonical...", flush=True)
    canonical = compute_canonical(A_all, "canon")
    if verbose:
        print(f"  Canonical done in {time.time()-t0:.1f}s", flush=True)
        print(f"  Computing comp_canonical...", flush=True)
    t0 = time.time()
    comp_canonical = compute_canonical(comp_A_all, "comp")
    if verbose:
        print(f"  Comp_canonical done in {time.time()-t0:.1f}s", flush=True)

    # ── Group by canonical form → iso classes ─────────────────────────────
    groups = defaultdict(list)
    for mask in range(num_t):
        groups[canonical[mask].tobytes()].append(mask)

    # ── For each class: compute data ──────────────────────────────────────
    # SC check: class is SC iff canonical[rep] == comp_canonical[rep]
    # Grid-sym count: for each mask in class, check if grid-symmetric
    gs_check_pairs = [(i, tm[i]) for i in range(m) if tm[i] > i]  # pairs (i,j), i<j

    def is_gs(mask):
        return all(((mask>>i)&1) == ((mask>>j)&1) for i,j in gs_check_pairs)

    all_gs = set(mask for mask in range(num_t) if is_gs(mask))

    classes = []
    sig_to_ci = {}
    for ci, (sig_bytes, members) in enumerate(sorted(groups.items())):
        rep = members[0]
        A_flat_rep = A_all[rep].tolist()
        H = H_count(A_flat_rep, n)
        # SC: canonical[rep] == comp_canonical[rep]
        sc = np.array_equal(canonical[rep], comp_canonical[rep])
        # comp sig
        comp_sig = comp_canonical[rep].tobytes()
        # grid-sym count
        gs_count = sum(1 for mk in members if mk in all_gs)
        aut_size = np_count * len(members) // num_t  # approx, use orbit-stab
        # Better: |Aut| = n! / orbit_size... but orbit ≠ class_size in tiling space!
        # Actually orbit = class_size in tiling space (each iso class = orbit of first representative)
        # Wait: #tilings in a class = H(C) / |Aut(C)| by orbit-stabilizer on tiling space.
        # We can verify: |Aut| from H and #tilings.
        # For now, compute |Aut| from perm check on representative.
        classes.append({
            'ci': ci, 'rep': rep, 'A_flat': A_flat_rep,
            'sig': sig_bytes, 'comp_sig': comp_sig,
            'members': members, 'nt': len(members),
            'H': H, 'sc': sc, 'gs': gs_count, 'au': None
        })
        sig_to_ci[sig_bytes] = ci

    # Compute |Aut| for each class via perm check on representative
    if verbose: print(f"  Computing |Aut| for {len(classes)} classes...", flush=True)
    t0 = time.time()
    A_reps = np.array([A_all[c['rep']] for c in classes], dtype=np.uint8)  # (#classes, n²)
    # For each perm p: check if A_reps[:, perm_flat[p]] == A_reps (automorphism test)
    aut_counts = np.zeros(len(classes), dtype=np.int32)
    for p_idx in range(np_count):
        perm_A = A_reps[:, perm_flat[p_idx]]  # (#classes, n²)
        is_aut = (perm_A == A_reps).all(axis=1)  # (#classes,) bool
        aut_counts += is_aut.astype(np.int32)
    for i, c in enumerate(classes):
        c['au'] = int(aut_counts[i])
    if verbose:
        print(f"  |Aut| done in {time.time()-t0:.1f}s", flush=True)

    # ── SC/NS classification ──────────────────────────────────────────────
    sc_classes = [c for c in classes if c['sc']]
    ns_classes = [c for c in classes if not c['sc']]

    # NS pairs: partner = class whose sig = comp_sig
    seen_ns = set(); ns_pairs = []
    for c in ns_classes:
        if c['ci'] in seen_ns: continue
        partner_ci = sig_to_ci.get(c['comp_sig'])
        if partner_ci is None:
            print(f"  WARNING: NS partner not found for ci={c['ci']}")
            continue
        partner = classes[partner_ci]
        seen_ns.add(c['ci']); seen_ns.add(partner['ci'])
        ns_pairs.append((c, partner))

    # ── Four-category data ────────────────────────────────────────────────
    sc_data = []
    for c in sc_classes:
        blue = c['gs']; black = c['nt'] - blue
        cat = 'pure_blue' if black == 0 else 'mixed'
        sc_data.append({'ci':c['ci'],'H':c['H'],'au':c['au'],'nt':c['nt'],
                        'blue':blue,'black':black,'cat':cat})

    ns_data = []
    for ca, cb in ns_pairs:
        nt = ca['nt'] + cb['nt']
        ns_data.append({'a':ca['ci'],'b':cb['ci'],'nt':nt,'black':nt,
                        'H_a':ca['H'],'H_b':cb['H'],'au_a':ca['au'],'au_b':cb['au']})

    total_t = num_t; total_blue = 1<<free_bits; total_black = total_t - total_blue
    elapsed = time.time() - t_start

    return {
        'n': n, 'm': m, 'free_bits': free_bits,
        'total_t': total_t, 'total_blue': total_blue, 'total_black': total_black,
        'classes': classes, 'sc_data': sc_data, 'ns_data': ns_data,
        'ns_pairs': ns_pairs, 'sc_classes': sc_classes, 'ns_classes': ns_classes,
        'num_classes': len(classes), 'num_sc': len(sc_classes), 'num_ns_pairs': len(ns_pairs),
        'elapsed': elapsed
    }

# ═══════════════════════════════════════════════════════════════════════════
# RUN
# ═══════════════════════════════════════════════════════════════════════════

results = {}
for n in range(3, 8):
    print(f"\n{'='*60}\nAnalyzing n={n}...", flush=True)
    data = analyze_n(n, verbose=True)
    results[n] = data
    d = data
    print(f"\n  RESULT n={n}: {d['num_classes']} classes, {d['num_sc']} SC, {d['num_ns_pairs']} NS-pairs  (total {d['elapsed']:.1f}s)")
    # Quick check
    sc = d['sc_data']; ns = d['ns_data']
    pb = [x for x in sc if x['cat']=='pure_blue']
    mx = [x for x in sc if x['cat']=='mixed']
    pb_v = sorted(x['blue'] for x in pb)
    mx_b = sorted(x['blue'] for x in mx)
    mx_bk = sorted(x['black'] for x in mx)
    ns_bk = sorted(x['black'] for x in ns)
    print(f"  Pure blue:  {pb_v}  sum={sum(pb_v)}")
    print(f"  Mixed blue: {mx_b}  sum={sum(mx_b)}")
    print(f"  Mixed blk:  {mx_bk}  sum={sum(mx_bk)}")
    print(f"  Pure blk count: {len(ns_bk)}, sum={sum(ns_bk)}")
    chk_b = sum(pb_v)+sum(mx_b); chk_bk = sum(mx_bk)+sum(ns_bk)
    print(f"  Check: blue {chk_b}/{d['total_blue']}{'✓' if chk_b==d['total_blue'] else '✗'}  "
          f"black {chk_bk}/{d['total_black']}{'✓' if chk_bk==d['total_black'] else '✗'}")

# ═══════════════════════════════════════════════════════════════════════════
# DETAILED OUTPUT
# ═══════════════════════════════════════════════════════════════════════════

print("\n" + "="*70)
print("CROSS-n EVOLUTION TABLE")
print("="*70)
print(f"{'n':>3} | {'m':>3} | {'2^m':>7} | {'fb':>3} | {'blue':>6} | {'#cls':>5} | {'#SC':>5} | {'#NS-pr':>7} | {'#mix':>5} | {'#pb':>4} | {'#pblk':>6}")
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    pb = [x for x in sc if x['cat']=='pure_blue']
    mx = [x for x in sc if x['cat']=='mixed']
    print(f"{n:>3} | {d['m']:>3} | {d['total_t']:>7} | {d['free_bits']:>3} | {d['total_blue']:>6} | "
          f"{d['num_classes']:>5} | {d['num_sc']:>5} | {d['num_ns_pairs']:>7} | {len(mx):>5} | {len(pb):>4} | {len(d['ns_data']):>6}")

print("\n" + "="*70)
print("FOUR CATEGORIES AT EACH n")
print("="*70)
for n in range(3,8):
    d = results[n]
    sc = d['sc_data']; ns = d['ns_data']
    pb = [x for x in sc if x['cat']=='pure_blue']
    mx = [x for x in sc if x['cat']=='mixed']
    pb_v = sorted(x['blue'] for x in pb)
    mx_b = sorted(x['blue'] for x in mx)
    mx_bk = sorted(x['black'] for x in mx)
    ns_bk = sorted(x['black'] for x in ns)

    print(f"\nn={n}:")
    print(f"  Pure blue:  {pb_v}")
    print(f"  Mixed blue: {mx_b}")
    print(f"  Mixed blk:  {mx_bk}")
    print(f"  Pure blk:   {ns_bk[:30]}{'...' if len(ns_bk)>30 else ''}")

print("\n" + "="*70)
print("MIXED PAIRS DETAIL FOR EACH n")
print("="*70)
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    mx = sorted([x for x in sc if x['cat']=='mixed'], key=lambda y:(y['H'],y['au'],y['black']))
    print(f"\nn={n} ({len(mx)} mixed):")
    print(f"  {'H':>6} | {'|Aut|':>5} | {'#t':>5} | {'black':>6} | {'blue':>4} | bk/bl")
    for x in mx[:40]:
        ratio = x['black']/x['blue'] if x['blue'] else float('inf')
        print(f"  {x['H']:>6} | {x['au']:>5} | {x['nt']:>5} | {x['black']:>6} | {x['blue']:>4} | {ratio:.2f}")
    if len(mx) > 40: print(f"  ... ({len(mx)} total)")

print("\n" + "="*70)
print("PURE BLUE CLASSES AT EACH n")
print("="*70)
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    pb = sorted([x for x in sc if x['cat']=='pure_blue'], key=lambda x:x['H'])
    print(f"n={n}: {len(pb)} pure-blue classes")
    for x in pb:
        print(f"  H={x['H']:>6}, |Aut|={x['au']}, #t={x['nt']}, blue={x['blue']}")

print("\n" + "="*70)
print("PURE BLACK NS PAIRS (sample) AT EACH n")
print("="*70)
for n in range(3,8):
    d = results[n]; ns = sorted(d['ns_data'], key=lambda x:x['black'])
    print(f"n={n}: {len(ns)} NS-pairs")
    for x in ns[:15]:
        print(f"  H_a={x['H_a']:>6}(au={x['au_a']}), H_b={x['H_b']:>6}(au={x['au_b']}): black={x['black']:>7}")
    if len(ns)>15: print(f"  ... ({len(ns)} total)")

print("\n" + "="*70)
print("UNIVERSAL (H, |Aut|) → (bk, bl) MAP")
print("="*70)
pair_map = {}
for n in range(3,8):
    d = results[n]; sc = d['sc_data']
    for x in sc:
        key = (x['H'], x['au'])
        val = (x['black'], x['blue'])
        if key not in pair_map: pair_map[key] = {}
        if val not in pair_map[key]: pair_map[key][val] = []
        pair_map[key][val].append(n)

univ = {k:v for k,v in pair_map.items() if len(v)==1}
multi = {k:v for k,v in pair_map.items() if len(v)>1}
print(f"Universal: {len(univ)}/{len(pair_map)}  Multiple: {len(multi)}")
print("\nMultiple-mapping entries:")
for k, v in sorted(multi.items()):
    print(f"  (H={k[0]:>7}, |Aut|={k[1]:>4}) → {sorted(v.items())}")
print("\nUniversal entries (sample):")
for k, v in sorted(univ.items())[:20]:
    (bk,bl), ns = list(v.items())[0]
    print(f"  (H={k[0]:>7}, |Aut|={k[1]:>4}) → (bk={bk:>6}, bl={bl:>4}) at n={ns}")
