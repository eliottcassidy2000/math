"""
scsf_small_n_analysis.py — Detailed SCSF kernel analysis at n=5,6,7,8 (small n).

For each n, compute:
- All SCSF kernel classes
- Their bl (GS tiling count), |Aut|, #t (total tilings)
- Whether blueself (GS flip-pair) or blackself (non-GS flip-pair)
- Score sequences
- Self-flip pair structure

oracle-2026-05-11
"""

import numpy as np
import sys, time
from itertools import permutations
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')

def analyze_scsf(n, verbose=False):
    t_start = time.time()
    tiles = [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
    m = len(tiles)
    verts = list(range(n, 0, -1))
    xi_arr = np.array([verts.index(x) for x,y in tiles], dtype=np.int32)
    yi_arr = np.array([verts.index(y) for x,y in tiles], dtype=np.int32)
    tile_idx = {t:i for i,t in enumerate(tiles)}
    tm = [tile_idx[(n-y+1, n-x+1)] for x,y in tiles]
    gs_pairs_check = [(i, tm[i]) for i in range(m) if tm[i] > i]
    self_paired = [i for i in range(m) if tm[i] == i]
    free_bits = len(self_paired) + len(gs_pairs_check)

    NT = 1 << m
    k_range = np.arange(m, dtype=np.int32)
    masks_all = np.arange(NT, dtype=np.uint32)
    flip_all = np.uint32(NT - 1)

    if verbose:
        print(f"  n={n}, m={m}, NT={NT}, free_bits={free_bits}", flush=True)

    # Build all adjacency matrices
    bits_all = ((masks_all[:, None] >> k_range[None, :]) & 1).astype(np.uint8)
    A_all = np.zeros((NT, n*n), dtype=np.uint8)
    for k in range(n-1):
        A_all[:, k*n+(k+1)] = 1
    for k in range(m):
        xk, yk = int(xi_arr[k]), int(yi_arr[k])
        A_all[:, xk*n+yk] += (1 - bits_all[:, k])
        A_all[:, yk*n+xk] += bits_all[:, k]

    # Complement
    trans_idx = np.array([j*n+i for i in range(n) for j in range(n)], dtype=np.int32)
    comp_A_all = A_all[:, trans_idx]

    # Flip (T̃): reverse all tile arcs, keep base path
    flip_masks = masks_all ^ flip_all
    A_flip = np.zeros((NT, n*n), dtype=np.uint8)
    flip_bits = 1 - bits_all
    for k in range(n-1):
        A_flip[:, k*n+(k+1)] = 1
    for k in range(m):
        xk, yk = int(xi_arr[k]), int(yi_arr[k])
        A_flip[:, xk*n+yk] += (1 - flip_bits[:, k])
        A_flip[:, yk*n+xk] += flip_bits[:, k]

    # Permutation indices
    perms_list = list(permutations(range(n)))
    NP = len(perms_list)
    perm_flat = np.zeros((NP, n*n), dtype=np.int32)
    for p_idx, sigma in enumerate(perms_list):
        sigma_arr = np.array(sigma, dtype=np.int32)
        for i in range(n):
            perm_flat[p_idx, i*n:(i+1)*n] = sigma_arr[i]*n + sigma_arr

    def compute_canonical(src):
        canon = np.full((NT, n*n), 255, dtype=np.uint8)
        for p_idx in range(NP):
            perm_A = src[:, perm_flat[p_idx]]
            diff = perm_A != canon
            has_diff = diff.any(axis=1)
            if not has_diff.any(): continue
            first = np.argmax(diff, axis=1)
            ridx = np.arange(NT)
            is_less = (perm_A[ridx, first] < canon[ridx, first]) & has_diff
            canon[is_less] = perm_A[is_less]
        return canon

    if verbose: print(f"  Canonicalizing T, T^op, T̃...", flush=True)
    canon = compute_canonical(A_all)
    comp_canon = compute_canonical(comp_A_all)
    flip_canon = compute_canonical(A_flip)

    # Group into iso classes
    groups = defaultdict(list)
    for mask in range(NT):
        groups[canon[mask].tobytes()].append(mask)

    sig_to_ci = {}
    classes = []
    for ci, (sig, members) in enumerate(groups.items()):
        rep = members[0]
        sc = np.array_equal(canon[rep], comp_canon[rep])
        # GS: mask is GS iff all GS-pair bits match
        gs = sum(1 for mk in members
                 if all(((mk>>i)&1)==((mk>>j)&1) for i,j in gs_pairs_check))
        # SF: any member has flip in same class
        sf = False
        for mk in members:
            flip_mk = int(mk ^ flip_all)
            if flip_mk < NT:
                fl_sig = canon[flip_mk].tobytes()
                if fl_sig == sig:
                    sf = True; break
        scores = sorted(A_all[rep].reshape(n,n).sum(axis=1).tolist(), reverse=True)
        classes.append({'ci':ci, 'sig':sig, 'rep':rep, 'members':members,
                        'nt':len(members), 'sc':sc, 'gs':gs, 'sf':sf,
                        'scores':tuple(scores)})
        sig_to_ci[sig] = ci

    # |Aut| via perm count
    A_reps = np.array([A_all[c['rep']] for c in classes], dtype=np.uint8)
    aut_counts = np.zeros(len(classes), dtype=np.int32)
    for p_idx in range(NP):
        perm_A = A_reps[:, perm_flat[p_idx]]
        aut_counts += (perm_A == A_reps).all(axis=1).astype(np.int32)
    for ci, c in enumerate(classes):
        c['au'] = int(aut_counts[ci])

    # Find SCSF kernel classes and their self-flip pair structure
    scsf = [c for c in classes if c['sc'] and c['sf']]

    for c in scsf:
        # Find self-flip pairs (pairs {M, flip(M)} with both in class)
        in_class = set(c['members'])
        flip_pairs = []
        seen = set()
        for mk in c['members']:
            if mk in seen: continue
            fmk = int(mk ^ flip_all)
            if fmk in in_class and fmk != mk:
                flip_pairs.append((mk, fmk))
                seen.add(mk); seen.add(fmk)
            elif fmk == mk:
                flip_pairs.append((mk, mk))  # self-flip (only if all bits flip to same: impossible)
                seen.add(mk)
        # Classify as blueself/blackself
        gs_set = set(mk for mk in c['members']
                     if all(((mk>>i)&1)==((mk>>j)&1) for i,j in gs_pairs_check))
        blueself_pairs = [(mk,fmk) for mk,fmk in flip_pairs
                          if mk in gs_set or fmk in gs_set]
        blackself_pairs = [(mk,fmk) for mk,fmk in flip_pairs
                           if mk not in gs_set and fmk not in gs_set]
        c['flip_pairs'] = flip_pairs
        c['blueself_pairs'] = blueself_pairs
        c['blackself_pairs'] = blackself_pairs
        c['type'] = 'BLUESELF' if len(blueself_pairs) > 0 else 'BLACKSELF'

    elapsed = time.time() - t_start
    return {
        'n': n, 'm': m, 'NT': NT, 'free_bits': free_bits,
        'classes': classes, 'scsf': scsf, 'elapsed': elapsed,
        'num_sc': sum(1 for c in classes if c['sc']),
        'num_sf': sum(1 for c in classes if c['sf']),
        'num_scsf': len(scsf),
    }

# ─── MAIN ─────────────────────────────────────────────────────────────────────
print("SCSF KERNEL ANALYSIS: n=4..7")
print("="*70)

all_data = {}
for n in [4, 5, 6, 7]:
    print(f"\nAnalyzing n={n}...", flush=True)
    data = analyze_scsf(n, verbose=True)
    all_data[n] = data
    print(f"  n={n}: {data['num_sc']} SC, {data['num_sf']} SF, {data['num_scsf']} SCSF  ({data['elapsed']:.1f}s)", flush=True)

print("\n" + "="*70)
print("KERNEL CLASSES BY n")
print("="*70)

for n in [4, 5, 6, 7]:
    data = all_data[n]
    scsf = data['scsf']
    print(f"\nn={n}: {len(scsf)} kernel classes")
    print(f"  {'#':>3} | {'type':>9} | {'bl':>4} | {'bsp':>4} | {'bksp':>5} | {'|Aut|':>6} | {'#t':>6} | scores")
    print("  " + "-"*65)
    for i, c in enumerate(scsf):
        bsp = len(c['blueself_pairs'])
        bksp = len(c['blackself_pairs'])
        print(f"  {i+1:>3} | {c['type']:>9} | {c['gs']:>4} | {bsp:>4} | {bksp:>5} | {c['au']:>6} | {c['nt']:>6} | {c['scores']}")

print("\n" + "="*70)
print("PATTERN SUMMARY")
print("="*70)
print(f"{'n':>3} | {'#SCSF':>6} | {'#blue':>6} | {'#black':>7} | {'bl_values'}")
print("-"*60)
for n in [4, 5, 6, 7]:
    scsf = all_data[n]['scsf']
    blue = [c for c in scsf if c['type']=='BLUESELF']
    black = [c for c in scsf if c['type']=='BLACKSELF']
    bl_vals = sorted([c['gs'] for c in scsf])
    print(f"{n:>3} | {len(scsf):>6} | {len(blue):>6} | {len(black):>7} | {bl_vals}")

print("\n" + "="*70)
print("SELF-FLIP PAIR STRUCTURE")
print("="*70)
for n in [4, 5, 6, 7]:
    data = all_data[n]
    scsf = data['scsf']
    flip_all_n = (1 << data['m']) - 1
    print(f"\nn={n}:")
    for c in scsf:
        print(f"  Class (bl={c['gs']}, {c['type']}, #fp={len(c['flip_pairs'])}):")
        for mk, fmk in c['flip_pairs'][:5]:
            is_gs_mk = all(((mk>>i)&1)==((mk>>j)&1) for i,j in
                           [(i,tm_j) for i,tm_j in
                            [(i,data['classes'][0]['ci']) for i in range(data['m'])]])
            # Actually just check membership directly
            gs_set_c = set(mk2 for mk2 in c['members']
                           if all(((mk2>>i)&1)==((mk2>>j)&1) for i,j in
                                  [(i, [tile_idx2[(n-y+1,n-x+1)] for x,y in
                                   [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]][i])
                                   for i in range(data['m'])
                                   if [tile_idx2[(n-y+1,n-x+1)] for x,y in
                                       [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]][i] > i]))
            pass
        # Simpler:
        tiles_n = [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
        tile_idx2 = {t:i for i,t in enumerate(tiles_n)}
        tm_n = [tile_idx2[(n-y+1,n-x+1)] for x,y in tiles_n]
        gp_n = [(i,tm_n[i]) for i in range(len(tiles_n)) if tm_n[i]>i]
        def is_gs_n(mk):
            return all(((mk>>i)&1)==((mk>>j)&1) for i,j in gp_n)
        for mk, fmk in c['flip_pairs'][:8]:
            gs1 = is_gs_n(mk); gs2 = is_gs_n(fmk)
            label = f"GS↔GS" if gs1 and gs2 else ("GS↔nonGS" if gs1 or gs2 else "nGS↔nGS")
            print(f"    {bin(mk)[2:].zfill(data['m'])} ↔ {bin(fmk)[2:].zfill(data['m'])}  [{label}]")

print(f"\nDone.")
