#!/usr/bin/env python3
"""
explorer_exact_s20de.py — EXACT replication of tournament-tiling-explorer logic
kind-pasteur-2026-03-23-S20de

This script replicates the EXACT JavaScript logic of the explorer:
- TILES in explorer order: (x,y) with y from 1 to N-2, x from N down to y+2
- VERTS = [N, N-1, ..., 1] (descending)
- Base path: A[k][k+1]=1 (vertex N → N-1 → ... → 1)
- bits[i]=0: x→y (forward). bits[i]=1: y→x (backward)
- All off = transitive
- TRANS_MAP: (x,y) → (N-y+1, N-x+1) [grid reflection = reversal anti-aut]
- isGridSym: bits[i] == bits[TRANS_MAP[i]] for all non-fixed tiles
- flipMask: invert ALL bits (complement tiling)
- BLUE line: BOTH endpoints grid-sym. BLACK: at least one not grid-sym.
- Since isGridSym(flip(t)) == isGridSym(t), each line is blue iff t is grid-sym.

KEY FACT: isGridSym(flip(bits)) = isGridSym(bits) always.
Proof: flip swaps 0↔1, but the grid-sym check only cares about equality
of paired positions, which is preserved under uniform bit-flip.
So: BLUE line ⟺ source tiling is grid-symmetric. Period.
"""

import sys
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EXACT EXPLORER REPLICATION")
print("  kind-pasteur-2026-03-23-S20de")
print("=" * 80)

def analyze(n):
    t0 = time.time()

    # === EXACT EXPLORER CONSTRUCTION ===
    VERTS = list(range(n, 0, -1))  # [n, n-1, ..., 1]
    TILES = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            TILES.append((x, y))
    m = len(TILES)

    tk = {}
    for i, (x, y) in enumerate(TILES):
        tk[f"{x},{y}"] = i
    TRANS_MAP = [tk.get(f"{n-y+1},{n-x+1}", -1) for x, y in TILES]

    def isGridSym(bits):
        for i in range(m):
            ti = TRANS_MAP[i]
            if ti != i and ti >= 0 and bits[i] != bits[ti]:
                return False
        return True

    def bitsToAdj(bits):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1):
            A[k][k+1] = 1
        for i in range(m):
            xL, yL = TILES[i]
            xi = VERTS.index(xL)
            yi = VERTS.index(yL)
            if bits[i] == 0:
                A[xi][yi] = 1
            else:
                A[yi][xi] = 1
        return A

    perms = list(permutations(range(n)))
    def canonicalize(A):
        best = None
        for p in perms:
            s = ''.join(str(A[p[i]][p[j]]) for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    # === ENUMERATE ALL TILINGS ===
    class_data = defaultdict(lambda: {
        'tilings': [], 'grid_sym_count': 0,
        'blue_lines': 0, 'black_lines': 0,
        'flip_targets': set()
    })
    tiling_info = {}

    for mask in range(1 << m):
        bits = [(mask >> k) & 1 for k in range(m)]
        A = bitsToAdj(bits)
        cn = canonicalize(A)
        gs = isGridSym(bits)
        flip_mask = sum((1 - b) << k for k, b in enumerate(bits))

        tiling_info[mask] = {'canon': cn, 'grid_sym': gs, 'flip_mask': flip_mask}
        class_data[cn]['tilings'].append(mask)
        if gs:
            class_data[cn]['grid_sym_count'] += 1

        if (mask + 1) % 5000 == 0:
            print(f"    {mask+1}/{1<<m} ({time.time()-t0:.0f}s)")

    # === COMPUTE FLIP-PAIR LINES ===
    # Each tiling pairs with its complement. We count per-class.
    drawn = set()
    total_blue_lines = 0
    total_black_lines = 0

    for mask in range(1 << m):
        ti = tiling_info[mask]
        fm = ti['flip_mask']
        pair_key = (min(mask, fm), max(mask, fm))
        if pair_key in drawn:
            continue
        drawn.add(pair_key)

        # Color: BLUE if source is grid-sym (which equals: both are grid-sym)
        is_blue = ti['grid_sym']
        if is_blue:
            total_blue_lines += 1
        else:
            total_black_lines += 1

        # Which class does flip partner belong to?
        src_class = ti['canon']
        tgt_class = tiling_info[fm]['canon']
        class_data[src_class]['flip_targets'].add(tgt_class)
        class_data[tgt_class]['flip_targets'].add(src_class)

        if is_blue:
            class_data[src_class]['blue_lines'] += 1
            if src_class != tgt_class:
                class_data[tgt_class]['blue_lines'] += 1
        else:
            class_data[src_class]['black_lines'] += 1
            if src_class != tgt_class:
                class_data[tgt_class]['black_lines'] += 1

    # === CLASSIFY CLASSES ===
    # Transpose (red) pairs
    class_list = sorted(class_data.keys())
    class_idx = {cn: i for i, cn in enumerate(class_list)}

    # Transpose: apply TRANS_MAP to each tiling
    trans_class = {}
    for cn in class_list:
        rep_mask = class_data[cn]['tilings'][0]
        rep_bits = [(rep_mask >> k) & 1 for k in range(m)]
        trans_bits = [0] * m
        for i in range(m):
            trans_bits[TRANS_MAP[i]] = rep_bits[i]
        trans_mask = sum(b << k for k, b in enumerate(trans_bits))
        trans_cn = tiling_info[trans_mask]['canon']
        trans_class[cn] = trans_cn

    # SC = class whose transpose is itself
    for cn in class_list:
        class_data[cn]['is_transpose_self'] = (trans_class[cn] == cn)

    # Merged classes: merge transpose pairs
    merge = {}
    mid = 0
    for cn in class_list:
        if cn in merge:
            continue
        merge[cn] = mid
        tc = trans_class[cn]
        if tc != cn and tc not in merge:
            merge[tc] = mid
        mid += 1
    V_merged = mid

    # === RESULTS ===
    num_classes = len(class_list)
    print(f"\n  n={n}: {num_classes} iso classes, {V_merged} merged classes")
    print(f"  {1<<m} tilings, {m} tiles")
    print(f"  Total flip-pair lines: {total_blue_lines + total_black_lines}")
    print(f"  BLUE lines: {total_blue_lines}")
    print(f"  BLACK lines: {total_black_lines}")
    print(f"  Grid-sym tilings: {sum(cd['grid_sym_count'] for cd in class_data.values())}/{1<<m}")

    # Per-class breakdown
    pure_blue = []
    pure_black = []
    mixed = []

    print(f"\n  PER-CLASS BREAKDOWN:")
    print(f"  {'idx':>4} {'fiber':>5} {'gs':>4} {'blue':>5} {'black':>5} "
          f"{'ts':>2} {'type':>12} {'merged':>4}")

    for cn in class_list:
        cd = class_data[cn]
        fiber = len(cd['tilings'])
        gs = cd['grid_sym_count']
        bl = cd['blue_lines']
        bk = cd['black_lines']
        ts = 'Y' if cd['is_transpose_self'] else 'N'
        mi = merge[cn]

        if gs == fiber:
            ctype = 'PURE BLUE'
            pure_blue.append(cn)
        elif gs == 0:
            ctype = 'PURE BLACK'
            pure_black.append(cn)
        else:
            ctype = f'MIXED {gs}/{fiber}'
            mixed.append(cn)

        print(f"  {class_idx[cn]:4d} {fiber:5d} {gs:4d} {bl:5d} {bk:5d} "
              f"{ts:>2} {ctype:>12} {mi:4d}")

    print(f"\n  SUMMARY:")
    print(f"    PURE BLUE classes: {len(pure_blue)}")
    print(f"    PURE BLACK classes: {len(pure_black)}")
    print(f"    MIXED classes: {len(mixed)}")
    print(f"    Transpose-self (SC analog): {sum(1 for cn in class_list if class_data[cn]['is_transpose_self'])}")

    # Check: is pure-blue ⊆ transpose-self?
    pb_ts = sum(1 for cn in pure_blue if class_data[cn]['is_transpose_self'])
    pb_nts = len(pure_blue) - pb_ts
    print(f"    PURE BLUE that are transpose-self: {pb_ts}")
    print(f"    PURE BLUE that are NOT transpose-self: {pb_nts}")

    # Red pairs (transpose pairs that are different classes)
    red_pairs = sum(1 for cn in class_list
                    if trans_class[cn] != cn and class_idx[cn] < class_idx[trans_class[cn]])
    print(f"    RED pairs (transpose, different class): {red_pairs}")

    print(f"    Time: {time.time()-t0:.1f}s")
    return {
        'n': n, 'classes': num_classes, 'merged': V_merged,
        'pure_blue': len(pure_blue), 'pure_black': len(pure_black),
        'mixed': len(mixed), 'blue_lines': total_blue_lines,
        'black_lines': total_black_lines
    }


# === MAIN ===
results = {}
for n in [3, 4, 5, 6, 7]:
    print(f"\n{'#'*70}")
    results[n] = analyze(n)

# === CROSS-n TABLE ===
print(f"\n\n{'='*70}")
print(f"  CROSS-n SUMMARY")
print(f"{'='*70}")
print(f"  {'n':>3} {'classes':>8} {'merged':>8} {'PB':>4} {'PBk':>4} {'Mix':>4} "
      f"{'blue_L':>7} {'black_L':>7} {'gs_frac':>8}")
for n in [3, 4, 5, 6, 7]:
    r = results[n]
    m = r['classes']  # not tiles
    gs_frac = r['blue_lines'] / (r['blue_lines'] + r['black_lines']) if r['blue_lines']+r['black_lines'] > 0 else 0
    print(f"  {n:3d} {r['classes']:8d} {r['merged']:8d} {r['pure_blue']:4d} {r['pure_black']:4d} "
          f"{r['mixed']:4d} {r['blue_lines']:7d} {r['black_lines']:7d} {gs_frac:8.4f}")

print(f"\n  PURE BLUE class count sequence: {[results[n]['pure_blue'] for n in [3,4,5,6,7]]}")
print(f"  PURE BLACK class count sequence: {[results[n]['pure_black'] for n in [3,4,5,6,7]]}")
print(f"  MIXED class count sequence: {[results[n]['mixed'] for n in [3,4,5,6,7]]}")

print(f"\n  DONE.")
print("=" * 80)
