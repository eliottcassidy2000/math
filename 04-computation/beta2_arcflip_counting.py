#!/usr/bin/env python3
"""
beta2_arcflip_counting.py - Counting formula for delta_|A_p| under arc flip

For a tournament T on n vertices, flipping arc (u->v) to (v->u):
- Lost p-paths: paths in T using u->v
- Gained p-paths: paths in T' using v->u

Hypothesis: |gained_3| - |lost_3| = 2 * (|gained_2| - |lost_2|)

Why? A path using arc u->v at position i in a (p+1)-vertex path has
the arc at one of p positions. The positions before u and after v are
independent "tails" in the tournament minus {u,v}.

For a tournament (complete), every pair has an arc, so the number of
paths through (u,v) depends only on the "reach" from/to u,v.

Let's define:
  in(u) = {w : w->u}, out(u) = {w : u->w}
  in(v) = {w : w->v}, out(v) = {w : v->w}

After flip: u's out-degree decreases by 1, v's increases by 1.

The key sets (excluding u,v):
  A = out(u) \ {v} intersect in(v) \ {u}  — common neighbors w with u->w, w->v
  B = in(u) \ {v} intersect out(v) \ {u}  — common neighbors w with w->u, v->w
  C = out(u) \ {v} intersect out(v) \ {u} — w with u->w, v->w
  D = in(u) \ {v} intersect in(v) \ {u}   — w with w->u, w->v

Since tournament: every w in V\{u,v} has exactly one direction to u and one to v.
So A + C = out(u)\{v}, A + D = in(v)\{u}, B + D = in(u)\{v}, B + C = out(v)\{u}.
And |A| + |B| + |C| + |D| = n-2.

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import enumerate_allowed_paths
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def get_ABCD(A, n, u, v):
    """Compute the partition of V\{u,v} into A,B,C,D sets."""
    A_set, B_set, C_set, D_set = [], [], [], []
    for w in range(n):
        if w == u or w == v:
            continue
        u_to_w = A[u][w]
        w_to_v = A[w][v]
        if u_to_w and w_to_v:
            A_set.append(w)   # u->w, w->v
        elif not u_to_w and not w_to_v:
            B_set.append(w)   # w->u, v->w
        elif u_to_w and not w_to_v:
            C_set.append(w)   # u->w, v->w
        else:
            D_set.append(w)   # w->u, w->v
    return A_set, B_set, C_set, D_set

def count_paths_using_arc(allowed_paths, u, v):
    """Count allowed paths using arc u->v, categorized by position."""
    by_position = defaultdict(int)
    total = 0
    for p in allowed_paths:
        for i in range(len(p) - 1):
            if p[i] == u and p[i+1] == v:
                by_position[i] += 1
                total += 1
                break  # Each path counted once even if it uses arc multiple times
    return total, dict(by_position)


n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

print("=" * 70)
print(f"ARC-FLIP COUNTING ANALYSIS AT n={n}")
print("=" * 70)

# Precompute all allowed paths for all tournaments
all_allowed = {}
for bits in range(total):
    A = build_adj(n, bits)
    allowed = {}
    for p in range(5):
        allowed[p] = enumerate_allowed_paths(A, n, p)
        if not allowed[p]:
            break
    all_allowed[bits] = allowed

# For each tournament and each arc, compute |A|, |B|, |C|, |D|
# and the exact path counts
formula_check = []

for bits in range(total):
    A = build_adj(n, bits)
    allowed = all_allowed[bits]

    for u in range(n):
        for v in range(n):
            if u == v or A[u][v] == 0:
                continue

            a, b, c, d = get_ABCD(A, n, u, v)
            la, lb, lc, ld = len(a), len(b), len(c), len(d)

            # Count paths through u->v at each dimension
            lost2, lost2_pos = count_paths_using_arc(allowed.get(2, []), u, v)
            lost3, lost3_pos = count_paths_using_arc(allowed.get(3, []), u, v)

            # After flip: compute paths through v->u
            B_adj = [row[:] for row in A]
            B_adj[u][v] = 0
            B_adj[v][u] = 1
            bits_flip = 0
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if B_adj[i][j] == 1:
                        bits_flip |= (1 << idx)
                    idx += 1

            allowed_flip = all_allowed[bits_flip]
            gained2, gained2_pos = count_paths_using_arc(allowed_flip.get(2, []), v, u)
            gained3, gained3_pos = count_paths_using_arc(allowed_flip.get(3, []), v, u)

            delta2 = gained2 - lost2
            delta3 = gained3 - lost3

            formula_check.append({
                'la': la, 'lb': lb, 'lc': lc, 'ld': ld,
                'lost2': lost2, 'gained2': gained2, 'delta2': delta2,
                'lost3': lost3, 'gained3': gained3, 'delta3': delta3,
                'ratio': delta3 / delta2 if delta2 != 0 else float('inf'),
                'out_u': sum(A[u]) - A[u][v],  # out-degree excluding v
                'in_v': sum(A[w][v] for w in range(n)) - A[u][v],  # in-degree excluding u
            })

# Check the 2:1 ratio
ratios = Counter(round(d['ratio'], 4) for d in formula_check if d['delta2'] != 0)
print(f"\n  delta3/delta2 ratios: {dict(sorted(ratios.items()))}")

# How many have delta2=0?
zero_d2 = sum(1 for d in formula_check if d['delta2'] == 0)
print(f"  delta2=0 cases: {zero_d2}/{len(formula_check)}")
if zero_d2 > 0:
    # What delta3 when delta2=0?
    d3_when_d2_zero = Counter(d['delta3'] for d in formula_check if d['delta2'] == 0)
    print(f"  delta3 when delta2=0: {dict(sorted(d3_when_d2_zero.items()))}")

# Try to find a formula for lost2 and gained2 in terms of |A|, |B|, |C|, |D|
print(f"\n  Trying formula: lost2 = f(|A|, |B|, |C|, |D|)")

# Group by (la, lb, lc, ld) and check if lost2 is determined
by_abcd = defaultdict(list)
for d in formula_check:
    key = (d['la'], d['lb'], d['lc'], d['ld'])
    by_abcd[key].append(d)

print(f"\n  # distinct (a,b,c,d) patterns: {len(by_abcd)}")
for key in sorted(by_abcd.keys()):
    vals = by_abcd[key]
    lost2_vals = set(d['lost2'] for d in vals)
    gained2_vals = set(d['gained2'] for d in vals)
    lost3_vals = set(d['lost3'] for d in vals)
    gained3_vals = set(d['gained3'] for d in vals)
    if len(lost2_vals) > 1 or len(gained2_vals) > 1:
        print(f"    a={key[0]}, b={key[1]}, c={key[2]}, d={key[3]}:")
        print(f"      lost2={sorted(lost2_vals)}, gained2={sorted(gained2_vals)}")
        print(f"      lost3={sorted(lost3_vals)}, gained3={sorted(gained3_vals)}")

# Try: is lost2 determined by just (|A|, |D|) or (|A|, |C|)?
print(f"\n  Is lost2 determined by (|A|, |D|) only?")
by_ad = defaultdict(set)
for d in formula_check:
    by_ad[(d['la'], d['ld'])].add(d['lost2'])
determined = all(len(v) == 1 for v in by_ad.values())
print(f"    Determined by (|A|,|D|): {determined}")
if not determined:
    for key, vals in sorted(by_ad.items()):
        if len(vals) > 1:
            print(f"      a={key[0]}, d={key[1]}: lost2 in {sorted(vals)}")

# Check if delta3 = 2*delta2 exactly
exact_2to1 = sum(1 for d in formula_check if d['delta3'] == 2 * d['delta2'])
print(f"\n  delta3 = 2*delta2 exactly: {exact_2to1}/{len(formula_check)} ({100*exact_2to1/len(formula_check):.1f}%)")

# If not exact, what's the distribution of delta3 - 2*delta2?
residual = Counter(d['delta3'] - 2*d['delta2'] for d in formula_check)
print(f"  delta3 - 2*delta2 distribution: {dict(sorted(residual.items()))}")

# Check: is the residual explained by some function of |A|,|B|,|C|,|D|?
print(f"\n  Residual delta3 - 2*delta2 by (a,b,c,d):")
by_abcd_resid = defaultdict(set)
for d in formula_check:
    key = (d['la'], d['lb'], d['lc'], d['ld'])
    by_abcd_resid[key].add(d['delta3'] - 2*d['delta2'])
for key in sorted(by_abcd_resid.keys()):
    vals = by_abcd_resid[key]
    if len(vals) > 1:
        print(f"    a={key[0]}, b={key[1]}, c={key[2]}, d={key[3]}: residual in {sorted(vals)}")

print("\nDone.")
