#!/usr/bin/env python3
"""
gn_blue_black_patterns_s211.py — Deep pattern analysis of blue/black degrees
opus-2026-03-22-S211

Focused on finding FORMULAS and LAWS governing the (blue, black) degree profile.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  DEEP PATTERN ANALYSIS: Blue/Black Degrees in G_n")
print("  opus-2026-03-22-S211")
print("=" * 80)

# ============================================================================
# HELPERS (same as before)
# ============================================================================

def tournament_adj(n, bits):
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
            idx += 1
    return adj

def H_dp(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            val = dp[S * n + v]
            if val == 0:
                continue
            for u in range(n):
                if S & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(S | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

def canonical_form(adj, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_sc(adj, n):
    comp = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    return canonical_form(adj, n) == canonical_form(comp, n)

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def count_3cycles(adj, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def aut_size(adj, n):
    count = 0
    for perm in permutations(range(n)):
        is_aut = True
        for i in range(n):
            for j in range(n):
                if adj[i][j] != adj[perm[i]][perm[j]]:
                    is_aut = False
                    break
            if not is_aut:
                break
        if is_aut:
            count += 1
    return count

# ============================================================================
# COMPUTE FOR n=3,4,5,6
# ============================================================================

all_data = {}

for n in range(3, 7):
    m = comb(n, 2)
    total = 1 << m
    t0 = time.time()

    iso_classes = defaultdict(list)
    adj_cache = {}

    for bits in range(total):
        adj = tournament_adj(n, bits)
        canon = canonical_form(adj, n)
        iso_classes[canon].append(bits)
        adj_cache[bits] = adj

    class_data = {}
    for canon, members in iso_classes.items():
        rep_adj = adj_cache[members[0]]
        h = H_dp(rep_adj, n)
        sc = is_sc(rep_adj, n)
        ss = score_seq(rep_adj, n)
        c3 = count_3cycles(rep_adj, n)
        aut = aut_size(rep_adj, n)
        class_data[canon] = {
            'h': h, 'sc': sc, 'scores': ss, 'c3': c3,
            'size': len(members), 'members': set(members),
            'aut': aut, 'id': None,
        }

    sorted_canons = sorted(class_data.keys(),
                           key=lambda c: (class_data[c]['h'], class_data[c]['scores'], class_data[c]['size']))
    for i, canon in enumerate(sorted_canons):
        class_data[canon]['id'] = i

    bits_to_id = {}
    for canon, data in class_data.items():
        for b in data['members']:
            bits_to_id[b] = data['id']

    # Build id_to_data
    id_to_data = {}
    for canon in sorted_canons:
        d = class_data[canon]
        id_to_data[d['id']] = d

    # Compute edges
    neighbor_set = defaultdict(set)
    edge_set = set()
    self_loop = Counter()

    for bits in range(total):
        id1 = bits_to_id[bits]
        for bit_pos in range(m):
            new_bits = bits ^ (1 << bit_pos)
            id2 = bits_to_id[new_bits]
            if id1 == id2:
                self_loop[id1] += 1
            else:
                neighbor_set[id1].add(id2)
                edge_set.add((min(id1, id2), max(id1, id2)))

    # Classify edges
    blue_edges = set()
    black_edges = set()
    for (i, j) in edge_set:
        if id_to_data[i]['sc'] == id_to_data[j]['sc']:
            blue_edges.add((i, j))
        else:
            black_edges.add((i, j))

    # Per-class degrees
    blue_deg = Counter()
    black_deg = Counter()
    for (i, j) in blue_edges:
        blue_deg[i] += 1
        blue_deg[j] += 1
    for (i, j) in black_edges:
        black_deg[i] += 1
        black_deg[j] += 1

    # Complement map
    full_mask = (1 << m) - 1
    complement_of = {}
    for i, d in id_to_data.items():
        rep = next(iter(d['members']))
        comp_rep = full_mask ^ rep
        complement_of[i] = bits_to_id[comp_rep]

    all_data[n] = {
        'id_to_data': id_to_data,
        'blue_deg': dict(blue_deg),
        'black_deg': dict(black_deg),
        'self_loop': dict(self_loop),
        'blue_edges': blue_edges,
        'black_edges': black_edges,
        'edge_set': edge_set,
        'complement_of': complement_of,
        'neighbor_set': dict(neighbor_set),
        'num_classes': len(sorted_canons),
    }

    print(f"\n  n={n}: {len(sorted_canons)} classes, "
          f"{len(blue_edges)} blue + {len(black_edges)} black = {len(edge_set)} edges "
          f"({time.time()-t0:.1f}s)")

# ============================================================================
# ANALYSIS 1: How SC status inverts blue/black
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 1: SC STATUS INVERTS BLUE/BLACK ROLES")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    sc_classes = [i for i, d in r['id_to_data'].items() if d['sc']]
    ns_classes = [i for i, d in r['id_to_data'].items() if not d['sc']]
    num_sc = len(sc_classes)
    num_ns = len(ns_classes)

    print(f"\n  n={n}: {num_sc} SC + {num_ns} NS = {num_sc + num_ns} classes")

    if sc_classes:
        for i in sorted(sc_classes):
            d = r['id_to_data'][i]
            b = r['blue_deg'].get(i, 0)
            k = r['black_deg'].get(i, 0)
            # SC neighbors are blue, NS neighbors are black
            print(f"    SC class {i:2d} (H={d['h']:3d}): "
                  f"connects to {b} SC class(es), {k} NS class(es)")

    if ns_classes and len(ns_classes) <= 20:
        for i in sorted(ns_classes):
            d = r['id_to_data'][i]
            b = r['blue_deg'].get(i, 0)
            k = r['black_deg'].get(i, 0)
            # NS neighbors are blue, SC neighbors are black
            print(f"    NS class {i:2d} (H={d['h']:3d}): "
                  f"connects to {b} NS class(es), {k} SC class(es)")

# ============================================================================
# ANALYSIS 2: Which SC classes does each SC class see?
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 2: THE SC↔SC SUBGRAPH (blue edges between SC classes)")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    sc_classes = sorted([i for i, d in r['id_to_data'].items() if d['sc']])

    print(f"\n  n={n}: SC classes = {sc_classes}")

    sc_sc_edges = [(i, j) for (i, j) in r['blue_edges']
                   if r['id_to_data'][i]['sc'] and r['id_to_data'][j]['sc']]

    for i in sc_classes:
        sc_neighbors = [j for (a, b) in sc_sc_edges for j in ([b] if a == i else ([a] if b == i else []))]
        d = r['id_to_data'][i]
        sc_nbr_str = ', '.join(f"{j}(H={r['id_to_data'][j]['h']})" for j in sorted(sc_neighbors))
        print(f"    SC {i:2d} (H={d['h']:3d}): SC-neighbors = [{sc_nbr_str}]")

    print(f"    Total SC↔SC edges: {len(sc_sc_edges)}")
    max_sc_edges = len(sc_classes) * (len(sc_classes) - 1) // 2
    print(f"    SC↔SC density: {len(sc_sc_edges)}/{max_sc_edges} = "
          f"{len(sc_sc_edges)/max_sc_edges:.3f}" if max_sc_edges > 0 else "")

# ============================================================================
# ANALYSIS 3: Black degree distribution — how many SC classes does each NS class touch?
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 3: NS→SC CONNECTIVITY (black degree of NS classes)")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    ns_classes = sorted([i for i, d in r['id_to_data'].items() if not d['sc']])

    if not ns_classes:
        print(f"\n  n={n}: no NS classes")
        continue

    black_degs = [r['black_deg'].get(i, 0) for i in ns_classes]
    print(f"\n  n={n}: NS classes' black_deg (= # SC neighbors):")
    print(f"    values: {sorted(black_degs)}")
    print(f"    distribution: {dict(Counter(black_degs))}")

    # Which NS classes are "SC-free"?
    sc_free = [i for i in ns_classes if r['black_deg'].get(i, 0) == 0]
    if sc_free:
        print(f"    SC-free NS classes (black=0): {sc_free}")
        for i in sc_free:
            d = r['id_to_data'][i]
            print(f"      class {i}: H={d['h']}, scores={d['scores']}, |Aut|={d['aut']}, c3={d['c3']}")

# ============================================================================
# ANALYSIS 4: blue+black = total_deg. What determines total_deg?
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 4: TOTAL DEGREE FORMULA SEARCH")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    m = comb(n, 2)
    nfact = factorial(n)

    print(f"\n  n={n}, m={m}, n!={nfact}:")
    print(f"    {'ID':>3} {'H':>4} {'SC':>3} {'|Aut|':>5} {'|C|':>5} "
          f"{'deg':>4} {'C(n,2)-1':>8} {'n-1':>4} "
          f"{'n!/(|Aut|*deg)':>14} {'|C|/deg':>8}")

    for i in sorted(r['id_to_data'].keys()):
        d = r['id_to_data'][i]
        deg = r['blue_deg'].get(i, 0) + r['black_deg'].get(i, 0)
        aut = d['aut']
        size = d['size']
        ratio1 = nfact / (aut * deg) if deg > 0 else float('inf')
        ratio2 = size / deg if deg > 0 else float('inf')
        print(f"    {i:3d} {d['h']:4d} {'SC' if d['sc'] else 'NS':>3} "
              f"{aut:5d} {size:5d} {deg:4d} {m-1:8d} {n-1:4d} "
              f"{ratio1:14.4f} {ratio2:8.4f}")

# ============================================================================
# ANALYSIS 5: blue / (blue+black) as a function of... what?
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 5: BLUE FRACTION b/(b+k)")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    sc_classes = sorted([i for i, d in r['id_to_data'].items() if d['sc']])
    ns_classes = sorted([i for i, d in r['id_to_data'].items() if not d['sc']])
    num_sc = len(sc_classes)
    num_ns = len(ns_classes)
    total_classes = num_sc + num_ns

    print(f"\n  n={n}: {num_sc} SC + {num_ns} NS = {total_classes}")

    # If connections were random, SC class would connect to:
    # (num_sc-1)/(total_classes-1) fraction SC neighbors (blue)
    # num_ns/(total_classes-1) fraction NS neighbors (black)
    expected_sc_blue_frac = (num_sc - 1) / (total_classes - 1) if total_classes > 1 else 1
    expected_ns_blue_frac = (num_ns - 1) / (total_classes - 1) if total_classes > 1 else 1

    print(f"    Random model prediction:")
    print(f"      SC blue fraction = (|SC|-1)/(|V|-1) = {num_sc-1}/{total_classes-1} = {expected_sc_blue_frac:.4f}")
    print(f"      NS blue fraction = (|NS|-1)/(|V|-1) = {num_ns-1}/{total_classes-1} = {expected_ns_blue_frac:.4f}")

    if sc_classes:
        actual_sc_blue_fracs = []
        for i in sc_classes:
            b = r['blue_deg'].get(i, 0)
            k = r['black_deg'].get(i, 0)
            if b + k > 0:
                actual_sc_blue_fracs.append(b / (b + k))
        if actual_sc_blue_fracs:
            print(f"    Actual SC blue fractions: {[f'{x:.3f}' for x in sorted(actual_sc_blue_fracs)]}")
            print(f"    Actual SC mean blue fraction: {sum(actual_sc_blue_fracs)/len(actual_sc_blue_fracs):.4f}")

    if ns_classes:
        actual_ns_blue_fracs = []
        for i in ns_classes:
            b = r['blue_deg'].get(i, 0)
            k = r['black_deg'].get(i, 0)
            if b + k > 0:
                actual_ns_blue_fracs.append(b / (b + k))
        if actual_ns_blue_fracs:
            print(f"    Actual NS blue fractions: min={min(actual_ns_blue_fracs):.3f}, "
                  f"max={max(actual_ns_blue_fracs):.3f}, "
                  f"mean={sum(actual_ns_blue_fracs)/len(actual_ns_blue_fracs):.4f}")

# ============================================================================
# ANALYSIS 6: Self-loop and blue/black degree interaction
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 6: SELF-LOOPS vs BLUE/BLACK")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    m = comb(n, 2)
    nfact = factorial(n)

    print(f"\n  n={n}:")
    has_sl = sorted([i for i in r['id_to_data'] if r['self_loop'].get(i, 0) > 0])
    no_sl = sorted([i for i in r['id_to_data'] if r['self_loop'].get(i, 0) == 0])

    print(f"    Classes WITH self-loop ({len(has_sl)}): {has_sl}")
    print(f"    Classes WITHOUT self-loop ({len(no_sl)}): {no_sl}")

    if has_sl:
        sl_blue = [r['blue_deg'].get(i, 0) for i in has_sl]
        sl_black = [r['black_deg'].get(i, 0) for i in has_sl]
        print(f"    Self-loop classes: avg_blue={sum(sl_blue)/len(sl_blue):.2f}, "
              f"avg_black={sum(sl_black)/len(sl_black):.2f}")

    if no_sl:
        nsl_blue = [r['blue_deg'].get(i, 0) for i in no_sl]
        nsl_black = [r['black_deg'].get(i, 0) for i in no_sl]
        print(f"    No-self-loop classes: avg_blue={sum(nsl_blue)/len(nsl_blue):.2f}, "
              f"avg_black={sum(nsl_black)/len(nsl_black):.2f}")

# ============================================================================
# ANALYSIS 7: H-level grouping — blue/black at each H level
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 7: BLUE/BLACK BY H-LEVEL")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]

    # Group by H
    h_groups = defaultdict(list)
    for i, d in r['id_to_data'].items():
        h_groups[d['h']].append(i)

    print(f"\n  n={n}:")
    for h_val in sorted(h_groups.keys()):
        classes = h_groups[h_val]
        entries = []
        for i in sorted(classes):
            b = r['blue_deg'].get(i, 0)
            k = r['black_deg'].get(i, 0)
            sc = 'SC' if r['id_to_data'][i]['sc'] else 'NS'
            entries.append(f"{sc}({b},{k})")
        print(f"    H={h_val:3d}: {', '.join(entries)}")

# ============================================================================
# ANALYSIS 8: Complement symmetry — show H↔H pairing
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 8: COMPLEMENT PAIRING WITH H SYMMETRY")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    comp = r['complement_of']

    print(f"\n  n={n}:")

    # Find the H axis of symmetry
    h_values = [r['id_to_data'][i]['h'] for i in r['id_to_data']]
    h_center = (min(h_values) + max(h_values)) / 2

    print(f"    H range: [{min(h_values)}, {max(h_values)}]")
    print(f"    H center: {h_center}")

    # Show complement pairs
    seen = set()
    for i in sorted(r['id_to_data'].keys()):
        if i in seen:
            continue
        ci = comp[i]
        seen.add(i)
        seen.add(ci)

        d_i = r['id_to_data'][i]
        d_ci = r['id_to_data'][ci]

        b_i = r['blue_deg'].get(i, 0)
        k_i = r['black_deg'].get(i, 0)

        if i == ci:
            marker = "SELF-COMP"
        else:
            marker = f"paired"

        print(f"    {i:2d}↔{ci:2d}: H=({d_i['h']:3d},{d_ci['h']:3d})  "
              f"blue=({b_i},{r['blue_deg'].get(ci, 0)})  "
              f"black=({k_i},{r['black_deg'].get(ci, 0)})  "
              f"|Aut|=({d_i['aut']},{d_ci['aut']})  "
              f"{'SC' if d_i['sc'] else 'NS'} {marker}")

# ============================================================================
# ANALYSIS 9: Look at degree sequences as a whole
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 9: DEGREE SEQUENCES")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    num = r['num_classes']

    all_degs = sorted([r['blue_deg'].get(i, 0) + r['black_deg'].get(i, 0)
                       for i in range(num)])
    all_blue = sorted([r['blue_deg'].get(i, 0) for i in range(num)])
    all_black = sorted([r['black_deg'].get(i, 0) for i in range(num)])

    print(f"\n  n={n} ({num} classes):")
    print(f"    Total deg: {all_degs}")
    print(f"    Blue deg:  {all_blue}")
    print(f"    Black deg: {all_black}")
    print(f"    Sum(total_deg) = 2×|E| = {sum(all_degs)} = 2×{sum(all_degs)//2}")
    print(f"    Sum(blue_deg)  = 2×|blue| = {sum(all_blue)} = 2×{sum(all_blue)//2}")
    print(f"    Sum(black_deg) = 2×|black| = {sum(all_black)} = 2×{sum(all_black)//2}")

# ============================================================================
# ANALYSIS 10: Which NS classes at n=6 connect to which SC classes?
# ============================================================================

print("\n" + "=" * 80)
print("  ANALYSIS 10: NS→SC NEIGHBORS AT n=6")
print("=" * 80)

if 6 in all_data:
    r = all_data[6]
    sc_classes = sorted([i for i, d in r['id_to_data'].items() if d['sc']])
    ns_classes = sorted([i for i, d in r['id_to_data'].items() if not d['sc']])

    # For each SC class, which NS classes connect to it?
    for sc_id in sc_classes:
        d = r['id_to_data'][sc_id]
        ns_nbrs = [j for j in r['neighbor_set'].get(sc_id, set())
                    if not r['id_to_data'][j]['sc']]
        print(f"  SC {sc_id:2d} (H={d['h']:3d}): {len(ns_nbrs)} NS neighbors = {sorted(ns_nbrs)}")

    # How many SC classes does each SC class see?
    print(f"\n  SC↔SC adjacency:")
    for sc_id in sc_classes:
        sc_nbrs = [j for j in r['neighbor_set'].get(sc_id, set())
                   if r['id_to_data'][j]['sc']]
        print(f"    SC {sc_id:2d} (H={d['h']:3d}): "
              f"SC neighbors = {sorted(sc_nbrs)}")

print("\n" + "=" * 80)
print("  DONE")
print("=" * 80)
