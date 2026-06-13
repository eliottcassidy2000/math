#!/usr/bin/env python3
"""
gn_black_structure_s211.py — What determines black degree at every n?
opus-2026-03-22-S211

KEY OBSERVATION: at n=6, most NS classes have black_deg = 2.
The SC-free classes (black=0) are special.
What property of a class determines its black degree?

Also: the SC↔SC subgraph has beautiful structure. Let's study it.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  BLACK DEGREE STRUCTURE IN G_n")
print("  opus-2026-03-22-S211")
print("=" * 80)

# ============================================================================
# HELPERS
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
# COMPUTE
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

    id_to_data = {}
    for canon in sorted_canons:
        d = class_data[canon]
        id_to_data[d['id']] = d

    neighbor_set = defaultdict(set)
    edge_set = set()
    self_loop = Counter()

    # Also track DIRECTED edge weights: out_weight[i][j] = # (T in class i, arc) -> class j
    out_weight = defaultdict(lambda: defaultdict(int))

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
                out_weight[id1][id2] += 1

    blue_edges = set()
    black_edges = set()
    for (i, j) in edge_set:
        if id_to_data[i]['sc'] == id_to_data[j]['sc']:
            blue_edges.add((i, j))
        else:
            black_edges.add((i, j))

    blue_deg = Counter()
    black_deg = Counter()
    for (i, j) in blue_edges:
        blue_deg[i] += 1
        blue_deg[j] += 1
    for (i, j) in black_edges:
        black_deg[i] += 1
        black_deg[j] += 1

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
        'out_weight': dict(out_weight),
    }

    print(f"  n={n}: done ({time.time()-t0:.1f}s)")

# ============================================================================
# ANALYSIS: What determines black_deg for NS classes?
# ============================================================================

print("\n" + "=" * 80)
print("  WHAT DETERMINES BLACK DEGREE?")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]
    m = comb(n, 2)
    nfact = factorial(n)

    sc_ids = sorted([i for i, d in r['id_to_data'].items() if d['sc']])
    ns_ids = sorted([i for i, d in r['id_to_data'].items() if not d['sc']])

    if not ns_ids:
        continue

    print(f"\n  n={n}: {len(sc_ids)} SC, {len(ns_ids)} NS classes")

    # For each NS class, look at: score_seq, H, c3, |Aut|, black_deg
    # And which specific SC classes it connects to
    print(f"\n    NS classes with their black (SC) neighbors:")
    print(f"    {'ID':>3} {'H':>4} {'scores':>20} {'|Aut|':>5} {'c3':>3} "
          f"{'black':>5} {'SC neighbors':>30}")

    for i in ns_ids:
        d = r['id_to_data'][i]
        sc_nbrs = sorted([j for j in r['neighbor_set'].get(i, set())
                         if r['id_to_data'][j]['sc']])
        sc_str = ', '.join(f"{j}(H={r['id_to_data'][j]['h']})" for j in sc_nbrs)
        print(f"    {i:3d} {d['h']:4d} {str(d['scores']):>20} {d['aut']:5d} "
              f"{d['c3']:3d} {r['black_deg'].get(i, 0):5d} {sc_str}")

    # Is there a formula? Check if black_deg correlates with score variance
    print(f"\n    Checking score-based predictors of black_deg:")
    for i in ns_ids:
        d = r['id_to_data'][i]
        ss = d['scores']
        mean_s = sum(ss) / n
        var_s = sum((s - mean_s)**2 for s in ss) / n  # score variance

        # Score "imbalance" = max(s) - min(s)
        imbalance = ss[-1] - ss[0]

        # Is score palindromic? s_i + s_{n-1-i} = n-1
        is_palin = all(ss[k] + ss[n-1-k] == n-1 for k in range(n))

        bd = r['black_deg'].get(i, 0)
        print(f"    class {i:2d}: var_s={var_s:.2f}  imbal={imbalance}  "
              f"palin={'Y' if is_palin else 'N'}  black={bd}")

# ============================================================================
# ANALYSIS: Score palindromicity and SC/black connectivity
# ============================================================================

print("\n" + "=" * 80)
print("  SCORE PALINDROMICITY AND BLACK DEGREE")
print("=" * 80)

for n in sorted(all_data.keys()):
    r = all_data[n]

    print(f"\n  n={n}:")
    for i in sorted(r['id_to_data'].keys()):
        d = r['id_to_data'][i]
        ss = d['scores']
        is_palin = all(ss[k] + ss[n-1-k] == n-1 for k in range(n))
        bd = r['black_deg'].get(i, 0)
        bu = r['blue_deg'].get(i, 0)

        print(f"    {i:2d} H={d['h']:3d} {'SC' if d['sc'] else 'NS':>3} "
              f"scores={str(ss):>20} {'PALIN' if is_palin else '     '} "
              f"b={bu:2d} k={bd:2d}")

# ============================================================================
# ANALYSIS: Directed weight from NS→SC
# ============================================================================

print("\n" + "=" * 80)
print("  DIRECTED WEIGHTS: NS→SC TRANSITIONS")
print("=" * 80)

for n in [5, 6]:
    r = all_data[n]
    m = comb(n, 2)

    print(f"\n  n={n}:")

    ns_ids = sorted([i for i, d in r['id_to_data'].items() if not d['sc']])
    sc_ids = sorted([i for i, d in r['id_to_data'].items() if d['sc']])

    for i in ns_ids:
        d = r['id_to_data'][i]
        total_out = d['size'] * m  # total transitions from this class
        sc_out = sum(r['out_weight'].get(i, {}).get(j, 0) for j in sc_ids)
        ns_out = sum(r['out_weight'].get(i, {}).get(j, 0) for j in ns_ids if j != i)
        sl = r['self_loop'].get(i, 0)

        print(f"    NS {i:2d} (H={d['h']:3d}, |C|={d['size']:4d}): "
              f"total={total_out:6d}  →SC={sc_out:5d}({sc_out/total_out:.3f})  "
              f"→NS={ns_out:5d}({ns_out/total_out:.3f})  "
              f"self={sl:4d}({sl/total_out:.3f})")

    print(f"\n  SC classes:")
    for i in sc_ids:
        d = r['id_to_data'][i]
        total_out = d['size'] * m
        sc_out = sum(r['out_weight'].get(i, {}).get(j, 0) for j in sc_ids if j != i)
        ns_out = sum(r['out_weight'].get(i, {}).get(j, 0) for j in ns_ids)
        sl = r['self_loop'].get(i, 0)

        print(f"    SC {i:2d} (H={d['h']:3d}, |C|={d['size']:4d}): "
              f"total={total_out:6d}  →SC={sc_out:5d}({sc_out/total_out:.3f})  "
              f"→NS={ns_out:5d}({ns_out/total_out:.3f})  "
              f"self={sl:4d}({sl/total_out:.3f})")

# ============================================================================
# ANALYSIS: Which pairs of SC classes does each NS class "bridge"?
# ============================================================================

print("\n" + "=" * 80)
print("  NS CLASSES AS BRIDGES BETWEEN SC CLASSES")
print("=" * 80)

for n in [5, 6]:
    r = all_data[n]
    ns_ids = sorted([i for i, d in r['id_to_data'].items() if not d['sc']])
    sc_ids = sorted([i for i, d in r['id_to_data'].items() if d['sc']])

    if n == 6 and len(ns_ids) > 30:
        print(f"\n  n={n}: (showing only NS classes with black_deg > 0)")
        ns_ids = [i for i in ns_ids if r['black_deg'].get(i, 0) > 0]

    print(f"\n  n={n}:")
    for i in ns_ids:
        d = r['id_to_data'][i]
        sc_nbrs = sorted([j for j in r['neighbor_set'].get(i, set())
                         if r['id_to_data'][j]['sc']])
        if len(sc_nbrs) >= 2:
            # This NS class bridges sc_nbrs[0] and sc_nbrs[1] (and more)
            bridge_str = f"bridges {sc_nbrs}"
        elif len(sc_nbrs) == 1:
            bridge_str = f"touches SC {sc_nbrs[0]} only"
        else:
            bridge_str = f"SC-free"

        print(f"    NS {i:2d} (H={d['h']:3d}): black={r['black_deg'].get(i, 0)}, {bridge_str}")

# ============================================================================
# GLOBAL SEQUENCES
# ============================================================================

print("\n" + "=" * 80)
print("  SEQUENCES FOR OEIS/PATTERN SEARCH")
print("=" * 80)

print("\n  Blue edges:   n=3,4,5,6 → ", [len(all_data[n]['blue_edges']) for n in [3,4,5,6]])
print("  Black edges:  n=3,4,5,6 → ", [len(all_data[n]['black_edges']) for n in [3,4,5,6]])
print("  Total edges:  n=3,4,5,6 → ", [len(all_data[n]['edge_set']) for n in [3,4,5,6]])
print("  SC-SC edges:  n=3,4,5,6 → ", end="")

sc_sc_edges = []
for n in [3,4,5,6]:
    r = all_data[n]
    count = len([e for e in r['blue_edges']
                if r['id_to_data'][e[0]]['sc'] and r['id_to_data'][e[1]]['sc']])
    sc_sc_edges.append(count)
print(sc_sc_edges)

ns_ns_edges = []
for n in [3,4,5,6]:
    r = all_data[n]
    count = len([e for e in r['blue_edges']
                if not r['id_to_data'][e[0]]['sc'] and not r['id_to_data'][e[1]]['sc']])
    ns_ns_edges.append(count)
print(f"  NS-NS edges:  n=3,4,5,6 → {ns_ns_edges}")

# Max blue degree, min blue degree
for label, deg_dict in [("blue", "blue_deg"), ("black", "black_deg")]:
    for typ in ["SC", "NS"]:
        degs = []
        for n in [3,4,5,6]:
            r = all_data[n]
            ids = [i for i, d in r['id_to_data'].items()
                   if (d['sc'] if typ == "SC" else not d['sc'])]
            if ids:
                dd = [r[deg_dict].get(i, 0) for i in ids]
                degs.append((min(dd), max(dd), sum(dd)/len(dd)))
            else:
                degs.append(None)
        print(f"  {typ} {label}_deg (min,max,avg): {degs}")

# SC backbone genus
print("\n  SC backbone: (vertices, edges, genus = E-V+1)")
for n in [3,4,5,6]:
    r = all_data[n]
    sc_ids = [i for i, d in r['id_to_data'].items() if d['sc']]
    sc_edges = [e for e in r['blue_edges']
               if r['id_to_data'][e[0]]['sc'] and r['id_to_data'][e[1]]['sc']]
    v_count = len(sc_ids)
    e_count = len(sc_edges)
    genus = e_count - v_count + 1
    print(f"    n={n}: V={v_count}, E={e_count}, genus={genus}")

print("\n" + "=" * 80)
print("  DONE")
print("=" * 80)
