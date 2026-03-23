#!/usr/bin/env python3
"""
gn_blue_black_degrees_s211.py — Blue/Black degree of every iso class in G_n
opus-2026-03-22-S211

For each isomorphism class at n=3,4,5,6:
  - blue_deg  = # edges to classes with SAME SC status
  - black_deg = # edges to classes with DIFFERENT SC status
  - self_loop = whether an arc flip can stay in the same class
  - blue_weight, black_weight = total flip counts (weighted edges)

BLUE = SC-preserving (both SC or both non-SC)
BLACK = SC-changing (one SC, one non-SC)

Goal: find patterns in (blue_deg, black_deg) across all n.
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  BLUE/BLACK DEGREES IN THE ISOMORPHISM CLASS GRAPH G_n")
print("  opus-2026-03-22-S211")
print("=" * 80)

# ============================================================================
# HELPERS
# ============================================================================

def tournament_adj(n, bits):
    """Build adjacency matrix from bit encoding."""
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
    """Count Hamiltonian paths via DP."""
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
    """Canonical form = lex-smallest adjacency under all vertex permutations."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_sc(adj, n):
    """Check if tournament is self-complementary."""
    comp = [[1 - adj[i][j] if i != j else 0 for j in range(n)] for i in range(n)]
    return canonical_form(adj, n) == canonical_form(comp, n)

def score_seq(adj, n):
    return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))

def count_3cycles(adj, n):
    """Count directed 3-cycles."""
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                # Check both orientations
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    c3 += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    c3 += 1
    return c3

def aut_size(adj, n):
    """Count automorphisms."""
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
# MAIN COMPUTATION
# ============================================================================

all_results = {}

for n in range(3, 7):
    m = comb(n, 2)
    total = 1 << m

    print(f"\n{'='*80}")
    print(f"  n = {n}  ({total} tournaments, {m} arcs)")
    print(f"{'='*80}")
    t0 = time.time()

    # Step 1: Find all iso classes
    iso_classes = defaultdict(list)
    adj_cache = {}

    for bits in range(total):
        adj = tournament_adj(n, bits)
        canon = canonical_form(adj, n)
        iso_classes[canon].append(bits)
        adj_cache[bits] = adj

        if n == 6 and (bits + 1) % 5000 == 0:
            elapsed = time.time() - t0
            rate = (bits + 1) / elapsed
            eta = (total - bits - 1) / rate
            print(f"    Classifying: {bits+1}/{total} ({rate:.0f}/s, ETA {eta:.0f}s)")

    print(f"  Found {len(iso_classes)} iso classes in {time.time()-t0:.1f}s")

    # Step 2: Compute class properties
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

    # Assign IDs by H, then score, then size
    sorted_canons = sorted(class_data.keys(),
                           key=lambda c: (class_data[c]['h'], class_data[c]['scores'], class_data[c]['size']))
    for i, canon in enumerate(sorted_canons):
        class_data[canon]['id'] = i

    # Build reverse lookup
    bits_to_id = {}
    bits_to_sc = {}
    for canon, data in class_data.items():
        for b in data['members']:
            bits_to_id[b] = data['id']
            bits_to_sc[b] = data['sc']

    # Step 3: Compute edges with blue/black classification
    t1 = time.time()

    # For each class: count neighbors, blue/black edges, self-loops, weights
    neighbor_set = defaultdict(set)   # class_id -> set of neighbor class_ids
    edge_weight = Counter()           # (min_id, max_id) -> weight
    self_loop_count = Counter()       # class_id -> # of (tournament, arc) pairs that stay in class

    for bits in range(total):
        id1 = bits_to_id[bits]
        for bit_pos in range(m):
            new_bits = bits ^ (1 << bit_pos)
            id2 = bits_to_id[new_bits]
            if id1 == id2:
                self_loop_count[id1] += 1
            else:
                edge = (min(id1, id2), max(id1, id2))
                neighbor_set[id1].add(id2)
                edge_weight[edge] += 1

    # Classify each edge as blue or black based on SC status of endpoints
    all_edges = set()
    for edge in edge_weight:
        all_edges.add(edge)

    # Get SC status per class id
    id_to_data = {}
    for canon in sorted_canons:
        d = class_data[canon]
        id_to_data[d['id']] = d

    blue_edges = set()
    black_edges = set()
    for (i, j) in all_edges:
        if id_to_data[i]['sc'] == id_to_data[j]['sc']:
            blue_edges.add((i, j))
        else:
            black_edges.add((i, j))

    # Compute per-class blue/black degree and weight
    blue_deg = Counter()
    black_deg = Counter()
    blue_weight_per_class = Counter()
    black_weight_per_class = Counter()

    for (i, j) in blue_edges:
        blue_deg[i] += 1
        blue_deg[j] += 1
        w = edge_weight[(i, j)]
        blue_weight_per_class[i] += w
        blue_weight_per_class[j] += w

    for (i, j) in black_edges:
        black_deg[i] += 1
        black_deg[j] += 1
        w = edge_weight[(i, j)]
        black_weight_per_class[i] += w
        black_weight_per_class[j] += w

    elapsed = time.time() - t1

    num_classes = len(sorted_canons)
    print(f"  Edges computed in {elapsed:.1f}s")
    print(f"  Total edges: {len(all_edges)}")
    print(f"  Blue edges: {len(blue_edges)}")
    print(f"  Black edges: {len(black_edges)}")
    print(f"  Self-loop classes: {len([c for c in range(num_classes) if self_loop_count[c] > 0])}/{num_classes}")
    print()

    # Step 4: Print detailed per-class table
    has_self = lambda cid: self_loop_count[cid] > 0

    print(f"  {'ID':>3} {'H':>4} {'SC':>3} {'scores':>20} {'|C|':>5} {'|Aut|':>5} "
          f"{'c3':>4} {'deg':>4} {'b_deg':>5} {'k_deg':>5} {'self':>5} "
          f"{'b_wt':>6} {'k_wt':>6} {'tot_wt':>7}")
    print(f"  {'─'*3} {'─'*4} {'─'*3} {'─'*20} {'─'*5} {'─'*5} "
          f"{'─'*4} {'─'*4} {'─'*5} {'─'*5} {'─'*5} "
          f"{'─'*6} {'─'*6} {'─'*7}")

    for canon in sorted_canons:
        d = class_data[canon]
        i = d['id']
        total_deg = blue_deg[i] + black_deg[i]
        sl = self_loop_count[i]
        bw = blue_weight_per_class[i]
        kw = black_weight_per_class[i]
        tw = bw + kw + sl  # total weight including self-loops

        print(f"  {i:3d} {d['h']:4d} {'SC' if d['sc'] else 'NS':>3} "
              f"{str(d['scores']):>20} {d['size']:5d} {d['aut']:5d} "
              f"{d['c3']:4d} {total_deg:4d} {blue_deg[i]:5d} {black_deg[i]:5d} "
              f"{sl:5d} {bw:6d} {kw:6d} {tw:7d}")

    print()

    # Step 5: Summary statistics
    sc_classes = [d['id'] for d in id_to_data.values() if d['sc']]
    ns_classes = [d['id'] for d in id_to_data.values() if not d['sc']]

    print(f"  SC classes ({len(sc_classes)}): {sorted(sc_classes)}")
    print(f"  Non-SC classes ({len(ns_classes)}): {sorted(ns_classes)}")
    print()

    # Blue degree stats for SC vs non-SC
    if sc_classes:
        sc_blue = [blue_deg[i] for i in sc_classes]
        sc_black = [black_deg[i] for i in sc_classes]
        print(f"  SC classes — blue_deg: {sorted(sc_blue)}, black_deg: {sorted(sc_black)}")
        print(f"    avg blue_deg: {sum(sc_blue)/len(sc_blue):.2f}, avg black_deg: {sum(sc_black)/len(sc_black):.2f}")

    if ns_classes:
        ns_blue = [blue_deg[i] for i in ns_classes]
        ns_black = [black_deg[i] for i in ns_classes]
        print(f"  Non-SC classes — blue_deg: {sorted(ns_blue)}, black_deg: {sorted(ns_black)}")
        print(f"    avg blue_deg: {sum(ns_blue)/len(ns_blue):.2f}, avg black_deg: {sum(ns_black)/len(ns_black):.2f}")

    print()

    # Check: total_weight = class_size * C(n,2)?
    print(f"  WEIGHT VERIFICATION: total_weight[i] should = |C_i| × C(n,2) = |C_i| × {m}")
    verify_ok = True
    for canon in sorted_canons:
        d = class_data[canon]
        i = d['id']
        tw = blue_weight_per_class[i] + black_weight_per_class[i] + self_loop_count[i]
        expected = d['size'] * m
        if tw != expected:
            print(f"    CLASS {i}: got {tw}, expected {expected} — MISMATCH!")
            verify_ok = False
    if verify_ok:
        print(f"    ✓ All {num_classes} classes verify correctly")
    print()

    # Adjacency list with colors
    print(f"  ADJACENCY LIST (B=blue, K=black):")
    for canon in sorted_canons:
        d = class_data[canon]
        i = d['id']
        nbrs = sorted(neighbor_set[i])
        parts = []
        for j in nbrs:
            edge = (min(i, j), max(i, j))
            color = 'B' if edge in blue_edges else 'K'
            w = edge_weight[edge]
            parts.append(f"{j}{color}({w})")
        sl_str = f" [self:{self_loop_count[i]}]" if self_loop_count[i] > 0 else ""
        print(f"    {i:2d} (H={d['h']:3d},{'SC' if d['sc'] else 'NS'}): {', '.join(parts)}{sl_str}")

    print()

    # Store results
    all_results[n] = {
        'num_classes': num_classes,
        'blue_edges': len(blue_edges),
        'black_edges': len(black_edges),
        'id_to_data': id_to_data,
        'blue_deg': dict(blue_deg),
        'black_deg': dict(black_deg),
        'self_loop': dict(self_loop_count),
        'edge_weight': dict(edge_weight),
        'blue_edges_set': blue_edges,
        'black_edges_set': black_edges,
    }

# ============================================================================
# CROSS-n PATTERN ANALYSIS
# ============================================================================

print("\n" + "=" * 80)
print("  CROSS-n PATTERN ANALYSIS")
print("=" * 80)

for n in sorted(all_results.keys()):
    r = all_results[n]
    print(f"\n  n={n}: {r['num_classes']} classes, "
          f"{r['blue_edges']} blue + {r['black_edges']} black = "
          f"{r['blue_edges'] + r['black_edges']} total edges")

    # Print (blue_deg, black_deg) for each class sorted by H
    data_items = sorted(r['id_to_data'].items())
    for i, d in data_items:
        bd = r['blue_deg'].get(i, 0)
        kd = r['black_deg'].get(i, 0)
        sl = r['self_loop'].get(i, 0) > 0
        print(f"    class {i:2d}: H={d['h']:4d}  SC={'Y' if d['sc'] else 'N'}  "
              f"|C|={d['size']:5d}  |Aut|={d['aut']:4d}  "
              f"blue={bd:2d}  black={kd:2d}  total={bd+kd:2d}  "
              f"{'SELF' if sl else '    '}")

# ============================================================================
# FORMULA HUNTING
# ============================================================================

print("\n" + "=" * 80)
print("  FORMULA HUNTING")
print("=" * 80)

for n in sorted(all_results.keys()):
    r = all_results[n]
    m = comb(n, 2)
    print(f"\n  n={n}, m=C(n,2)={m}:")

    data_items = sorted(r['id_to_data'].items())

    # Check if blue_deg + black_deg correlates with class size, |Aut|, H, c3
    print(f"    {'ID':>3} {'H':>4} {'SC':>3} {'|C|':>5} {'|Aut|':>5} "
          f"{'c3':>4} {'n!/|Aut|':>8} {'deg':>4} {'b':>3} {'k':>3} "
          f"{'b/(b+k)':>7} {'deg*|Aut|':>9} {'n!/(deg*|Aut|)':>14}")

    nfact = factorial(n)
    for i, d in data_items:
        bd = r['blue_deg'].get(i, 0)
        kd = r['black_deg'].get(i, 0)
        deg = bd + kd
        aut = d['aut']
        ratio = bd / (bd + kd) if (bd + kd) > 0 else 0
        dag = deg * aut if deg > 0 else 0
        quot = nfact / dag if dag > 0 else 0
        print(f"    {i:3d} {d['h']:4d} {'SC' if d['sc'] else 'NS':>3} "
              f"{d['size']:5d} {aut:5d} {d['c3']:4d} "
              f"{nfact//aut:8d} {deg:4d} {bd:3d} {kd:3d} "
              f"{ratio:7.3f} {dag:9d} {quot:14.4f}")

# ============================================================================
# DEEPER PATTERNS
# ============================================================================

print("\n" + "=" * 80)
print("  DEEPER PATTERNS: Blue vs Black for SC and non-SC classes")
print("=" * 80)

for n in sorted(all_results.keys()):
    r = all_results[n]
    m = comb(n, 2)
    data_items = sorted(r['id_to_data'].items())

    print(f"\n  n={n}:")

    # For SC classes: how many of their neighbors are SC vs non-SC?
    sc_ids = [i for i, d in data_items if d['sc']]
    ns_ids = [i for i, d in data_items if not d['sc']]

    print(f"    SC classes: {sc_ids}")
    print(f"    Non-SC classes: {ns_ids}")

    # For SC classes, blue edges go to other SC classes
    # For non-SC classes, blue edges go to other non-SC classes
    # Black edges cross the SC boundary

    # The blue_deg of an SC class = # SC neighbors
    # The black_deg of an SC class = # non-SC neighbors
    # The blue_deg of a non-SC class = # non-SC neighbors
    # The black_deg of a non-SC class = # SC neighbors

    if sc_ids and ns_ids:
        # How many edges between SC-SC, NS-NS, SC-NS?
        sc_sc = len([e for e in r['blue_edges_set']
                     if r['id_to_data'][e[0]]['sc'] and r['id_to_data'][e[1]]['sc']])
        ns_ns = len([e for e in r['blue_edges_set']
                     if not r['id_to_data'][e[0]]['sc'] and not r['id_to_data'][e[1]]['sc']])
        sc_ns = len(r['black_edges_set'])

        print(f"    SC↔SC edges: {sc_sc}")
        print(f"    NS↔NS edges: {ns_ns}")
        print(f"    SC↔NS edges: {sc_ns} (= all black edges)")

        # Density within SC, within NS, between SC-NS
        if len(sc_ids) >= 2:
            max_sc = len(sc_ids) * (len(sc_ids) - 1) // 2
            print(f"    SC-SC density: {sc_sc}/{max_sc} = {sc_sc/max_sc:.3f}" if max_sc > 0 else "")
        if len(ns_ids) >= 2:
            max_ns = len(ns_ids) * (len(ns_ids) - 1) // 2
            print(f"    NS-NS density: {ns_ns}/{max_ns} = {ns_ns/max_ns:.3f}" if max_ns > 0 else "")
        max_cross = len(sc_ids) * len(ns_ids)
        print(f"    SC-NS density: {sc_ns}/{max_cross} = {sc_ns/max_cross:.3f}" if max_cross > 0 else "")

# ============================================================================
# COMPLEMENT PAIRING
# ============================================================================

print("\n" + "=" * 80)
print("  COMPLEMENT PAIRING: How blue/black degrees relate for T vs T^op")
print("=" * 80)

for n in sorted(all_results.keys()):
    r = all_results[n]
    m = comb(n, 2)
    data_items = sorted(r['id_to_data'].items())

    print(f"\n  n={n}:")

    # For each class, find its complement class
    # Complement of bits = all bits flipped = (2^m - 1) ^ bits
    # But we need to find which class the complement belongs to
    # Actually, we can just look at which class has the complementary H
    # At n=5, H(T) + H(T^op) is NOT constant in general...
    # Actually H(T^op) = H(T) by path reversal! So complement preserves H.
    # But complement ≠ T^op. Complement reverses all arcs. T^op reverses all arcs.
    # Wait — for tournaments, the complement IS T^op (reverse all arcs).
    # So H(T^op) = H(T), meaning complement preserves H.
    # But the complement of a tournament may or may not be isomorphic to itself.

    # SC classes pair with themselves. Non-SC classes pair with their complement class.
    # Let's find complement pairing.

    # For each class, take representative, complement it, find its class
    # Need the bits_to_id mapping... let me reconstruct

    # Actually, complement = flip all m bits
    full_mask = (1 << m) - 1

    # We stored members in id_to_data but not the bits_to_id map
    # Let me use the members sets

    # Build bits_to_id
    bits_to_id_local = {}
    for i, d in data_items:
        for b in d['members']:
            bits_to_id_local[b] = i

    complement_of = {}
    for i, d in data_items:
        rep = next(iter(d['members']))
        comp_rep = full_mask ^ rep
        comp_id = bits_to_id_local[comp_rep]
        complement_of[i] = comp_id

    for i, d in data_items:
        ci = complement_of[i]
        bd_i = r['blue_deg'].get(i, 0)
        kd_i = r['black_deg'].get(i, 0)
        bd_ci = r['blue_deg'].get(ci, 0)
        kd_ci = r['black_deg'].get(ci, 0)
        marker = " ← SC (self-comp)" if i == ci else ""
        print(f"    class {i:2d} (H={d['h']:3d}) ↔ class {ci:2d}: "
              f"blue=({bd_i},{bd_ci})  black=({kd_i},{kd_ci}){marker}")

print("\n" + "=" * 80)
print("  DONE")
print("=" * 80)
