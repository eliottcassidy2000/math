#!/usr/bin/env python3
"""
n7_path_analysis_s20cr.py — Principal path analysis at n=7 using cached BFS tree
kind-pasteur-2026-03-23-S20cr

Uses the BFS tree from principal_line_n7_s20cp.out to analyze:
1. Shortest path from transitive to Paley (H=189)
2. H-monotonicity along the path
3. All paths to H=189 through the SC blue spine
4. Branch structure and local optima
"""

import sys
sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  N=7 PRINCIPAL PATH ANALYSIS (FROM CACHED BFS TREE)")
print("  kind-pasteur-2026-03-23-S20cr")
print("=" * 80)

# Parse the BFS tree from previous computation output
# Format from principal_line_n7_s20cp.out:
# [mid] H=... score=(...) c3=... |Aut|=... ch=[children]

# Reconstructed parent-child relationships from the output:
tree = {
    0: {'H': 1, 'score': (0,1,2,3,4,5,6), 'c3': 0, 'aut': 1, 'children': [10, 19, 7]},
    10: {'H': 9, 'score': (0,2,2,3,4,4,6), 'c3': 3, 'aut': 1, 'children': [101, 14, 17, 88]},
    19: {'H': 33, 'score': (1,1,2,3,4,5,5), 'c3': 5, 'aut': 1, 'children': [90]},
    7: {'H': 3, 'score': (0,1,3,3,3,5,6), 'c3': 1, 'aut': 3, 'children': []},
    101: {'H': 93, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': [132, 102, 152, 190]},
    14: {'H': 9, 'score': (0,2,2,3,4,4,6), 'c3': 3, 'aut': 1, 'children': [89]},
    17: {'H': 11, 'score': (0,2,3,3,3,4,6), 'c3': 4, 'aut': 1, 'children': [120]},
    88: {'H': 15, 'score': (0,2,3,3,3,4,6), 'c3': 4, 'aut': 3, 'children': []},
    90: {'H': 51, 'score': (1,1,3,3,3,5,5), 'c3': 6, 'aut': 3, 'children': []},
    132: {'H': 93, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': [201]},
    102: {'H': 85, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': [169, 191, 111]},
    152: {'H': 111, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': [238]},
    190: {'H': 123, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 3, 'children': []},
    89: {'H': 13, 'score': (0,2,3,3,3,4,6), 'c3': 4, 'aut': 1, 'children': []},
    120: {'H': 15, 'score': (0,3,3,3,3,3,6), 'c3': 5, 'aut': 5, 'children': []},
    201: {'H': 117, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': []},
    169: {'H': 105, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 1, 'children': [170, 171, 187, 219]},
    191: {'H': 111, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 3, 'children': [197]},
    111: {'H': 81, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': []},
    238: {'H': 135, 'score': (1,3,3,3,3,3,5), 'c3': 10, 'aut': 5, 'children': []},
    170: {'H': 105, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 1, 'children': [242, 183, 27]},
    171: {'H': 123, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [28, 243, 188]},
    187: {'H': 137, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [233, 189]},
    219: {'H': 135, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 3, 'children': []},
    197: {'H': 99, 'score': (1,3,3,3,3,3,5), 'c3': 10, 'aut': 3, 'children': []},
    242: {'H': 141, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [128, 259, 262, 241, 246]},
    183: {'H': 125, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [260, 268, 184]},
    27: {'H': 45, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': []},
    28: {'H': 51, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': [129]},
    243: {'H': 155, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [264, 247]},
    188: {'H': 151, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [265, 269]},
    233: {'H': 159, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 3, 'children': []},
    189: {'H': 145, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [263]},
    128: {'H': 77, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': [125]},
    259: {'H': 159, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [256, 230, 271]},
    262: {'H': 145, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [163, 231]},
    241: {'H': 143, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [226, 150]},
    246: {'H': 137, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': []},
    260: {'H': 153, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [261, 270, 253]},
    268: {'H': 147, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [258, 266, 245]},
    184: {'H': 117, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 1, 'children': []},
    129: {'H': 83, 'score': (1,3,3,3,3,3,5), 'c3': 10, 'aut': 1, 'children': []},
    264: {'H': 157, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [257]},
    247: {'H': 151, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [244]},
    265: {'H': 135, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': []},
    269: {'H': 175, 'score': (3,3,3,3,3,3,3), 'c3': 14, 'aut': 7, 'children': []},
    263: {'H': 129, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 1, 'children': [254]},
    125: {'H': 79, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': [24, 123]},
    256: {'H': 159, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [162, 248]},
    230: {'H': 129, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [232]},
    271: {'H': 189, 'score': (3,3,3,3,3,3,3), 'c3': 14, 'aut': 21, 'children': []},
    163: {'H': 75, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': [211]},
    231: {'H': 115, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 1, 'children': []},
    226: {'H': 117, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 1, 'children': [80]},
    150: {'H': 67, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': []},
    261: {'H': 139, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [267, 174]},
    270: {'H': 171, 'score': (3,3,3,3,3,3,3), 'c3': 14, 'aut': 3, 'children': []},
    253: {'H': 133, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': []},
    258: {'H': 151, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': []},
    266: {'H': 133, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [138]},
    245: {'H': 139, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': []},
    257: {'H': 159, 'score': (2,3,3,3,3,3,4), 'c3': 13, 'aut': 1, 'children': [214]},
    244: {'H': 141, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': [203, 117]},
    254: {'H': 127, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 1, 'children': [209, 147]},
    24: {'H': 31, 'score': (1,1,3,3,3,5,5), 'c3': 6, 'aut': 1, 'children': [22]},
    123: {'H': 69, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': []},
    162: {'H': 93, 'score': (1,3,3,3,3,3,5), 'c3': 10, 'aut': 1, 'children': [206]},
    248: {'H': 141, 'score': (2,2,3,3,3,4,4), 'c3': 12, 'aut': 1, 'children': []},
    232: {'H': 99, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 9, 'children': []},
    211: {'H': 61, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': []},
    80: {'H': 49, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': []},
    267: {'H': 129, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 3, 'children': [139]},
    174: {'H': 71, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': []},
    138: {'H': 77, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': []},
    214: {'H': 95, 'score': (1,3,3,3,3,3,5), 'c3': 10, 'aut': 1, 'children': []},
    203: {'H': 91, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': [65]},
    117: {'H': 67, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': []},
    209: {'H': 69, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': [72]},
    147: {'H': 65, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': []},
    22: {'H': 25, 'score': (1,1,2,3,4,5,5), 'c3': 5, 'aut': 1, 'children': []},
    206: {'H': 87, 'score': (1,2,3,3,3,4,5), 'c3': 9, 'aut': 1, 'children': []},
    139: {'H': 65, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': [44]},
    65: {'H': 41, 'score': (1,1,3,3,3,5,5), 'c3': 6, 'aut': 1, 'children': []},
    72: {'H': 31, 'score': (1,1,2,3,4,5,5), 'c3': 5, 'aut': 1, 'children': []},
    44: {'H': 27, 'score': (1,1,2,3,4,5,5), 'c3': 5, 'aut': 1, 'children': [45, 141]},
    45: {'H': 9, 'score': (1,1,1,3,5,5,5), 'c3': 2, 'aut': 9, 'children': []},
    141: {'H': 65, 'score': (1,2,2,3,4,4,5), 'c3': 8, 'aut': 1, 'children': [250]},
    250: {'H': 135, 'score': (2,2,2,3,4,4,4), 'c3': 11, 'aut': 3, 'children': []},
}

# Build parent map
parent = {}
for mid, data in tree.items():
    for ch in data['children']:
        parent[ch] = mid

# ============================================================================
# 1. PATH FROM TRANSITIVE TO PALEY (H=189)
# ============================================================================
print(f"\n{'='*70}")
print("  PATH FROM TRANSITIVE (H=1) TO PALEY (H=189)")
print(f"{'='*70}")

path = []
v = 271  # Paley
while v in parent:
    path.append(v)
    v = parent[v]
path.append(0)  # transitive (root)
path.reverse()

print(f"\n  Path length: {len(path)-1} blue edges")
print(f"  Path (mid): {path}")
print(f"\n  Detailed path:")
for i, mid in enumerate(path):
    d = tree[mid]
    delta = ""
    if i > 0:
        dh = d['H'] - tree[path[i-1]]['H']
        delta = f"  Delta_H={dh:+d}"
        if dh < 0:
            delta += " <-- DOWNHILL!"
    print(f"    [{i}] mid={mid:3d}  H={d['H']:3d}  c3={d['c3']:2d}  score={d['score']}  |Aut|={d['aut']}{delta}")

h_path = [tree[m]['H'] for m in path]
c3_path = [tree[m]['c3'] for m in path]

print(f"\n  H sequence: {h_path}")
print(f"  c3 sequence: {c3_path}")

drops = [(i, h_path[i], h_path[i+1]) for i in range(len(h_path)-1) if h_path[i+1] < h_path[i]]
print(f"\n  H-drops (downhill steps): {len(drops)}")
for i, h_before, h_after in drops:
    print(f"    Step {i}->{i+1}: H {h_before} -> {h_after} (Delta = {h_after-h_before})")

# ============================================================================
# 2. ALL PATHS TO REGULAR TOURNAMENTS
# ============================================================================
print(f"\n{'='*70}")
print("  ALL PATHS TO REGULAR (c3=14) TOURNAMENTS")
print(f"{'='*70}")

# There are 3 regular classes at n=7: 269 (H=175), 270 (H=171), 271 (H=189=Paley)
regular = [mid for mid, d in tree.items() if d['score'] == (3,3,3,3,3,3,3)]
print(f"\n  Regular SC classes: {regular}")
for mid in regular:
    d = tree[mid]
    print(f"    mid={mid}: H={d['H']}, |Aut|={d['aut']}")

for target in regular:
    p = []
    v = target
    while v in parent:
        p.append(v); v = parent[v]
    p.append(0); p.reverse()
    hp = [tree[m]['H'] for m in p]
    drops_t = sum(1 for i in range(len(hp)-1) if hp[i+1]<hp[i])
    print(f"\n  Path to mid={target} (H={tree[target]['H']}):")
    print(f"    Length: {len(p)-1}, drops: {drops_t}")
    print(f"    H: {hp}")

# ============================================================================
# 3. BRANCH ANALYSIS AT CRITICAL JUNCTION
# ============================================================================
print(f"\n{'='*70}")
print("  CRITICAL JUNCTION: NODE 101 (H=93)")
print(f"{'='*70}")

# Node 101 is the gateway — all paths to high-H go through it
node = 101
d = tree[node]
print(f"\n  Node {node}: H={d['H']}, score={d['score']}, c3={d['c3']}")
print(f"  Children: {d['children']}")
for ch in d['children']:
    cd = tree[ch]
    # Find max H reachable from this child
    def max_h_subtree(mid):
        best = tree[mid]['H']
        for c in tree[mid].get('children', []):
            best = max(best, max_h_subtree(c))
        return best
    def subtree_size(mid):
        return 1 + sum(subtree_size(c) for c in tree[mid].get('children', []))
    mh = max_h_subtree(ch)
    sz = subtree_size(ch)
    print(f"    [{ch}] H={cd['H']:3d} -> max_H_reachable={mh}, subtree_size={sz}")
    print(f"         score={cd['score']}, c3={cd['c3']}, |Aut|={cd['aut']}")

print("""
  KEY INSIGHT: From node 101 (H=93), the branch structure is:
    - [132] H=93 -> leads to max H=117 (dead end)
    - [102] H=85 -> leads to max H=189 (PALEY!) but requires going DOWN
    - [152] H=111 -> leads to max H=135 (dead end)
    - [190] H=123 -> leaf (dead end)

  The ONLY path to Paley goes through node 102 (H=85), which is a DOWNHILL
  step from 101 (H=93). A greedy H-maximizing walk would NEVER take this path.

  This is the CRITICAL BOTTLENECK: the path to the global maximum (H=189)
  passes through a local minimum relative to the junction node.
""")

# ============================================================================
# 4. DEPTH PROFILE
# ============================================================================
print(f"{'='*70}")
print("  DEPTH PROFILE OF SC BLUE SPINE")
print(f"{'='*70}")

# Compute depth for each node (iterative BFS)
depth = {0: 0}
queue = [0]
while queue:
    mid = queue.pop(0)
    for ch in tree[mid].get('children', []):
        if ch not in depth:
            depth[ch] = depth[mid] + 1
            queue.append(ch)

# Max depth
max_d = max(depth.values())
print(f"\n  Max depth: {max_d}")

# Depth of key nodes
for mid in [271, 269, 270]:
    d = tree[mid]
    dep = depth[mid]
    print(f"  Regular node {mid} (H={d['H']}): depth={dep}")

# Nodes per depth level
from collections import Counter
depth_count = Counter(depth.values())
print(f"\n  Nodes per depth level:")
for d in range(max_d+1):
    nodes_at_d = [mid for mid, dp in depth.items() if dp == d]
    h_vals = sorted([tree[mid]['H'] for mid in nodes_at_d])
    avg_h = sum(h_vals)/len(h_vals)
    max_h = max(h_vals)
    print(f"    d={d:2d}: {depth_count[d]:3d} nodes, H range [{min(h_vals)},{max(h_vals)}], avg_H={avg_h:.1f}")

# Leaves (nodes with no children)
leaves = [mid for mid, d in tree.items() if not d['children']]
print(f"\n  Leaves: {len(leaves)} of {len(tree)} nodes ({100*len(leaves)/len(tree):.1f}%)")
leaf_depths = [depth[mid] for mid in leaves]
print(f"  Leaf depths: min={min(leaf_depths)}, max={max(leaf_depths)}, avg={sum(leaf_depths)/len(leaf_depths):.1f}")
leaf_H = [tree[mid]['H'] for mid in leaves]
print(f"  Leaf H: min={min(leaf_H)}, max={max(leaf_H)}, avg={sum(leaf_H)/len(leaf_H):.1f}")

# ============================================================================
# 5. H-GRADIENT STATISTICS
# ============================================================================
print(f"\n{'='*70}")
print("  H-GRADIENT ON SC BLUE TREE")
print(f"{'='*70}")

uphill = 0; downhill = 0; level = 0
for mid, d in tree.items():
    for ch in d['children']:
        dh = tree[ch]['H'] - d['H']
        if dh > 0: uphill += 1
        elif dh < 0: downhill += 1
        else: level += 1

total_edges = uphill + downhill + level
print(f"  Total tree edges: {total_edges}")
print(f"  Uphill (child H > parent H): {uphill} ({100*uphill/total_edges:.1f}%)")
print(f"  Downhill (child H < parent H): {downhill} ({100*downhill/total_edges:.1f}%)")
print(f"  Level (same H): {level} ({100*level/total_edges:.1f}%)")

# List all downhill edges
print(f"\n  Downhill edges (child H < parent H):")
for mid, d in tree.items():
    for ch in d['children']:
        dh = tree[ch]['H'] - d['H']
        if dh < 0:
            print(f"    [{mid}] H={d['H']} -> [{ch}] H={tree[ch]['H']} (Delta={dh})")

# ============================================================================
# 6. OPUS S223 COMPARISON: PRINCIPAL PATH SEQUENCES
# ============================================================================
print(f"\n{'='*70}")
print("  PRINCIPAL PATH SEQUENCES (cross-n)")
print(f"{'='*70}")

print("""
  Opus S223 found the principal path (greedy H-increasing along SC backbone)
  splits by PARITY:
    Odd n: 1 -> 3 -> 15 -> 123 -> ...
    Even n: 1 -> 5 -> 37 -> 389 -> ...

  At n=7, the greedy path reaches H=123 (node 190), NOT H=189 (Paley).
  Node 190 is a LEAF with |Aut|=3 and score (1,2,3,3,3,4,5).

  The path from transitive to Paley requires going THROUGH the bottleneck:
    0(H=1) -> 10(H=9) -> 101(H=93) -> 102(H=85) -> ... -> 271(H=189)
                                         ^^^^^^^^
                              DOWNHILL from 93 to 85

  This means the greedy/principal path and the path to global max DIVERGE
  at node 101. The principal path takes [190] (H=123, leaf, uphill) while
  the Paley path takes [102] (H=85, deep subtree, initially downhill).

  STRUCTURAL INSIGHT:
  The SC blue spine at n=7 has a RUGGED FITNESS LANDSCAPE.
  The global H-maximum (Paley, H=189) is hidden behind a local barrier.
  A greedy climber gets trapped at H=123 (a leaf dead-end).

  The "valley" from 93 to 85 has depth Delta_H = -8.
  After crossing this valley, H climbs monotonically: 85 -> 105 -> 105 -> 141 -> 159 -> 189.
""")

# Delta H = 2^(n-2) from transitive to first blue neighbor
print(f"  FIRST BLUE NEIGHBOR Delta_H:")
first_neighbors = {
    3: 3, 4: 5, 5: 9, 6: 5, 7: 33
}
# Actually from the data:
# n=3: H=1->3, delta=2
# n=4: H=1->5, delta=4
# n=5: H=1->9, delta=8 (but also H=1->3)
# n=6: H=1->5, delta=4 (but also H=1->17)
# n=7: H=1->33, delta=32 (but also H=1->9 and H=1->3)

# The Delta_H = 2^(n-2) pattern:
print(f"  n=3: Delta = 2 = 2^1 (highest first neighbor)")
print(f"  n=4: Delta = 4 = 2^2 (only blue neighbor)")
print(f"  n=5: Delta = 8 = 2^3 (highest: mid=4, H=9)")
print(f"  n=6: Delta = 16 = 2^4 (highest: mid=11, H=17)")
print(f"  n=7: Delta = 32 = 2^5 (highest: mid=19, H=33)")
print(f"  Pattern: 2^(n-2) is the H-value of the highest first blue neighbor")
print(f"  Equivalently: max Delta_H from transitive = 2^(n-2)")

print(f"\n  DONE.")
print("=" * 80)
