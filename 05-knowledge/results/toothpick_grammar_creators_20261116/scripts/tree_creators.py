#!/usr/bin/env python3
"""tree_creators.py (R-7): first creators in the canonical-parent tree.  counts_z.txt (order h per class, digraph6 stream
order) aligns with cparent_n.txt.  For each order: first creators, how many have a first creator as tree parent, how many
are pure chains (every tree ancestor back to the 3-cycle or the single vertex a first creator), and the tree-children
statistics of first creators versus other classes (unequal growth)."""
from collections import defaultdict, Counter
H = defaultdict(list)
for line in open('counts_z.txt'):
    n, h = map(int, line.split()); H[n].append(h)
par = {2: [0]}
for n in range(3, 10): par[n] = [int(x) for x in open(f'cparent_{n}.txt')]
seen = set(); first = {}
for n in range(1, 10):
    new = set(H[n]) - seen; first[n] = [h in new for h in H[n]]; seen |= set(H[n])
nchild = {n: Counter(par[n + 1]) for n in range(1, 9)}
pure = {1: [True], 2: [False]}
print("order: first creators | tree parent a first creator | pure chains | mean tree-children of first creators vs others | share of next-order first creators whose tree parent is a first creator")
for n in range(3, 10):
    fc = [i for i in range(len(H[n])) if first[n][i]]
    pf = sum(1 for i in fc if first[n - 1][par[n][i]])
    pure[n] = [first[n][i] and (pure[n - 1][par[n][i]] or H[n - 1][par[n][i]] == 1 and n - 1 <= 2) for i in range(len(H[n]))]
    npure = sum(pure[n])
    if n < 9:
        mc_f = sum(nchild[n][i] for i in fc) / max(1, len(fc)); mc_o = sum(nchild[n][i] for i in range(len(H[n])) if not first[n][i]) / max(1, len(H[n]) - len(fc))
        nxt = [i for i in range(len(H[n + 1])) if first[n + 1][i]]; share = sum(1 for i in nxt if first[n][par[n + 1][i]]) / max(1, len(nxt))
        print(f"  {n}: {len(fc)} | {pf} | {npure} | {mc_f:.2f} vs {mc_o:.2f} | {share:.3f}")
    else:
        print(f"  {n}: {len(fc)} | {pf} | {npure} | - | -")
