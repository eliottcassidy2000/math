#!/usr/bin/env python3
"""seeds2.py W: enumerate caterpillar words over {L, S} up to length W by the one-port automaton alone (fast), keep the
shortest word per ready count (b = s - r - g = 0) among the targets 13, 23, 41, 59, 75, 95, then build the tournaments
and verify their exact rooted profiles."""
import sys
sys.setrecursionlimit(10000)
W = int(sys.argv[1]); targets = {13, 23, 41, 59, 75, 95}; found = {}
def rec(word, state):
    h, g, r, s = state
    if s - r - g == 0 and h in targets and (h not in found or len(word) < len(found[h])): found[h] = word
    if len(word) >= W: return
    rec(word + 'L', (h + 2*r, g, r + s - g, r + s)); rec(word + 'S', (h + 2*r, h, s, h + r))
rec('', (1, 1, 0, 1))
print("shortest ready words by count:", {h: found[h] for h in sorted(found)}, "; missing:", sorted(targets - set(found)))
exec(open('seeds.py').read().split("targets = {")[0].split("def append")[0])  # reuse profile()
def append(masks, root, move):
    n = len(masks); z = n; masks = masks[:] + [0]
    for v in range(n):
        if v == root: masks[z] |= 1 << v
        else: masks[v] |= 1 << z
    return masks, (z if move == 'S' else root)
for h in sorted(found):
    masks, root = [0], 0
    for mv in found[h]: masks, root = append(masks, root, mv)
    H, E, O, g = profile(masks, root)
    print(f"seed h={h}: word {found[h]} order {len(masks)} root {root} profile (H,E,O,g)=({H},{E},{O},{g}) charge a={2*E+1-H} defect b={O-E-g} masks {masks}")
