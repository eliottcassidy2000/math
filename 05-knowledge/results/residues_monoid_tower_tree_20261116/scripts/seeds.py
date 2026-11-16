#!/usr/bin/env python3
"""seeds.py: wave 45's one-port caterpillar automaton on (h, g, r, s) at the end x: L = hang a leaf on x (new vertex z
beats only x; root stays x), S = extend the spine (z beats only x; z becomes the root).  Build the actual tournaments
for all words up to length W, verify the maps against the exact rooted profile (E = r, O = s), and list the ready
seeds (b = s - r - g = 0) with counts 13, 23, 41, 59, 75, 95 as mask lists for the port-move search."""
import sys, subprocess
W = int(sys.argv[1])
def profile(masks, root):
    n = len(masks); N = 1 << n
    dp = [[0]*n for _ in range(N)]; dpr = [[0]*n for _ in range(N)]
    for v in range(n): dp[1<<v][v] = 1
    dpr[1<<root][root] = 1
    for m in range(1, N):
        for l in range(n):
            a, b = dp[m][l], dpr[m][l]
            if not a and not b: continue
            nx = masks[l] & ~m
            while nx:
                v = (nx & -nx).bit_length() - 1; nx &= nx - 1
                dp[m|(1<<v)][v] += a; dpr[m|(1<<v)][v] += b
    Hm = [sum(dp[m]) for m in range(N)]; Hm[0] = 1; rs = [sum(dpr[m]) for m in range(N)]
    H = Hm[N-1]; g = Hm[(N-1) & ~(1<<root)]
    E = sum(rs[m]*Hm[(N-1) & ~m] for m in range(1, N) if (m >> root) & 1 and bin(m).count('1') % 2 == 0)
    O = sum(rs[m]*Hm[(N-1) & ~m] for m in range(1, N) if (m >> root) & 1 and bin(m).count('1') % 2 == 1)
    return H, E, O, g
def append(masks, root, move):
    n = len(masks); z = n; masks = masks[:] + [0]
    for v in range(n):
        if v == root: masks[z] |= 1 << v       # z beats only the root
        else: masks[v] |= 1 << z               # every other old vertex beats z
    return masks, (z if move == 'S' else root)
targets = {13, 23, 41, 59, 75, 95}
found = {}; checked = 0; mism = 0
def rec(masks, root, word, state):
    global checked, mism
    H, E, O, g = profile(masks, root); h, gg, r, s = state
    checked += 1
    if (H, g, E, O) != (h, gg, r, s): mism += 1; print("MISMATCH", word, (H, g, E, O), state)
    b = s - r - g
    if b == 0 and h in targets and h not in found: found[h] = (word, masks, root, (h, gg, r, s))
    if len(word) >= W: return
    for mv in 'LS':
        nm, nr = append(masks, root, mv)
        h2 = h + 2*r
        st = (h2, gg, r + s - gg, r + s) if mv == 'L' else (h2, h, s, h + r)
        rec(nm, nr, word + mv, st)
rec([0], 0, '', (1, 1, 0, 1))
print(f"words up to length {W}: {checked} checked, automaton mismatches {mism}")
for h in sorted(found): print(f"ready seed h={h}: word {found[h][0]} order {len(found[h][1])} root {found[h][2]} masks {found[h][1]} (h,g,r,s)={found[h][3]}")
print("ready counts found among targets:", sorted(found), "; missing:", sorted(targets - set(found)))
