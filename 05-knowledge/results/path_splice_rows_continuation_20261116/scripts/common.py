#!/usr/bin/env python3
"""common.py <n> <root> <gentourng-string>: convert to masks, run ./profile search2, and list the continuations (X', g')
and (E', O', g') shared by the three digits d = 0, 2, 4, with the realising extension masks."""
import sys, subprocess, re
n = int(sys.argv[1]); root = int(sys.argv[2]); s = sys.argv[3]
masks = [0]*n; k = 0
for i in range(n):
    for j in range(i+1, n):
        if s[k] == '1': masks[i] |= 1 << j
        else: masks[j] |= 1 << i
        k += 1
print("masks:", masks)
out = subprocess.run(['./profile', 'search2', str(n), str(root)] + [str(m) for m in masks], capture_output=True, text=True).stdout
print(out.splitlines()[0])
hits = {0: [], 2: [], 4: []}
lines = out.splitlines()
for i, line in enumerate(lines):
    m = re.match(r'd=(\d) strong=(\d) masks=(\d+),? \(H,E,O,g\)=\((\d+),(\d+),(\d+),(\d+)\)', line)
    if m:
        d, st, cx, H, E, O, g = map(int, m.groups())
        cy = lines[i+1].strip() if i+1 < len(lines) and lines[i+1].startswith('   cy=') else ''
        hits[d].append((E+O, g, E, O, st, cx, cy))
for key, name in [((0, 1), '(X,g)'), ((2, 3, 1), '(E,O,g)')]:
    sets = [set(tuple(h[i] for i in key) for h in hits[d]) for d in (0, 2, 4)]
    common = sets[0] & sets[1] & sets[2]
    print(f"continuations {name} shared by all three digits:", sorted(common))
    for c in sorted(common):
        for d in (0, 2, 4):
            ex = [h for h in hits[d] if tuple(h[i] for i in key) == c]
            print(f"   d={d}: {len(ex)} extensions, e.g. masks cx={ex[0][5]} {ex[0][6]} (E,O,g)=({ex[0][2]},{ex[0][3]},{ex[0][1]}) strong={ex[0][4]}")
