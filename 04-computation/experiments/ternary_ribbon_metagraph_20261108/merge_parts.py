#!/usr/bin/env python3
"""merge_parts.py n -- merge the part files of lines_census.c for order n, compare the k- and cut-set tables with
lines_formula.py's output, and compute the value-colour and pair statistics (THM-4930's eternal colours, exact
below 1131; values >= 1131 are 'unknown')."""
import sys, glob, math
from collections import Counter
n = int(sys.argv[1])
BLACK = set([7, 21, 19, 23, 35, 39, 43, 47, 53, 55, 59, 63, 101, 103, 107, 109, 113, 119, 131, 163, 167, 181, 183, 185, 203, 205, 207, 209, 215,
             231, 235, 239, 251, 267, 271, 291, 299, 301, 303, 311, 317, 321, 327, 331, 335, 347, 359, 377, 415, 427, 431, 439, 451, 463, 483, 487, 503, 515, 519, 523, 527, 563,
             671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949, 983, 1003, 1009, 1067, 1087, 1115])
def colour(h): return 'unknown' if h >= 1131 else ('black' if h in BLACK else 'blue')
k = Counter(); kb = Counter(); cuts = Counter(); cutsb = Counter(); pairs = Counter(); meta = Counter()
for path in sorted(glob.glob(f"lines_n{n}_p*.txt")):
    sec = None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('# n='):
            toks = line.replace(':', '').split()
            for key in ('lines', 'blue', 'both-ends-strong'):
                i = toks.index(key); meta[key] += int(toks[i + 1])
            continue
        if line.startswith('#k'): sec = 'k'; continue
        if line.startswith('#cuts'): sec = 'cuts'; continue
        if line.startswith('#pairs'): sec = 'pairs'; continue
        if line.startswith('#'): sec = None; continue
        if not line: continue
        p = line.split()
        if sec == 'k': k[int(p[0])] += int(p[1]); kb[int(p[0])] += int(p[2])
        elif sec == 'cuts': cuts[int(p[0])] += int(p[1]); cutsb[int(p[0])] += int(p[2])
        elif sec == 'pairs': pairs[(int(p[0]), int(p[1]), int(p[2]), int(p[3]))] += int(p[4])
print(f"n={n}: parts {len(glob.glob(f'lines_n{n}_p*.txt'))}, lines {meta['lines']}, blue {meta['blue']}, both ends strong {meta['both-ends-strong']} ({meta['both-ends-strong']/meta['lines']:.4f})")
# formula comparison
fk = {}; fc = {}; sec = None
for line in open(f"formula_n{n}.txt"):
    line = line.strip()
    if line.startswith('#k'): sec = 'k'; continue
    if line.startswith('#cuts'): sec = 'cuts'; continue
    if line.startswith('#'): sec = None; continue
    if sec == 'k': a, b, c = line.split(); fk[int(a)] = (int(b), int(c))
    elif sec == 'cuts': a, b, c = line.split(); fc[int(a)] = (int(b), int(c))
ck = {kk: (k[kk], kb[kk]) for kk in range(1, n + 1)}
cc = {c: (cuts[c], cutsb[c]) for c in cuts}
print("  k-table census:", ck)
print("  k-table == formula:", ck == {kk: fk.get(kk, (0, 0)) for kk in range(1, n + 1)}, "; cut-set table == formula:", cc == fc)
if pairs:
    total = sum(pairs.values())
    V = Counter()
    for (a, b, bl, st), c in pairs.items():
        V[(tuple(sorted((colour(a), colour(b)))), 'blue' if bl else 'black')] += c
    print("  value colours (lines):", dict(V))
    bb = sum(c for (a, b, bl, st), c in pairs.items() if colour(a) == colour(b) == 'blue'); kk2 = sum(c for (a, b, bl, st), c in pairs.items() if colour(a) == colour(b) == 'black')
    unk = sum(c for (a, b, bl, st), c in pairs.items() if 'unknown' in (colour(a), colour(b)))
    print(f"  both eternally blue {bb} ({bb/total:.4f}), both eternally black {kk2} ({kk2/total:.4f}), with a value >= 1131 (undecided) {unk} ({unk/total:.4f})")
    same = sum(c for (a, b, bl, st), c in pairs.items() if a == b); lt = sum(c for (a, b, bl, st), c in pairs.items() if a < b); gt = sum(c for (a, b, bl, st), c in pairs.items() if a > b)
    EHc = sum(a * c for (a, b, bl, st), c in pairs.items()) / total; EHo = sum(b * c for (a, b, bl, st), c in pairs.items()) / total
    cov = sum((a - EHc) * (b - EHo) * c for (a, b, bl, st), c in pairs.items()) / total
    va = sum((a - EHc) ** 2 * c for (a, b, bl, st), c in pairs.items()) / total; vb = sum((b - EHo) ** 2 * c for (a, b, bl, st), c in pairs.items()) / total
    print(f"  same-value {same} ({same/total:.4f}); closed larger {gt} ({gt/total:.4f}), smaller {lt} ({lt/total:.4f}); means {EHc:.1f} / {EHo:.1f}; correlation {cov/math.sqrt(va*vb):.4f}; max values {max(a for (a,b,bl,st) in pairs)} / {max(b for (a,b,bl,st) in pairs)}")
    ho = Counter()
    for (a, b, bl, st), c in pairs.items(): ho[b] += c
    print("  open end H = 1, 3, 5, 9, 15:", [ho.get(h, 0) for h in (1, 3, 5, 9, 15)], f"(n-2 = {n-2}, 5(n-3) = {5*(n-3)})")
