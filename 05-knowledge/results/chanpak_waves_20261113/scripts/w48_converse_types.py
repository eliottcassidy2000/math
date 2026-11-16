"""Isomorphism classes of tournaments of order n = 5, 6, 7 grouped by H, and their converse structure
(s self-converse classes, t converse pairs).  'Platonic' type = (1, 2) (five solids: one self-dual + two dual pairs);
'tiling' type = (1, 1) (three flat tilings: one self-dual + one pair)."""
import itertools, sys
sys.path.insert(0, '.')
from w45_hk import H_of
def run(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]; m = len(pairs)
    idx = {pr: t for t, pr in enumerate(pairs)}
    perms = list(itertools.permutations(range(n)))
    seen = bytearray(1 << m); res = {}
    def relabel(bits, p):
        b = 0
        for t, (i, j) in enumerate(pairs):
            x = bits >> t & 1          # 1 means i -> j
            a, c = p[i], p[j]
            if a < c: b |= x << idx[(a, c)]
            else: b |= (1 - x) << idx[(c, a)]
        return b
    full = (1 << m) - 1
    classes = []
    for bits in range(1 << m):
        if seen[bits]: continue
        orb = set(relabel(bits, p) for p in perms)
        for o in orb: seen[o] = 1
        classes.append((bits, orb))
    for bits, orb in classes:
        out = [0] * n
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        h = H_of(n, out)
        conv = full ^ bits                       # reversing every arc
        selfc = conv in orb
        res.setdefault(h, [0, 0])
        if selfc: res[h][0] += 1
        else: res[h][1] += 1
    print(f"n={n}: {len(classes)} classes")
    plat = [h for h, (s, t2) in sorted(res.items()) if s == 1 and t2 == 4]
    til = [h for h, (s, t2) in sorted(res.items()) if s == 1 and t2 == 2]
    print("   counts with 1 self-converse class + 2 converse pairs (Platonic type):", plat)
    print("   counts with 1 self-converse class + 1 converse pair (tiling type):", til)
    print("   full table H: (self-converse, pairs):", {h: (s, t2 // 2) for h, (s, t2) in sorted(res.items())})
for n in (5, 6, 7):
    run(n)
