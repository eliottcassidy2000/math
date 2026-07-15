# opus-2026-07-15-S316 -- HYP-6935: parity is the tamer -- three exact probes.
# (i)  THE AXIS LEVEL-POPULATION LAW: populated x-values (score sequences,
#      Landau) vs mod-8/mod-16 residue laws and sum-of-squares representability.
# (ii) H mod 4 / mod 8 vs invariants (c3 = C(n,3) - sum C(s_v,2), x, z).
# (iii) THE STAIRCASE LOCKER: interval-containment zeta mod 2; Z(all-ones) =
#      gap = 2,3 mod 4 stripes (subtile count of gap-g tile = T_{g-1});
#      + is Z class-respecting? (transport census at n = 4, 5).
import sys, itertools
from collections import defaultdict
sys.path.insert(0, '04-computation')
from smith_diagram_of_the_metagraph_opus_S307 import build

# ---------- (i) level population
print("(i) AXIS LEVEL-POPULATION (score sequences via Landau):")
for n in range(4, 9):
    seqs = []
    def rec(prefix, rem_sum):
        k = len(prefix)
        if k == n:
            if rem_sum == 0: seqs.append(tuple(prefix))
            return
        lo = prefix[-1] if prefix else 0
        for s in range(lo, n):
            # Landau partial: sum of first k+1 >= C(k+1,2); prune by remaining
            pre = prefix + [s]
            ps = sum(pre)
            if ps < (k+1)*k//2: continue
            rest = n - k - 1
            if ps + rest*(n-1) < n*(n-1)//2: continue
            if ps > n*(n-1)//2: continue
            rec(pre, rem_sum)
    rec([], 0)
    xs = sorted({sum((2*s - (n-1))**2 for s in q) for q in seqs
                 if sum(q) == n*(n-1)//2 and
                 all(sum(q[:k+1]) >= (k+1)*k//2 for k in range(n))})
    mods8 = sorted({x % 8 for x in xs})
    mods16 = sorted({x % 16 for x in xs})
    # candidate residue-law gaps: which x with the right mod-8 residue in range
    # are NOT populated?
    cand = [x for x in range(min(xs), max(xs)+1)
            if x % 8 == (n % 8 if n % 2 == 0 else min(mods8))]
    holes = [x for x in cand if x not in xs] if n % 2 == 0 else []
    print(f"   n={n}: levels={xs if len(xs)<20 else str(xs[:12])+'...'}")
    print(f"        mod 8: {mods8}   mod 16: {mods16}"
          + (f"   holes-in-mod8-class: {holes}" if n % 2 == 0 else ""))

# ---------- (ii) H mod 4/8
print("\n(ii) H mod 4 / mod 8 vs c3 (n = 4..7):")
for n in range(4, 8):
    B = build(n)
    H_of, x_of, cls_of = B['H_of'], B['x_of'], B['cls_of']
    # c3 from a representative's scores: c3 = C(n,3) - sum C(s_v,2)
    # scores from x?? need scores: rebuild from a representative tiling
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
    rep = {}
    for t in range(1 << len(tiles)):
        c = cls_of[t]
        if c not in rep: rep[c] = t
    from math import comb
    rows = defaultdict(set)
    for c, t in rep.items():
        s = [0]*(n+1)
        for k in range(2, n+1): s[k] += 1
        for i, (xx, yy) in enumerate(tiles):
            if (t >> i) & 1: s[xx] += 1
            else: s[yy] += 1
        sc = [s[v] - 1 + 0 for v in range(1, n+1)]
        # careful: base path arcs k->k-1 give s[k]+=1 for k=2..n; vertex 1 gets 0
        sc = []
        wins = [0]*(n+1)
        for k in range(2, n+1): wins[k] += 1
        for i, (xx, yy) in enumerate(tiles):
            if (t >> i) & 1: wins[xx] += 1
            else: wins[yy] += 1
        sc = [wins[v] for v in range(1, n+1)]
        c3 = comb(n, 3) - sum(comb(w, 2) for w in sc)
        rows[(c3 % 2, c3 % 4)].add((H_of[c] % 4, H_of[c] % 8))
    print(f"   n={n}: (c3 mod 2, c3 mod 4) -> set of (H mod 4, H mod 8):")
    for k in sorted(rows): print(f"      c3={k}: {sorted(rows[k])}")

# ---------- (iii) the staircase locker
print("\n(iii) THE STAIRCASE LOCKER (interval zeta mod 2):")
for n in (4, 5, 6):
    tiles = [(x, y) for y in range(1, n-1) for x in range(n, y+1, -1)]
    m = len(tiles)
    tidx = {t: i for i, t in enumerate(tiles)}
    # subtile relation: [y',x'] subseteq [y,x] with gap >= 2
    sub = [[j for j, (x2, y2) in enumerate(tiles)
            if y1 <= y2 and x2 <= x1] for (x1, y1) in tiles]
    # zeta: (Zt)(tau) = xor of t over subtiles of tau (including tau)
    def Z(t):
        out = 0
        for i in range(m):
            b = 0
            for j in sub[i]: b ^= (t >> j) & 1
            out |= b << i
        return out
    ones = (1 << m) - 1
    zo = Z(ones)
    stripes = 0
    for i, (x1, y1) in enumerate(tiles):
        g = x1 - y1
        # subtile count = T_{g-1} = g(g-1)/2; odd iff g = 2,3 mod 4
        if (g*(g-1)//2) % 2 == 1: stripes |= 1 << i
    print(f"   n={n}: Z(all-ones) == gap-(2,3 mod 4)-stripes: {zo == stripes}")
    if n <= 5:
        B = build(n)
        cls_of = B['cls_of']
        preserved = sum(1 for t in range(1 << m) if cls_of[Z(t)] == cls_of[t])
        # class-functional test: does cls(t)=cls(t') => cls(Zt)=cls(Zt')?
        img = defaultdict(set)
        for t in range(1 << m): img[cls_of[t]].add(cls_of[Z(t)])
        func = all(len(v) == 1 for v in img.values())
        print(f"        Z fixes class for {preserved}/{1<<m} tilings; "
              f"class-functional: {func}; image-class-count histogram: "
              f"{sorted(set(len(v) for v in img.values()))}")
