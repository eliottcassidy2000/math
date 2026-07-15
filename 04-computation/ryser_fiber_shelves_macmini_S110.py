#!/usr/bin/env python3
"""ryser_fiber_shelves_macmini_S110.py -- mac-mini-2026-07-15-S110.
(A) THM-869 probes: the Morse shelves of the axis.
    stuck (no descending flip) <=> d in {0,+-2} with every -2 beating every +2
    (odd n; even n: none above floor). Verify the characterization exhaustively n=5,6,7
    (against the raw census), get fiber sizes & upset fractions, and CONSTRUCT an
    x=16 (a=2) stuck tournament at n=9 (predicted first shelf-2).
(B) Ryser-fiber monodromy: score fibers are connected under directed-3-cycle reversals
    (Ryser); alpha_1 mod 2 (= digit 1 of H) is fiber-constant at n=4, escapes at n=5.
    Per-reversal local-law hunt: Delta-alpha_1 mod 2 vs local statistics of the triangle.
"""
import sys
from itertools import combinations, permutations
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def all_tournaments(n):
    pairs = list(combinations(range(n), 2))
    m = len(pairs)
    for code in range(1 << m):
        adj = [[0] * n for _ in range(n)]
        for i, (u, v) in enumerate(pairs):
            if (code >> i) & 1: adj[u][v] = 1
            else: adj[v][u] = 1
        yield code, adj, pairs

def scores(adj, n): return [sum(row) for row in adj]

def alpha1(adj, n):
    """total number of directed odd cycles (simple), via DFS from each min vertex."""
    total = 0
    for a in range(n):
        # cycles with min vertex a
        stack = [(a, 1 << a, a, 0)]
        # iterative DFS counting closed walks that are simple cycles with min a
        # (path stored implicitly via bitmask; length parity via count)
        st = [(a, 1 << a, 0)]
        while st:
            v, mask, ln = st.pop()
            for w in range(a, n):
                if adj[v][w]:
                    if w == a and ln >= 2:
                        if (ln + 1) % 2 == 1: total += 1
                    elif not (mask >> w) & 1:
                        st.append((w, mask | (1 << w), ln + 1))
    return total

def is_stuck(adj, n, d):
    for u in range(n):
        for v in range(n):
            if u != v and adj[u][v] and d[u] - d[v] > 2:
                return False
    return True

print("(A) THM-869 -- the Morse shelves")
print("  characterization check: stuck <=> d in {0,+-2} and all (-2)->(+2) upsets:")
for n in (5, 6, 7):
    n_stuck = 0; n_char = 0; mism = 0
    fib_size = Counter(); stuck_by_score = Counter()
    for code, adj, pairs in all_tournaments(n):
        s = scores(adj, n)
        d = [2 * x - (n - 1) for x in s]
        x = sum(v * v for v in d)
        floor = 0 if n % 2 else n
        fib_size[tuple(sorted(s))] += 1
        if x == floor: continue
        stuck = is_stuck(adj, n, d)
        char = set(d) <= {0, 2, -2} and all(
            adj[u][v] for u in range(n) for v in range(n)
            if u != v and d[u] == -2 and d[v] == 2) and any(v != 0 for v in d)
        if stuck: n_stuck += 1; stuck_by_score[tuple(sorted(s))] += 1
        if char: n_char += 1
        if stuck != char: mism += 1
    print(f"  n={n}: stuck={n_stuck}, characterized={n_char}, mismatches={mism}")
    if n_stuck:
        for sc, c in stuck_by_score.items():
            print(f"        score {sc}: stuck {c} / fiber {fib_size[sc]}  (upset fraction {c}/{fib_size[sc]})")

print("\n  n=9, a=2 construction (predicted first x=16 shelf):")
n = 9
# lows L={0,1} (score 3), highs H={2,3} (score 5), mids M={4..8} (score 4)
adj = [[0] * n for _ in range(n)]
def setarc(a, b):
    adj[a][b] = 1; adj[b][a] = 0
# lows beat highs (upset saturation)
for l in (0, 1):
    for h in (2, 3): setarc(l, h)
# low 0 beats low 1; low 0 loses to all mids; low 1 beats mid 4, loses to mids 5..8
setarc(0, 1)
for mid in range(4, 9): setarc(mid, 0)
setarc(1, 4)
for mid in range(5, 9): setarc(mid, 1)
# highs: 2 beats 3; each high must reach score 5: 2 has {0?no} wins so far: beats 3 -> beat 4 mids
setarc(2, 3)
for mid in range(4, 8): setarc(2, mid)   # H0=2 beats mids 4..7 (score 1+4 = 5)
setarc(8, 2)                             # loses to mid 8
for mid in range(4, 9): setarc(3, mid)   # H1=3 beats ALL mids (score 5)
# mids among themselves: need each mid score 4 total; current mid scores:
cur = [sum(adj[v][w] for w in range(n) if w != v) for v in range(n)]
# mids 4..8 currently have wins from beating lows/highs; give each mid the right number
# of internal (mid-vs-mid) wins: target 4 - current
need = {mid: 4 - cur[mid] for mid in range(4, 9)}
# internal round-robin on 5 mids: score sequence must be need[4..8] summing to C(5,2)=10
mids = list(range(4, 9))
assert sum(need[m] for m in mids) == 10, (need, cur)
# realize a tournament on mids with given scores (small: brute force)
import itertools
mp = list(combinations(mids, 2))
done = False
for mcode in range(1 << len(mp)):
    tmp = {m: 0 for m in mids}
    arcs = []
    for i, (u, v) in enumerate(mp):
        if (mcode >> i) & 1: tmp[u] += 1; arcs.append((u, v))
        else: tmp[v] += 1; arcs.append((v, u))
    if all(tmp[m] == need[m] for m in mids):
        for a, b in arcs: setarc(a, b)
        done = True; break
assert done, "no internal mid tournament with needed scores"
s = scores(adj, n)
d = [2 * x - (n - 1) for x in s]
x = sum(v * v for v in d)
ok_t = all(adj[u][v] + adj[v][u] == 1 for u in range(n) for v in range(n) if u != v)
print(f"    constructed: scores={sorted(s)}  x={x}  tournament-valid={ok_t}  "
      f"STUCK={is_stuck(adj, n, d)}   (shelf level a=2, x=16: {'CONFIRMED' if x == 16 and is_stuck(adj, n, d) else 'FAIL'})")

print("\n(B) RYSER-FIBER MONODROMY (fibers = score VECTORS, per Ryser; n=4 baseline + n=5)")
for nn in (4, 5):
    tt = {}
    for code, adj, pairs in all_tournaments(nn):
        tt[code] = alpha1(adj, nn) % 2
    fl = ke = 0
    pp = list(combinations(range(nn), 2))
    pix = {p: i for i, p in enumerate(pp)}
    for code, adj, pairs in all_tournaments(nn):
        for a, b, c in combinations(range(nn), 3):
            tri = None
            if adj[a][b] and adj[b][c] and adj[c][a]: tri = ((a, b), (b, c), (c, a))
            elif adj[a][c] and adj[c][b] and adj[b][a]: tri = ((a, c), (c, b), (b, a))
            if tri:
                nc = code
                for (u, v) in tri: nc ^= (1 << pix[(min(u, v), max(u, v))])
                if tt[code] != tt[nc]: fl += 1
                else: ke += 1
    print(f"   n={nn}: reversal edges: parity-flipping {fl//2}, keeping {ke//2}"
          f"  -> digit 1 {'MONODROMY-CARRIED' if fl else 'fiber-constant'}")
n = 5
tourns = {}
for code, adj, pairs in all_tournaments(n):
    tourns[code] = (adj, tuple(scores(adj, n)), alpha1(adj, n))
fibers = defaultdict(list)
for code, (adj, sc, a1) in tourns.items(): fibers[sc].append(code)
pairs5 = list(combinations(range(n), 2))
pidx = {p: i for i, p in enumerate(pairs5)}
def triangles(adj):
    out = []
    for a, b, c in combinations(range(n), 3):
        # directed 3-cycle on {a,b,c}?
        if adj[a][b] and adj[b][c] and adj[c][a]: out.append((a, b, c))
        elif adj[a][c] and adj[c][b] and adj[b][a]: out.append((a, c, b))
    return out
nonconst = 0
edge_stats = Counter()
for sc, codes in fibers.items():
    pars = {c: tourns[c][2] % 2 for c in codes}
    if len(set(pars.values())) > 1: nonconst += 1
    # Ryser connectivity + parity-flip census on reversal edges
    seen = set(); comp = 0
    for c0 in codes:
        if c0 in seen: continue
        comp += 1
        stack = [c0]; seen.add(c0)
        while stack:
            c = stack.pop()
            adj = tourns[c][0]
            for tri in triangles(adj):
                a, b, cc = tri
                nc = c
                for (u, v) in ((a, b), (b, cc), (cc, a)):
                    i = pidx[(min(u, v), max(u, v))]
                    nc ^= (1 << i)
                # parity flip stat
                d_par = (tourns[c][2] + tourns[nc][2]) % 2
                # local statistic: number of odd cycles through at least one triangle vertex
                edge_stats[(d_par,)] += 1
                if nc not in seen and nc in pars:
                    seen.add(nc); stack.append(nc)
    if comp > 1: print(f"   fiber {sc}: NOT connected under 3-cycle reversals ({comp} comps)")
nonconst_sorted = len({tuple(sorted(sc)) for sc, codes in fibers.items()
                       if len({tourns[c][2] % 2 for c in codes}) > 1})
print(f"   score-vector fibers: {len(fibers)}; vector-fibers with parity split: {nonconst};"
      f" sorted-classes affected: {nonconst_sorted} (THM-466: 1 of 9)")
print(f"   Ryser connectivity: every score-vector fiber connected (no disconnection lines above)")
flips = edge_stats[(1,)] // 2; keeps = edge_stats[(0,)] // 2
print(f"   reversal edges: parity-flipping {flips}, parity-keeping {keeps} "
      f"(digit 1 is monodromy-carried iff flips > 0: {'YES' if flips else 'no'})")

# local-law hunt: Delta-alpha1 mod 2 vs simple local counts
print("\n   local-law hunt (n=5, all reversal edges):")
rows = []
for code, (adj, sc, a1) in tourns.items():
    for tri in triangles(adj):
        a, b, c = tri
        nc = code
        for (u, v) in ((a, b), (b, c), (c, a)):
            i = pidx[(min(u, v), max(u, v))]
            nc ^= (1 << i)
        a1b = tourns[nc][2]
        dpar = (a1 + a1b) % 2
        tv = {a, b, c}
        # candidate local statistics in T:
        s5 = 0; s5_1arc = 0; s5_2arc = 0
        for perm in permutations(range(n)):
            if perm[0] != 0: continue
            if all(adj[perm[i]][perm[(i + 1) % n]] for i in range(n)):
                arcs = {(perm[i], perm[(i + 1) % n]) for i in range(n)}
                tri_arcs = {(a, b), (b, c), (c, a)}
                k = len(arcs & tri_arcs)
                if k >= 1: s5 += 1
                if k == 1: s5_1arc += 1
                if k == 2: s5_2arc += 1
        deg_in_tri = sum(1 for v in tv for w in range(n) if w not in tv and adj[v][w])
        rows.append((dpar, s5 % 2, s5_1arc % 2, s5_2arc % 2, deg_in_tri % 2))
labels = ["#5cyc >=1 tri-arc", "#5cyc exactly-1 arc", "#5cyc exactly-2 arcs", "tri out-degree"]
for col, lab in enumerate(labels, start=1):
    tab = Counter((r[0], r[col]) for r in rows)
    law = all(k[0] == k[1] for k in tab) or all(k[0] != k[1] for k in tab)
    print(f"   Delta-par vs ({lab}) mod 2: {dict(tab)}  LAW: {law}")
print("\nDONE")
