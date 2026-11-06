#!/usr/bin/env python3
"""level_tables.py (klein-2026-11-06): the Hamiltonian-path spectrum resolved by the level k = c3 (number of cyclic
triangles) of the tournament metagraph G_n, from the exact censuses census_<n>.txt.
usage: level_tables.py census_dir nmax
Prints, for each order n:
  (1) per level k: classes, strong classes, distinct values, min h, max h, and whether the level range is solid;
  (2) the level windows K(n,h) = {k : a class at level k has h paths}: interval or not, lengths;
  (3) where the new values N(n) = S(n) - S(n-1) sit: the least level carrying a new value, and the fraction of new
      values at each level;
  (4) the pigeonhole frontier: the least x such that fewer classes have h >= x than there are odd numbers in [x, max],
      against the onset of the top band of unattained values;
  (5) the (n, k) support of each small value h: for h < 160, the levels at each order (the staircase).
Also the per-order sequence of the minimum h at each level, and of the maximum h at each level."""
import sys
from collections import defaultdict
D = sys.argv[1]; NMAX = int(sys.argv[2])
cnt = {}
for n in range(3, NMAX + 1):
    T = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    with open(f"{D}/census_{n}.txt") as fh:
        for line in fh:
            a = line.split()
            if a and a[0] == 'T': T[int(a[1])][int(a[2])][int(a[3])] += int(a[4])
    cnt[n] = T
S = {n: set(cnt[n]) for n in cnt}; S[2] = {1}
def ranges(vals):
    vals = sorted(vals); out = []; i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[j + 1] == vals[j] + 2: j += 1
        out.append(f"{vals[i]}" if i == j else f"{vals[i]}..{vals[j]}({(j - i) + 1})"); i = j + 1
    return " ".join(out)
for n in range(5, NMAX + 1):
    T = cnt[n]; maxH = max(T)
    lev = defaultdict(lambda: {'cls': 0, 'str': 0, 'vals': set()})
    for h in T:
        for k, (a, b) in T[h].items():
            lev[k]['cls'] += a + b; lev[k]['str'] += b; lev[k]['vals'].add(h)
    kmax = max(lev)
    print(f"\n=== order {n}: levels 0..{kmax}, max H = {maxH} ===")
    print("(1) k | classes | strong | distinct | min h | max h | odd values in [min,max] | solid?")
    for k in range(kmax + 1):
        L = lev[k]; v = L['vals']; lo, hi = min(v), max(v); odd = (hi - lo) // 2 + 1
        solid = all((x in v) for x in range(lo, hi + 1, 2) if x not in (7, 21))
        print(f"  {k:2d} | {L['cls']:8d} | {L['str']:8d} | {len(v):6d} | {lo:6d} | {hi:6d} | {odd:6d} | {'solid' if solid else 'gaps %d' % (odd - len(v) - sum(1 for x in (7, 21) if lo <= x <= hi))}")
    mins = [min(lev[k]['vals']) for k in range(kmax + 1)]; maxs = [max(lev[k]['vals']) for k in range(kmax + 1)]
    print(f"  min h by level: {mins}")
    print(f"  max h by level: {maxs}")
    print(f"  min h nondecreasing in k: {all(mins[i] <= mins[i + 1] for i in range(kmax))}; max h nondecreasing in k: {all(maxs[i] <= maxs[i + 1] for i in range(kmax))}")
    K = {h: sorted(T[h]) for h in T}
    nonint = [(h, K[h]) for h in K if K[h] != list(range(K[h][0], K[h][-1] + 1))]
    lens = defaultdict(int)
    for h in K: lens[len(K[h])] += 1
    print(f"(2) level windows K(n,h): {len(K)} values; not intervals: {len(nonint)}" + (f" e.g. {nonint[:6]}" if nonint else "") + f"; window-length distribution {dict(sorted(lens.items()))}")
    new = S[n] - S[n - 1]
    bylev = defaultdict(int); least_new = min(min(K[h]) for h in new)
    for h in new:
        for k in K[h]: bylev[k] += 1
    old_top = defaultdict(int)
    for h in S[n] & S[n - 1]:
        for k in K[h]: old_top[k] += 1
    print(f"(3) new values: {len(new)}; least level carrying a new value: {least_new}; values per level (new / inherited): " +
          " ".join(f"{k}:{bylev[k]}/{old_top[k]}" for k in range(kmax + 1)))
    # (4) pigeonhole frontier
    classes_ge = {}; tot = 0
    for h in sorted(T, reverse=True):
        tot += sum(a + b for a, b in T[h].values()); classes_ge[h] = tot
    x_front = None
    for h in sorted(T):
        if classes_ge[h] < (maxH - h) // 2 + 1: x_front = h; break
    un = [h for h in range(1, maxH, 2) if h not in S[n] and h not in (7, 21)]
    print(f"(4) pigeonhole frontier x = {x_front} (classes with h >= x: {classes_ge.get(x_front)}, odd values in [x, max]: {(maxH - x_front) // 2 + 1 if x_front else None}); " +
          f"first unattained value above 21: {un[0] if un else None}; last: {un[-1] if un else None}; strong classes above the frontier: {sum(b for h in T if h >= (x_front or 0) for a, b in T[h].values())}")
    # top band versus levels: which levels have classes above the first gap
    if un:
        above = {k: sum(a + b for h in T if h > un[0] for kk, (a, b) in T[h].items() if kk == k) for k in range(kmax + 1)}
        print(f"    classes with h above the first gap, by level: {{{', '.join(f'{k}:{c}' for k, c in above.items() if c)}}}")
print("\n(5) the (n, k) support of each value h < 160: for each order n, the levels carrying h  [n: k-list]")
for h in range(3, 160, 2):
    row = []
    for n in range(3, NMAX + 1):
        if h in cnt[n]: row.append(f"{n}:{min(cnt[n][h])}-{max(cnt[n][h])}")
    if row: print(f"  {h}: " + " ".join(row))
print("\n(6) for each order, the values attained only at the top level and only at the bottom level")
for n in range(5, NMAX + 1):
    T = cnt[n]; kmax = max(k for h in T for k in T[h])
    onlytop = sorted(h for h in T if set(T[h]) == {kmax}); onlybot = sorted(h for h in T if set(T[h]) == {0})
    print(f"  n = {n}: only at level {kmax}: {ranges(onlytop)}; only at level 0: {ranges(onlybot)}")

print("\n(7) zone laws: the level floor phi_n(k) = min h at level k and ceiling psi_n(k) = max h at level k")
def Lam(m): return (m**3 - m) // 24 if m % 2 else (m**3 - 4 * m) // 24   # the maximum number of cyclic triangles at order m
maxH = {n: max(cnt[n]) for n in cnt}
for n in range(5, NMAX + 1):
    T = cnt[n]; lev = defaultdict(set)
    for h in T:
        for k in T[h]: lev[k].add(h)
    kmax = max(lev); phi = {k: min(lev[k]) for k in lev}; psi = {k: max(lev[k]) for k in lev}
    L = Lam(n - 1)
    inh_levels = sorted(k for k in lev if any(h in S[n - 1] for h in lev[k]))
    new_levels = sorted(k for k in lev if any(h not in S[n - 1] for h in lev[k]))
    un_prev = [h for h in range(1, maxH[n - 1], 2) if h not in S[n - 1] and h not in (7, 21)] if n - 1 >= 5 else []
    u_prev = un_prev[0] if un_prev else maxH[n - 1] + 2
    print(f"  n = {n}: Lambda(n-1) = {L}, kmax = {kmax}; top inherited level = {max(inh_levels)}, least new level = {min(new_levels)}; "
          f"phi_n(Lambda(n-1)+1) = {phi.get(L + 1)} vs max H(n-1) = {maxH[n - 1]}; psi_n(least new level - 1) = {psi.get(min(new_levels) - 1)} vs first gap of n-1 = {u_prev}; psi_n(least new level) = {psi[min(new_levels)]}")
    print(f"    phi_n(k): {[phi[k] for k in range(kmax + 1)]}")
    print(f"    psi_n(k): {[psi[k] for k in range(kmax + 1)]}")
print("  stabilisation of the level floor in n (phi_n(k) for n = 6..nmax, k <= 12):")
for k in range(0, 13):
    row = []
    for n in range(6, NMAX + 1):
        lev = defaultdict(set)
        for h in cnt[n]:
            if k in cnt[n][h]: lev[k].add(h)
        row.append(min(lev[k]) if lev[k] else None)
    print(f"    k = {k:2d}: {row}")

print("\n(8) the top band against the inter-level gap: the ceiling of level Lambda(n)-1, the floor of level Lambda(n), the band's blocks")
for n in range(6, NMAX + 1):
    T = cnt[n]; lev = defaultdict(set)
    for h in T:
        for k in T[h]: lev[k].add(h)
    kmax = max(lev); un = [h for h in range(1, maxH[n], 2) if h not in S[n] and h not in (7, 21)]
    psi1 = max(lev[kmax - 1]); phi0 = min(lev[kmax]); cls_top = sum(sum(T[h][kmax]) for h in lev[kmax])
    gap_block = [h for h in un if psi1 < h < phi0]
    print(f"  n = {n}: psi(Lambda-1) = {psi1}, phi(Lambda) = {phi0}, inter-level gap holds {len(gap_block)} of the {len(un)} unattained values; top level: {cls_top} classes, {len(lev[kmax])} distinct values over {(maxH[n] - phi0) // 2 + 1} odd values; "
          f"unattained values inside the top level's range: {sum(1 for h in un if h >= phi0)}; unattained below psi(Lambda-1): {sum(1 for h in un if h <= psi1)} (first {un[0] if un else None})")
print("\n(9) the level floor minimisers: for k <= 17, phi(k) with the orders (and class counts, strong counts) attaining it at orders <= nmax")
for k in range(0, 18):
    rows = []
    for n in range(3, NMAX + 1):
        vals = {h: cnt[n][h][k] for h in cnt[n] if k in cnt[n][h]}
        if vals:
            m = min(vals); rows.append((n, m, vals[m][0] + vals[m][1], vals[m][1]))
    if rows:
        best = min(r[1] for r in rows)
        print(f"  k = {k:2d}: phi = {best}; by order (n: min h [classes, strong]): " + " ".join(f"{n}:{m}[{c},{s}]" for n, m, c, s in rows))
print("\n(10) Moon's bound: the least level with a strong class at each order (expected n - 2)")
print("  " + ", ".join(f"n={n}: {min(k for h in cnt[n] for k, v in cnt[n][h].items() if v[1])}" for n in range(3, NMAX + 1)))
