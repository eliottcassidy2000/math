"""S51 audit c3: THM-4748 (6), three reversed arcs, by the auditor's own code (alib51).
(a) every triple of back arcs with spans >= 2 on m = 5..13 vertices (exact DP): configurations, distinct H, values = 7
    (mod 8), the largest odd part of H + 1 among them, first order at which each 2^k - 1 occurs;
(b) the same with ALL spans (all tournaments with fa <= 3 on m vertices), to see whether span-1 arcs add repunits;
(c) a one-line family proving that with three reversed arcs the values = 7 (mod 8) are NOT on finitely many lines
    c*2^j - 1: a crossing pair with p = q = 0 (H = 2^(r+2) + 2*3^r - 1 = 5 mod 8 for odd r) disjoint from a C_3 block:
    H = 3(2^(r+2) + 2*3^r - 1) = 7 (mod 8), and the odd part of H + 1 is at least 3^(r+1)/(4(r+1)) by LTE.
Prints ALL CHECKS PASSED or the failures."""
import itertools, time
import numpy as np
import numba as nb
import alib51 as L

FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


@nb.njit(parallel=True, cache=True)
def triples_H(m, arcs, combos):
    base = np.zeros(m, dtype=np.int64)
    for i in range(m):
        for j in range(i + 1, m):
            base[i] |= 1 << j
    out = np.zeros(combos.shape[0], dtype=np.int64)
    for t in nb.prange(combos.shape[0]):
        adj = base.copy()
        for k in range(combos.shape[1]):
            y = arcs[combos[t, k], 0]; x = arcs[combos[t, k], 1]
            adj[x] &= ~(1 << y)
            adj[y] |= 1 << x
        out[t] = L.ham(m, adj)
    return out


def oddpart(x):
    while x % 2 == 0:
        x //= 2
    return x


mers = {2 ** k - 1: k for k in range(2, 64)}
sess = {5: (20, 5, 1, 1), 6: (120, 14, 3, 3), 7: (455, 44, 8, 13), 8: (1330, 113, 24, 37), 9: (3276, 256, 66, 97),
        10: (7140, 527, 128, 343), 11: (14190, 1034, 239, 829), 12: (26235, 1907, 413, 3737), 13: (45760, 3320, 658, 7267)}
first = {}; firstall = {}
mine = {}
print('(a) three back arcs with spans >= 2; (b) all spans (every tournament with fa <= 3)')
for m in range(5, 14):
    arcs2 = np.array([(y, x) for x in range(m) for y in range(x + 2, m)], dtype=np.int64)
    combos = np.array(list(itertools.combinations(range(len(arcs2)), 3)), dtype=np.int64)
    H = triples_H(m, arcs2, combos)
    vals = set(int(h) for h in H)
    v7 = [h for h in vals if h % 8 == 7]
    mine[m] = (len(combos), len(vals), len(v7), max(oddpart(h + 1) for h in v7))
    for h in sorted(vals):
        if h in mers and h not in first:
            first[h] = m
    arcs1 = np.array([(y, x) for x in range(m) for y in range(x + 1, m)], dtype=np.int64)
    allv = set()
    for r in (1, 2, 3):
        cb = np.array(list(itertools.combinations(range(len(arcs1)), r)), dtype=np.int64)
        allv |= set(int(h) for h in triples_H(m, arcs1, cb))
    for h in sorted(allv):
        if h in mers and h not in firstall:
            firstall[h] = m
    print(f'     m = {m:2d}: spans >= 2: {mine[m]} (configs, distinct H, # = 7 mod 8, max odd part of H+1); session {sess[m]}; '
          f'fa <= 3 world: {len(allv)} distinct H [{time.time() - T0:.0f}s]', flush=True)
check(all(mine[m] == sess[m] for m in sess), 'three_arc_census.out reproduced at every m = 5..13 (configurations, distinct H, count = 7 mod 8, largest odd part)')
print(f'     first order of 2^k - 1 (three arcs, spans >= 2): {dict(sorted(first.items()))}')
print(f'     first order of 2^k - 1 (fa <= 3, all spans):      {dict(sorted(firstall.items()))}')
check(first == {15: 5, 31: 6, 63: 8, 127: 8, 255: 9, 511: 11, 1023: 10, 2047: 13, 32767: 13}, 'repunit first orders match the session table')
check(all(2 ** k - 1 in first for k in range(4, 12)) and 2 ** 15 - 1 in first and not any(2 ** k - 1 in first for k in (12, 13, 14)),
      'three arcs reach 2^k - 1 for 4 <= k <= 11 and k = 15 at m <= 13, and NOT k = 12, 13, 14 (the list in (6) is exact for m <= 13)')
check({h: v for h, v in firstall.items() if h >= 15} == first, 'allowing span-1 arcs and fewer arcs (all fa <= 3 tournaments) changes no first-occurrence order of 2^k - 1, k >= 4')

print('(c) an explicit three-arc family with unbounded odd parts of H + 1 (H = 7 mod 8)')
ok = True; ods = []
for r in range(1, 41, 2):
    hc = 2 ** (r + 2) + 2 * 3 ** r - 1
    h = 3 * hc
    ok &= hc % 8 == 5 and h % 8 == 7
    od = oddpart(h + 1)
    ods.append(od)
    if r >= 5:                                 # then v_2(3*2^(r+1) + 3^(r+1) - 1) = v_2(r+1) + 2 exactly (LTE)
        ok &= od * 4 * (r + 1) >= 3 ** (r + 1)
for r in (1, 3, 5, 7):
    m = r + 7                                  # crossing block on {0..r+3} (arcs (r+2 -> 0), (r+3 -> 1)) + C_3 on the last 3
    adj = L.reverse_arcs(m, [(r + 2, 0), (r + 3, 1), (m - 1, m - 3)])
    ok &= int(L.ham(m, adj)) == 3 * (2 ** (r + 2) + 2 * 3 ** r - 1)
check(ok, f'family H = 3(2^(r+2) + 2*3^r - 1), r odd: = 7 mod 8, DP-confirmed for r = 1,3,5,7, odd parts of H+1 {ods[:8]}... '
          f'are >= 3^(r+1)/(4(r+1)) for r >= 5 (so "no finitely many lines" is PROVABLE, not only FINITE-EXACT)')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
