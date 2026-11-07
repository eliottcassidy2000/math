#!/usr/bin/env python3
"""value_tables.py (klein-2026-11-06): the per-value tables of the Hamiltonian-path spectrum.
usage: value_tables.py census_dir nmax > least_order.txt
For every odd h up to max H(nmax) prints one line
  h  n(h)  m(h)  d  strong_window  level_windows
where n(h) = least order attaining h ('>nmax' if none up to nmax), m(h) = min{n : max H(n) >= h} (max H(11) = 95095 from
OEIS A038375 is used when needed), d = n(h) - m(h), strong_window = the orders n <= nmax with a strong realiser with the
class count in brackets, marked 'complete' when h < f(nmax+1) (no strong realiser of a larger order can exist, Busch), and
level_windows = for each order n <= nmax attaining h the range kmin-kmax of the levels (numbers of cyclic triangles) of its
realisers, as n:kmin-kmax.  Also prints the level-width table (classes, strong classes and distinct values per level)."""
import sys
from collections import defaultdict
import gzip, os
def open_census(path):
    """open census_<n>.txt, or its gzip-compressed twin census_<n>.txt.gz, as text"""
    if os.path.exists(path): return open(path)
    return gzip.open(path + ".gz", "rt")
D = sys.argv[1]; NMAX = int(sys.argv[2])
def busch(n):
    if n == 3: return 3
    if n == 4: return 5
    k, r = divmod(n, 3); return {1: 5**k, 2: 9 * 5**(k - 1), 0: 3 * 5**(k - 1)}[r]
cnt = {}
for n in range(3, NMAX + 1):
    T = defaultdict(lambda: defaultdict(lambda: [0, 0]))
    with open_census(f"{D}/census_{n}.txt") as fh:
        for line in fh:
            a = line.split()
            if a and a[0] == 'T': T[int(a[1])][int(a[2])][int(a[3])] += int(a[4])
    cnt[n] = T
maxH = {n: max(cnt[n]) for n in cnt}; maxH[1] = maxH[2] = 1; maxH[11] = max(maxH.get(11, 0), 95095)
def m_of(h):
    for n in range(1, 12):
        if maxH[n] >= h: return n
    return None
print(f"# least orders, strong windows and level windows of every odd h <= {maxH[NMAX]} from the exact censuses of orders 3..{NMAX}")
print("# h n(h) m(h) n(h)-m(h) strong_window level_windows")
lim = busch(NMAX + 1)
for h in range(1, maxH[NMAX] + 1, 2):
    ords = [n for n in range(3, NMAX + 1) if h in cnt[n]]
    nh = ords[0] if ords else None
    if h == 1: nh = 1
    mh = m_of(h)
    sw = [(n, sum(v[1] for v in cnt[n][h].values())) for n in ords]
    sw = [(n, c) for n, c in sw if c]
    swt = " ".join(f"{n}[{c}]" for n, c in sw) if sw else "-"
    if h < lim: swt += " complete"
    lw = " ".join(f"{n}:{min(cnt[n][h])}-{max(cnt[n][h])}" for n in ords)
    print(f"{h} {nh if nh else '>' + str(NMAX)} {mh} {(nh - mh) if nh else '?'} {swt} {lw}")
print("# level widths: for each order n, the number of classes at each level k = 0..Lambda(n), then the strong classes, then the distinct values")
for n in range(3, NMAX + 1):
    lev = defaultdict(lambda: [0, 0, set()])
    for h in cnt[n]:
        for k, (a, b) in cnt[n][h].items(): lev[k][0] += a + b; lev[k][1] += b; lev[k][2].add(h)
    K = max(lev)
    print(f"# n = {n} classes:  " + " ".join(str(lev[k][0]) for k in range(K + 1)))
    print(f"# n = {n} strong:   " + " ".join(str(lev[k][1]) for k in range(K + 1)))
    print(f"# n = {n} distinct: " + " ".join(str(len(lev[k][2])) for k in range(K + 1)))
