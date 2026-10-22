# Wave-1 analysis of the exhaustive per-class OCF data (ocf_n5..9.txt): realized small cycle-packing profiles,
# H mod 3 against Aut and the Euler-characteristic term, and the Steiner signature.  mac-mini-2026-10-22-sierp.
import re, collections, gzip, os
def opn(n):
    f = f'ocf_n{n}.txt'
    return gzip.open(f + '.gz', 'rt') if os.path.exists(f + '.gz') else open(f)
rows = {}
for n in range(5, 10):
    rows[n] = []
    for line in opn(n):
        m = re.match(r'(\S+) H=(\d+) a=(\d+),(\d+),(\d+),(\d+) Im1=(-?\d+) J=(-?\d+) S=(\d+),(\d+),(\d+) aut=(\d+) sc=(\d+)', line)
        g = m.groups()
        rows[n].append(dict(code=g[0], H=int(g[1]), a=tuple(map(int, g[2:6])), Im1=int(g[6]), J=int(g[7]),
                            S=tuple(map(int, g[8:11])), aut=int(g[11]), sc=g[12]))
print("== realized profiles (alpha1, alpha2, alpha3) with H <= 45, all n <= 9 ==")
prof = collections.defaultdict(set)
for n, R in rows.items():
    for r in R:
        if r['H'] <= 45: prof[r['H']].add(r['a'][:3])
for H in sorted(prof): print(H, sorted(prof[H]))
print("\n== all candidate profiles for H = 7 and 21 (alpha_k >= 0, alpha_2 <= C(alpha_1,2), triples need alpha_2 >= 3) ==")
from math import comb
for target in (7, 21):
    c = []
    for a3 in range(0, 3):
        for a2 in range(0, 11):
            for a1 in range(0, 11):
                if 1 + 2*a1 + 4*a2 + 8*a3 != target: continue
                if a2 > comb(a1, 2): continue
                if a3 > 0 and (a2 < 3 * a3 or a1 < 3): continue
                c.append((a1, a2, a3))
    print(target, c)
print("\n== H mod 3 against |Aut| (n <= 8) ==")
for n in range(5, 9):
    tab = collections.Counter()
    for r in rows[n]: tab[(r['H'] % 3, r['aut'] % 3 == 0)] += 1
    print(n, dict(tab))
print("\n== distribution of Im1 = I(Omega,-1) ==")
for n in range(5, 10):
    c = collections.Counter(r['Im1'] for r in rows[n])
    print(n, "min", min(c), "max", max(c), "values with counts (first 25):", sorted(c.items())[:25])
