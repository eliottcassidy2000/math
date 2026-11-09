"""Aggregate the merged line multigraph (pairs_<n>.txt from line_census.c) by node classifications.
Classifications (all functions of the isomorphism class, from alpha_<n>.txt):
  nu       packing number = max number of vertex-disjoint odd cycles (degree of the packing polynomial)
  carry    carry-free (nu <= 1, intercyclic) vs carried (nu >= 2)
  comps    number of strong components
  H4       H mod 4 (= 1 + 2 alpha_1 mod 4, THM-466)
  sc       self-converse or not
  fc       first creator (H new at order n) or not
usage: python line_classes.py n CENSUS_DIR PAIRS_FILE
"""
import sys
import os
import collections


def load(d, n):
    inv = {}
    for line in open(os.path.join(d, f"alpha_{n}.txt")):
        f = line.split()
        code = int(f[0]); H = int(f[1]); strong = int(f[2]); sc = int(f[3])
        a = [int(x) for x in f[4:8]]
        nu = max([k + 1 for k in range(4) if a[k] > 0], default=0)
        comps = int(f[10])
        inv[code] = dict(H=H, strong=strong, sc=sc, a=a, nu=nu, comps=comps)
    return inv


def main():
    n = int(sys.argv[1]); d = sys.argv[2]; pf = sys.argv[3]
    inv = load(d, n)
    prev = set()
    if n > 1:
        for line in open(os.path.join(d, f"alpha_{n - 1}.txt")):
            prev.add(int(line.split()[1]))
    for c in inv.values():
        c['fc'] = int(c['H'] not in prev)
        c['carry'] = 'free' if c['nu'] <= 1 else 'carried'
        c['H4'] = c['H'] % 4
    tables = {k: collections.Counter() for k in ('nu', 'carry', 'comps', 'H4', 'sc', 'fc')}
    total = collections.Counter()
    nu_sum_min = {}
    for line in open(pf):
        a, b, blue, cnt = line.split()
        A, B, blue, cnt = inv[int(a)], inv[int(b)], int(blue), int(cnt)
        total[blue] += cnt
        for k in tables:
            x, y = sorted((A[k], B[k]), key=str)
            tables[k][(x, y, 'blue' if blue else 'black')] += cnt
        s = A['nu'] + B['nu']
        nu_sum_min[s] = nu_sum_min.get(s, 0) + cnt
    print(f"n = {n}: lines {sum(total.values())} (blue {total[1]}, black {total[0]}); floor(n/3) = {n // 3}")
    print(f"  distribution of nu(t) + nu(t+1) over lines: {dict(sorted(nu_sum_min.items()))}")
    for k, tab in tables.items():
        print(f"  by {k}:")
        for key in sorted(tab, key=str):
            print(f"    {key}: {tab[key]}")


if __name__ == "__main__":
    main()
