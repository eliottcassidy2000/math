"""Independent verification of the band escapers band_3(13) and band_4(16) (genesis lane, mac-mini-2026-11-09).

Pure Python; shares no code with tprop.c / wit.c.  For each band B = band_s(n):
  * H(B) and H(B - v) for every v by a subset DP written here;
  * membership H(B) in S(n-1) and H(B - v) in S(n-2): either via klein's census (orders <= 11, read from the
    census T-lines) or via an explicit witness tournament from the witness files (its H recomputed here);
  * lambda(B) by max-flow and primality by module closure (band_lemmas.py);
  * twin distance (fewest distinguishers of a vertex pair).
An escaper of THM-4942's grammar: not a first creator, no first-creator deletion, strong, prime, lambda >= 3.
"""
import os
import sys
from band_lemmas import band_out, lam, is_prime
from verify_escapers import census_values

HERE = os.path.dirname(os.path.abspath(__file__))


def ham_from_out(out_sets, verts):
    vs = sorted(verts)
    k = len(vs)
    if k <= 1:
        return 1
    pos = {v: i for i, v in enumerate(vs)}
    loc = [0] * k
    for v in vs:
        for w in out_sets[v]:
            if w in pos:
                loc[pos[v]] |= 1 << pos[w]
    full = (1 << k) - 1
    dp = [None] * (1 << k)
    for a in range(k):
        dp[1 << a] = {a: 1}
    total = 0
    for m in range(1, 1 << k):
        row = dp[m]
        if not row:
            continue
        if m == full:
            total = sum(row.values())
            break
        for a, c in row.items():
            nx = loc[a] & ~m
            while nx:
                b = (nx & -nx).bit_length() - 1
                nx &= nx - 1
                mm = m | (1 << b)
                d = dp[mm]
                if d is None:
                    d = dp[mm] = {}
                d[b] = d.get(b, 0) + c
        dp[m] = None
    return total


def code_out(code, n):
    out = [set() for _ in range(n)]
    p = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[p] == '1':
                out[i].add(j)
            else:
                out[j].add(i)
            p += 1
    return out


def load_witnesses(files):
    wit = {}
    for fn in files:
        for line in open(os.path.join(HERE, fn)):
            p = line.split()
            if len(p) == 3 and p[0] == 'W':
                h, code = int(p[1]), p[2]
                n = int(round((1 + (1 + 8 * len(code)) ** 0.5) / 2))
                assert n * (n - 1) // 2 == len(code)
                wit[(n, h)] = code
    return wit


def member(m, h, census, wit):
    if m <= 11:
        return h in census[m], 'census'
    for mm in range(m, 1, -1):          # a witness of smaller order also certifies (pad with a sink: S(mm) <= S(m))
        if (mm, h) in wit:
            hv = ham_from_out(code_out(wit[(mm, h)], mm), range(mm))
            return hv == h, 'witness order %d (H recomputed %d)' % (mm, hv)
        if mm <= 11 and h in census[mm]:
            return True, 'census order %d' % mm
    return None, 'no certificate'


def twin(out, n):
    return min(sum(1 for x in range(n) if x not in (u, w) and ((u in out[x]) != (w in out[x])))
               for u in range(n) for w in range(u + 1, n))


def main():
    census = {m: census_values(m) for m in (10, 11)}
    wit = load_witnesses(['wit_band3_13_H.out', 'wit_band4_16_dels.out', 'wit_band4_16_H.out'])
    all_ok = True
    for s, n in ((3, 13), (4, 16)):
        out = band_out(n, s)
        H = ham_from_out(out, range(n))
        dels = [ham_from_out(out, [x for x in range(n) if x != v]) for v in range(n)]
        okH, whyH = member(n - 1, H, census, wit)
        dres = [member(n - 2, d, census, wit) for d in dels]
        L = lam(n, out)
        pr = is_prime(n, out)
        tw = twin(out, n)
        esc = okH and all(r[0] for r in dres) and L >= 3 and pr
        all_ok &= bool(esc)
        print('band_%d(%d): H=%d in S(%d): %s [%s]' % (s, n, H, n - 1, okH, whyH))
        print('   deletions:', sorted(set(dels)))
        print('   all deletions in S(%d): %s  [%s]' % (n - 2, all(r[0] for r in dres), sorted(set(r[1] for r in dres))))
        print('   lambda=%d prime=%s twin=%d  ESCAPER=%s' % (L, pr, tw, esc))
    print('# all band escapers verified:', all_ok)


if __name__ == '__main__':
    main()
