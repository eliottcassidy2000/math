"""Independent verification of order-12 escapers (pure Python; shares no code with tprop.c / grammar.py).

For each escaper code (12 vertices):
  (1) H(T) by a forward subset DP written here, and H(T - v) for all v;
  (2) H(T) in S(11) and every H(T - v) in S(10), with S(10), S(11) read from klein's census T-lines
      (census_10.txt, census_11.txt.gz), not from the unattained lists;
  (3) strong connectivity and arc-strong connectivity lambda by enumerating all cuts;
  (4) primality by module closure: for every pair {u, w} the smallest module containing it is V.
Usage: python3 verify_escapers.py ext12_escaper_canon.txt
"""
import gzip
import os
import sys

KLEIN = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '05-knowledge', 'results',
                     'path_spectrum_two_angles_metagraph_20261106', 'data')


def census_values(n):
    fn = os.path.join(KLEIN, 'census_%d.txt' % n)
    op = open
    if not os.path.exists(fn):
        fn += '.gz'
        op = gzip.open
    vals = set()
    with op(fn, 'rt') as f:
        for line in f:
            if line.startswith('T '):
                vals.add(int(line.split()[1]))
    return vals


def adj_from(code, n):
    out = [0] * n
    p = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code[p] == '1':
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            p += 1
    return out


def ham(out, verts):
    vs = [v for v in range(len(out)) if verts >> v & 1]
    k = len(vs)
    if k <= 1:
        return 1
    loc = [0] * k
    for a in range(k):
        for b in range(k):
            if out[vs[a]] >> vs[b] & 1:
                loc[a] |= 1 << b
    dp = [[0] * k for _ in range(1 << k)]
    for a in range(k):
        dp[1 << a][a] = 1
    for m in range(1, 1 << k):
        row = dp[m]
        for a in range(k):
            c = row[a]
            if c:
                nx = loc[a] & ~m
                b = 0
                while nx:
                    if nx & 1:
                        dp[m | (1 << b)][b] += c
                    nx >>= 1
                    b += 1
    return sum(dp[(1 << k) - 1])


def lam(out):
    n = len(out)
    full = (1 << n) - 1
    best = None
    for S in range(1, full):
        c = 0
        for v in range(n):
            if S >> v & 1:
                c += bin(out[v] & ~S & full).count('1')
        if best is None or c < best:
            best = c
    return best


def prime(out):
    n = len(out)
    full = (1 << n) - 1
    for u in range(n):
        for w in range(u + 1, n):
            M = (1 << u) | (1 << w)
            changed = True
            while changed:
                changed = False
                for x in range(n):
                    if M >> x & 1:
                        continue
                    o = out[x] & M
                    if o != 0 and o != M:
                        M |= 1 << x
                        changed = True
            if M != full:
                return False
    return True


def main():
    S10 = census_values(10)
    S11 = census_values(11)
    print('# |S(10)| =', len(S10), ' |S(11)| =', len(S11))
    ok_all = True
    seen = 0
    for line in open(sys.argv[1]):
        n, code = line.split()[:2]
        n = int(n)
        out = adj_from(code, n)
        full = (1 << n) - 1
        H = ham(out, full)
        dels = [ham(out, full & ~(1 << v)) for v in range(n)]
        L = lam(out)
        pr = prime(out)
        c1 = H in S11
        c2 = all(d in S10 for d in dels)
        ok = c1 and c2 and L >= 3 and pr
        ok_all &= ok
        seen += 1
        print('H=%d inS11=%s dels_all_inS10=%s lambda=%d prime=%s ESCAPER=%s' % (H, c1, c2, L, pr, ok))
    print('# checked', seen, 'all escapers:', ok_all)


if __name__ == '__main__':
    main()
