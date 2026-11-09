"""AUDIT G1: misc exact checks with own DP (aud_core props).
 - THM-329's 13-vertex tournament (adjacency bitmask list) claimed H = 3,719,831 (> circulant 3,711,175);
 - the Z_13 circulant S = {1..6} (claimed 3,711,175);
 - band_4(14), band_4(15) (claimed 4,013,709 and 15,927,743);
 - C3[R5, R5, TT_2] (claimed H = 229,725, lambda = 4) and S17's C3[R5, R5, 1] (11,475)."""
import os, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
core = os.path.join(HERE, 'aud_core')
def props(n, out):
    code = ''.join('1' if out[i] >> j & 1 else '0' for i in range(n) for j in range(i + 1, n))
    # sanity: tournament
    for i in range(n):
        assert not out[i] >> i & 1
        for j in range(i + 1, n):
            assert (out[i] >> j & 1) ^ (out[j] >> i & 1), (i, j)
    line = subprocess.run([core, 'props'], input='%d %s\n' % (n, code), capture_output=True, text=True, check=True).stdout
    a = line.split('|')[0].split()
    return int(a[1]), a[2:], line.split('|')[2].split()[1:]
# THM-329 n=13 best tournament: adj = bitmask of out-neighbours (bit j = 1 iff i -> j)
adj = [3666, 5352, 6411, 3009, 686, 3213, 6708, 5701, 243, 3366, 4444, 5522, 825]
print('THM-329 13-vertex:', props(13, adj)[:2])
circ = [sum(1 << ((i + d) % 13) for d in range(1, 7)) for i in range(13)]
print('circulant Z13 {1..6}:', props(13, circ)[:2])
def band(s, n):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if j - i <= s: out[j] |= 1 << i
            else: out[i] |= 1 << j
    return out
for n in (14, 15):
    print('band_4(%d):' % n, props(n, band(4, n))[:2])
def subst(Q, parts):
    """Q: quotient out-masks on k vertices; parts: list of (size, out-masks) modules; returns substituted tournament"""
    offs = []; o = 0
    for sz, _ in parts: offs.append(o); o += sz
    n = o; out = [0] * n
    for a, (sa, oa) in enumerate(parts):
        for x in range(sa):
            v = offs[a] + x
            out[v] |= sum(1 << (offs[a] + y) for y in range(sa) if oa[x] >> y & 1)
            for b, (sb, ob) in enumerate(parts):
                if b != a and Q[a] >> b & 1:
                    out[v] |= sum(1 << (offs[b] + y) for y in range(sb))
    return n, out
C3 = [0b010, 0b100, 0b001]          # 0->1->2->0
R5 = [sum(1 << ((i + d) % 5) for d in (1, 2)) for i in range(5)]
TT2 = [0b10, 0b00]
ONE = [0]
for name, parts in (('C3[R5,R5,TT2]', [(5, R5), (5, R5), (2, TT2)]), ('C3[R5,R5,1]', [(5, R5), (5, R5), (1, ONE)]),
                    ('C3[R5,TT2,R5]', [(5, R5), (2, TT2), (5, R5)])):
    n, out = subst(C3, parts)
    H, kv, sc = props(n, out)
    print(name, 'n=%d H=%d' % (n, H), kv, 'scores', sorted(map(int, sc)))
