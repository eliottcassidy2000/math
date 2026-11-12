# Verifier (lane architecture): independent recomputation of the base-2 port-chain block matrices of order <= 4
# (THM-5070 (1)) straight from the definition, by fitting the 3x3 matrix on base states and checking it on further
# bases; comparison with the S93 files used by the lane's det_density.py; exact sums of 1/|det|.
# State of a ported tournament (T, u): (h, E, O), h = H(T), E / O = sum of H(T - V(P)) over the directed paths P
# starting at the port with an even / odd number of vertices.  Gluing (B, w, u'): every old vertex beats every block
# vertex, except w -> u; the new port is u'.
# Also checks the identity End_T(u) = O - E (so that (h, e, u) of the lane's P4 is the state in other coordinates).
import itertools, sys
from fractions import Fraction

def tables(n, out):
    full = 1 << n
    inn = [0] * n
    for v in range(n):
        for w in range(n):
            if out[v] >> w & 1:
                inn[w] |= 1 << v
    en = [[0] * n for _ in range(full)]
    st = [[0] * n for _ in range(full)]
    for v in range(n):
        en[1 << v][v] = 1
        st[1 << v][v] = 1
    for mask in range(1, full):
        re = en[mask]; rs = st[mask]
        for v in range(n):
            c = re[v]
            if c:
                nxt = out[v] & ~mask
                while nxt:
                    b = nxt & -nxt
                    en[mask | b][b.bit_length() - 1] += c
                    nxt ^= b
            c = rs[v]
            if c:
                nxt = inn[v] & ~mask
                while nxt:
                    b = nxt & -nxt
                    st[mask | b][b.bit_length() - 1] += c
                    nxt ^= b
    hp = [sum(en[mask]) for mask in range(full)]
    hp[0] = 1
    return st, en, hp

def state(n, out, port):
    st, en, hp = tables(n, out)
    full = (1 << n) - 1
    E = O = 0
    for S in range(1, full + 1):
        if S >> port & 1:
            c = st[S][port] * hp[full ^ S]
            if bin(S).count("1") % 2 == 0: E += c
            else: O += c
    return (hp[full], E, O), en[full][port]

def tour(n, bits):
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits >> k & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
            k += 1
    return out

def glue(nT, outT, u, nB, outB, w, up):
    # old vertices 0..nT-1, block vertices nT..nT+nB-1
    n = nT + nB
    out = outT[:] + [o << nT for o in outB]
    for a in range(nT):
        for b in range(nB):
            if a == u and b == w: out[nT + b] |= 1 << a
            else: out[a] |= 1 << (nT + b)
    return n, out, nT + up

# bases: all ported tournaments of order <= 3 (labelled), plus a few of order 4
bases = []
for n in (1, 2, 3):
    for bits in range(1 << (n * (n - 1) // 2)):
        o = tour(n, bits)
        for u in range(n):
            bases.append((n, o, u))
for bits in (0, 21, 38, 63, 10):
    o = tour(4, bits)
    for u in (0, 2):
        bases.append((4, o, u))
bstates = []
bad_id = 0
for (n, o, u) in bases:
    s, e = state(n, o, u)
    bstates.append(s)
    if e != s[2] - s[1]: bad_id += 1
print("bases:", len(bases), "distinct base states:", sorted(set(bstates)))
print("identity End_T(port) = O - E on the bases: failures", bad_id)

def solve3(S, R):
    # S: three independent states (rows), R: their images; return matrix N (tuple of 9 Fractions, row-major)
    a = [[Fraction(x) for x in row] for row in S]
    det = (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
           + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))
    assert det != 0
    inv = [[0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            m = [[a[r][c] for c in range(3) if c != i] for r in range(3) if r != j]
            inv[i][j] = (-1) ** (i + j) * (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det
    # N * S^T = R^T  ->  N = R^T * (S^T)^-1 ; (S^T)^-1 = (S^-1)^T
    N = [[sum(Fraction(R[k][i]) * inv[j][k] for k in range(3)) for j in range(3)] for i in range(3)]
    # check: N s_k = r_k
    return N

# pick three independent base states
idx = None
for c in itertools.combinations(range(len(bases)), 3):
    S = [bstates[i] for i in c]
    d = (S[0][0] * (S[1][1] * S[2][2] - S[1][2] * S[2][1]) - S[0][1] * (S[1][0] * S[2][2] - S[1][2] * S[2][0])
         + S[0][2] * (S[1][0] * S[2][1] - S[1][1] * S[2][0]))
    if d != 0:
        idx = c; break
print("fitting bases:", idx, [bstates[i] for i in idx])

def det3(m):
    a, b, c, d, e, f, g, h, i = m
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

mats_by_order = {}
lin_fail = 0; nonint = 0; neg = 0; nblocks = 0; id_fail = 0
for k in (1, 2, 3, 4):
    mats = set()
    for bits in range(1 << (k * (k - 1) // 2)):
        oB = tour(k, bits)
        for w in range(k):
            for up in range(k):
                nblocks += 1
                imgs = []
                for (n, o, u) in bases:
                    N_, out_, port_ = glue(n, o, u, k, oB, w, up)
                    s2, e2 = state(N_, out_, port_)
                    if e2 != s2[2] - s2[1]: id_fail += 1
                    imgs.append(s2)
                Nm = solve3([bstates[i] for i in idx], [imgs[i] for i in idx])
                flat = tuple(x for row in Nm for x in row)
                if any(x.denominator != 1 for x in flat): nonint += 1
                if any(x < 0 for x in flat): neg += 1
                for s, r in zip(bstates, imgs):
                    if tuple(sum(Nm[i][j] * s[j] for j in range(3)) for i in range(3)) != tuple(r):
                        lin_fail += 1
                mats.add(tuple(int(x) for x in flat))
    mats_by_order[k] = mats
print("marked blocks of order <= 4:", nblocks, "| linearity failures on all bases:", lin_fail,
      "| non-integer matrices:", nonint, "| matrices with a negative entry:", neg,
      "| End = O - E failures on glued tournaments:", id_fail)
D = "C:/tmp/wt-s93/04-computation/experiments/s93_chanpak_wave_20261111/thinorbit/"
def load(k):
    s = set()
    for line in open(D + "blocks2_k%d.txt" % k):
        t = line.split()
        if len(t) >= 9: s.add(tuple(int(x) for x in t[:9]))
    return s
def transpose(m):
    return (m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8])
for k in (3, 4):
    mine = set()
    for j in range(1, k + 1): mine |= mats_by_order[j]
    theirs = load(k)
    inv = [abs(det3(m)) for m in mine if det3(m) != 0]
    print("order <= %d: mine %d distinct, %d invertible, %d unimodular, sum 1/|det| = %s ; S93 file %d distinct ; "
          "equal as sets (row-major): %s ; equal after transposition: %s"
          % (k, len(mine), len(inv), sum(1 for d in inv if d == 1), sum(Fraction(1, d) for d in inv), len(theirs),
             mine == theirs, set(transpose(m) for m in mine) == theirs))
for k in (5, 6):
    theirs = load(k)
    inv = [abs(det3(m)) for m in theirs if det3(m) != 0]
    print("order <= %d (S93 file only, not recomputed): %d distinct, %d invertible, %d unimodular, sum 1/|det| = %.6f"
          % (k, len(theirs), len(inv), sum(1 for d in inv if d == 1), float(sum(Fraction(1, d) for d in inv))))
tower = (1, 2, 0, 0, 0, 1, 1, 1, 0)
print("one-vertex block matrix:", sorted(mats_by_order[1]), "equals the tower matrix [[1,2,0],[0,0,1],[1,1,0]]:",
      mats_by_order[1] == {tower})
