"""Sandwich modules at a fixed root x: G_u above T-x, G_d below T-x, G_u => G_d, arbitrary arcs to x.
Registers W[r1][r2] over ordered pairs (R1, R2) of disjoint paths covering T-x (R1 = before x, R2 = after x):
r1 in {e, i (last -> x), o (x -> last)},  r2 in {e, i (first -> x), o (x -> first)} (same constants E_, I_, O_).  H = W[ee]+W[eo]+W[ie]+W[io]."""
import itertools, random
E_, I_, O_ = 0, 1, 2
HIDX = (E_ * 3 + E_, E_ * 3 + O_, I_ * 3 + E_, I_ * 3 + O_)   # R1 empty or ends at an in-neighbour of x; R2 empty or starts at an out-neighbour

def path_table(n, out, verts):
    """dict (mask, first, last) -> count of directed paths on exactly that vertex set"""
    tab = {}
    def ext(p, used):
        key = (used, p[0], p[-1]); tab[key] = tab.get(key, 0) + 1
        for w in verts:
            if not used >> w & 1 and out[p[-1]] >> w & 1:
                ext(p + [w], used | 1 << w)
    for v in verts:
        ext([v], 1 << v)
    return tab

def regs(n, out, x=0):
    verts = [v for v in range(n) if v != x]
    full = 0
    for v in verts: full |= 1 << v
    tab = path_table(n, out, verts)
    bymask = {}
    for (m, f, l), c in tab.items():
        bymask.setdefault(m, []).append((f, l, c))
    W = [0] * 9
    def t1(l):  return I_ if out[l] >> x & 1 else O_
    def t2(f):  return O_ if out[x] >> f & 1 else I_
    sub_masks = [0] + list(bymask.keys())
    for m1 in sub_masks:
        if m1 & ~full: continue
        m2 = full & ~m1
        if m2 != 0 and m2 not in bymask: continue
        L1 = [(E_, 1)] if m1 == 0 else [(t1(l), c) for (f, l, c) in bymask[m1]]
        L2 = [(E_, 1)] if m2 == 0 else [(t2(f), c) for (f, l, c) in bymask[m2]]
        for a, ca in L1:
            for b, cb in L2:
                W[a * 3 + b] += ca * cb
    return W

def H_of(W):
    return sum(W[i] for i in HIDX)

def splits(k, out):
    """ordered pairs (P1, P2) of disjoint paths covering all k vertices (either may be empty): list of
    (P1 first, P1 last, P2 first, P2 last) with None for empty"""
    verts = list(range(k)); full = (1 << k) - 1
    tab = path_table(k, out, verts) if k else {}
    res = []
    if k == 0:
        return [(None, None, None, None)]
    bym = {}
    for (m, f, l), c in tab.items():
        bym.setdefault(m, []).append((f, l, c))
    for m1 in [0] + list(bym.keys()):
        m2 = full & ~m1
        if m2 != 0 and m2 not in bym: continue
        A = [(None, None, 1)] if m1 == 0 else bym[m1]
        B = [(None, None, 1)] if m2 == 0 else bym[m2]
        for f1, l1, c1 in A:
            for f2, l2, c2 in B:
                res.extend([(f1, l1, f2, l2)] * (c1 * c2))
    return res

def module_matrix(gu, ou, xu, gd, od, xd):
    """gu = |G_u|, ou = out masks of G_u, xu = bitmask of G_u vertices that x beats; same for G_d.
    returns 9x9 C with W' = C W"""
    C = [[0] * 9 for _ in range(9)]
    SU, SD = splits(gu, ou), splits(gd, od)
    def tend_u(v):  return O_ if xu >> v & 1 else I_      # last vertex v of a run, in G_u
    def tend_d(v):  return O_ if xd >> v & 1 else I_
    def tst_u(v):   return O_ if xu >> v & 1 else I_      # first vertex of R2 in G_u: x -> v gives 'o'
    def tst_d(v):   return O_ if xd >> v & 1 else I_
    for r1 in range(3):
        for r2 in range(3):
            for (a1f, a1l, a2f, a2l) in SU:
                for (d1f, d1l, d2f, d2l) in SD:
                    # R1' = alpha1 R1 delta1 : end type
                    if d1l is not None: n1 = tend_d(d1l)
                    elif r1 != E_: n1 = r1
                    elif a1l is not None: n1 = tend_u(a1l)
                    else: n1 = E_
                    # R2' = alpha2 R2 delta2 : start type
                    if a2f is not None: n2 = tst_u(a2f)
                    elif r2 != E_: n2 = r2
                    elif d2f is not None: n2 = tst_d(d2f)
                    else: n2 = E_
                    C[n1 * 3 + n2][r1 * 3 + r2] += 1
    return C

def apply(C, W):
    return [sum(C[i][j] * W[j] for j in range(9)) for i in range(9)]

def build(n, out, gu, ou, xu, gd, od, xd, x=0):
    """the literal tournament T' (root x stays vertex x)"""
    N = n + gu + gd; o = out[:] + [0] * (gu + gd)
    U = list(range(n, n + gu)); D = list(range(n + gu, N)); rest = [v for v in range(n) if v != x]
    for a in range(gu):
        o[U[a]] |= sum(1 << U[b] for b in range(gu) if ou[a] >> b & 1)
        for v in rest: o[U[a]] |= 1 << v
        for d in D: o[U[a]] |= 1 << d
        if xu >> a & 1: o[x] |= 1 << U[a]
        else: o[U[a]] |= 1 << x
    for a in range(gd):
        o[D[a]] |= sum(1 << D[b] for b in range(gd) if od[a] >> b & 1)
        for v in rest: o[v] |= 1 << D[a]
        if xd >> a & 1: o[x] |= 1 << D[a]
        else: o[D[a]] |= 1 << x
    return N, o

def randtour(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < .5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out

def tours(k):
    pairs = [(i, j) for i in range(k) for j in range(i + 1, k)]
    for bits in range(1 << len(pairs)):
        out = [0] * k
        for t, (i, j) in enumerate(pairs):
            if bits >> t & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out

if __name__ == '__main__':
    from w45_hk import H_of as Hbrute
    rng = random.Random(1); bad = badH = 0
    for trial in range(120):
        n = rng.randint(1, 5); T = randtour(n, rng)
        gu = rng.randint(0, 2); gd = rng.randint(0, 2)
        ou = randtour(gu, rng); od = randtour(gd, rng)
        xu = rng.randrange(1 << gu) if gu else 0; xd = rng.randrange(1 << gd) if gd else 0
        W = regs(n, T); C = module_matrix(gu, ou, xu, gd, od, xd)
        N, o = build(n, T, gu, ou, xu, gd, od, xd)
        W2 = regs(N, o)
        if apply(C, W) != W2: bad += 1
        if H_of(W2) != Hbrute(N, o): badH += 1
    print("sandwich register machine: W' = C W failures", bad, "; H readout failures", badH, "(120 random)")
