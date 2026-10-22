"""Audit S50, second pass - shared helpers for b1/b2 (own code): which (A, B, i, j, orientation) configurations the
written proof of THM-4745 (5) leaves uncovered, and the gluing builder."""
import alib as L


def bounds(parts):
    ls = [0]
    for p in parts:
        ls.append(ls[-1] + p)
    return ls


def pos_info(parts, i):
    ls = bounds(parts)
    if i in ls[1:-1]:
        return ('J',)
    j = max(t for t in range(len(parts)) if ls[t] <= i)
    if ls[j] < i < ls[j + 1] and parts[j] >= 3:
        return ('I', j)
    return ('E', j, ls[j], ls[j + 1])


def uncovered_cases(nmax):
    comps = [c for N in range(2, nmax + 1) for c in L.compositions(N, 2, 4)]
    out = []
    for pa in comps:
        for pb in comps:
            NA, NB = sum(pa), sum(pb)
            for i in range(NA + 1):
                for j in range(NB + 1):
                    ia, ib = pos_info(pa, i), pos_info(pb, j)
                    if ia[0] != 'E' or ib[0] != 'E':
                        continue
                    _, ja, sa, ea = ia
                    _, jb, sb, eb = ib
                    c3a, c3b = pa[ja] == 2, pb[jb] == 2
                    wholeA, wholeB = pa == (2,), pb == (2,)
                    for o in ('fwd', 'bwd'):
                        if o == 'fwd':
                            allowed = (c3a or i == ea) and (c3b or j == sb)
                            covered = (wholeA or i == NA) and (wholeB or j == 0)
                        else:
                            allowed = (c3a or i == sa) and (c3b or j == eb)
                            covered = (wholeA or i == 0) and (wholeB or j == NB)
                        if not allowed or covered:
                            continue
                        types = []
                        if not wholeA:
                            if i not in (0, NA):
                                types.append('MID_A')
                            elif (o == 'fwd' and i != NA) or (o == 'bwd' and i != 0):
                                types.append('END_A')
                        if not wholeB:
                            if j not in (0, NB):
                                types.append('MID_B')
                            elif (o == 'fwd' and j != 0) or (o == 'bwd' and j != NB):
                                types.append('END_B')
                        out.append((pa, pb, i, j, o, tuple(types)))
    return out


def glue(pa, pb, i, j, orient):
    """T on A' = 0..a-2, v = a-1, B' = a..n-1, with A, B arcs and the LOCAL cross arcs oriented by `orient`;
    returns (adjacency, free cross pairs, n, A vertex set, B vertex set, H(A)H(B), info)."""
    A, B = L.chain(pa), L.chain(pb)
    a, b = len(A), len(B)
    n = a + b - 1
    mA = {}; k = 0
    for x in range(a):
        if x == i:
            mA[x] = a - 1
        else:
            mA[x] = k; k += 1
    mB = {}; k = a
    for y in range(b):
        if y == j:
            mB[y] = a - 1
        else:
            mB[y] = k; k += 1
    base = [0] * n
    for x in range(a):
        for x2 in range(a):
            if (A[x] >> x2) & 1:
                base[mA[x]] |= 1 << mA[x2]
    for y in range(b):
        for y2 in range(b):
            if (B[y] >> y2) & 1:
                base[mB[y]] |= 1 << mB[y2]
    ia, ib = pos_info(pa, i), pos_info(pb, j)
    locA = [mA[x] for x in range(ia[2], ia[3] + 1) if x != i]
    locB = [mB[y] for y in range(ib[2], ib[3] + 1) if y != j]
    for x in locA:
        for y in locB:
            if orient == 'fwd':
                base[x] |= 1 << y
            else:
                base[y] |= 1 << x
    free = [(x, y) for x in range(a - 1) for y in range(a, n) if not (x in locA and y in locB)]
    Aset = set(range(a)); Bset = set([a - 1] + list(range(a, n)))
    info = dict(v=a - 1, mA=mA, mB=mB, locA=locA, locB=locB)
    return base, free, n, Aset, Bset, L.H(A) * L.H(B), info
