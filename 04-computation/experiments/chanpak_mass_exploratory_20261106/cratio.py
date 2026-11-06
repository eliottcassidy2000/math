import itertools
from fractions import Fraction
def paths_by_subset(n, out):
    end = [[0]*n for _ in range(1<<n)]
    for v in range(n): end[1<<v][v] = 1
    for S in range(1, 1<<n):
        for v in range(n):
            c = end[S][v]
            if not c: continue
            nb = out[v] & ~S
            while nb:
                w = (nb & -nb).bit_length()-1; nb &= nb-1
                end[S|(1<<w)][w] += c
    return [sum(end[S]) if S else 1 for S in range(1<<n)]
ratios = {}
for n in range(1, 7):
    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()
    for bits in range(1 << len(arcs)):
        out = [0]*n
        for k,(i,j) in enumerate(arcs):
            if bits>>k & 1: out[i] |= 1<<j
            else: out[j] |= 1<<i
        hs = paths_by_subset(n, out)
        full = (1<<n)-1
        H = hs[full]
        C = sum(hs[A]*hs[full & ~A] for A in range(1<<n))
        r = Fraction(C, H)
        ratios.setdefault(n, set()).add((H, C, r))
for n in sorted(ratios):
    rs = sorted(set(x[2] for x in ratios[n]))
    print(n, 'distinct C/H ratios:', [str(r) for r in rs][:30], ('...' if len(rs) > 30 else ''))
