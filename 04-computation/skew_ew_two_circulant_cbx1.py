#!/usr/bin/env python3
"""THM-476 (claudebox-2026-06-11-S1): the skew Ehlich-Wojtas square law (2n-3 must be a
square) + explicit witnesses at n = 14, 26, 62 + the n=10 golden-ratio best-known (HYP-2405)
+ the multiplier-symmetric two-circulant search (negative at the open order n=86).
Pure python, exact (Bareiss / integer PAF)."""

def det_bareiss(mat):
    a = [row[:] for row in mat]; n = len(a); sign = 1; prev = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            for i in range(k + 1, n):
                if a[i][k] != 0: a[k], a[i] = a[i], a[k]; sign = -sign; break
            else: return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * a[k][k] - a[i][k] * a[k][j]) // prev
        prev = a[k][k]
    return sign * a[n - 1][n - 1]

def EW(n): return 2 * (n - 1) * (n - 2) ** ((n - 2) // 2)

def two_circulant_M(ap, bp):
    """M = [[A, B], [-B^T, A^T]] from circulant first rows (strings over +-)."""
    m = len(ap); n = 2 * m
    a = [1 if c == '+' else -1 for c in ap]; b = [1 if c == '+' else -1 for c in bp]
    assert a[0] == 1 and all(a[k] == -a[(m - k) % m] for k in range(1, m)), "A not skew-type"
    M = [[0] * n for _ in range(n)]
    for i in range(m):
        for j in range(m):
            M[i][j] = a[(j - i) % m]; M[i][j + m] = b[(j - i) % m]
            M[i + m][j] = -b[(i - j) % m]; M[i + m][j + m] = a[(i - j) % m]
    return M

def S_from_rows(rows):
    return [[{'+': 1, '-': -1, '0': 0}[c] for c in r] for r in rows]

def M_of_S(S):
    return [[(1 if i == j else 0) + S[i][j] for j in range(len(S))] for i in range(len(S))]

def charpoly(A):
    """Faddeev-LeVerrier, exact integers."""
    nn = len(A); M = None; coeffs = [1]
    for k in range(1, nn + 1):
        if k == 1: M = [row[:] for row in A]
        else:
            T = [[M[i][j] + (coeffs[-1] if i == j else 0) for j in range(nn)] for i in range(nn)]
            M = [[sum(A[i][x] * T[x][j] for x in range(nn)) for j in range(nn)] for i in range(nn)]
        coeffs.append(-sum(M[i][i] for i in range(nn)) // k)
    return coeffs

# ---- stored witnesses ----
W14 = ['0+--+-+++++---','-0-+--+-+-++++','++0+-+++-++-++','+--0---++---+-','-+++0++++-----',
       '++-+-0--++-+-+','---+-+0--++---','-+---++0-+-++-','--+---++0+---+','-+-++----0--++',
       '---+++-+++0+++','+-+++-+-++-0+-','+---+++-+---0+','+--++-++---+-0']   # anneal witness
A26, B26 = '++-+---+++-+-', '++--+--------'      # float+PAF search witness (m=13, H of order 3)
A62 = '+-+-+-++--+-----+++++-++--+-+-+'           # m=31, H of order 3
B62 = '+++--+--+++---+----+-----+-----'
W10 = ['0+----+++-','-0+++-++++','+-0-+---+-','+-+0+++++-','+---0-++-+',
       '+++-+0-+-+','--+--+0-++','--+---+0--','----++-+0-','+-++---++0']   # HYP-2405 best known

def subgroup_orbits(m, h):
    def order(g):
        x, k = g % m, 1
        while x != 1: x = x * g % m; k += 1
        return k
    g = next(x for x in range(2, m) if order(x) == m - 1)
    H = sorted({pow(g, ((m - 1) // h) * i, m) for i in range(h)})
    seen, orbs = set(), []
    for k in range(1, m):
        if k in seen: continue
        o = sorted({k * x % m for x in H})
        orbs.append(o); seen.update(o)
    return orbs

def paf(row, m):
    return tuple(sum(row[k] * row[(k + j) % m] for k in range(m)) for j in range(1, (m - 1) // 2 + 1))

def search(m, ha, hb, verbose=True):
    """Multiplier-symmetric two-circulant skew-EW search at n = 2m (m prime).
    EW <=> PAF_a'(j) + PAF_b(j) = 2 for all j != 0 (integer conditions, no floats)."""
    orbA = subgroup_orbits(m, ha)
    apairs, used = [], set()
    for i, o in enumerate(orbA):
        if i in used: continue
        neg = sorted({(-x) % m for x in o})
        hit = next((j for j, o2 in enumerate(orbA) if j != i and j not in used and o2 == neg), None)
        if hit is None: return None  # -1 in H
        apairs.append((o, orbA[hit])); used.update({i, hit})
    orbB = subgroup_orbits(m, hb)
    amap = {}
    for mask in range(1 << len(apairs)):
        row = [0] * m; row[0] = 1
        for i, (o, oneg) in enumerate(apairs):
            s = 1 if (mask >> i) & 1 else -1
            for k in o: row[k] = s
            for k in oneg: row[k] = -s
        amap.setdefault(paf(row, m), []).append(tuple(row))
    for b0 in (1, -1):
        for mask in range(1 << len(orbB)):
            row = [b0] * m
            for i, o in enumerate(orbB):
                s = 1 if (mask >> i) & 1 else -1
                for k in o: row[k] = s
            need = tuple(2 - v for v in paf(row, m))
            if need in amap: return (amap[need][0], tuple(row))
    return None

if __name__ == '__main__':
    import math
    print("== THM-476 witnesses (exact Bareiss) ==")
    for n, M in [(14, M_of_S(S_from_rows(W14))),
                 (26, two_circulant_M(A26, B26)),
                 (62, two_circulant_M(A62, B62))]:
        k2 = 2 * n - 3; k = math.isqrt(k2)
        d = det_bareiss(M)
        print(f"n={n}: 2n-3 = {k2} = {k}^2; det = {d}; EW = {EW(n)}; ATTAINED: {d == EW(n)}")
        assert d == EW(n) and k * k == k2
    print("\n== HYP-2405: n=10 best known (EW impossible: 17 not a square) ==")
    M10 = M_of_S(S_from_rows(W10))
    d10 = det_bareiss(M10)
    S10 = S_from_rows(W10)
    SST = [[sum(S10[i][k] * S10[j][k] for k in range(10)) for j in range(10)] for i in range(10)]
    cp = charpoly(SST)
    # (x^2-18x+61)^4 (x-9)^2 expanded check
    import functools
    def polymul(p, q):
        r = [0] * (len(p) + len(q) - 1)
        for i, pi in enumerate(p):
            for j, qj in enumerate(q): r[i + j] += pi * qj
        return r
    target = functools.reduce(polymul, [[1, -18, 61]] * 4 + [[1, -9]] * 2)
    print(f"det = {d10} = 2^9*125; EW(10) = {EW(10)} (unattainable by THM-476)")
    print(f"SS^T char poly == (x^2-18x+61)^4 (x-9)^2: {cp == target}  (eigenvalues 9 +- 2*sqrt(5) in Q(sqrt 5))")
    assert d10 == 64000 and cp == target
    print("\n== search reproduction (multiplier-symmetric two-circulant) ==")
    for m, hs in [(13, [3]), (31, [3]), (43, [21, 7, 3])]:
        got = False
        for ha in hs:
            for hb in hs:
                r = search(m, ha, hb)
                if r:
                    print(f"m={m} (n={2*m}) H_a={ha} H_b={hb}: FOUND")
                    got = True; break
            if got: break
        if not got:
            print(f"m={m} (n={2*m}): NONE under multiplier symmetry orders {hs} (negative result; first open order n=86 resists this ansatz)")
    print("\nall checks passed")
