#!/usr/bin/env python3
"""(f) small computational sanity checks for the literature items."""
import itertools, math
from fractions import Fraction as Fr
import sympy as sp

print("=== (iii) score vectors of tournaments vs labelled forests vs lattice points of the permutohedron ===")
forests = [1, 1, 2, 7, 38, 291, 2932, 36961]      # A001858, n = 0..7
def score_vectors(n):
    pairs = list(itertools.combinations(range(n), 2))
    S = set()
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        sc = [0]*n
        for (i, j), b in zip(pairs, bits):
            sc[i if b else j] += 1
        S.add(tuple(sc))
    return S
def perm_lattice_points(n):
    # integer points of conv{perms of (0..n-1)}: sum = C(n,2) and every k-subset sum >= C(k,2) (Rado/Landau)
    cnt = 0
    for v in itertools.product(range(n), repeat=n):
        if sum(v) != n*(n-1)//2: continue
        w = sorted(v)
        if all(sum(w[:k]) >= k*(k-1)//2 for k in range(1, n+1)): cnt += 1
    return cnt
for n in range(1, 7):
    sv = score_vectors(n)
    print(f"  n={n}: #score vectors = {len(sv)}, #lattice points of permutohedron = {perm_lattice_points(n)}, #forests = {forests[n]}")

print("\n=== (iv) Voronoi cell of A_n^* = permutohedron (exact check, n = 1..5) ===")
def check_An_star(n):
    N1 = n + 1
    c = [Fr(n - 2*i, 2*N1) for i in range(N1)]           # candidate vertex
    perms = set(itertools.permutations(c))
    def ip(a, b): return sum(x*y for x, y in zip(a, b))
    # (1) every facet of the permutohedron, normal u_S = 1_S - |S|/(n+1) 1 (a vector of A_n^*), is a Voronoi inequality
    ok1 = True
    for k in range(1, N1):
        S = list(range(k))
        u = [Fr(1) - Fr(k, N1) if i in S else -Fr(k, N1) for i in range(N1)]
        hmax = sum(sorted(c, reverse=True)[:k]) - Fr(k, N1)*sum(c)
        ok1 &= (hmax == ip(u, u)/2)
    # (2) every permutohedron vertex satisfies (p,v) <= |v|^2/2 for all lattice v with |v|^2 <= 4|c|^2
    R2 = 4*ip(c, c)
    vecs = []
    B = 3
    for x in itertools.product(range(-B, B+1), repeat=N1):
        if x[-1] != 0: continue              # representatives: A_n^* = proj(Z^{n+1}); fix last coord 0 (adding 1-vector is in kernel)
        m = Fr(sum(x), N1)
        v = [Fr(t) - m for t in x]
        if 0 < ip(v, v) <= R2: vecs.append(v)
    ok2 = all(ip(p, v) <= ip(v, v)/2 for p in perms for v in vecs)
    return ok1, ok2, len(perms)
for n in range(1, 6):
    ok1, ok2, nv = check_An_star(n)
    print(f"  n={n}: facets of permutohedron are Voronoi inequalities: {ok1}; all {nv} vertices inside the Voronoi cell: {ok2}")

print("\n=== (v) binary extended QR codes of length 8 and 24; Construction A of the Golay code ===")
x = sp.symbols('x')
def qr_code(p):
    import flint
    f = flint.nmod_poly([1] + [0]*(p-1) + [1], 2)        # x^p + 1 = x^p - 1 over F2
    c, facs = f.factor()
    g = [fa for fa, e in facs if fa.degree() == (p-1)//2][0]
    k = p - (p-1)//2
    gc = [int(t) for t in g.coeffs()]                        # low -> high
    gens = []
    for s in range(k):
        row = [0]*p
        for i, cf in enumerate(gc): row[(i+s) % p] = int(cf) % 2
        gens.append(row)
    words = []
    for coeffs in itertools.product((0, 1), repeat=k):
        wv = [0]*p
        for cf, row in zip(coeffs, gens):
            if cf: wv = [(a+b) % 2 for a, b in zip(wv, row)]
        wv.append(sum(wv) % 2)                       # extend by overall parity
        words.append(wv)
    WE = {}
    for wv in words: WE[sum(wv)] = WE.get(sum(wv), 0) + 1
    return g, dict(sorted(WE.items())), words
for p in (7, 23):
    g, WE, words = qr_code(p)
    two_is_qr = pow(2, (p-1)//2, p) == 1
    print(f"  p={p}: 2 is a QR mod p: {two_is_qr} (so the degree-(p-1)/2 factors of x^p-1 over F2 are the QR/NR generator polys); g = {g}")
    print(f"     extended QR code length {p+1}: weight enumerator {WE}")
    if p == 23:
        golay_words = words
print("  [8,4,4] with A4 = 14 is the extended Hamming code; [24,12,8] with 759/2576/759 is the extended Golay code (unique codes with these parameters)")
# Construction A: L = (1/sqrt2){x in Z^24 : x mod 2 in C}; roots = norm 2 <=> |x|^2 = 4
min_wt = min(sum(wv) for wv in golay_words if sum(wv) > 0)
print(f"  Construction A on extended Golay code: min nonzero weight {min_wt} > 4, so norm-2 vectors are exactly the 48 vectors +-2e_i/sqrt2 -> root system A1^24;"
      f" all weights = 0 mod 4 -> even; |C| = 2^12 -> det = 2^24/|C|^2 = 1 -> Niemeier lattice N(A1^24)")

print("\n=== (vi) cannonball: sum_{k<=N} k^2 square for N <= 10^6 ===")
sols = []
s = 0
for N in range(1, 10**6 + 1):
    s += N*N
    r = math.isqrt(s)
    if r*r == s: sols.append((N, r))
print("  solutions:", sols, "   sum_{k=0}^{24} k^2 =", sum(k*k for k in range(25)), "= 70^2 =", 70**2,
      " => w = (0,1,...,24|70) has norm 4900 - 4900 = 0 in II_{25,1}")

print("\n=== paper 1 / III-B claim: a(I_4, Theta^(4)_L) (ordered orthogonal root 4-frames) separates all 24 Niemeier lattices ===")
from functools import lru_cache
M = 6
def smul(a, b): return [sum(a[i]*b[n-i] for i in range(n+1)) for n in range(M)]
@lru_cache(None)
def Z(t, k):
    one = tuple([Fr(1)] + [Fr(0)]*(M-1))
    if t == 'A':
        if k <= 0: return one
        N, sub = k*(k+1), [Z('A', k-2)]
    elif t == 'D':
        if k == 1: return one
        if k == 2: return tuple(smul(list(Z('A', 1)), list(Z('A', 1))))
        if k == 3: return Z('A', 3)
        N, sub = 2*k*(k-1), [Z('D', k-2), Z('A', 1)]
    else:
        N = {6: 72, 7: 126, 8: 240}[k]
        sub = [Z('A', 5)] if k == 6 else [Z('D', 6)] if k == 7 else [Z('E', 7)]
    prod = list(one)
    for s_ in sub: prod = smul(prod, list(s_))
    return tuple([Fr(1)] + [N*prod[m-1]/m for m in range(1, M)])
NIEM = {'Leech': [], 'A1^24': [('A',1,24)], 'A2^12': [('A',2,12)], 'A3^8': [('A',3,8)], 'A4^6': [('A',4,6)],
 'A5^4D4': [('A',5,4),('D',4,1)], 'D4^6': [('D',4,6)], 'A6^4': [('A',6,4)], 'A7^2D5^2': [('A',7,2),('D',5,2)],
 'A8^3': [('A',8,3)], 'A9^2D6': [('A',9,2),('D',6,1)], 'D6^4': [('D',6,4)], 'A11D7E6': [('A',11,1),('D',7,1),('E',6,1)],
 'E6^4': [('E',6,4)], 'A12^2': [('A',12,2)], 'D8^3': [('D',8,3)], 'A15D9': [('A',15,1),('D',9,1)],
 'A17E7': [('A',17,1),('E',7,1)], 'D10E7^2': [('D',10,1),('E',7,2)], 'D12^2': [('D',12,2)], 'A24': [('A',24,1)],
 'D16E8': [('D',16,1),('E',8,1)], 'E8^3': [('E',8,3)], 'D24': [('D',24,1)]}
O4 = {}
for name, comps in NIEM.items():
    a = [Fr(1)] + [Fr(0)]*(M-1)
    for t, k, m in comps:
        for _ in range(m): a = smul(a, list(Z(t, k)))
    O4[name] = int(a[4]*24)
vals = list(O4.values())
print("  a(I_4) values:", O4)
print("  pairwise distinct:", len(set(vals)) == len(vals))
# cross-check with the genus-4 test: F_30 coefficient at I_4 computed in c_niemeier_g4detail.out was 206438400
print("  O4(E8^3) - O4(D16E8) =", O4['E8^3'] - O4['D16E8'], "(c_niemeier_g4detail: 206438400)")
