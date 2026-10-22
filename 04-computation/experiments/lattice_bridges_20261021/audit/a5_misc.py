#!/usr/bin/env python3
"""AUDIT of the results note lattices_codes_tournaments_20261021.md, sections 2-10: independent numerical checks.
 (1) binary QR codes spanned by the out-neighbourhoods of P_q (q = 7, 23, 31, 47, 71-dim only); parity extension;
     weight enumerators (meet-in-the-middle over 2^((q+1)/2) codewords); doubly even / self-dual; GF(2) kernels of A(P_7),
     A(P_23); perfectness of [7,4,3] and [23,12,7].
 (2) Construction A theta series from the weight enumerators vs E_4, E_4^3 - 672 Delta; Delta coefficients; Leech theta;
     THM-489's W_g24 = W_e8^3 - 42 P_24 and P_24 -> 16 Delta.
 (3) 196560 decomposition; E_8 = A_8 + glue [3], [6]: explicit norm-2 count and unimodularity.
 (4) score vectors (DP over arcs) vs labelled forests (EGF); c_3 formula; cannonball.
 (5) H(P_7 with one cyclic triangle reversed); cube pattern H(P_q)/|Aut P_q| - 1; |M_12|; bP_(4m).
 (6) w(2^k -+ 1) for k <= 40 under V = {odd >= 3} minus {7, 21}; which k have every prime factor below 80,265.
Output: a5_misc.out
"""
import os, sys, math, itertools, random
from fractions import Fraction
from math import comb, isqrt
import numpy as np
import sympy
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'a5_misc.out'), 'w')
def say(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

# ---------------- (1) codes ----------------
def gf2_rank_basis(rows):
    """rows: python ints; returns a reduced basis (list of ints)."""
    basis = []
    for r in rows:
        for b in basis:
            r = min(r, r ^ b)
        if r: basis.append(r)
    return basis
def gf2_kernel(rows, ncols):
    """kernel of the matrix with given rows (ints over ncols bits): returns basis of {x : r.x = 0 for all r}."""
    # Gaussian elimination on the matrix, then back-substitute
    M = [r for r in rows]; piv = []; R = []
    for c in range(ncols):
        p = next((i for i, r in enumerate(M) if (r >> c) & 1), None)
        if p is None: continue
        r = M.pop(p)
        M = [x ^ r if (x >> c) & 1 else x for x in M]
        R = [x ^ r if (x >> c) & 1 else x for x in R]
        R.append(r); piv.append(c)
    free = [c for c in range(ncols) if c not in piv]
    ker = []
    for f in free:
        x = 1 << f
        for r, c in zip(R, piv):
            if (r >> f) & 1: x |= 1 << c
        ker.append(x)
    return ker
POP16 = np.array([bin(i).count('1') for i in range(1 << 16)], dtype=np.int64)
def popcount64(a):
    a = a.astype(np.uint64)
    return (POP16[(a & np.uint64(0xFFFF)).astype(np.int64)] + POP16[((a >> np.uint64(16)) & np.uint64(0xFFFF)).astype(np.int64)] +
            POP16[((a >> np.uint64(32)) & np.uint64(0xFFFF)).astype(np.int64)] + POP16[((a >> np.uint64(48)) & np.uint64(0xFFFF)).astype(np.int64)])
def span_all(basis):
    vals = np.zeros(1 << len(basis), dtype=np.uint64)
    for i, b in enumerate(basis):
        vals[1 << i:1 << (i + 1)] = vals[:1 << i] ^ np.uint64(b)
    return vals
def weight_enumerator(basis, n):
    k = len(basis); h = k // 2
    A = span_all(basis[:h]); B = span_all(basis[h:])
    W = np.zeros(n + 1, dtype=np.int64)
    for a in A:
        W += np.bincount(popcount64(B ^ a), minlength=n + 1)
    return W
def qr(q): return sorted({x * x % q for x in range(1, q)})
for q in (7, 23, 31, 47):
    Q = qr(q)
    rows = [sum(1 << ((i + r) % q) for r in Q) for i in range(q)]        # out-neighbourhood of i
    basis = gf2_rank_basis(rows)
    ext = [b | ((bin(b).count('1') & 1) << q) for b in basis]           # parity extension
    W = weight_enumerator(ext, q + 1)
    nz = {w: int(c) for w, c in enumerate(W) if c}
    de = all(w % 4 == 0 for w in nz)
    # self-duality: k = n/2 and all inner products of basis vectors even
    sd = (2 * len(ext) == q + 1) and all(bin(a & b).count('1') % 2 == 0 for a in ext for b in ext)
    dmin = min(w for w in nz if w > 0)
    say(f"q={q}: rank of A(P_q) over GF(2) = {len(basis)} (= (q+1)/2: {len(basis) == (q + 1)//2}); extended code [{q+1},{len(ext)},{dmin}], "
        f"doubly even {de}, self-dual {sd}; weight enumerator {nz}")
    if q in (7, 23):
        Wq = weight_enumerator(basis, q); nzq = {w: int(c) for w, c in enumerate(Wq) if c}
        t = (min(w for w in nzq if w > 0) - 1) // 2
        perfect = (2 ** (q - len(basis))) == sum(comb(q, i) for i in range(t + 1))
        ker = gf2_kernel(rows, q); Wk = weight_enumerator(ker, q) if ker else None
        nzk = {w: int(c) for w, c in enumerate(Wk) if c}
        # kernel = even-weight subcode of the row space?
        rowspace = set(int(x) for x in span_all(basis)); evensub = {x for x in rowspace if bin(x).count('1') % 2 == 0}
        kerspace = set(int(x) for x in span_all(ker))
        say(f"   row space [{q},{len(basis)},{min(w for w in nzq if w>0)}] perfect (t={t}): {perfect}; weights {nzq}")
        say(f"   kernel of A(P_{q}) over GF(2): dim {len(ker)}, weights {nzk}; kernel == even-weight subcode of row space: {kerspace == evensub}")
say("Hamming [7,4,3] and Golay [23,12,7] sphere-packing equalities:", 2 ** 3 == 1 + 7, 2 ** 11 == 1 + 23 + 253 + 1771,
    "; extended [24,12,8] perfect?", 2 ** 12 == sum(comb(24, i) for i in range(4)))
# the P_71 gauge: only the dimension (d = 12 is quoted as known)
q = 71; Q = qr(q); rows = [sum(1 << ((i + r) % q) for r in Q) for i in range(q)]
say(f"q=71: rank of A(P_71) over GF(2) = {len(gf2_rank_basis(rows))} (extended QR code [72,36]; d = 12 is classical, not recomputed)")

# ---------------- (2) modular forms ----------------
N = 8
def sigma(k, n): return sum(d ** k for d in range(1, n + 1) if n % d == 0)
E4 = [1] + [240 * sigma(3, n) for n in range(1, N)]
E6 = [1] + [-504 * sigma(5, n) for n in range(1, N)]
def mul(a, b):
    c = [0] * N
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if i + j < N: c[i + j] += x * y
    return c
E43 = mul(mul(E4, E4), E4); E62 = mul(E6, E6)
Delta = [Fraction(x - y, 1728) for x, y in zip(E43, E62)]
say("E_4 =", E4[:5], "; Delta =", [int(x) for x in Delta[:5]])
say("E_8 = E_4^2:", mul(E4, E4)[:4] == [1] + [480 * sigma(7, n) for n in range(1, 4)],
    "; E_10 = E_4 E_6:", mul(E4, E6)[:4] == [1] + [-264 * sigma(9, n) for n in range(1, 4)],
    "; E_14 = E_4^2 E_6:", mul(mul(E4, E4), E6)[:4] == [1] + [-24 * sigma(13, n) for n in range(1, 4)])
say("E_4^3 - 672 Delta =", [int(a - 672 * b) for a, b in zip(E43, Delta)][:4], "; E_4^3 - 720 Delta =", [int(a - 720 * b) for a, b in zip(E43, Delta)][:4])
# Construction A theta from weight enumerators, exponents in units of q^(1/4)
M4 = 4 * N
def theta3():   # sum_z q^(z^2)  -> exponent 4 z^2
    t = [0] * M4
    for z in range(-N, N + 1):
        e = 4 * z * z
        if e < M4: t[e] += 1
    return t
def theta2():   # sum_z q^((z+1/2)^2) -> exponent (2z+1)^2
    t = [0] * M4
    for z in range(-N, N + 1):
        e = (2 * z + 1) ** 2
        if e < M4: t[e] += 1
    return t
def pmul(a, b):
    c = [0] * M4
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y and i + j < M4: c[i + j] += x * y
    return c
def ppow(a, k):
    r = [1] + [0] * (M4 - 1)
    for _ in range(k): r = pmul(r, a)
    return r
def theta_from_WE(we, n):
    t3, t2 = theta3(), theta2(); tot = [0] * M4
    for w, c in we.items():
        term = pmul(ppow(t3, n - w), ppow(t2, w))
        tot = [x + c * y for x, y in zip(tot, term)]
    return tot
th8 = theta_from_WE({0: 1, 4: 14, 8: 1}, 8)
th24 = theta_from_WE({0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, 24)
say("Construction A of [8,4,4]: theta (integer powers) =", th8[0:M4:4][:5], "; only integer powers:", all(th8[i] == 0 for i in range(M4) if i % 4))
say("Construction A of [24,12,8]: theta =", th24[0:M4:4][:4], "; equals E_4^3 - 672 Delta:",
    th24[0:M4:4][:N] == [int(a - 672 * b) for a, b in zip(E43, Delta)][:N])
# THM-489: W_g24 = W_e8^3 - 42 P_24, P_24 = x^4 y^4 (x^4 - y^4)^4; and P_24 -> 16 Delta under Construction A
x, y = sympy.symbols('x y')
We8 = x**8 + 14*x**4*y**4 + y**8; P24 = x**4*y**4*(x**4 - y**4)**4
Wg24 = x**24 + 759*x**16*y**8 + 2576*x**12*y**12 + 759*x**8*y**16 + y**24
say("THM-489 identity W_g24 = W_e8^3 - 42 P_24:", sympy.expand(We8**3 - 42*P24 - Wg24) == 0)
P24poly = sympy.Poly(sympy.expand(P24), x, y)
weP24 = {}
for (ex, ey), c in P24poly.terms(): weP24[ey] = int(c)
t3, t2 = theta3(), theta2(); tot = [0] * M4
for w, c in weP24.items():
    term = pmul(ppow(t3, 24 - w), ppow(t2, w)); tot = [a + c * b for a, b in zip(tot, term)]
say("P_24 under Construction A =", tot[0:M4:4][:5], "; 16 Delta =", [int(16 * d) for d in Delta[:5]])
say("Leech minimal vectors 4*C(24,2) + 759*2^7 + 24*2^12 =", 4 * comb(24, 2) + 759 * 2 ** 7 + 24 * 2 ** 12)
jm720 = None
# j - 720 = theta_Leech / Delta : first terms q^-1 + 24 + ...
thL = [int(a - 720 * b) for a, b in zip(E43, Delta)]
# divide power series thL by Delta = q * (1 - 24 q + ...): thL/Delta = q^-1 * thL / (Delta/q)
D1 = [Delta[i + 1] for i in range(N - 1)]
inv = [Fraction(0)] * (N - 1); inv[0] = 1 / D1[0]
for k in range(1, N - 1):
    inv[k] = -sum(D1[i] * inv[k - i] for i in range(1, k + 1)) / D1[0]
quo = [sum(Fraction(thL[i]) * inv[k - i] for i in range(k + 1)) for k in range(N - 1)]
say("theta_Leech / Delta = q^-1 *", [int(c) for c in quo[:4]], " (j - 720 = q^-1 + 24 + 196884 q + ...)")

# E_8 via A_8 + glue
import numpy.linalg as la
glue3 = np.array([3] * 6 + [-6] * 3, dtype=float) / 9
basisA8 = [np.eye(9)[i] - np.eye(9)[i + 1] for i in range(8)]
# generate E_8 = A_8 + Z*glue3 lattice vectors of norm 2 by brute force on small coefficient ranges: count norm-2 vectors as
# (i) e_i - e_j, (ii) all coordinate permutations of glue classes [3], [6] of norm 2
roots_A8 = {tuple(np.eye(9)[i] - np.eye(9)[j]) for i in range(9) for j in range(9) if i != j}
cls3 = set(); cls6 = set()
for T in itertools.combinations(range(9), 3):
    v = np.full(9, 1 / 3); v[list(T)] = -2 / 3; cls3.add(tuple(np.round(v, 9)))      # (1/3)sum e - (e_i+e_j+e_k)
    w = np.full(9, -1 / 3); w[list(T)] = 2 / 3; cls6.add(tuple(np.round(w, 9)))       # e_i+e_j+e_k - (1/3)sum e
norm = lambda v: float(np.dot(v, v))
ok = all(abs(norm(np.array(v)) - 2) < 1e-9 for v in list(cls3) + list(cls6))
# check membership: v - glue3 in A_8 (integer coordinates summing to 0)?  class [3] = A_8 + glue3, glue3 = (1/3)(1,..,1,-2,-2,-2)
def in_coset(v, g):
    d = np.array(v) - g; return np.allclose(d, np.round(d), atol=1e-6) and abs(d.sum()) < 1e-6
say(f"E_8 via A_8: |roots A_8| = {len(roots_A8)}; class [3] norm-2 vectors = {len(cls3)} (all in A_8+[3]: {all(in_coset(v, glue3) for v in cls3)}); "
    f"class [6] = {len(cls6)} (in A_8+[6]: {all(in_coset(v, 2*glue3) for v in cls6)}); norms 2: {ok}; total {len(roots_A8)+len(cls3)+len(cls6)}")
# completeness: shortest vectors of A_8 + [3] have norm 2 and are exactly cls3: enumerate v = glue3 + sum c_i e_i with sum c = 0, |c_i| <= 2
cnt3 = 0
for c in itertools.product(range(-2, 3), repeat=9):
    if sum(c) != 0: continue
    v = glue3 + np.array(c)
    if abs(norm(v) - 2) < 1e-9: cnt3 += 1
say(f"   brute-force norm-2 vectors in A_8+[3] with |c_i| <= 2: {cnt3}")
# determinant of L = A_8 + Z*glue3: integer HNF of the generators scaled by 3, then Gram determinant / 3^16
gens = [[int(round(3 * x)) for x in v] for v in basisA8 + [glue3]]
Mh = sympy.Matrix(gens)
from sympy.matrices.normalforms import hermite_normal_form
Hn = hermite_normal_form(Mh.T).T          # rows span the same Z-module
rowsH = [Hn.row(i) for i in range(Hn.rows) if any(Hn.row(i))]
B = sympy.Matrix.vstack(*rowsH)
detL = (B * B.T).det() / sympy.Integer(3) ** (2 * B.rows)
say(f"   rank of A_8 + Z[3] = {B.rows}; determinant of its Gram matrix = {detL} (E_8 is unimodular; det A_8 = {sympy.Matrix([[np.dot(a,b) for b in basisA8] for a in basisA8]).det()})")
say("248 = 80 + 84 + 84:", 80 + 84 + 84 == 248, "; dim sl_9 =", 81 - 1, "; C(9,3) =", comb(9, 3), "; 240 = 72 + 168:", 72 + 2 * comb(9, 3) == 240)

# ---------------- (4) scores ----------------
def score_vectors(n):
    S = {tuple([0] * n)}
    for i in range(n):
        for j in range(i + 1, n):
            new = set()
            for s in S:
                a = list(s); a[i] += 1; new.add(tuple(a))
                b = list(s); b[j] += 1; new.add(tuple(b))
            S = new
    return S
def forests(n):
    # labelled forests via f(n) = sum_k C(n-1, k-1) k^(k-2) f(n-k)  (tree containing vertex 1 has size k)
    f = [1]
    for m in range(1, n + 1):
        f.append(sum(comb(m - 1, k - 1) * (k ** (k - 2) if k >= 2 else 1) * f[m - k] for k in range(1, m + 1)))
    return f
fo = forests(8)
sv = [len(score_vectors(n)) for n in range(1, 9)]
say("labelled score vectors n=1..8:", sv, "; labelled forests:", fo[1:9], "; equal:", sv == fo[1:9])
# c_3 formula
def c3_direct(n, adj):
    return sum(1 for a, b, c in itertools.combinations(range(n), 3)
               if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[b][a] and adj[c][b] and adj[a][c]))
okc3 = True
for n in range(3, 31):
    for _ in range(20):
        adj = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                if random.random() < 0.5: adj[i][j] = 1
                else: adj[j][i] = 1
        s = [sum(r) for r in adj]; sb = Fraction(n - 1, 2)
        if Fraction(n ** 3 - n, 24) - sum((Fraction(x) - sb) ** 2 for x in s) / 2 != c3_direct(n, adj): okc3 = False
say("c_3 = (n^3-n)/24 - |s - sbar|^2/2 on 560 random tournaments n = 3..30:", okc3)
nn = sympy.symbols('n')
say("c_3 identity symbolic: C(n,3) - (sum s^2 - C(n,2))/2 with sum s^2 = |s-sbar|^2 + n(n-1)^2/4 gives",
    sympy.factor(sympy.binomial(nn, 3).expand(func=True) + nn * (nn - 1) / 4 - nn * (nn - 1) ** 2 / 8))
cb = [N_ for N_ in range(1, 10 ** 6) if isqrt(N_ * (N_ + 1) * (2 * N_ + 1) // 6) ** 2 == N_ * (N_ + 1) * (2 * N_ + 1) // 6]
say("N < 10^6 with 1^2 + ... + N^2 a square:", cb)

# ---------------- (5) small tournament facts ----------------
def ham_paths(n, adj):
    f = {}
    for v in range(n): f[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = f.get((mask, v));
            if not c: continue
            for w in range(n):
                if not (mask >> w) & 1 and adj[v][w]: f[(mask | (1 << w), w)] = f.get((mask | (1 << w), w), 0) + c
    full = (1 << n) - 1
    return sum(f.get((full, v), 0) for v in range(n))
P7 = [[1 if (j - i) % 7 in (1, 2, 4) else 0 for j in range(7)] for i in range(7)]
cyc = [T for T in itertools.combinations(range(7), 3) if (P7[T[0]][T[1]] and P7[T[1]][T[2]] and P7[T[2]][T[0]]) or (P7[T[1]][T[0]] and P7[T[2]][T[1]] and P7[T[0]][T[2]])]
a, b, c = cyc[0]
Q7 = [r[:] for r in P7]
for u, v in ((a, b), (b, c), (a, c)):
    Q7[u][v], Q7[v][u] = Q7[v][u], Q7[u][v]
say(f"P_7: H = {ham_paths(7, P7)}; reversing the cyclic triangle {cyc[0]}: H = {ham_paths(7, Q7)}, scores unchanged: {[sum(r) for r in Q7] == [sum(r) for r in P7]}")
Hs = {3: 3, 7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349, 27: 418425945688463432031, 31: 19876729650849384218015865}
aut = {q: q * (q - 1) // 2 for q in (3, 7, 11, 19, 23, 31)}; aut[27] = 27 * 13 * 3
for q in sorted(Hs):
    v = Hs[q] // aut[q]; r = round((v - 1) ** (1 / 3)); cube = any((r + d) ** 3 == v - 1 for d in (-1, 0, 1))
    say(f"q={q}: H/|Aut| = {v} (exact: {v * aut[q] == Hs[q]}), minus 1 is a cube: {cube}")
say("|M_12| = 95040 = 1728*55:", 95040 == 1728 * 55, "; 95040 + 55 =", 95040 + 55, "; 1729 = 7*13*19 Korselt:", all((1729 - 1) % (p - 1) == 0 for p in (7, 13, 19)))
def bP(m):
    # Kervaire-Milnor with modern Bernoulli numbers: |bP_4m| = 2^(2m-2) (2^(2m-1) - 1) numer(4 B_2m / m)
    # (gives 28, 992, 8128, 261632, 1448424448; |Theta_19| = 523264 = 2 |bP_20|)
    B = sympy.bernoulli(2 * m)
    return 2 ** (2 * m - 2) * (2 ** (2 * m - 1) - 1) * sympy.numer(sympy.Rational(4) * B / m)
say("|bP_4m| m=2..6:", [int(abs(bP(m))) for m in range(2, 7)], "; divisible by 2^(2m-1)-1:", [int(abs(bP(m))) % (2 ** (2 * m - 1) - 1) == 0 for m in range(2, 7)])

# ---------------- (6) repunit degrees ----------------
from functools import lru_cache
def inV(d): return d >= 3 and d % 2 == 1 and d not in (7, 21)
@lru_cache(maxsize=None)
def w(h):
    if h == 1: return 0
    best = -10 ** 9
    for d in sympy.divisors(h):
        if d > 1 and inV(d):
            r = w(h // d)
            if r >= 0: best = max(best, 1 + r)
    return best if best > -10 ** 8 else -1      # -1: no factorization (h not attained)
rows = []
for k in range(2, 41):
    wm, wp = w(2 ** k - 1), w(2 ** k + 1)
    pf = max(max(sympy.factorint(2 ** k - 1)), max(sympy.factorint(2 ** k + 1)))
    rows.append((k, wm, wp, pf < 80265))
say("k: w(2^k-1), w(2^k+1), all prime factors < 80265 ->", ' '.join(f"{k}:{a},{b}{'' if c else '*'}" for k, a, b, c in rows))
ties_odd = [k for k, a, b, c in rows if k % 2 == 1 and a == b]
mers_larger = [k for k, a, b, c in rows if k % 2 == 1 and a > b]
say("odd k <= 40 with tie:", ties_odd, "; odd k with Mersenne side larger:", mers_larger,
    "; certified k (<= 40):", [k for k, a, b, c in rows if c])
say("rarer (w_- < w_+) for k <= 16:", [k for k, a, b, c in rows if k <= 16 and a < b], "; commoner:", [k for k, a, b, c in rows if k <= 16 and a > b],
    "; ties:", [k for k, a, b, c in rows if k <= 16 and a == b])
out.close()
