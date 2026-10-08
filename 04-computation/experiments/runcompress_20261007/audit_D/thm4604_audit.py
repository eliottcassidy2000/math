#!/usr/bin/env python3
"""Audit D of THM-4604 (Collatz words are the reducible locus).  Independent of fricke_check.py.

Standard library only (fractions, decimal, random, itertools).  Exact rational arithmetic everywhere except the
explicit 60-digit Decimal checks of the cosh formula.  Light: well under a minute.

Parts
  A  matrix conventions, projective action, the repo's carry (THM-4600 / THM-4555 / frieze note (1)), actual U-iterates
  B  cocycle identity, positivity / oddness / 3-adic unit, coboundaries = conjugation by translations, normalised cocycle,
     non-splitness (no global coboundary), anti-homomorphism remark
  C  cycle points: B_(w^j) = c_w (2^(jS) - 3^(jm)), G_w (c_w,1) = 2^S (c_w,1), c_w is the 2-adic periodic point with U-word w,
     centraliser in the Borel group (rank computation), anchored states, the unipotent exception in the GROUP
  D  traces: tr = 3^m + 2^S, t = 2cosh(delta/2) (60 digits), traces of group words, t injective on count pairs
  E  Cayley cubic: exact Fricke identity in GL2 form, tr[G_u,G_v] = 2, grid proof of the (a+1/a, b+1/b, ab+1/ab)
     parametrisation, the (O,E) example, Markov level -2 (Cohn matrices, x^2+y^2+z^2 = xyz solutions = 3 x Markov)
  F  mutation (x,y,z) -> (x, z, xz - y) = (u,v) -> (u,uv); Christoffel tree from (O,E) walks Farey mediants (m, L);
     HYP-9230 tuning errors theta_f = delta/ln2 when Sigma w = round(|w| log2 3)
  G  friezes: M(a) = [[a,-1],[1,0]]; Collatz valuation words that ARE quiddities with product -I (frieze note (15));
     G_w never +-I / scalar; ear identity is an M-relation, changes (|w|, Sigma w); minors of (3^i, B_i) = carries
     (cluster/Pluecker relations among carries); faithfulness of w -> G_w (no relations at all); an equal-count merge
     decided purely by carries (trace-invisible)
"""
import random, itertools, math
from fractions import Fraction as Fr
from decimal import Decimal, getcontext

getcontext().prec = 60
rnd = random.Random(20261007)
CHECKS = {}
def ok(name, cond, n=1):
    if not cond: raise AssertionError(name)   # explicit (survives python -O)
    CHECKS[name] = CHECKS.get(name, 0) + n

# ---------- 2x2 exact matrices ----------
def mul(X, Y):
    return ((X[0][0]*Y[0][0] + X[0][1]*Y[1][0], X[0][0]*Y[0][1] + X[0][1]*Y[1][1]),
            (X[1][0]*Y[0][0] + X[1][1]*Y[1][0], X[1][0]*Y[0][1] + X[1][1]*Y[1][1]))
def det(X): return X[0][0]*X[1][1] - X[0][1]*X[1][0]
def tr(X): return X[0][0] + X[1][1]
def inv(X):
    d = Fr(det(X))
    return ((X[1][1]/d, -X[0][1]/d), (-X[1][0]/d, X[0][0]/d))
I2 = ((1, 0), (0, 1))
def mpow(X, k):
    R = I2
    for _ in range(k): R = mul(R, X)
    return R

O = ((3, 1), (0, 2))          # x -> (3x+1)/2
E = ((1, 0), (0, 2))          # x -> x/2
def Gletter(a): return mul(mpow(E, a - 1), O)
def G(w):                      # chronological: G_(uv) = G_v G_u
    R = I2
    for a in w: R = mul(Gletter(a), R)
    return R
def carry_rec(w):              # frieze note (1): B_(i+1) = 3 B_i + Q_i, Q_i = 2^(a_1+...+a_i)
    B, A = 0, 0
    for a in w:
        B, A = 3*B + 2**A, A + a
    return B
def carry_closed(w):           # B_w = sum_i 3^(m-i) 2^(A_(i-1))
    m, A, B = len(w), 0, 0
    for i, a in enumerate(w, start=1):
        B += 3**(m - i) * 2**A
        A += a
    return B
def v2(n):
    n = abs(n); k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k
def Uword(n, r):               # actual U-letters of an odd integer (or 2-adic rational with odd denominator)
    x = Fr(n); out = []
    for _ in range(r):
        if not (x.numerator % 2 == 1 and x.denominator % 2 == 1): raise AssertionError('Uword: not a 2-adic unit')
        y = 3*x + 1
        a = v2(y.numerator)
        out.append(a); x = y / 2**a
    return out, x
def rand_word(maxlen=8, maxa=6): return [rnd.randint(1, maxa) for _ in range(rnd.randint(1, maxlen))]

# ======================= A. conventions =======================
for a in range(1, 13):
    ok("A1 G_a = E^(a-1) O = [[3,1],[0,2^a]]", Gletter(a) == ((3, 1), (0, 2**a)))
for _ in range(3000):
    w = rand_word(10, 7); m, S = len(w), sum(w); M = G(w)
    ok("A2 G_w = [[3^|w|, B_w],[0, 2^Sw]] with B_w = repo carry (recursion = closed form)",
       M == ((3**m, carry_rec(w)), (0, 2**S)) and carry_rec(w) == carry_closed(w))
for _ in range(3000):
    n = 2*rnd.randint(0, 10**12) + 1
    r = rnd.randint(1, 12)
    w, y = Uword(n, r); M = G(w)
    ok("A3 projective action of G_w on (n,1) = actual U^r(n) (F_w(x) = (3^m x + B_w)/2^S)",
       Fr(M[0][0]*n + M[0][1], M[1][1]) == y)
    ok("A4 THM-4555 inverse I_w(z) = (2^A_w z - B_w)/3^|w| recovers the source",
       Fr(2**sum(w) * y - carry_rec(w), 3**len(w)) == n)
# Terras-step (O/E) version of the same action
for _ in range(1000):
    n = rnd.randint(1, 10**9); x = n; R = I2
    for _ in range(rnd.randint(1, 30)):
        if x % 2: x = (3*x + 1)//2; R = mul(O, R)
        else: x //= 2; R = mul(E, R)
    ok("A5 O/E Terras matrices act chronologically (R = T_last ... T_first)", Fr(R[0][0]*n + R[0][1], R[1][1]) == x)

# ======================= B. cocycle =======================
for _ in range(4000):
    u, v = rand_word(), rand_word()
    Bu, Bv, Buv = carry_rec(u), carry_rec(v), carry_rec(u + v)
    ok("B1 cocycle B_uv = 3^|v| B_u + 2^Su B_v", Buv == 3**len(v)*Bu + 2**sum(u)*Bv)
    ok("B2 G_(uv) = G_v G_u (chronological composition; see H1 for non-commutativity)",
       G(u + v) == mul(G(v), G(u)))
    # normalised cocycle z(w) = B_w/3^|w| is an honest left 1-cocycle for concatenation, twist 2^S/3^m
    zu, zv, zuv = Fr(Bu, 3**len(u)), Fr(Bv, 3**len(v)), Fr(Buv, 3**len(u + v))
    ok("B3 z(uv) = z(u) + (2^Su/3^|u|) z(v) for z = B/3^|w| (frieze note's marked coordinate)",
       zuv == zu + Fr(2**sum(u), 3**len(u))*zv)
for _ in range(4000):
    w = rand_word(12, 8); B = carry_rec(w)
    ok("B4 B_w > 0, odd, prime to 3 (nonempty w)", B > 0 and B % 2 == 1 and B % 3 != 0)
# coboundary = conjugation by a translation T_c = [[1,c],[0,1]]
for _ in range(500):
    w = rand_word(); M = G(w); c = Fr(rnd.randint(-50, 50), rnd.randint(1, 30))
    T, Ti = ((1, c), (0, 1)), ((1, -c), (0, 1))
    C1 = mul(Ti, mul(M, T)); C2 = mul(T, mul(M, Ti))
    m, S = len(w), sum(w)
    ok("B5 T_c^-1 G T_c changes B by c(chi1 - chi2); T_c G T_c^-1 by c(chi2 - chi1) (same set of coboundaries)",
       C1[0][1] == M[0][1] + c*(3**m - 2**S) and C2[0][1] == M[0][1] + c*(2**S - 3**m))
# non-split: no single c with B_w = c(2^S - 3^m) for both letters 1 and 2
c1 = Fr(carry_rec([1]), 2**1 - 3); c2 = Fr(carry_rec([2]), 2**2 - 3)
ok("B6 cocycle is not a global coboundary (c_1 = -1 != c_2 = +1): the extension is non-split", (c1, c2) == (-1, 1))

# ======================= C. cycle points, centraliser =======================
def cw(w): return Fr(carry_rec(w), 2**sum(w) - 3**len(w))
for _ in range(1500):
    w = rand_word(7, 6); m, S = len(w), sum(w); c = cw(w); M = G(w)
    ok("C1 2^S != 3^m and c_w != 0 for nonempty w", 2**S != 3**m and c != 0)
    ok("C2 G_w (c_w,1)^T = 2^S (c_w,1)^T", (M[0][0]*c + M[0][1], M[1][1]) == (2**S * c, 2**S))
    for j in range(0, 5):
        ok("C3 B_(w^j) = c_w (2^(jS) - 3^(jm)), j = 0..4", carry_rec(w*j) == c*(2**(j*S) - 3**(j*m)))
    # c_w is a 2-adic unit whose actual (2-adic) U-word is w, periodically
    word3, back = Uword(c, 3*m)
    ok("C4 c_w is a 2-adic unit with U-word w repeated, returning to c_w (rational cycle point)",
       word3 == w*3 and back == c)
    # centraliser in the Borel group: X = [[p,q],[0,r]], XG = GX  <=>  q(2^S - 3^m) = B (r - p)
    p, r = Fr(rnd.randint(1, 40)), Fr(rnd.randint(1, 40))
    X = ((p, c*(r - p)), (0, r))
    ok("C5 [[p, c_w(r-p)],[0,r]] commutes with G_w", mul(X, M) == mul(M, X))
    # it fixes c_w and infinity projectively; it is T_c diag(p,r) T_c^-1
    ok("C6 the centraliser element fixes c_w and infinity", Fr(p*c + c*(r - p), r) == c and X[1][0] == 0)
    # any commuting Borel X has this form: the linear system in (p,q,r) has rank 1 (solution space dim 2)
    # coefficients of XG - GX in the (1,2) slot: p*B + q*2^S - 3^m*q - B*r = 0, other slots vanish identically
    row = (M[0][1], 2**S - 3**m, -M[0][1])
    ok("C7 commuting condition is one nonzero linear equation (centraliser = 2-dim torus)", any(row))
    # anchored states of THM-4600: A(x) = c + 3^k (x - c) commute with F_w; general torus element has any ratio
    for k in (-2, -1, 1, 3):
        Ak = ((Fr(3)**k, c*(1 - Fr(3)**k)), (0, 1))
        ok("C8 anchored states x -> c_w + 3^k (x - c_w) lie in the centraliser torus", mul(Ak, M) == mul(M, Ak))
# the hypothesis 2^S != 3^m only bites in the GROUP: G_(12) G_(21)^-1 is a translation (unipotent)
U12 = mul(G([1, 2]), inv(G([2, 1])))
ok("C9 group element G_(1,2) G_(2,1)^-1 = translation x -> x - 1/4 (unipotent; centraliser not a torus)",
   U12 == ((1, Fr(-1, 4)), (0, 1)))

# ======================= D. traces =======================
def tnorm_dec(m, S):  # (3^m + 2^S)/sqrt(3^m 2^S) in Decimal
    return (Decimal(3)**m + Decimal(2)**S) / (Decimal(3)**m * Decimal(2)**S).sqrt()
def twocosh_half(m, S):
    d = Decimal(m) * Decimal(3).ln() - Decimal(S) * Decimal(2).ln()
    return (d/2).exp() + (-d/2).exp()
for _ in range(1500):
    w = rand_word(12, 8); m, S = len(w), sum(w); M = G(w)
    ok("D1 tr G_w = 3^m + 2^S, det = 3^m 2^S", tr(M) == 3**m + 2**S and det(M) == 3**m * 2**S)
    ok("D2 t(w) = tr/sqrt(det) = 2cosh(delta/2), delta = m ln3 - S ln2 (60 digits)",
       abs(tnorm_dec(m, S) - twocosh_half(m, S)) < Decimal(10)**-45)
# traces of arbitrary GROUP words in the letters (with inverses) depend only on the summed counts
for _ in range(800):
    R = I2; m = 0; S = 0
    for _ in range(rnd.randint(1, 7)):
        w = rand_word(4, 5); e = rnd.choice((1, -1))
        R = mul(G(w) if e == 1 else inv(G(w)), R); m += e*len(w); S += e*sum(w)
    ok("D3 trace of any group word = 3^(sum e|w|) + 2^(sum e Sw) (character of the semisimplification)",
       tr(R) == Fr(3)**m + Fr(2)**S)
# the normalised trace is injective on count pairs of nonempty words (t^2 = (3^m+2^S)^2/(3^m 2^S))
seen = {}
for m in range(1, 61):
    for S in range(m, 2*m + 40):
        t2 = Fr((3**m + 2**S)**2, 3**m * 2**S)
        ok("D4 t(w) determines (|w|, Sw) (no two count pairs share a trace; t > 2)", t2 not in seen and t2 > 4)
        seen[t2] = (m, S)

# ======================= E. Cayley cubic, Fricke =======================
def kappa_gl2(A, B):   # x^2+y^2+z^2-xyz with SL2-normalised traces, exact (sqrt(detA detB) = sqrt(det AB))
    AB = mul(A, B)
    x2 = Fr(tr(A)**2, det(A)); y2 = Fr(tr(B)**2, det(B)); z2 = Fr(tr(AB)**2, det(AB))
    xyz = Fr(tr(A)*tr(B)*tr(AB), det(AB))
    return x2 + y2 + z2 - xyz
def commut_tr(A, B): return tr(mul(mul(A, B), mul(inv(A), inv(B))))
# Fricke identity tr[A,B] = kappa - 2 for random GL2(Q) with positive determinant (normalisation-invariant)
for _ in range(500):
    while True:
        A = tuple(tuple(Fr(rnd.randint(-9, 9)) for _ in range(2)) for _ in range(2))
        B = tuple(tuple(Fr(rnd.randint(-9, 9)) for _ in range(2)) for _ in range(2))
        if det(A) > 0 and det(B) > 0: break
    ok("E1 Fricke identity tr[A,B] = x^2+y^2+z^2-xyz-2 (normalised traces, random GL2+ matrices)",
       commut_tr(A, B) == kappa_gl2(A, B) - 2)
for _ in range(2000):
    u, v = rand_word(), rand_word()
    Gu, Gv = G(u), G(v)
    ok("E2 Collatz pair: kappa = 4 and tr[G_u, G_v] = 2 (commutator unipotent)",
       kappa_gl2(Gu, Gv) == 4 and commut_tr(Gu, Gv) == 2)
    x3 = rand_word()
    ok("E3 products of several words: (t(uv), t(x), t(uvx)) also on the cubic", kappa_gl2(G(u + v), G(x3)) == 4)
# rigorous: a^2 b^2 (x^2+y^2+z^2-xyz-4) is a polynomial of degree <= 4 in a and in b; vanishing on a 6x6 grid proves it
pts = [Fr(k, 7) for k in (1, 2, 3, 5, 11, 13)]
for a in pts:
    for b in pts:
        x, y, z = a + 1/a, b + 1/b, a*b + 1/(a*b)
        ok("E4 (a+1/a, b+1/b, ab+1/(ab)) on the Cayley cubic: 6x6 grid => polynomial identity", x*x + y*y + z*z - x*y*z == 4)
ok("E5 example (O,E): (5/sqrt6, 3/sqrt2, 7/sqrt12), 25/6 + 9/2 + 49/12 - 105/12 = 4",
   Fr(25, 6) + Fr(9, 2) + Fr(49, 12) - Fr(105, 12) == 4 and kappa_gl2(O, E) == 4
   and (tr(O), det(O), tr(E), det(E), tr(mul(O, E)), det(mul(O, E))) == (5, 6, 3, 2, 7, 12))
# Markov at kappa = -2: Cohn matrices and the scaled equation
A_ = ((1, 1), (1, 2)); B_ = ((1, -1), (-1, 2))
ok("E6 Cohn pair: traces (3,3,3), tr[A,B] = -2", (tr(A_), tr(B_), tr(mul(A_, B_))) == (3, 3, 3) and commut_tr(A_, B_) == -2)
sols = [(x, y, z) for x in range(1, 200) for y in range(x, 200) for z in range(y, 200) if x*x + y*y + z*z == x*y*z]
markov = {(a, b, c) for a in range(1, 70) for b in range(a, 70) for c in range(b, 70) if a*a + b*b + c*c == 3*a*b*c}
ok("E7 positive solutions of x^2+y^2+z^2 = xyz (<200) are exactly 3 x Markov triples",
   set(sols) == {(3*a, 3*b, 3*c) for (a, b, c) in markov if 3*c < 200} and len(sols) > 3)

# ======================= F. mutation, Christoffel / Farey =======================
for _ in range(1000):
    u, v = rand_word(5, 5), rand_word(5, 5)
    Gu, Guv = G(u), G(u + v); Gv = G(v)
    # Cayley-Hamilton: tr(A^2 B) = tr(A) tr(AB) - det(A) tr(B); normalised: t(u.uv) = t(u) t(uv) - t(v)
    ok("F1 tr(G_(u u v)) = tr(G_u) tr(G_(uv)) - det(G_u) tr(G_v) (=> t(u.uv) = t(u)t(uv) - t(v))",
       tr(G(u + u + v)) == tr(Gu)*tr(Guv) - det(Gu)*tr(Gv))
    x, y, z = (float(tnorm_dec(len(w), sum(w))) for w in (u, v, u + v))
    znew = float(tnorm_dec(len(u + u + v), sum(u + u + v)))
    ok("F2 mutation (x,y,z) -> (x, z, xz - y) is (u,v) -> (u,uv) (float check)", abs(znew - (x*z - y)) < 1e-9*max(1, znew))
# Christoffel tree from (O,E) as O/E words; counts (m,L): m = #O, L = #letters (Terras steps)
def counts(word): return (word.count('O'), len(word))
def tOE(word):
    m, L = counts(word); return tnorm_dec(m, L)
frontier = [('O', 'E')]; nodes = 0
for depth in range(10):
    nxt = []
    for (u, v) in frontier:
        (mu, Lu), (mv, Lv) = counts(u), counts(v)
        ok("F3 Christoffel basis pairs from (O,E) have Farey-neighbour slopes m/L (|mu Lv - mv Lu| = 1)", abs(mu*Lv - mv*Lu) == 1)
        x, y, z = tOE(u), tOE(v), tOE(u + v)
        # left child (u, uv): (x, z, xz - y); right child (uv, v): (z, y, yz - x)
        ok("F4 left mutation (x, z, xz - y) gives t(u.uv)", abs(tOE(u + u + v) - (x*z - y)) < Decimal(10)**-40)
        ok("F5 right mutation (z, y, yz - x) gives t(uv.v)", abs(tOE(u + v + v) - (y*z - x)) < Decimal(10)**-40)
        mm, LL = counts(u + v)
        d = Decimal(mm)*Decimal(3).ln() - Decimal(LL)*Decimal(2).ln()
        ok("F6 values are 2cosh(delta/2) at the mediant slope (mu+mv)/(Lu+Lv)",
           (mm, LL) == (mu + mv, Lu + Lv) and abs(z - ((d/2).exp() + (-d/2).exp())) < Decimal(10)**-40)
        nxt += [(u, u + v), (u + v, v)]; nodes += 1
    frontier = nxt
CHECKS["F7 Christoffel-tree nodes checked (depth 10)"] = nodes
# HYP-9230: theta_f = f log2 3 - round(f log2 3) equals delta/ln2 for words with (|w|, Sw) = (f, round(f log2 3))
log23 = Decimal(3).ln()/Decimal(2).ln()
for f in (5, 12, 41, 53, 94, 147, 306, 665):
    L = int((Decimal(f)*log23).to_integral_value())
    delta = Decimal(f)*Decimal(3).ln() - Decimal(L)*Decimal(2).ln()
    theta = Decimal(f)*log23 - L
    ok("F8 HYP-9230 theta_f = delta/ln2 exactly when Sw = round(f log2 3) (nearest-integer words only)",
       abs(theta - delta/Decimal(2).ln()) < Decimal(10)**-40)

# ======================= G. friezes, minors, faithfulness, merges =======================
def M(a): return ((a, -1), (1, 0))
def Mprod(q):
    R = I2
    for a in q: R = mul(R, M(a))
    return R
mI = ((-1, 0), (0, -1))
ok("G1 quiddity (1,1,1) = triangle: M(1)^3 = -I", Mprod([1, 1, 1]) == mI)
w15 = [1, 2, 2, 2, 2, 2, 2, 1, 7]
ok("G2 Collatz valuation word (1,2,2,2,2,2,2,1,7) is a fan quiddity with M-product -I (frieze note (15))",
   Mprod(w15) == mI and sum(w15) == 3*len(w15) - 6)
# and these words are ACTUAL Collatz words at the frieze note's sources, merging at 24245
wn, yn = Uword(2583211, 9); wh, yh = Uword(7183, 3)
ok("G3 2583211 has U-word (1,2,2,2,2,2,2,1,7), 7183 has (1,1,1), both reach 24245 (a merge of frieze-closing words)",
   wn == w15 and wh == [1, 1, 1] and yn == yh == 24245)
for _ in range(2000):
    w = rand_word(10, 8); Mw = G(w)
    ok("G4 G_w is never +-I nor scalar (diagonal 3^m != 2^S)", Mw[0][0] != Mw[1][1])
for a in range(1, 9):
    for b in range(1, 9):
        ok("G5 ear identity M(a)M(b) = M(a+1)M(1)M(b+1)", mul(M(a), M(b)) == Mprod([a + 1, 1, b + 1]))
        g1, g2 = G([a, b]), G([a + 1, 1, b + 1])
        ok("G6 the ear move is not a G-relation: it changes (|w|, Sw) by (+1, +3), so traces already differ",
           (g1[0][0], g1[1][1]) == (9, 2**(a + b)) and (g2[0][0], g2[1][1]) == (27, 2**(a + b + 3))
           and Fr(tr(g1)**2, det(g1)) != Fr(tr(g2)**2, det(g2)))
# minors of the marked configuration v_i = (3^i, B_i) are carries of subwords: cluster coordinates SEE the cocycle
for _ in range(600):
    w = rand_word(9, 6); r = len(w)
    Bp = [carry_rec(w[:i]) for i in range(r + 1)]; Ap = [sum(w[:i]) for i in range(r + 1)]
    D = {}
    for i in range(r + 1):
        for j in range(i + 1, r + 1):
            D[i, j] = 3**i * Bp[j] - 3**j * Bp[i]
            ok("G7 det((3^i,B_i),(3^j,B_j)) = 3^i 2^(A_i) B_(w[i+1..j]) > 0", D[i, j] == 3**i * 2**Ap[i] * carry_rec(w[i:j]) > 0)
    for i, j, k, l in itertools.combinations(range(r + 1), 4):
        ok("G8 three-term Pluecker (type-A cluster exchange) relation holds among carries of overlapping subwords",
           D[i, k]*D[j, l] == D[i, j]*D[k, l] + D[i, l]*D[j, k])
# faithfulness: w -> G_w injective (also projectively) on all words with Sw <= 16; explicit 2-adic decoder
def compositions(S):
    for mask in range(2**(S - 1)):
        parts, cur = [], 1
        for bit in range(S - 1):
            if mask >> bit & 1: parts.append(cur); cur = 1
            else: cur += 1
        parts.append(cur); yield parts
def decode(m, S, B):
    w, A = [], 0
    for i in range(1, m):
        # B = sum_{i'} 3^(m-i') 2^(A_(i'-1)); peel the known head and read the next prefix sum as a 2-adic valuation
        B -= 3**(m - i) * 2**A
        nxtA = v2(B)
        w.append(nxtA - A); A = nxtA
    B -= 2**A
    if B != 0: raise AssertionError('decode: residue')
    w.append(S - A); return w
keys = set(); nwords = 0
for S in range(1, 17):
    for w in compositions(S):
        g = G(w); key = (g[0][0], g[0][1], g[1][1])
        ok("G9 2-adic decoder recovers w from (|w|, Sw, B_w)", decode(len(w), S, g[0][1]) == w)
        keys.add(key); nwords += 1
ok("G10 w -> G_w injective on all words with Sw <= 16 (projectively too: scalar would need 3^dm = 2^dS)",
   len(keys) == nwords == 2**16 - 1)
CHECKS["G10 words enumerated"] = nwords
# an equal-count merge decided purely by carries (identical traces): u = (1,7), v = (7,1), (m,S) = (2,8)
u, v = [1, 7], [7, 1]
ok("G11 equal counts => equal traces, different carries; F_u(n) = F_v(h) iff 9(n - h) = B_v - B_u = 126",
   tr(G(u)) == tr(G(v)) and carry_rec(v) - carry_rec(u) == 126)
found = None
for n in range(15, 2*10**6, 2):
    if Uword(n, 2)[0] == u and Uword(n - 14, 2)[0] == v:
        found = (n, n - 14, Uword(n, 2)[1]); break
ok("G12 an actual equal-count merge exists (sources n, n-14 with U-words (1,7), (7,1))", found is not None)
CHECKS["G12 first equal-count merge (n, h, common value)"] = found

# ======================= H. extra checks used in the report =======================
# order matters: G_1 G_2 != G_2 G_1, so "representation" must be read with chronological (reversed) composition
ok("H1 G_(1,2) = G_2 G_1 = [[9,5],[0,8]] differs from G_1 G_2 = [[9,7],[0,8]] (anti-homomorphism is not cosmetic)",
   G([1, 2]) == ((9, 5), (0, 8)) and mul(G([1]), G([2])) == ((9, 7), (0, 8)))
# the frieze transfer representation a -> M(a) is irreducible: tr[M(a),M(b)] = 2 + (a-b)^2 (never -2, never 2 for a != b)
for a in range(-6, 12):
    for b in range(-6, 12):
        ok("H2 tr[M(a),M(b)] = 2 + (a-b)^2: frieze transfer pairs are irreducible for a != b, never at the Markov level -2",
           commut_tr(M(a), M(b)) == 2 + (a - b)**2)
# U-letter example on the cubic: letters 1 and 2 give (5/sqrt6, 7/sqrt12, 17/sqrt72)
g1, g2 = G([1]), G([2]); g12 = G([1, 2])
ok("H3 U-letters (1),(2): traces/dets (5,6), (7,12), (17,72) and kappa = 4",
   (tr(g1), det(g1), tr(g2), det(g2), tr(g12), det(g12)) == (5, 6, 7, 12, 17, 72) and kappa_gl2(g1, g2) == 4)
# THM-4600's integral cycle points
for w, c in (([1], -1), ([2], 1), ([1, 2], -5), ([1, 1, 1, 2, 1, 1, 4], -17)):
    ok("H4 THM-4600 integral anchors: c_1 = -1, c_2 = 1, c_12 = -5, c_1112114 = -17", cw(w) == c)
# c_w = c_w' iff w, w' are powers of one primitive word (maximal split submonoids are <p>, p primitive), Sw <= 12
def prim_root(w):
    n = len(w)
    for d in range(1, n + 1):
        if n % d == 0 and w[:d] * (n // d) == w: return tuple(w[:d])
groups = {}
for S in range(1, 13):
    for w in compositions(S):
        groups.setdefault(cw(w), set()).add(prim_root(w))
ok("H5 c_w = c_w' iff same primitive root (Sw <= 12): the cocycle splits exactly on submonoids of <p>",
   all(len(s) == 1 for s in groups.values()))
CHECKS["H5 distinct cycle points (Sw <= 12)"] = len(groups)
# endpoint label mod 3 = (-1)^(last letter): not a function of the counts (trace-invisible), decides the ear obstruction
for _ in range(2000):
    w = rand_word(8, 7); m, S = len(w), sum(w)
    lab = (carry_rec(w) * pow(2, -S, 3**m)) % 3**m          # frieze note (12): c(w) = B_w Q_w^-1 mod 3^|w|
    ok("H6 endpoint label mod 3 = (-1)^(a_last) mod 3", lab % 3 == (1 if w[-1] % 2 == 0 else 2))
ok("H7 equal counts, different labels: (1,2) and (2,1) have equal traces but labels 1 and 2 mod 3",
   tr(G([1, 2])) == tr(G([2, 1])) and (carry_rec([1, 2])*pow(2, -3, 9)) % 3 != (carry_rec([2, 1])*pow(2, -3, 9)) % 3)

# carry exchange relation (the 3-term Pluecker relation of (3^i, B_i) divided by common factors):
#   B_(xy) B_(yz) = B_y B_(xyz) + 3^|y| 2^(Sy) B_x B_z   for all words x, y, z (empty allowed)
for _ in range(5000):
    x, y, z = ([rnd.randint(1, 6) for _ in range(rnd.randint(0, 6))] for _ in range(3))
    B = carry_rec
    ok("H8 carry exchange relation B_xy B_yz = B_y B_xyz + 3^|y| 2^(Sy) B_x B_z (a cluster identity ON the cocycle side)",
       B(x + y)*B(y + z) == B(y)*B(x + y + z) + 3**len(y) * 2**sum(y) * B(x)*B(z))

print("THM-4604 audit D: all assertions passed")
for k, val in CHECKS.items():
    print(f"  {k}: {val}")
