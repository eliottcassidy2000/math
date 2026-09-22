#!/usr/bin/env python3
"""Pythagorean semicircle lane of session collatz-mod6-20260917.

Exact controls for the user's geometric claims about primitive Pythagorean
triples (PPTs): hypotenuse-plus-leg odd squares, the "hypotenuse = square+1"
family, the Thales semicircle normalization, the one-parameter (e/d, l, theta)
variation, angle doubling (Gaussian squaring, x^2-2), the near-isosceles Pell
family, and the 29 = 5^2+2^2 coincidence test against THM-4146.  Section S6
carries the 2026-09-21 audit-driven additions (THM-2142 half-angle map, parity
split, Fermat closure of reading (b1), THM-3341 citations).

Run from the repository root:
  python3 04-computation/experiments/collatz_mod6_20260917_pythagorean_semicircle.py
Requires sympy (symbolic identities only); every load-bearing count uses exact
integer / Fraction arithmetic. All checks stay active under python -O.
"""
from fractions import Fraction
from math import gcd, isqrt
import sys

import sympy as sp


def require(predicate, message):
    if not predicate:
        raise RuntimeError(message)


def is_square(n):
    return n >= 0 and isqrt(n) ** 2 == n


def totient(n):
    result, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result


def omega(n):
    count, p = 0, 2
    while p * p <= n:
        if n % p == 0:
            count += 1
            while n % p == 0:
                n //= p
        p += 1
    return count + (1 if n > 1 else 0)


def ppts(cmax):
    """All PPTs (a,b,c) with a odd leg, b even leg, c<=cmax, with Euclid (m,n)."""
    out = []
    m = 2
    while m * m + 1 <= cmax:
        for n in range(1, m):
            if (m - n) % 2 == 1 and gcd(m, n) == 1:
                c = m * m + n * n
                if c <= cmax:
                    out.append((m * m - n * n, 2 * m * n, c, m, n))
        m += 1
    return out


def ppts_by_hyp(cmax):
    table = {}
    for a, b, c, m, n in ppts(cmax):
        table.setdefault(c, []).append((a, b, c, m, n))
    return table


def banner(text):
    print()
    print("=" * 78)
    print(text)
    print("=" * 78)


# ---------------------------------------------------------------------------
banner("S1. Hypotenuse +- even leg are odd squares; (s,t)=(m+n,m-n) is the identity")
CMAX1 = 1000
P1 = ppts(CMAX1)
print(f"universe: all PPTs with c <= {CMAX1}; count = {len(P1)}")
fibre_count = {}
for a, b, c, m, n in P1:
    s, t = m + n, m - n
    require(a * a + b * b == c * c, "not Pythagorean")
    require(gcd(a, b) == 1, "not primitive")
    require(c + b == s * s and c - b == t * t, "c+-b not (m+-n)^2")
    require(s % 2 == 1 and t % 2 == 1 and gcd(s, t) == 1, "(s,t) not odd coprime")
    require(c + a == 2 * m * m and c - a == 2 * n * n, "c+-a not twice squares")
    require(not is_square(c + a) and not is_square(c - a), "c+-odd leg is a square")
    require(a == s * t and 2 * b == s * s - t * t and 2 * c == s * s + t * t,
            "(s,t) chart fails")
    # half-angle readings, exact: tan(2u) = 2u/(1-u^2)
    u = Fraction(n, m)
    require(2 * u / (1 - u * u) == Fraction(b, a), "tan(phi/2)=n/m fails")
    v = Fraction(t, s)
    require(2 * v / (1 - v * v) == Fraction(a, b), "tan(psi/2)=t/s fails")
    fibre_count[s] = fibre_count.get(s, 0) + 1
print("PROVED and checked on the universe: c+b=(m+n)^2, c-b=(m-n)^2 (odd squares),")
print("c+a=2m^2, c-a=2n^2 (never squares), a=st, b=(s^2-t^2)/2, c=(s^2+t^2)/2,")
print("tan(phi/2)=n/m for phi opposite the even leg, tan(psi/2)=(m-n)/(m+n) for psi")
print("opposite the odd leg.  Exactly ONE leg (the even one) satisfies the claim.")
# fibre of s: number of PPTs with c+b = s^2 equals phi(s)/2 (all t odd coprime, t<s)
SMAX = 201
fib_ok = 0
for s in range(3, SMAX + 1, 2):
    cnt = sum(1 for t in range(1, s, 2) if gcd(s, t) == 1)
    require(cnt == totient(s) // 2, f"fibre size at s={s} is not phi(s)/2")
    fib_ok += 1
print(f"'half of their identity' made exact: for odd s>=3 the number of PPTs with")
print(f"c+b=s^2 is phi(s)/2 (checked s=3..{SMAX}, {fib_ok} values); s alone is one of")
print("two coordinates, t=m-n is the other.  First fibres:")
for s in (3, 5, 7, 9, 15):
    members = [(s * t, (s * s - t * t) // 2, (s * s + t * t) // 2)
               for t in range(1, s, 2) if gcd(s, t) == 1]
    print(f"  s={s}: c+b={s*s}: {members}")

# ---------------------------------------------------------------------------
banner("S2. Readings of 'hypotenuse = square + 1, altitude = square root'")
k, u_, j_ = sp.symbols('k u j', positive=True)
print("(b1) literal: hypotenuse k^2+1, altitude-to-hypotenuse sqrt(k), k a POSITIVE INTEGER")
print("    (for rational k the rational-leg argument below needs k a rational square).")
print("    A right triangle with hypotenuse C and altitude H exists iff C >= 2H")
print("    (legs are roots of X^2 - S X + C H with S^2 = C^2 + 2CH, (a-b)^2 = C^2-2CH).")
poly = sp.expand(u_ ** 4 - 2 * u_ + 1)  # u = sqrt(k): k^2+1-2sqrt(k)
fac = sp.factor(poly)
print(f"    with u=sqrt(k): C-2H = u^4-2u+1 = {fac}")
require(sp.expand(fac - poly) == 0, "factorization wrong")
cubic = sp.Poly(u_ ** 3 + u_ ** 2 + u_ - 1, u_)
root = sp.nroots(cubic.as_expr())[0]
print(f"    real root of u^3+u^2+u-1: u0 = {sp.N(root, 12)}, k0 = u0^2 = {sp.N(root**2, 12)}")
print("    so a real right triangle exists iff k >= 1 or k <= k0; k=1 gives legs sqrt2, sqrt2, hyp 2.")
# integer/rational legs impossible: altitude ab/c = sqrt(k) with integer a,b,c=k^2+1
# => k(k^2+1)^2 = (ab)^2 => k = j^2 => c = j^4+1, ab = j(j^4+1)
# (a+b)^2 = (j^4+1)(j^4+1+2j), (a-b)^2 = (j^4+1)(j^4+1-2j)
JMAX = 2000
rational_leg_hits = []
for j in range(1, JMAX + 1):
    c = j ** 4 + 1
    p, q = c * (c + 2 * j), c * (c - 2 * j)
    if is_square(p) and is_square(q):
        rational_leg_hits.append(j)
print(f"    rational-leg instances need k=j^2 and both (j^4+1)(j^4+1+-2j) squares;")
print(f"    j<=({JMAX}): hits = {rational_leg_hits} (none).  j even: gcd(j^4+1,2j)=1 forces")
print("    j^4+1 square, impossible (PROVED).  j odd: forces (j^4+1)/2 square, i.e.")
print("    j^4+1=2w^2, equivalently a near-isosceles PPT whose leg sum is a square;")
# Pell check of leg sums: x^2-2y^2=-1 numerators 1,7,41,239,...
x, y = 1, 1
pell_square_hits = []
for idx in range(400):
    if is_square(x):
        pell_square_hits.append((idx, x))
    x, y = 3 * x + 4 * y, 2 * x + 3 * y
print(f"    squares among x^2-2y^2=-1 numerators through index 400: {pell_square_hits}")
print("    (only x=1, i.e. j=1, and j=1 is NOT a rational-leg member: (j^4+1+2j)/2 = 2")
print("    is not a square, its legs are sqrt2, sqrt2).")
require(not is_square((1 + 1 + 2) // 2), "j=1 would have rational legs")
# Fermat closure of the j-odd branch: 2w^2 = j^4+1  =>  w^4 - j^4 = ((j^4-1)/2)^2
w_ = sp.symbols('w', positive=True)
lhs = sp.expand((w_ ** 4 - j_ ** 4).subs(w_ ** 4, (j_ ** 4 + 1) ** 2 / 4))
require(sp.expand(lhs - ((j_ ** 4 - 1) / 2) ** 2) == 0, "Fermat identity fails")
require(sp.expand((j_ ** 4 + 1) ** 2 - 4 * j_ ** 4 - (j_ ** 4 - 1) ** 2) == 0, "square identity")
JODD = 10 ** 5
odd_hits = [j for j in range(1, JODD + 1, 2) if is_square((j ** 4 + 1) // 2)]
require(odd_hits == [1], f"j^4+1=2w^2 hits {odd_hits}")
print("    Fermat closure: if j^4+1 = 2w^2 then w^4 - j^4 = ((j^4-1)/2)^2 (identity")
print("    (j^4+1)^2 - 4j^4 = (j^4-1)^2, checked symbolically), a solution of X^4-Y^4=Z^2")
print("    with Z != 0 whenever j > 1.  CITED: Fermat's descent theorem (X^4-Y^4=Z^2 has")
print("    no solution with XYZ != 0; classical, e.g. Hardy-Wright ch. XIII, Mordell")
print("    'Diophantine Equations' ch. 4; theorem numbers not verified here) forces j=1.")
print(f"    Corroboration: odd j<={JODD} with (j^4+1)/2 square: {odd_hits}.")
print("    Hence reading (b1) has NO rational-leg member for any positive integer k:")
print("    PROVED modulo the classical citation (the earlier 'Ljunggren-type' boundary")
print("    is retired; Ljunggren's equation x^2-2y^4=-1 is a different equation).")
# 3-4-5 under (b1): hyp 5 -> k=2, altitude would be sqrt2; actual altitude 12/5
alt345 = Fraction(3 * 4, 5)
require(alt345 == Fraction(12, 5), "3-4-5 altitude")
print(f"    3-4-5: hypotenuse 5 = 2^2+1 but altitude ab/c = {alt345}, not sqrt(2);")
print("    REFUTED as a description of 3-4-5 (minimal witness: the triple itself).")

print("(b2) family (k^2-1, 2k, k^2+1) = Euclid (m,n)=(k,1):")
rows = []
for kk in range(2, 12):
    a, b, c = kk * kk - 1, 2 * kk, kk * kk + 1
    g = gcd(gcd(a, b), c)
    alt = Fraction(a * b, c)
    require(alt == Fraction(2 * kk * (kk * kk - 1), kk * kk + 1), "altitude formula")
    rows.append((kk, (a, b, c), "primitive" if g == 1 else f"gcd {g}", alt))
for r in rows:
    print(f"    k={r[0]:2d}: {r[1]}  {r[2]:>10}  altitude={r[3]}")
print("    hypotenuse = k^2+1 EXACTLY; primitive iff k even; k>=2 (k=1 is the degenerate")
print("    (0,2,2)); first member k=2 is 3-4-5;")
print("    the half-angle triangle (k,1,sqrt(k^2+1)) has leg k = sqrt(k^2): this is the")
print("    only length equal to 'the square root of that same number'.  Angle opposite")
print("    the even leg = 2*arctan(1/k) -> 0, vertex -> diameter endpoint.  HEURISTIC")
print("    reading match: (b2) fits the words 'hypotenuse one greater than a square,")
print("    3-4-5 first'; the 'altitude sqrt' clause fits only the half-angle leg k.")
print("(b6) THM-3335 family: hypotenuse = (even square) + 1, even leg = that square:")
n_, q_ = 0, 0
x, y = 3, 1  # x^2-8y^2=1 -> n=(x-1)/2, q=y
hits = []
for _ in range(4):
    n_, q_ = (x - 1) // 2, y
    if n_ > 0:
        s2 = (2 * q_) ** 2
        a, b, c = 2 * n_ + 1, s2, s2 + 1
        require(a * a + b * b == c * c, "THM-3335 member not Pythagorean")
        hits.append((a, b, c))
    x, y = 3 * x + 8 * y, x + 3 * y
print(f"    first members {hits}: (3,4,5) is first here too, but 'altitude sqrt'")
print("    matches nothing in the triangle (sqrt(4)=2 is not a side or altitude).")

# ---------------------------------------------------------------------------
banner("S3. Thales semicircle normalization and the (e/d, l, theta) chart")
th = sp.symbols('theta', positive=True)
ed = sp.tan(th) ** 2
l_expr = sp.sin(2 * th) / 2
print("Thales: the right-angle vertex lies on the circle with the hypotenuse as")
print("DIAMETER (radius 1/2 after normalization), not as radius; the central angle")
print("is 2*theta (inscribed-angle theorem).  With a<b, theta=arcsin(a/c) in (0,pi/4).")
print("e=a^2/c^2=sin^2, d=b^2/c^2=cos^2, l=ab/c^2=sin(2theta)/2, e/d=tan^2 theta.")
d_ed = sp.simplify(sp.diff(ed, th))
d_l = sp.simplify(sp.diff(l_expr, th))
print(f"  d/dtheta (e/d) = {d_ed}  > 0 on (0,pi/2)")
print(f"  d/dtheta l     = {d_l}   > 0 on (0,pi/4)")
require(sp.simplify(d_l - sp.cos(2 * th)) == 0, "l derivative")
require(sp.simplify(sp.tan(th) ** 2 - sp.sin(th) ** 2 / sp.cos(th) ** 2) == 0, "e/d")
lim = (sp.limit(ed, th, sp.pi / 4), sp.limit(l_expr, th, sp.pi / 4))
print(f"  limits at theta=pi/4: (e/d, l) = {lim}; isosceles endpoint is not a PPT")
print("  PROVED: (e/d, l, theta) is a monotone one-parameter chart in theta.")
CMAX3 = 1000
P3 = ppts(CMAX3)
count_ok = 0
for a, b, c, m, n in P3:
    a1, b1 = (a, b) if a < b else (b, a)   # a1 = shorter leg
    e = Fraction(a1 * a1, c * c)
    d = Fraction(b1 * b1, c * c)
    l = Fraction(a1 * b1, c * c)
    require(e + d == 1 and l * l == e * d and e / d == Fraction(a1 * a1, b1 * b1), "chart")
    # semicircle vertex (centre at origin, radius 1/2) is half the doubled triple
    vx = (d - e) / 2
    vy = l
    require(vx == Fraction(b1 * b1 - a1 * a1, 2 * c * c) and vy == Fraction(2 * a1 * b1, 2 * c * c),
            "vertex is not the half-scaled doubled triple")
    require(vx * vx + vy * vy == Fraction(1, 4), "vertex not on the semicircle")
    # tan(theta/2): which chart?  shorter leg even -> n/m ; shorter leg odd -> t/s
    if a1 == b:
        half = Fraction(n, m)
    else:
        half = Fraction(m - n, m + n)
    require(2 * half / (1 - half * half) == Fraction(a1, b1), "tan(theta/2) chart")
    require(half > 0 and half * half + 2 * half - 1 < 0, "tan(theta/2) not in (0, sqrt2-1)")
    count_ok += 1
print(f"  FINITE-EXACT on {count_ok} PPTs with c<={CMAX3}: e+d=1, l^2=ed, e/d=(a/b)^2,")
print("  vertex ((d-e)/2, l) = (1/2)*(b^2-a^2, 2ab)/c^2 = half-scaled ANGLE-DOUBLED triple,")
print("  and tan(theta/2) is n/m (even shorter leg) or (m-n)/(m+n) (odd shorter leg).")
print("  Parity refinement: tan(theta/2)=p/q reduced with p,q opposite parity <-> theta")
print("  opposite the even leg; p,q both odd <-> theta opposite the odd leg.  Both occur.")
# rational tan(theta) without PPT: legs (1,2)
e12, d12, l12 = Fraction(1, 5), Fraction(4, 5), Fraction(2, 5)
require(l12 * l12 == e12 * d12, "legs 1,2 chart")
print("  Lost-information witness: legs (1,2) give rational e/d=1/4, l=2/5 but hypotenuse")
print("  sqrt5; (e/d,l) depend on tan(theta) only; PPT <=> tan(theta/2) also rational.")

print("Extremal family A (theta->0): (k^2-1,2k,k^2+1), k even, and (s,(s^2-1)/2,(s^2+1)/2), s odd.")
print("Extremal family B (theta->pi/4): near-isosceles |a-b|=1, Pell (m,n)->(2m+n,m):")
m, n = 2, 1
pell = []
for _ in range(8):
    a, b, c = m * m - n * n, 2 * m * n, m * m + n * n
    require(abs(a - b) == 1, "not near-isosceles")
    s, t = m + n, m - n
    require((m - n) ** 2 - 2 * n * n in (1, -1), "t^2-2n^2 != +-1")
    require(s * s - 2 * m * m in (1, -1), "s^2-2m^2 != +-1")
    pell.append((a, b, c, m, n))
    m, n = 2 * m + n, m
for i in range(2, len(pell)):
    c0, c1, c2 = pell[i - 2][2], pell[i - 1][2], pell[i][2]
    require(c2 == 6 * c1 - c0, "hypotenuse recursion fails")
    x0, x1, x2 = (min(p[0], p[1]) for p in pell[i - 2:i + 1])
    require(x2 == 6 * x1 - x0 + 2, "short-leg recursion fails")
# convergents h_i/k_i of sqrt2-1 = [0;2,2,2,...] and their mediants (intermediate fractions)
conv = [Fraction(0, 1), Fraction(1, 2)]
while len(conv) < len(pell) + 1:
    h = 2 * conv[-1].numerator + conv[-2].numerator
    kk_ = 2 * conv[-1].denominator + conv[-2].denominator
    conv.append(Fraction(h, kk_))
mediants = [Fraction(conv[i - 1].numerator + conv[i].numerator,
                     conv[i - 1].denominator + conv[i].denominator) for i in range(1, len(conv))]
pell_numbers = [0, 1]
while len(pell_numbers) < 2 * len(pell) + 2:
    pell_numbers.append(2 * pell_numbers[-1] + pell_numbers[-2])
small_seq = []
for i, p in enumerate(pell):
    a, b, c, m, n = p
    phi_half, psi_half = Fraction(n, m), Fraction(m - n, m + n)
    require(phi_half == conv[i + 1], f"n/m is not convergent {i+1} of sqrt2-1")
    require(psi_half == mediants[i], f"t/s is not the mediant of convergents {i},{i+1}")
    require(psi_half not in conv, "t/s unexpectedly a convergent")
    require(c == pell_numbers[2 * i + 3], "hypotenuse is not the odd-index Pell number pell_(2i+3)")
    theta_half = psi_half if a < b else phi_half     # theta = smaller angle (S3 convention)
    small_seq.append(theta_half)
    which = "psi (odd leg shorter)" if a < b else "phi (even leg shorter)"
    print(f"    (a,b,c)=({a},{b},{c}), (m,n)=({m},{n}), (s,t)=({m+n},{m-n}), "
          f"tan(phi/2)=n/m={phi_half}, tan(psi/2)=t/s={psi_half}, smaller angle = {which}")
require(all(small_seq[i] < small_seq[i + 1] for i in range(len(small_seq) - 1)), "smaller-angle sequence not monotone")
require(all(h * h + 2 * h - 1 < 0 for h in small_seq), "smaller-angle half-tangent not < sqrt2-1")
print("    PROVED: |a-b|=1 <=> (m-n)^2-2n^2=+-1 (Pell); c_{i+1}=6c_i-c_{i-1}, x_{i+1}=6x_i-x_{i-1}+2.")
print(f"    tan(phi/2)=n/m = {', '.join(str(x) for x in conv[1:len(pell)+1])} are the CONVERGENTS of sqrt2-1=[0;2,2,...];")
print(f"    tan(psi/2)=t/s = {', '.join(str(x) for x in mediants[:len(pell)])} are the INTERMEDIATE fractions")
print("    (mediants of consecutive convergents; ratios of consecutive companion Pell numbers 1,3,7,17,41,99,...),")
print("    never convergents.  With theta the SMALLER angle (S3 convention) tan(theta/2) alternates:")
print(f"    {', '.join(str(x) for x in small_seq)} (monotone, < sqrt2-1).")
print(f"    hypotenuses = odd-index Pell numbers pell_3, pell_5, ... = {[pell_numbers[2*i+3] for i in range(len(pell))]}")
print("    (pell_0=0, pell_1=1, pell_(r+2)=2pell_(r+1)+pell_r; THM-3341 (12) also lists the degenerate m_0=pell_1=1):")
print("    INHERITED: this family, its recursions and the odd-index-Pell hypotenuses are")
print("    THM-3341 (6)-(12) and (24) (the Berggren A_+ ray); THM-3357 (22) lists the same ray.")
CMAX_NI = 10 ** 6
ni_all = sorted(c for a, b, c, m, n in ppts(CMAX_NI) if abs(a - b) == 1)
ni_pell = sorted(p[2] for p in pell if p[2] <= CMAX_NI)
require(ni_all == ni_pell, f"near-isosceles census mismatch {ni_all} vs {ni_pell}")
print(f"    FINITE-EXACT: near-isosceles hypotenuses c<=10^6 are exactly {ni_all}.")

print("29 test against THM-4146 (D=a^2+2b^2-ab, cycle -(a+b)->h->-(b-a)):")


def G_orbit(a, b, h):
    D = a * a + 2 * b * b - a * b
    pts = [Fraction(-(a + b))]
    for _ in range(3):
        pts.append((pts[-1] ** 2 - D) / b)
    return D, pts


D345, orb345 = G_orbit(3, 4, 5)
require(D345 == 29 and orb345[0] == orb345[3] and orb345[1] == 5 and orb345[2] == -1, "3-4-5 cycle")
D2029, orb2029 = G_orbit(20, 21, 29)
D2120, orb2120 = G_orbit(21, 20, 29)
print(f"    (3,4,5): D={D345}, orbit {[str(v) for v in orb345]}  (closes, THM-4146 (33))")
print(f"    (20,21,29): D={D2029}, orbit {[str(v) for v in orb2029]}  (does NOT close)")
print(f"    (21,20,29): D={D2120}, orbit {[str(v) for v in orb2120]}  (does NOT close either)")
require(orb2029[0] != orb2029[3] and orb2029[1] != 29, "(20,21,29) unexpectedly cycles")
require(D2029 == 862 and D2120 == 821, "D values of the two labelings")
require(orb2120[0] != orb2120[3] and orb2120[1] != 29 and 29 not in orb2120, "(21,20,29) unexpectedly cycles")
require(3 * 20 - 21 != 29 and 3 * 21 - 20 != 29, "THM-4146 (31) 3a-b=h holds for a (20,21,29) labeling?")
mm, nn = sp.symbols('m n', integer=True, positive=True)
A_, B_, C_ = mm ** 2 - nn ** 2, 2 * mm * nn, mm ** 2 + nn ** 2
P_D = sp.expand(A_ ** 2 + 2 * B_ ** 2 - A_ * B_)          # THM-4146 constant at Euclid (m,n)
P_pell = sp.expand((2 * mm + nn) ** 2 + mm ** 2)           # Pell successor hypotenuse
diff = sp.factor(P_D - P_pell)
print(f"    D(m,n) = {P_D}")
print(f"    c_next(m,n) = {P_pell}")
print(f"    D - c_next = {diff}")
agree = [(m, n) for m in range(1, 301) for n in range(1, m) if P_D.subs({mm: m, nn: n}) == P_pell.subs({mm: m, nn: n})]
print(f"    integer agreement with 300>=m>n>=1: {agree}")
require(agree == [(2, 1)], "unexpected agreement locus")
# all-(m,n) proof of the agreement locus
at_n1 = sp.factor((P_D - P_pell).subs(nn, 1))
require(sp.expand(at_n1 - mm * (mm - 2) * (mm ** 2 + 1)) == 0, "n=1 factorization")
decomp = (mm ** 2 * ((mm - nn) ** 2 - 1) + mm ** 2 * (5 * nn ** 2 - 4)
          + 2 * mm * nn * (nn ** 2 - 2) + nn ** 2 * (nn ** 2 - 1))
require(sp.expand(P_D - P_pell - decomp) == 0, "positivity decomposition")
require(sp.expand(P_D - (C_ ** 2 + B_ * (B_ - A_))) == 0, "D = c^2 + b(b-a)")
print(f"    at n=1: D - c_next = {at_n1}, zero only at m=2;")
print("    for n>=2: D - c_next = m^2((m-n)^2-1) + m^2(5n^2-4) + 2mn(n^2-2) + n^2(n^2-1) > 0")
print("    (every term >= 0 for m>n>=2, the last > 0).  PROVED: D(m,n)=c_next(m,n) only at")
print("    (m,n)=(2,1) for ALL m>n>=1; the 300-search is corroboration only.")
print("    REFUTED as a map: 29 is 5^2+2^2 in THM-4146 (c^2 + b(b-a), b(b-a)=4 at 3-4-5)")
print("    and c=m^2+n^2 for the Pell pair (5,2); the two polynomials in (m,n) agree only")
print("    at (2,1).  Neither labeling of (20,21,29) is in the 3:4:5 class (THM-4146 (31)")
print("    needs 3a-b=h) and it plays no role in the cycle.  THM-4146 (33a) already types")
print("    29=5^2+2^2=2^2+3^2+4^2 as a non-analogy; new content here is the negative test only.")

# ---------------------------------------------------------------------------
banner("S4. Angle doubling = Gaussian squaring; hypotenuse squaring forest; x^2-2")
X = sp.symbols('x')
require(sp.expand(2 * (2 * X ** 2 - 1) - ((2 * X) ** 2 - 2)) == 0, "conjugacy 2x^2-1 ~ y^2-2")
print("PROVED: cos(2theta)=2cos^2-1; with y=2cos(theta), y -> y^2-2 (Chebyshev T_2).")
a_, b_ = sp.symbols('a b', real=True)
z2 = sp.expand((b_ + sp.I * a_) ** 2)
require(sp.re(z2) == b_ ** 2 - a_ ** 2 and sp.im(z2) == 2 * a_ * b_, "Gaussian square")
print("PROVED: (b+ai)^2 = (b^2-a^2) + 2ab i, |.|=c^2: the doubling (a,b,c)->(2ab,|b^2-a^2|,c^2)")
print("is Gaussian squaring; THM-3333's lift Phi(m,n) = (m+ni)^2 is the SAME map, so the")
print("half-angle (m,n) triangle -> PPT is one doubling step.  INHERITED: this map is")
print("THM-3341 (23), Gamma(a,b,c)=(|a^2-b^2|,2ab,c^2), with (3,4,5)->(7,24,25) and")
print("(21,20,29)->(41,840,841) at THM-3341 (26).")
CMAX4 = 10 ** 4
T4 = ppts_by_hyp(CMAX4)
all4 = [t for lst in T4.values() for t in lst]
print(f"universe: PPTs with c<={CMAX4}: {len(all4)}")
# doubles of PPTs are PPTs (primitive)
for a, b, c, m, n in ppts(100):
    A2, B2, C2 = 2 * a * b, abs(b * b - a * a), c * c
    require(A2 * A2 + B2 * B2 == C2 * C2 and gcd(A2, B2) == 1, "double not a PPT")
    require(any((t[0], t[1]) in ((B2, A2), (A2, B2)) for t in T4[C2]), "double missing")
square_hyp = [t for t in all4 if is_square(t[2])]
print(f"PPTs with square hypotenuse, c<={CMAX4}: {len(square_hyp)}")
halves = 0
for a, b, c, m, n in square_hyp:
    r = isqrt(c)
    require(m * m + n * n == r * r, "Euclid parameters not a Pythagorean pair")
    half = [t for t in T4[r] if (t[0], t[1]) in ((m, n), (n, m))]
    require(len(half) == 1, "half not unique")
    halves += 1
    print(f"  ({a},{b},{c}) = double of {half[0][:3]}   (m,n)=({m},{n}), sqrt(c)={r}")
require(halves == len(square_hyp), "halves count")
n_small = sum(len(T4[c]) for c in T4 if c <= isqrt(CMAX4))
require(n_small == len(square_hyp), "doubling not a bijection onto square hypotenuses")
for c in sorted(T4):
    if c <= isqrt(CMAX4):
        require(len(T4[c]) == 2 ** (omega(c) - 1) == len(T4[c * c]), "2^(omega-1) count")
print(f"PROVED: (a,b,c) is an angle-double of a PPT <=> c is a perfect square <=> (n,m,sqrt c)")
print(f"is a PPT; doubling is a bijection PPT(c) -> PPT(c^2), both of size 2^(omega(c)-1)")
print("(for c a PPT hypotenuse, i.e. every prime factor 1 mod 4: THM-3334 (36), Tripathi;")
print("THM-3341 (31) uses the same count).  The all-PPT bijection is the only new statement;")
print("THM-3341 section 2 has square hypotenuse <=> Gamma-image on the U-spine only.")
print(f"FINITE-EXACT: #PPT(c<={isqrt(CMAX4)}) = {n_small} = #square-hypotenuse PPTs (c<={CMAX4}).")
print("Hypotenuse orbit c, c^2, c^4, ... is a ray of the squaring forest of")
print("arithmetic_braids_20260917_summand.md section 2; on hypotenuses the square-root")
print("ancestry lifts bijectively to PPT halving (the summand note's hostile concerned")
print("arbitrary integers, not hypotenuses).")
chain = (3, 4, 5)
print("  doubling chain from 3-4-5:", end=" ")
for _ in range(4):
    print(chain, end=" -> ")
    a, b, c = chain
    chain = (2 * a * b, abs(b * b - a * a), c * c)
print("...")

# integer / rational preperiodic set of y^2-2
pre = {}
for y0 in range(-3, 4):
    y, seen = y0, []
    while abs(y) <= 2 and y not in seen:
        seen.append(y)
        y = y * y - 2
    pre[y0] = (seen, "escapes" if abs(y) > 2 else f"returns to {y}")
for y0 in range(-2, 3):
    require(pre[y0][1].startswith("returns"), f"{y0} not preperiodic")
for y0 in (-3, 3):
    require(pre[y0][1] == "escapes", f"{y0} preperiodic?")
print("PROVED: |y|>=3 => |y^2-2| > |y| (escape); odd-prime denominators double their")
print("valuation and 2-adic ones too, so PrePer(y^2-2, Q) = {-2,-1,0,1,2}.  Graph:")
angle = {2: "0", 1: "pi/3", 0: "pi/2", -1: "2pi/3", -2: "pi"}
for y0 in (2, -2, 0, 1, -1):
    print(f"  {y0:2d} (2cos {angle[y0]:>5}) -> {y0*y0-2:2d}")
print("  = angle doubling on {0, pi/3, pi/2, 2pi/3, pi}: pi->0, pi/2->pi, pi/3->2pi/3->4pi/3~2pi/3.")
# 3-cycles of y^2-2 at 2cos(2pi k/9) and 2cos(2pi k/7)
Y = sp.symbols('y')
for N, minpoly in ((9, Y ** 3 - 3 * Y + 1), (7, Y ** 3 + Y ** 2 - 2 * Y - 1)):
    f3 = Y
    for _ in range(3):
        f3 = sp.expand(f3 ** 2 - 2)
    r3 = sp.rem(f3 - Y, minpoly, Y)
    r1 = sp.rem(sp.expand(Y ** 2 - 2 - Y), minpoly, Y)
    require(r3 == 0 and r1 != 0, f"2cos(2pi/{N}) not exact period 3")
    ordN = next(e for e in range(1, 20) if pow(2, e, N) == 1)
    print(f"  N={N}: ord_{N}(2)={ordN}; 2cos(2pi k/{N}), gcd(k,{N})=1 (min poly {minpoly}) has exact period 3.")
require(pow(2, 3, 7) == 1 and pow(2, 3, 9) != 1, "ord_7(2)=3 / ord_9(2)!=3")
require((-1) ** 2 - 2 == -1, "2cos(2pi*3/9) = -1 is fixed")
print("  (k=3 at N=9 gives 2cos(2pi/3) = -1, a fixed point, hence the gcd condition.)")
print("  N=7 has ord_7(2)=3 already, so its 3-cycle needs no -1 fold: the 'order six over")
print("  three via a central -1' shape below is N=9-specific.")
print("  N=9: the angle orbit {1,2,4,8,7,5} has order 6 = ord_9(2) (the same fact as the")
print("  inherited 'exponents 1,5,3 mod 6' row law); its fold by -1 is the 3-cycle {1,2,4}.")
rows_res = {2, 5, 8}
require(all((-r) % 9 not in rows_res for r in rows_res), "row residues not a -1 transversal")
print("  The Collatz row residues {2,5,8} mod 9 meet each +-pair {1,8},{2,7},{4,5} exactly once")
print("  (PROVED trivially: -1 swaps the classes 1,2 mod 3).  Order-six-over-three via a")
print("  central -1 is the shape shared with THM-4139's B^3=-I; no map between the two")
print("  3-cycles exists (cubic irrationals versus quarter-integers).")

# ---------------------------------------------------------------------------
banner("S5. The user's 'x^2+{0,1,2}' graphs are x^2-{0,1,2}; sign audit")


def int_graph(cc, box=3):
    edges, cyc = [], set()
    for y0 in range(-box, box + 1):
        y1 = y0 * y0 + cc
        if abs(y1) <= box:
            edges.append((y0, y1))
    for y0 in range(-box, box + 1):
        y, k = y0, 0
        while abs(y) <= box and k <= 2 * box + 2:
            y = y * y + cc
            k += 1
            if y == y0:
                cyc.add(y0)
                break
    return edges, sorted(cyc)


for cc in (0, -1, -2, 1, 2):
    e, cyc = int_graph(cc)
    print(f"  x^2{cc:+d}: edges inside [-3,3]: {e}; periodic points: {cyc}")
e0, c0 = int_graph(0)
e1, c1 = int_graph(-1)
e2, c2 = int_graph(-2)
require(c0 == [0, 1] and (-1, 1) in e0, "x^2 graph")
require(c1 == [-1, 0] and (1, 0) in e1, "x^2-1 graph")
require(c2 == [-1, 2] and (1, -1) in e2 and (0, -2) in e2 and (-2, 2) in e2, "x^2-2 graph")
require(int_graph(1)[1] == [] and int_graph(2)[1] == [], "x^2+1/x^2+2 have integer cycles?")
print("  User's descriptions (self-loops 0,1 with -1->1; cycle -1<->0 with 1->0; self-loops")
print("  -1,2 fed by 1 and by 0->-2) match x^2-0, x^2-1, x^2-2 exactly and x^2+1, x^2+2 not")
print("  at all (no integer cycles).  REFUTED as written, survivor: the sign is minus,")
print("  consistent with the user's own x^2-7/4 and x^2-29/16.  Typing: x^2 is z->z^2 on")
print("  the unit circle (angle doubling on the circle), x^2-2 its trace projection")
print("  y=z+1/z; x^2-1 has no circle/angle reading (not conjugate to a power or Chebyshev map).")

# ---------------------------------------------------------------------------
banner("S6. Audit-driven additions (2026-09-21): THM-2142 map, parity split, half-angle triangles")
n_odd_short = n_even_short = 0
for a, b, c, m, n in P1:
    s, t = m + n, m - n
    # THM-2142: b(a(cos phi)) = (1+cos phi)/2 = cos^2(phi/2), evaluated on a PPT
    require(Fraction(c + a, 2 * c) == Fraction(m * m, c), "cos^2(phi/2) != m^2/c")
    require(Fraction(c - a, 2 * c) == Fraction(n * n, c), "sin^2(phi/2) != n^2/c")
    require(Fraction(c + b, 2 * c) == Fraction(s * s, 2 * c), "cos^2(psi/2) != s^2/(2c)")
    require(Fraction(c - b, 2 * c) == Fraction(t * t, 2 * c), "sin^2(psi/2) != t^2/(2c)")
    # consistency with tan(phi/2)=n/m: cos^2 = 1/(1+tan^2)
    require(1 / (1 + Fraction(n * n, m * m)) == Fraction(m * m, c), "half-angle cosine vs tangent")
    require(1 / (1 + Fraction(t * t, s * s)) == Fraction(s * s, 2 * c), "half-angle cosine vs tangent (psi)")
    # the (n,m,sqrt c) triangle halves phi ONLY: its other acute angle is pi/2-phi/2, not psi/2
    require(Fraction(m, n) != Fraction(t, s), "m/n = t/s?")
    if a < b:
        n_odd_short += 1
    else:
        n_even_short += 1
print("PROVED (map to THM-2142 section 1, b(a(cos theta)) = (1+cos theta)/2 = cos^2(theta/2)):")
print("  on a PPT, cos^2(phi/2) = (c+a)/(2c) = m^2/c, sin^2(phi/2) = (c-a)/(2c) = n^2/c,")
print("  cos^2(psi/2) = (c+b)/(2c) = s^2/(2c), sin^2(psi/2) = (c-b)/(2c) = t^2/(2c);")
print("  i.e. 'hypotenuse +- leg' is 2c times the half-angle cosine/sine squared, and")
print(f"  Theorem 1 is THM-2142's half-angle functional evaluated at PPT angles (checked on {len(P1)} PPTs).")
print("  The (n,m,sqrt c) triangle is the half-angle triangle of phi only (its other acute")
print("  angle is pi/2-phi/2, and m/n != t/s always); psi is halved by the (t,s,sqrt(2c)) triangle.")
require(n_odd_short + n_even_short == len(P1), "parity split total")
print(f"FINITE-EXACT parity split of the smaller angle over the {len(P1)} PPTs with c<={CMAX1}:")
print(f"  opposite the odd leg (tan(theta/2)=t/s, both odd): {n_odd_short};")
print(f"  opposite the even leg (tan(theta/2)=n/m, opposite parity): {n_even_short}.  Both charts occur.")
print("SCOPE: the arithmetic_braids2_20260917_* notes contain no Pythagorean, Pell or")
print("Gaussian-squaring content (grep: no match), so no braids2 inheritance applies here;")
print("the braids geometry note's Gaussian squaring (its section 5, via THM-3336) concerns")
print("fresh leg primes and content, not the hypotenuse bijection of S4.")

print()
print("ALL CHECKS PASSED")
