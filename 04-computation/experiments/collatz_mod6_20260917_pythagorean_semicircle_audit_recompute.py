#!/usr/bin/env python3
"""Independent recomputation for the pythagorean_semicircle lane (adversarial audit).

Does NOT import the explorer's script.  PPTs are enumerated by brute force
a^2+b^2=c^2 (not via Euclid parameters); half-angle tangents use
tan(x/2)=sin x/(1+cos x); the near-isosceles census uses a^2+(a+1)^2=c^2
directly; continued-fraction convergents of sqrt2-1 are computed and compared
with the explorer's t/s sequence; both (a,b) assignments are tried in the
THM-4146 cycle test; the j-odd branch of reading (b1) is closed by an explicit
reduction to Fermat's x^4-y^4=z^2.

Run:  python3 04-computation/experiments/collatz_mod6_20260917_pythagorean_semicircle_audit_recompute.py
"""
from fractions import Fraction
from math import gcd, isqrt
import numpy as np


def need(ok, msg):
    if not ok:
        raise RuntimeError("AUDIT FAIL: " + msg)


def is_sq(n):
    return n >= 0 and isqrt(n) ** 2 == n


def phi(n):
    r, m, p = n, n, 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            r -= r // p
        p += 1
    if m > 1:
        r -= r // m
    return r


def omega(n):
    k, p = 0, 2
    while p * p <= n:
        if n % p == 0:
            k += 1
            while n % p == 0:
                n //= p
        p += 1
    return k + (n > 1)


def brute_ppts(cmax):
    """All PPTs (odd, even, c) with c<=cmax by brute force over a<b, numpy assisted."""
    out = []
    bmax = cmax
    for b in range(2, bmax):
        a = np.arange(1, b, dtype=np.int64)
        s = a * a + b * b
        r = np.floor(np.sqrt(s.astype(np.float64))).astype(np.int64)
        for rr in (r - 1, r, r + 1):
            pass
        ok = (r * r == s) | ((r + 1) * (r + 1) == s) | ((r - 1) * (r - 1) == s)
        for aa in a[ok]:
            aa = int(aa)
            c2 = aa * aa + b * b
            c = isqrt(c2)
            if c * c == c2 and c <= cmax and gcd(aa, b) == 1:
                odd, even = (aa, b) if aa % 2 else (b, aa)
                out.append((odd, even, c))
    return sorted(set(out), key=lambda t: (t[2], t[0]))


def sec(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


# --------------------------------------------------------------------------
sec("A1. PS1/PS2/PS3: brute-force PPTs c<=1000, c+-leg squares, fibre phi(s)/2, half angles")
P = brute_ppts(1000)
print("brute-force #PPT(c<=1000) =", len(P))
need(len(P) == 158, "count 158")
both_odd_cases = 0
opp_par_cases = 0
for a, b, c in P:
    need(a % 2 == 1 and b % 2 == 0 and gcd(a, b) == 1 and a * a + b * b == c * c, "PPT")
    need(is_sq(c + b) and is_sq(c - b) and isqrt(c + b) % 2 == 1 and isqrt(c - b) % 2 == 1, "c+-b odd squares")
    need(not is_sq(c + a) and not is_sq(c - a), "c+-a never square")
    need((c + a) % 2 == 0 and is_sq((c + a) // 2) and is_sq((c - a) // 2), "c+-a twice squares")
    s, t = isqrt(c + b), isqrt(c - b)
    need(gcd(s, t) == 1 and a == s * t and 2 * b == s * s - t * t and 2 * c == s * s + t * t, "(s,t) chart")
    # half-angle tangents via sin/(1+cos): angle opposite even leg b: tan(phi/2)=b/(c+a)
    hb = Fraction(b, c + a)          # = n/m
    ha = Fraction(a, c + b)          # = t/s
    need(ha == Fraction(t, s), "tan(psi/2)=t/s")
    m, n = (s + t) // 2, (s - t) // 2
    need(hb == Fraction(n, m) and (m - n) % 2 == 1 and gcd(m, n) == 1, "tan(phi/2)=n/m")
    # parity of the reduced half-angle tangent of the SMALLER angle
    if a < b:
        h = ha
        need(h.numerator % 2 == 1 and h.denominator % 2 == 1, "both odd <-> odd shorter leg")
        both_odd_cases += 1
    else:
        h = hb
        need((h.numerator + h.denominator) % 2 == 1, "opposite parity <-> even shorter leg")
        opp_par_cases += 1
    need(h * h + 2 * h - 1 < 0, "tan(theta/2) < sqrt2-1")
print("smaller angle opposite odd leg (both-odd tan(theta/2)):", both_odd_cases,
      "; opposite even leg (opposite-parity):", opp_par_cases)
need(both_odd_cases > 0 and opp_par_cases > 0, "both parity cases occur")
# fibre over s by brute force: PPTs with c+b=s^2, from a bigger brute-force universe
P2 = brute_ppts(6000)
fib = {}
for a, b, c in P2:
    fib[isqrt(c + b)] = fib.get(isqrt(c + b), 0) + 1
for s in range(3, 76, 2):
    # need c up to (s^2+(s-1)^2)/2 <= 6000  -> s<=77
    need(fib.get(s, 0) == phi(s) // 2, f"fibre s={s}: {fib.get(s,0)} vs {phi(s)//2}")
print("fibre |{PPT: c+b=s^2}| = phi(s)/2 confirmed by brute force for odd s=3..75")
need(1 not in fib and fib.get(3) == 1 and fib.get(5) == 2 and fib.get(7) == 3 and fib.get(9) == 3 and fib.get(15) == 4, "first fibres")
# the (m,n) 'half-angle triangle' claim: half-angle identity c+-a = 2c cos^2, sin^2 of half of phi
for a, b, c in P[:40]:
    # cos^2(phi/2) = (1+cos phi)/2 = (c+a)/(2c) where phi opposite b
    need(Fraction(c + a, 2 * c) == Fraction(((isqrt(c + b) + isqrt(c - b)) // 2) ** 2, c), "cos^2(phi/2)=m^2/c")
print("half-angle map to THM-2142: cos^2(phi/2)=(c+a)/(2c)=m^2/c, cos^2(psi/2)=(c+b)/(2c)=s^2/(2c)")

# --------------------------------------------------------------------------
sec("A2. PS4/PS5: reading (b1) hypotenuse k^2+1, altitude sqrt k")
# factorization check numerically at many u
for u in [Fraction(p, 7) for p in range(0, 30)]:
    need(u ** 4 - 2 * u + 1 == (u - 1) * (u ** 3 + u ** 2 + u - 1), "factorization")
# real root by bisection
lo, hi = 0.0, 1.0
for _ in range(80):
    mid = (lo + hi) / 2
    if mid ** 3 + mid ** 2 + mid - 1 < 0:
        lo = mid
    else:
        hi = mid
u0 = lo
print(f"u0 = {u0:.12f}, k0 = {u0*u0:.12f}")
need(abs(u0 - 0.543689012692) < 1e-9 and abs(u0 * u0 - 0.295597742522) < 1e-9, "root values")
need(Fraction(12, 5) == Fraction(3 * 4, 5), "3-4-5 altitude 12/5")
# rational legs: k=j^2; need (j^4+1)(j^4+1+-2j) both squares
hits = [j for j in range(1, 2001) if is_sq((j ** 4 + 1) * (j ** 4 + 1 + 2 * j)) and is_sq((j ** 4 + 1) * (j ** 4 + 1 - 2 * j))]
print("j<=2000 with both (j^4+1)(j^4+1+-2j) square:", hits)
need(hits == [], "no rational-leg member j<=2000")
for j in range(2, 2001, 2):
    need(gcd(j ** 4 + 1, 2 * j) == 1, "j even gcd 1")
for j in range(1, 2001, 2):
    need(gcd(j ** 4 + 1, 2 * j) == 2, "j odd gcd 2")
# j odd: j^4+1 = 2 w^2 -> Pell x^2-2y^2=-1 numerators; independent generation via (1+sqrt2)^(2k+1)
x, y = 1, 1
sq = []
for idx in range(400):
    if is_sq(x):
        sq.append((idx, x))
    need(x * x - 2 * y * y == -1, "pell -1")
    x, y = 3 * x + 4 * y, 2 * x + 3 * y
print("squares among x with x^2-2y^2=-1, first 400:", sq)
need(sq == [(0, 1)], "only x=1")
# j=1 gives legs sqrt2,sqrt2: (j^4+1+2j)/2 = 2 is NOT a square -> j=1 is NOT a rational-leg member
need(not is_sq((1 + 1 + 2) // 2), "j=1 has irrational legs")
print("j=1: (j^4+1+2j)/2 = 2 not a square -> legs sqrt2: NOT a rational-leg member (explorer's 'only j=1 occurs' is about the Pell equation, not about a member)")
# Fermat closure: j^4+1=2w^2  =>  w^4 - j^4 = ((j^4-1)/2)^2 ; check identity on Pell pairs
x, y = 1, 1
for _ in range(50):
    # x = j^2 would be needed; identity holds for any x odd with x^2+1=2y^2:  y^4 - x^2... general: (x^2+1)^2-(x^2-1)^2=4x^2
    need((x * x + 1) ** 2 - (x * x - 1) ** 2 == 4 * x * x, "identity")
    need(4 * y ** 4 - (x * x - 1) ** 2 == 4 * x * x, "identity with 2y^2=x^2+1")
    x, y = 3 * x + 4 * y, 2 * x + 3 * y
print("Fermat reduction: j^4+1=2w^2 => w^4 - j^4 = ((j^4-1)/2)^2, nonzero square for j>1,")
print("contradicting Fermat's theorem that X^4-Y^4=Z^2 has no solution with XYZ!=0 -> j=1 only.")
print("Hence reading (b1) has NO rational-leg member at all (classical, not Ljunggren).")

# --------------------------------------------------------------------------
sec("A3. PS6: (k^2-1,2k,k^2+1) and THM-3335 family")
for k in range(2, 200):
    a, b, c = k * k - 1, 2 * k, k * k + 1
    need(a * a + b * b == c * c, "pyth")
    need((gcd(a, b) == 1) == (k % 2 == 0), "primitive iff k even")
    need(Fraction(a * b, c) == Fraction(2 * k * (k * k - 1), k * k + 1), "altitude")
    need(Fraction(b, c + a) == Fraction(1, k), "tan(phi/2)=1/k")
print("(k^2-1,2k,k^2+1): primitive iff k even, k=2..199; tan(phi/2)=1/k; k=2 is (3,4,5)")
fam = [(3, 4, 5), (17, 144, 145), (99, 4900, 4901), (577, 166464, 166465)]
for a, b, c in fam:
    need(a * a + b * b == c * c and gcd(a, b) == 1 and is_sq(b) and isqrt(b) % 2 == 0 and c == b + 1, "THM-3335 member")
    need(2 not in (Fraction(a * b, c), a, b, c), "sqrt(4)=2 matches nothing in 3-4-5")
print("THM-3335 members check; they are k=(2,12,70,408) rows of THM-3335 table (20)")

# --------------------------------------------------------------------------
sec("A4. PS7/PS8: Thales normalization and the chart")
for a, b, c in P:
    a1, b1 = min(a, b), max(a, b)
    e, d, l = Fraction(a1 * a1, c * c), Fraction(b1 * b1, c * c), Fraction(a1 * b1, c * c)
    need(e + d == 1 and l * l == e * d and e / d == Fraction(a1 * a1, b1 * b1), "chart identities")
    # vertex: foot of altitude at signed distance from centre; height l
    vx, vy = (d - e) / 2, l
    need(vx * vx + vy * vy == Fraction(1, 4), "vertex on diameter circle")
    # central angle 2theta: cos2theta = 1-2sin^2 = (b1^2-a1^2)/c^2 ; vertex/(1/2) = (cos2theta, sin2theta)
    need(2 * vx == 1 - 2 * e and 2 * vy == Fraction(2 * a1 * b1, c * c), "central angle 2theta")
    need(not (vx * vx + vy * vy == 1), "vertex NOT on radius-1 circle")
    # e/d and l depend only on tan theta
    tt = Fraction(a1, b1)
    need(e / d == tt * tt and l == tt / (1 + tt * tt), "functions of tan theta")
print("Thales: vertex ((d-e)/2,l) on circle radius 1/2, = (cos2theta, sin2theta)/2; on all 158 PPTs")
# monotonicity: numerical sample
import math
prev = (-1, -1)
for i in range(1, 1000):
    th = i * (math.pi / 4) / 1000
    cur = (math.tan(th) ** 2, math.sin(2 * th) / 2)
    need(cur[0] > prev[0] and cur[1] > prev[1], "monotone")
    prev = cur
print("(tan^2 theta, sin2theta/2) strictly increasing on a 1000-point grid of (0,pi/4); endpoint (1,1/2)")
need(Fraction(1, 4) == Fraction(1, 2) ** 2 and Fraction(2, 5) == Fraction(1 * 2, 5), "legs 1,2 witness")

# --------------------------------------------------------------------------
sec("A5. PS9: near-isosceles census by direct a^2+(a+1)^2=c^2, recursions, convergents")
ni = []
a = 1
while 2 * a * a + 2 * a + 1 <= 10 ** 12:
    c2 = a * a + (a + 1) ** 2
    if c2 > 10 ** 12:
        break
    c = isqrt(c2)
    if c * c == c2:
        ni.append((a, a + 1, c))
    a += 1
    if a > 710000:
        break
ni_c = [c for _, _, c in ni if c <= 10 ** 6]
print("near-isosceles hypotenuses c<=10^6 (direct):", ni_c)
need(ni_c == [5, 29, 169, 985, 5741, 33461, 195025], "census")
for i in range(2, len(ni_c)):
    need(ni_c[i] == 6 * ni_c[i - 1] - ni_c[i - 2], "c recursion")
shorts = [x for x, _, c in ni if c <= 10 ** 6]
for i in range(2, len(shorts)):
    need(shorts[i] == 6 * shorts[i - 1] - shorts[i - 2] + 2, "short-leg recursion")
# Euclid pairs recovered from (s,t)
pairs = []
for x, y, c in ni[:8]:
    odd, even = (x, y) if x % 2 else (y, x)
    s, t = isqrt(c + even), isqrt(c - even)
    m, n = (s + t) // 2, (s - t) // 2
    pairs.append((m, n, s, t, odd, even))
    need((m - n) ** 2 - 2 * n * n in (1, -1) and s * s - 2 * m * m in (1, -1), "Pell forms")
for i in range(1, len(pairs)):
    need(pairs[i][0] == 2 * pairs[i - 1][0] + pairs[i - 1][1] and pairs[i][1] == pairs[i - 1][0], "(m,n)->(2m+n,m)")
print("(m,n) chain:", [(p[0], p[1]) for p in pairs])
# convergents of sqrt2-1 = [0;2,2,2,...]
h0, k0, h1, k1 = 1, 0, 0, 1   # p_{-2}/q_{-2}=1/0? standard: p_{-1}=1,q_{-1}=0; p_0=0,q_0=1
conv = []
p_prev, q_prev, p_cur, q_cur = 1, 0, 0, 1
for _ in range(10):
    p_prev, q_prev, p_cur, q_cur = p_cur, q_cur, 2 * p_cur + p_prev, 2 * q_cur + q_prev
    conv.append(Fraction(p_cur, q_cur))
print("convergents of sqrt2-1:", [str(f) for f in conv[:7]])
ts = [Fraction(p[3], p[2]) for p in pairs]
nm = [Fraction(p[1], p[0]) for p in pairs]
print("explorer's t/s sequence:", [str(f) for f in ts[:6]])
print("n/m sequence:          ", [str(f) for f in nm[:6]])
need(all(f in conv for f in nm), "n/m ARE convergents")
need(sum(1 for f in ts if f in conv) == 0, "t/s are NOT convergents")
# t/s are mediants of consecutive convergents (semiconvergents)
for i in range(1, 6):
    med = Fraction(conv[i - 1].numerator + conv[i].numerator, conv[i - 1].denominator + conv[i].denominator)
    need(med == ts[i], f"t/s[{i}] is the mediant of convergents {i-1},{i}")
need(ts[0] == Fraction(1, 3) == Fraction(0 + 1, 1 + 2), "1/3 mediant of 0/1 and 1/2")
print("REFUTED detail: t/s = 1/3,3/7,7/17,17/41,... are the intermediate fractions (mediants of")
print("consecutive convergents) of sqrt2-1, not its convergents; n/m = 1/2,2/5,5/12,12/29,... are.")
# theta = smaller angle: tan(theta/2) = short/(c+long)
print("tan(theta/2) for theta = angle opposite the SHORTER leg (note's Theorem 3 convention):")
for x, y, c in ni[:5]:
    print(f"   ({x},{y},{c}): tan(theta/2) = {Fraction(x, c + y)}   (t/s would be {Fraction(isqrt(c-(y if y%2==0 else x)), isqrt(c+(y if y%2==0 else x)))})")
need(Fraction(20, 29 + 21) == Fraction(2, 5), "(20,21,29): tan(theta/2)=2/5 not 3/7")
need(Fraction(696, 985 + 697) == Fraction(12, 29), "(696,697,985): tan(theta/2)=12/29 not 17/41")

# --------------------------------------------------------------------------
sec("A6. PS10: the 29 test, both assignments, and the agreement locus for all (m,n)")


def orbit(a, b, steps=3):
    D = a * a + 2 * b * b - a * b
    y = [Fraction(-(a + b))]
    for _ in range(steps):
        y.append((y[-1] ** 2 - D) / b)
    return D, y


D, y = orbit(3, 4)
need(D == 29 and y == [-7, 5, -1, -7], "3-4-5 cycle")
for (a, b) in ((20, 21), (21, 20)):
    D, y = orbit(a, b)
    print(f"   (a,b)=({a},{b}): D={D}, orbit {[str(v) for v in y]}")
    need(y[3] != y[0] and y[1] != 29, "no cycle for (20,21,29)")
need(orbit(20, 21)[0] == 862 and orbit(20, 21)[1][3] == Fraction(54139, 9261), "explorer's numbers")


def Dmn(m, n):
    a, b = m * m - n * n, 2 * m * n
    return a * a + 2 * b * b - a * b


agree = [(m, n) for m in range(2, 301) for n in range(1, m) if Dmn(m, n) == (2 * m + n) ** 2 + m * m]
need(agree == [(2, 1)], "agreement locus")
# all (m,n): D = c^2 + b(b-a) >= c^2 - ab >= c^2/2 ; c_next < 10 m^2 <= 10 c ; c^2/2 > 10c iff c > 20
for m in range(2, 6):
    for n in range(1, m):
        if (m, n) != (2, 1):
            need(Dmn(m, n) != (2 * m + n) ** 2 + m * m, "small locus")
print("agreement locus for 300>=m>n>=1 is [(2,1)]; and PROVABLE for all (m,n): D>=c^2-ab>=c^2/2>10c>c_next once c>20")
need(all(2 * m * n * (2 * m * n - m * m + n * n) != 4 for m in range(2, 200) for n in range(1, m) if (m, n) != (2, 1)), "b(b-a)=4 only at (2,1)")

# --------------------------------------------------------------------------
sec("A7. PS11: doubling forest, square hypotenuses c<=10^4, 2^(omega-1)")
P4 = brute_ppts(10 ** 4)
print("brute-force #PPT(c<=10^4) =", len(P4))
need(len(P4) == 1593, "1593")
byc = {}
for a, b, c in P4:
    byc.setdefault(c, []).append((a, b))
sqh = [(a, b, c) for a, b, c in P4 if is_sq(c)]
print("square-hypotenuse PPTs c<=10^4:", len(sqh))
need(len(sqh) == 16, "16")
small = [(a, b, c) for a, b, c in P4 if c <= 100]
need(len(small) == 16, "#PPT(c<=100)=16")
doubles = {}
for a, b, c in small:
    A, B, C = 2 * a * b, abs(b * b - a * a), c * c
    need(A * A + B * B == C * C and gcd(A, B) == 1, "double is PPT")
    key = (B, A, C)  # odd, even, c^2
    need(key in [(x, y, z) for x, y, z in sqh], "double present among square-hyp PPTs")
    need(key not in doubles, "injective")
    doubles[key] = (a, b, c)
need(len(doubles) == 16, "bijection onto square hypotenuse PPTs")
for c in byc:
    if c <= 100:
        need(len(byc[c]) == 2 ** (omega(c) - 1) == len(byc[c * c]), "2^(omega-1)")
print("doubling PPT(c<=100) -> square-hyp PPT(c<=10^4) is a bijection (16 each); |PPT(c)|=|PPT(c^2)|=2^(omega(c)-1)")
need(doubles[(7, 24, 25)] == (3, 4, 5) and doubles[(119, 120, 169)] == (5, 12, 13) and doubles[(41, 840, 841)] == (21, 20, 29), "listed pairs")
chain = (3, 4, 5)
for expect in [(24, 7, 25), (336, 527, 625), (354144, 164833, 390625)]:
    a, b, c = chain
    chain = (2 * a * b, abs(b * b - a * a), c * c)
    need(chain == expect, "chain")
print("chain (3,4,5)->(24,7,25)->(336,527,625)->(354144,164833,390625) confirmed")
# conjugacy 2(2x^2-1) = (2x)^2-2
for xq in [Fraction(p, 11) for p in range(-20, 21)]:
    need(2 * (2 * xq * xq - 1) == (2 * xq) ** 2 - 2, "conjugacy")

# --------------------------------------------------------------------------
sec("A8. PS12: y^2-2 preperiodic, 3-cycles, ord_9(2), row transversal")
graph = {y: y * y - 2 for y in range(-2, 3)}
need(graph == {2: 2, -2: 2, 0: -2, 1: -1, -1: -1}, "graph")
for y0 in range(-50, 51):
    if abs(y0) >= 3:
        need(abs(y0 * y0 - 2) > abs(y0), "escape")
# rational non-integers: denominator squares
for p in range(-30, 31):
    for q in range(2, 12):
        if gcd(p, q) == 1:
            r = Fraction(p, q) ** 2 - 2
            need(r.denominator == q * q, "denominator squares")
# angles
for y0, ang in ((2, 0), (1, math.pi / 3), (0, math.pi / 2), (-1, 2 * math.pi / 3), (-2, math.pi)):
    need(abs(2 * math.cos(2 * ang) - graph[y0]) < 1e-12, "angle doubling")
for N in (9, 7):
    y0 = 2 * math.cos(2 * math.pi / N)
    y1 = y0 * y0 - 2
    y2 = y1 * y1 - 2
    y3 = y2 * y2 - 2
    need(abs(y3 - y0) < 1e-9 and abs(y1 - y0) > 1e-3, f"period 3 at N={N}")
    need(abs(y0 ** 3 - 3 * y0 + 1) < 1e-9 if N == 9 else abs(y0 ** 3 + y0 ** 2 - 2 * y0 - 1) < 1e-9, "min poly")
ordn = min(e for e in range(1, 20) if pow(2, e, 9) == 1)
need(ordn == 6, "ord_9(2)=6")
orb = [pow(2, e, 9) for e in range(6)]
need(sorted(orb) == [1, 2, 4, 5, 7, 8], "orbit")
need(all((-r) % 9 not in {2, 5, 8} for r in {2, 5, 8}), "transversal")
print("PrePer graph, escape, 3-cycles at N=7,9, ord_9(2)=6, {2,5,8} transversal all confirmed")
# NOTE: cycle at N=7 exists because ord_7(2)=3 directly, no -1 fold needed
need(min(e for e in range(1, 20) if pow(2, e, 7) == 1) == 3, "ord_7(2)=3")

# --------------------------------------------------------------------------
sec("A9. PS13: integer graphs of x^2+c")
for c in (1, 2):
    for x0 in range(-100, 101):
        need(x0 * x0 + c > x0, "x^2+c>x")
need({x: x * x for x in (-1, 0, 1)} == {-1: 1, 0: 0, 1: 1}, "x^2")
need({x: x * x - 1 for x in (-1, 0, 1)} == {-1: 0, 0: -1, 1: 0}, "x^2-1")
need({x: x * x - 2 for x in (-2, -1, 0, 1, 2)} == {-2: 2, -1: -1, 0: -2, 1: -1, 2: 2}, "x^2-2")
print("x^2, x^2-1, x^2-2 integer graphs confirmed; x^2+1, x^2+2 have no integer cycles (x^2+c>x)")

print("\nAUDIT RECOMPUTE: ALL CHECKS PASSED")
