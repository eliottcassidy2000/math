"""Audit THM-4824 (F) and the results note's numerics (G), with methods independent of counts.py."""
from math import isqrt, gcd
import sympy as sp
from sympy import divisor_sigma, isprime, primerange

print("(F) zero-sum quartic identity, symbolic")
a, b = sp.symbols('a b'); c = -a - b
lhs = sp.expand(a**4 + b**4 + c**4)
print("   a^4+b^4+c^4 - 2(a^2+ab+b^2)^2 =", sp.expand(lhs - 2*(a*a + a*b + b*b)**2),
      "; - (a^2+b^2+c^2)^2/2 =", sp.expand(lhs - (a*a + b*b + c*c)**2 / 2))
Qw = sp.expand(a*a + b*b + c*c - a*b - b*c - c*a)
print("   Q(a+bw+cw^2) on a+b+c=0 equals 3(a^2+ab+b^2):", sp.expand(Qw - 3*(a*a + a*b + b*b)) == 0,
      "; (2/9)Q^2 - lhs =", sp.expand(sp.Rational(2, 9)*Qw**2 - lhs))
# Q(x + yw) = x^2 - xy + y^2 with a + bw + cw^2 = (a-c) + (b-c)w
x_, y_ = a - c, b - c
print("   norm of (a-c)+(b-c)w equals a^2+b^2+c^2-ab-bc-ca (general a,b,c):",
      sp.expand((lambda A, B, C: (A - C)**2 - (A - C)*(B - C) + (B - C)**2 - (A*A + B*B + C*C - A*B - B*C - C*A))(*sp.symbols('A B C'))) == 0)
eis_norms = sorted({p*p - p*q + q*q for p in range(-20, 21) for q in range(-20, 21)})
print("   2N^2 (N an Eisenstein norm) below 51:", [2*n*n for n in eis_norms if 0 < 2*n*n < 51],
      "; a^2+ab+b^2 takes every Eisenstein norm (|a|,|b|<=20 sample):",
      {p*p + p*q + q*q for p in range(-20, 21) for q in range(-20, 21)} >= {n for n in eis_norms if n <= 300})

print("(G1) sums of four fourth powers (zeros allowed) = n needs at most 4 positive fourth powers; DP on min count")
N = 10**4
fp = [k**4 for k in range(1, 11) if k**4 <= N]
INFc = 99
mc = [0] + [INFc]*N
for n in range(1, N + 1):
    mc[n] = 1 + min(mc[n - f] for f in fp if f <= n)
rep = [mc[n] <= 4 for n in range(N + 1)]
for X in (5, 10, 25, 50, 100, 200, 500, 1000, 10**4):
    r = sum(rep[1:X + 1]); print(f"   N = {X}: representable {r}, not {X - r}")
print("   representable in 1..50:", [n for n in range(1, 51) if rep[n]])
gaps = []; n = 1
while n <= 80:
    if not rep[n]:
        s0 = n
        while n <= 80 and not rep[n]: n += 1
        gaps.append((s0, n - 1, n - s0))
    else: n += 1
print("   all gaps in 1..80 (start, end, length):", gaps, "; 81 representable:", rep[81])
print("   set below 81 equals {16a+b : a+b <= 4}:", {n for n in range(81) if rep[n]} == {16*a_ + b_ for a_ in range(5) for b_ in range(5 - a_)},
      "; size", len({16*a_ + b_ for a_ in range(5) for b_ in range(5 - a_)}))
print("   residues mod 16 of representable n <= 10^4:", sorted({n % 16 for n in range(N + 1) if rep[n]}), "; mc(79) =", mc[79],
      "; n <= 10^4 needing 19:", [n for n in range(N + 1) if mc[n] == 19])

print("(G2) square-later: n = 2^a u, a >= 1, u odd, qualifying iff u > 2^a")
def v2(n):
    a_ = (n & -n).bit_length() - 1; return a_, n >> a_
for X in (5, 10, 25, 50, 100, 200, 1000, 10**4):
    ev = list(range(2, X + 1, 2)); qq = sum(1 for n in ev if v2(n)[1] > 2**v2(n)[0])
    print(f"   X = {X}: evens {len(ev)}, qualifying {qq}, failures {len(ev) - qq}")
print("   failures <= 100:", [n for n in range(2, 101, 2) if v2(n)[1] <= 2**v2(n)[0]])

print("(G3) lattice-ball counts (closed-form row sums)")
def ball_Z2(t): return sum(2*isqrt(t - x*x) + 1 for x in range(-isqrt(t), isqrt(t) + 1))
def ball_Zw(t):  # a^2 - ab + b^2 <= t  <=>  (2a - b)^2 + 3b^2 <= 4t
    tot = 0
    for b_ in range(-isqrt(4*t // 3) - 1, isqrt(4*t // 3) + 2):
        rem = 4*t - 3*b_*b_
        if rem < 0: continue
        r = isqrt(rem)          # |2a - b| <= r
        lo = -(( r + b_) // 2) if True else 0
        # count integers a with -r <= 2a - b <= r
        lo = -((r - b_) // 2); hi = (r + b_) // 2
        lo = (-r + b_ + 1) // 2 if (-r + b_) % 2 else (-r + b_) // 2
        hi = (r + b_) // 2
        tot += max(0, hi - lo + 1)
    return tot
brute_w = lambda t: sum(1 for p in range(-2*isqrt(t)-2, 2*isqrt(t)+3) for q in range(-2*isqrt(t)-2, 2*isqrt(t)+3) if p*p - p*q + q*q <= t)
assert all(ball_Zw(t) == brute_w(t) for t in range(0, 60))
S2 = sorted({ball_Z2(t) for t in range(0, 400)}); Sw = sorted({ball_Zw(t) for t in range(0, 400)})
common = sorted(set(S2) & set(Sw))
print("   common values <= 400:", [v for v in common if v <= 400], "; all = 1 mod 12:", all(v % 12 == 1 for v in common))
print("   #Z2(x^2+y^2<=4) =", ball_Z2(4), " #Z2(<=10) =", ball_Z2(10), " #Zw(Q<=3) =", ball_Zw(3), " #Zw(Q<=9) =", ball_Zw(9),
      "; S_2, S_3 =", 6*2*1 + 1, 6*3*2 + 1, "; H_4 =", 3*4*3 + 1)

print("(G4) primes = 1 mod 12 below 200:", [p for p in primerange(2, 200) if p % 12 == 1])
print("     13, 37 as x^2+y^2 and x^2+3y^2:", {p: ([(x, y) for x in range(1, 8) for y in range(x, 8) if x*x + y*y == p],
      [(x, y) for x in range(0, 8) for y in range(1, 5) if x*x + 3*y*y == p]) for p in (13, 37)})

print("(G5) triangles")
leg12 = sorted({(min(b_, 12), max(b_, 12), cc) for b_ in range(1, 200) for cc in range(b_, 400) if b_*b_ + 144 == cc*cc and gcd(gcd(b_, 12), cc) == 1})
print("   primitive Pythagorean triangles with a leg 12:", leg12)
side7_short = sorted({(7, b_, cc) for b_ in range(1, 500) for cc in range(1, 600) if 49 + 7*b_ + b_*b_ == cc*cc and gcd(gcd(7, b_), cc) == 1})
side7_long = sorted({(x, y, 7) for x in range(1, 7) for y in range(x, 7) if x*x + x*y + y*y == 49 and gcd(gcd(x, y), 7) == 1})
print("   primitive 120-degree triangles with 7 adjacent to the 120-degree angle:", side7_short,
      "; with 7 as the LONGEST side:", side7_long)
def pyth_legs(cc): return {s for x in range(1, cc) for y in range(x, cc) if x*x + y*y == cc*cc and gcd(gcd(x, y), cc) == 1 for s in (x, y)}
def eis_sides(cc): return {s for x in range(1, cc) for y in range(x, cc) if x*x + x*y + y*y == cc*cc and gcd(gcd(x, y), cc) == 1 for s in (x, y)}
ps = [p for p in primerange(13, 200) if p % 12 == 1]
L = {p: (pyth_legs(p), eis_sides(p)) for p in ps}
both = [(p, q) for i, p in enumerate(ps) for q in ps[i+1:] if (L[p][0] & L[q][0]) and (L[p][1] & L[q][1])]
pyonly = [(p, q, sorted(L[p][0] & L[q][0])) for i, p in enumerate(ps) for q in ps[i+1:] if (L[p][0] & L[q][0]) and not (L[p][1] & L[q][1])]
eonly = [(p, q, sorted(L[p][1] & L[q][1])) for i, p in enumerate(ps) for q in ps[i+1:] if (L[p][1] & L[q][1]) and not (L[p][0] & L[q][0])]
print("   pairs sharing BOTH a Pythagorean leg and a 120-degree side:", both, "; Pythagorean only:", pyonly, "; 120 only:", eonly)
print("   61:", sorted(L[61][0]), sorted(L[61][1]), " 109:", sorted(L[109][0]), sorted(L[109][1]))
print("   13^2 = 11^2 + 3*4^2:", 169 == 121 + 48, "; 37^2 = 13^2 + 3*20^2:", 1369 == 169 + 1200,
      "; reps x^2+3y^2 of 169:", [(x, y) for x in range(0, 14) for y in range(0, 8) if x*x + 3*y*y == 169],
      " of 1369:", [(x, y) for x in range(0, 38) for y in range(0, 22) if x*x + 3*y*y == 1369])

print("(G6) norm / aliquot table")
gauss = lambda n: any(isqrt(n - x*x)**2 == n - x*x for x in range(0, isqrt(n) + 1))
eisn = lambda n: any((lambda r: isqrt(r)**2 == r)(4*n - 3*y*y) for y in range(0, isqrt(4*n // 3) + 1))   # n = Q iff 4n = s^2 + 3y^2
def sql(n):
    if n % 2: return 'odd'
    a_, u = v2(n); return 'qualifies' if u > 2**a_ else 'fails'
s = lambda n: int(divisor_sigma(n)) - n
for n in (5, 10, 25, 50, 100, 11, 12, 13):
    print(f"   {n:4d}: Gaussian {gauss(n)}, Eisenstein {eisn(n)}, square-later {sql(n)}, s(n) = {s(n)}")
print("   s(25) = 6 is perfect:", s(6) == 6)
fib = {m: [n for n in range(2, 20000) if s(n) == m] for m in (6, 7, 8, 11, 12, 13)}
print("   aliquot fibres (n < 20000; complete since composite n has s(n) > sqrt n):", fib)
two = [n for n in range(1, 200) if len({(x, y) for x in range(1, isqrt(n) + 1) for y in range(x, isqrt(n) + 1) if x*x + y*y == n}) >= 2]
two_distinct = [n for n in range(1, 200) if len({(x, y) for x in range(1, isqrt(n) + 1) for y in range(x + 1, isqrt(n) + 1) if x*x + y*y == n}) >= 2]
print("   two representations as a sum of two positive squares:", two[:4], "; with distinct squares:", two_distinct[:3])
