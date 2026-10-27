"""Owner's 2026-10-27 grid / 13-37 directive (mac-mini-2026-10-27-slab): exact counts behind the observations.
(1) sums of four fourth powers (zeros allowed): representable n <= N, the gap structure below 81, residues mod 16;
(2) session69's square-later integers (n = 2^a u, a >= 1, u odd, u > 2^a) among the evens <= X;
(3) lattice-point counts of the square grid Z^2 and the Eisenstein grid Z[omega] in norm balls, and their common values;
(4) the zero-sum quartic identity a + b + c = 0  =>  a^4 + b^4 + c^4 = 2(a^2 + ab + b^2)^2 = (2/9) Q(a + b w + c w^2)^2."""
from itertools import combinations_with_replacement
from math import isqrt
def fourth_power_sums(N):
    p = [x**4 for x in range(0, int(N**0.25) + 2) if x**4 <= N]
    S = set()
    for c in combinations_with_replacement(p, 4):
        s = sum(c)
        if s <= N: S.add(s)
    return S
S = fourth_power_sums(10**4)
print("(1) sums of four fourth powers")
for N in (5, 10, 25, 50, 100, 200, 500, 1000, 10**4):
    rep = sum(1 for n in S if 1 <= n <= N); print(f"   N={N}: representable positives {rep}, non-representable {N - rep}")
rep50 = sorted(n for n in S if 1 <= n <= 50); print("   representable in 1..50:", rep50)
gaps = []; prev = 0
for n in sorted(x for x in S if x <= 80):
    if n - prev > 1: gaps.append((prev + 1, n - 1, n - prev - 1))
    prev = n
print("   gaps below 81 (start, end, length):", gaps)
print("   below 81 the representable set is {16a + b : a + b <= 4}:", sorted(x for x in S if x <= 80) == sorted({16*a + b for a in range(5) for b in range(5 - a)} - {x for x in range(81, 10**5)}))
print("   residues mod 16 of all sums <= 10^4:", sorted({n % 16 for n in S}))
print("(2) square-later integers (session69): n = 2^a u, a >= 1, u odd; qualifying iff u > 2^a")
def split2(n):
    a = 0
    while n % 2 == 0: n //= 2; a += 1
    return a, n
for X in (5, 10, 25, 50, 100, 200, 1000, 10**4):
    q = f = 0
    for n in range(2, X + 1, 2):
        a, u = split2(n)
        if u > 2**a: q += 1
        else: f += 1
    print(f"   X={X}: evens {X//2}, qualifying {q}, failures {f}")
print("   failures <= 100:", [n for n in range(2, 101, 2) if split2(n)[1] <= 2**split2(n)[0]])
print("(3) lattice-ball counts")
def ball_sq(t): return sum(1 for x in range(-isqrt(t), isqrt(t) + 1) for y in range(-isqrt(t), isqrt(t) + 1) if x*x + y*y <= t)
def ball_eis(t):
    R = isqrt(4*t) + 2
    return sum(1 for a in range(-R, R + 1) for b in range(-R, R + 1) if a*a - a*b + b*b <= t)
sq = sorted({ball_sq(t) for t in range(0, 200)}); eis = sorted({ball_eis(t) for t in range(0, 200)})
common = sorted(set(sq) & set(eis))
print("   Z^2 counts:", sq[:24]); print("   Z[w] counts:", eis[:24]); print("   common values:", [c for c in common if c <= 400])
print("   13 = #Z^2 in x^2+y^2 <= 4:", ball_sq(4), " 37 = #Z^2 in x^2+y^2 <= 10:", ball_sq(10))
print("   13 = #Z[w] in Q <= 3:", ball_eis(3), " 37 = #Z[w] in Q <= 9:", ball_eis(9))
star = [6*n*(n-1) + 1 for n in range(1, 8)]; hexc = [3*n*(n-1) + 1 for n in range(1, 8)]
print("   star numbers:", star, " centred hexagonal:", hexc)
def isprime(n): return n > 1 and all(n % d for d in range(2, isqrt(n) + 1))
print("   primes = 1 mod 12 below 200 (split in both Z[i] and Z[w]):", [p for p in range(2, 200) if isprime(p) and p % 12 == 1])
for p in (13, 37):
    r1 = [(x, y) for x in range(0, 7) for y in range(x, 7) if x*x + y*y == p]
    r3 = [(x, y) for x in range(0, 7) for y in range(0, 4) if x*x + 3*y*y == p]
    print(f"   {p} = x^2+y^2 at {r1}; = x^2+3y^2 at {r3}")
print("(4) zero-sum quartic identity")
ok = all(a**4 + b**4 + (a+b)**4 == 2*(a*a + a*b + b*b)**2 for a in range(-60, 61) for b in range(-60, 61))
def Qw(p, q, r): return p*p + q*q + r*r - p*q - q*r - r*p   # Q(p + q w + r w^2)
ok2 = all(9*(a**4 + b**4 + c**4) == 2*Qw(a, b, c)**2 for a in range(-40, 41) for b in range(-40, 41) for c in [-(a+b)])
print("   a^4 + b^4 + (a+b)^4 = 2(a^2+ab+b^2)^2 for |a|,|b| <= 60:", ok, "; a+b+c=0 => 9 sum a^4 = 2 Q(a+bw+cw^2)^2:", ok2)
print("   values 2(a^2+ab+b^2)^2 <= 50:", sorted({2*(a*a + a*b + b*b)**2 for a in range(-6, 7) for b in range(-6, 7) if 0 < 2*(a*a+a*b+b*b)**2 <= 50}))
print("(5) 13 and 37 in both grids")
from math import gcd
pyth12 = sorted({c for a in range(1, 200) for c in range(a + 1, 400) if a*a + 144 == c*c and gcd(gcd(a, 12), c) == 1})
print("   primitive Pythagorean triangles with a leg 12: hypotenuses", pyth12)
eis7 = sorted({c for b in range(1, 400) for c in range(b, 800) if 49 + 7*b + b*b == c*c and gcd(gcd(7, b), c) == 1})
print("   primitive 120-degree triangles (c^2 = a^2 + ab + b^2) with a side 7: longest sides", eis7)
print("   13^2 = 11^2 + 3*4^2:", 13**2 == 11**2 + 3*4**2, "; 37^2 = 13^2 + 3*20^2:", 37**2 == 13**2 + 3*20**2)
two = [n for n in range(1, 200) if sum(1 for a in range(1, isqrt(n) + 1) for b in range(a, isqrt(n) + 1) if a*a + b*b == n) >= 2]
print("   smallest integers with two representations as a sum of two positive squares:", two[:5])
for n in (5, 10, 25, 50, 100, 11, 12, 13, 37):
    gauss = any(a*a + b*b == n for a in range(0, isqrt(n) + 1) for b in range(0, isqrt(n) + 1))
    eis = any(a*a - a*b + b*b == n for a in range(-2*isqrt(n) - 2, 2*isqrt(n) + 3) for b in range(-2*isqrt(n) - 2, 2*isqrt(n) + 3))
    print(f"   {n}: Gaussian norm {gauss}, Eisenstein norm {eis}")
print("(6) do the next primes = 1 mod 12 share sides in both grids?")
def pyth_sides(c): return sorted({x for a in range(1, c) for b in range(a, c) if a*a + b*b == c*c and gcd(gcd(a, b), c) == 1 for x in (a, b)})
def eis120_sides(c): return sorted({x for a in range(1, c) for b in range(a, c) if a*a + a*b + b*b == c*c and gcd(gcd(a, b), c) == 1 for x in (a, b)})
ps = [p for p in range(13, 200) if isprime(p) and p % 12 == 1]
for i, p in enumerate(ps):
    for q in ps[i+1:]:
        sp = set(pyth_sides(p)) & set(pyth_sides(q)); se = set(eis120_sides(p)) & set(eis120_sides(q))
        if sp and se: print(f"   {p}, {q}: shared Pythagorean side(s) {sorted(sp)} AND shared 120-degree side(s) {sorted(se)}")
print("   61:", pyth_sides(61), eis120_sides(61), " 109:", pyth_sides(109), eis120_sides(109))
