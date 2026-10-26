# Fixed-twist rational lifts of hourglass boundary points over Q.
# A base point y with vanishing set Z(y) lifts to X_d(Q) iff the nonzero d_i L_i(y) share one squareclass.
from fractions import Fraction
from math import gcd, isqrt
import itertools
forms = {'E':(1,0,0),'E+x':(1,1,0),'E-x':(1,-1,0),'E+z':(1,0,1),'E-z':(1,0,-1),'E+x+z':(1,1,1),'E-x-z':(1,-1,-1)}
names = list(forms)
def sqclass(q):
    # squarefree part with sign of a nonzero rational
    q = Fraction(q); n, d = q.numerator, q.denominator
    n *= d; sgn = -1 if n < 0 else 1; n = abs(n)
    f = 1; p = 2
    while p*p <= n:
        while n % (p*p) == 0: n //= p*p
        if n % p == 0: f *= p; n //= p
        p += 1
    return sgn * f * n
def lifts(d, y):
    vals = [Fraction(a*y[0]+b*y[1]+c*y[2]) for (a,b,c) in forms.values()]
    nz = [sqclass(di*vi) for di, vi in zip(d, vals) if vi != 0]
    return len(set(nz)) == 1, [n for n,vi in zip(names, vals) if vi == 0]
def central_line_points(d, H=60):
    # points (0:x:z) with x,z coprime integers of height <= H that lift
    pts = []
    for x in range(-H, H+1):
        for z in range(-H, H+1):
            if (x, z) == (0, 0) or gcd(abs(x), abs(z)) != 1: continue
            if z < 0 or (z == 0 and x < 0): continue
            ok, van = lifts(d, (0, x, z))
            if ok: pts.append((x, z, tuple(van)))
    return pts
untwisted = (1,1,1,1,1,1,1)
pyth = (1,1,-1,1,-1,1,-1)
print("untwisted: liftable points on E=0 up to height 60:", len(central_line_points(untwisted)))
P = central_line_points(pyth)
print("twist (1,1,-1,1,-1,1,-1): liftable points on E=0 up to height 60:", len(P), "examples:", P[:8])
# all should come from Pythagorean data: x z square and x(x+z) square (up to the twist signs) -> x=a^2, z=b^2, a^2+b^2=c^2 (or sign variants)
def is_sq(n): return n >= 0 and isqrt(n)**2 == n
print("all liftable central points have x*z and x*(x+z) perfect squares (times sign bookkeeping):",
      all(is_sq(abs(x*z)) and is_sq(abs(x*(x+z))) for x,z,_ in P))
# other strata for the untwisted surface: triple points and double points
tri = [(0,0,1),(0,1,0),(0,1,-1)]
print("untwisted triple points lift?", [lifts(untwisted, y)[0] for y in tri])
# a twist making the triple point (0:0:1) lift: need d4*1, d5*(-1), d6*1, d7*(-1) in one class: d = (1,1,1,1,-1,1,-1)
print("twist (1,1,1,1,-1,1,-1): triple point (0:0:1) lifts:", lifts((1,1,1,1,-1,1,-1),(0,0,1)))
# generic point lift locus for the Pythagorean twist: count points of P^2(Q) of height<=12 lifting
cnt = 0; tot = 0
for E in range(0,13):
    for x in range(-12,13):
        for z in range(-12,13):
            if (E,x,z)==(0,0,0) or gcd(gcd(E,abs(x)),abs(z))!=1: continue
            tot += 1
            if lifts(pyth,(E,x,z))[0]: cnt += 1
print("Pythagorean twist: rational points of P^2 (height<=12) lifting to X_d:", cnt, "of", tot)
