"""Fast exact test over Q: at random integer circle triples compute the four Apollonius discriminants Delta_e (one per sign
class) as rationals, and test which products are rational squares.  If Delta_1 Delta_2 Delta_3 Delta_4 is a square at every
random point, the product is a square as a rational function (the monodromy group is contained in the parity subgroup of (Z/2)^4)."""
from fractions import Fraction as F
from math import isqrt
import random
random.seed(11)
def is_sq(fr):
    if fr == 0: return True
    if fr < 0: return False
    n, d = fr.numerator, fr.denominator
    return isqrt(n)**2 == n and isqrt(d)**2 == d
def discs(circ):
    out = []
    (a1,b1,r1),(a2,b2,r2),(a3,b3,r3) = circ
    for e in [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1)]:
        rows = []
        for (a,b,rr), ei in ((circ[1], e[1]), (circ[2], e[2])):
            rows.append((F(2*(a-a1)), F(2*(b-b1)), F(2*(ei*rr - e[0]*r1)), F(a1*a1 + b1*b1 - r1*r1 - a*a - b*b + rr*rr)))
        (m11, m12, c1, d1), (m21, m22, c2, d2) = rows
        det = m11*m22 - m12*m21
        if det == 0: return None
        xr = (-(c1*m22) + c2*m12)/det; x0 = (-(d1*m22) + d2*m12)/det
        yr = (-(c2*m11) + c1*m21)/det; y0 = (-(d2*m11) + d1*m21)/det
        A = xr*xr + yr*yr - 1; B = 2*xr*(x0 - a1) + 2*yr*(y0 - b1) - 2*e[0]*r1; Cc = (x0 - a1)**2 + (y0 - b1)**2 - r1*r1
        out.append(B*B - 4*A*Cc)
    return out
import itertools
stats = {"all4": 0, "pairs": {pr: 0 for pr in itertools.combinations(range(4), 2)}, "triples": {t: 0 for t in itertools.combinations(range(4), 3)}, "n": 0}
examples = []
while stats["n"] < 400:
    circ = [(random.randint(-30, 30), random.randint(-30, 30), random.randint(1, 30)) for _ in range(3)]
    D = discs(circ)
    if D is None or any(d == 0 for d in D): continue
    stats["n"] += 1
    if is_sq(D[0]*D[1]*D[2]*D[3]): stats["all4"] += 1
    for pr in stats["pairs"]:
        if is_sq(D[pr[0]]*D[pr[1]]): stats["pairs"][pr] += 1
    for t in stats["triples"]:
        if is_sq(D[t[0]]*D[t[1]]*D[t[2]]): stats["triples"][t] += 1
    if len(examples) < 3: examples.append((circ, D))
print("random integer triples:", stats["n"])
print("product of all four discriminants a rational square:", stats["all4"], "of", stats["n"])
print("pairwise products squares:", stats["pairs"])
print("triple products squares:", stats["triples"])
for circ, D in examples:
    print("example", circ, "discriminants", [str(d) for d in D])
