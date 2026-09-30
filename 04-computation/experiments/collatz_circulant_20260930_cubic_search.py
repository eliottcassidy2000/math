"""x^3+2x+y^3+2y+z^3+2z = xyz+1: small solutions, structure (e3 determined by e1,e2), the x=y family, mod obstructions."""
import math, sys
from sympy import integer_nthroot
def f(t): return t**3 + 2*t
def int_roots_cubic(p, q):
    """integer roots of z^3 + p z + q = 0 (p,q ints): rational-root test via divisors of q is slow; use float roots + check."""
    import numpy as np
    if q == 0: cands = {0}
    else: cands = set()
    r = np.roots([1, 0, p, q])
    for c in r:
        if abs(c.imag) < 1e-6*max(1,abs(c.real)):
            for z in (int(round(c.real)),):
                for d in (-1,0,1):
                    zz = z+d
                    if zz**3 + p*zz + q == 0: cands.add(zz)
    return sorted(cands)
B = int(sys.argv[1]) if len(sys.argv)>1 else 400
sols = set()
for x in range(-B, B+1):
    for y in range(x, B+1):
        p = 2 - x*y; q = f(x)+f(y)-1
        for z in int_roots_cubic(p, q):
            sols.add(tuple(sorted((x,y,z))))
print(f"solutions with |x|,|y| <= {B} (z any):", len(sols))
for s in sorted(sols, key=lambda t: max(abs(v) for v in t))[:60]: print("  ", s, " e1 =", sum(s))
# x = y family
fam = []
for x in range(-10**6, 10**6+1):
    p = 2 - x*x; q = 2*f(x) - 1
    for z in int_roots_cubic(p, q): fam.append((x,x,z))
print("x=y family, |x| <= 1e6:", fam[:20], "count", len(fam))
