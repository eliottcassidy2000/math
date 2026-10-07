# Hoffman lower bound for kappa(p) = chi(Cay(F_{p^2}, mu_{p+1})), p odd prime: eigenvalues are character sums over the norm-1 circle.
import cmath, math
from sympy import primerange, legendre_symbol
for p in primerange(3, 60):
    nu = next(a for a in range(2, p) if legendre_symbol(a, p) == -1)
    circ = [(a, b) for a in range(p) for b in range(p) if (a*a - nu*b*b) % p == 1]
    assert len(circ) == p + 1
    lmin = min(sum(cmath.exp(2j*math.pi*2*(c*a + nu*d*b)/p) for a, b in circ).real
               for c in range(p) for d in range(p) if (c, d) != (0, 0))
    print(f"p={p:2d}  degree={p+1:2d}  lambda_min={lmin:8.4f}  (Weil -2sqrt(p)={-2*math.sqrt(p):7.3f})  Hoffman kappa >= {math.ceil(1+(p+1)/(-lmin)-1e-9)}")
