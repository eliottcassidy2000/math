# audit_H1: explicit counterexample to the multiplicative form of HYP-9490's Consequence when c3(T) = c3(T'): order-5 classes with c3 = 4 have different sw_2, so H(Gamma^r T)/H(Gamma^r T') != 1 for every r, while 1 + 2*0*rho*(1+o(1)) = 1.
from fractions import Fraction
from audit_tournament_lib import classes, c3, cover_poly, gate

cl = [A for A in classes(5) if c3(A) == 4]
cps = [cover_poly(A) for A in cl]
print("order-5 classes with c3 = 4: cover polynomials (low -> high) =", cps)
W = {1: 1}
Wl = [dict(W)]
for r in range(60):
    new = {}
    for k, c in W.items():
        if k == 0:
            continue
        for kk, wt in ((k + 1, k), (k, k * k), (k - 1, -k)):
            new[kk] = new.get(kk, 0) + wt * c
    W = new
    Wl.append(W)
U = lambda d, r: sum(c * (k ** d if d else 1) for k, c in Wl[r].items())
H = lambda co, r: sum(co[e] * U(e, r) for e in range(len(co)))
a, b = cps[0], cps[-1]
for r in (0, 1, 5, 20, 60):
    ratio = Fraction(H(a, r), H(b, r))
    rho = Fraction(U(3, r), U(5, r))
    print(f"r = {r:2d}: H ratio - 1 = {float(ratio - 1):+.3e}  (rho_(5,r) = {float(rho):.3e}; "
          f"(ratio-1)/rho = {float((ratio - 1) / rho):+.3e}; predicted second-order term "
          f"(sw2 - sw2')U_(1,r)/U_(5,r) = {float(Fraction((a[1] - b[1]) * U(1, r), U(5, r))):+.3e})")
