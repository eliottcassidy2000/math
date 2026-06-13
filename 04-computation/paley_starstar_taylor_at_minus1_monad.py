#!/usr/bin/env python3
"""
paley_starstar_taylor_at_minus1_monad.py
monad-explorer-2026-06-07 (deep-research, 11th session)

The TWO THM-438 handoffs are the first two Taylor coefficients of V(s,y) at s=-1:
   a_0(y) = V(-1,y)         = -y                          (handoff #1: deg P_m=m-2)
   a_1(y) = V_s(-1,y)       = sum (2^m-1) y^m = y/((1-y)(1-2y))   (handoff #2: lead=2^m-1)
where a_n(y) = sum_m [Q_m^{(n)}(-1)/n!] y^m.

QUESTION: is a_n(y) rational with a CLEAN denominator (a product of (1-j y))?  If the
pattern is denom_n = prod_{j=1}^{?}(1-j y), that is a NEW structural law and PREDICTS
all P_m coefficients top-down -- in particular c_2 of P_5 (testable by the k=7 run).
"""
import sympy as sp

s, y = sp.symbols('s y')
Q = {
 1: s,
 2: 3*s**2 + 3*s**3,
 3: 13*s**3 + 33*s**4 + 20*s**5,
 4: 69*s**4 + 304*s**5 + 416*s**6 + 181*s**7,
}
# also a partial Q5/Q6 lower coeffs (independent of c2,c3,c4 at low n) for cross-check
c2,c3,c4 = sp.symbols('c2 c3 c4')
Q[5] = 421*s**5 + 2740*s**6 + (5694+c2)*s**7 + (4852+2*c2+c3)*s**8 + (1477+c2+c3+c4)*s**9

NMAX = 4
print("Taylor coefficients a_n(y) = sum_m Q_m^{(n)}(-1)/n!  y^m  (m=1..4 exact)")
A = {}
for n in range(0, NMAX+1):
    coeffs = {}
    for m in range(1, 5):
        val = sp.diff(Q[m], s, n).subs(s, -1) / sp.factorial(n)
        coeffs[m] = sp.nsimplify(val)
    A[n] = coeffs
    seq = [coeffs[m] for m in range(1,5)]
    print(f"  a_{n}: [m=1..4] = {seq}")

print("\nFit each a_n(y) to rational num/prod(1-j y):")
def fit_rational(seq, nmax_num=4):
    """seq = [a_1,a_2,a_3,a_4] coeffs of y^1..y^4. Try denominators
       prod_{j in J}(1-j y) for small J; find smallest with polynomial numerator."""
    f = sum(seq[m-1]*y**m for m in range(1,5))  # series y^1..y^4
    best = None
    import itertools
    cand_factors = [1,2,3,4,5]
    for r in range(0, 5):
        for combo in itertools.combinations(cand_factors, r):
            denom = sp.prod([(1-j*y) for j in combo]) if combo else sp.Integer(1)
            num = sp.series(f*denom, y, 0, 5).removeO()
            num_poly = sp.Poly(sp.expand(num), y)
            # check if num is "low degree" (a polynomial of degree <= len(combo)+1, i.e. truncation-stable)
            deg = num_poly.degree()
            # heuristic: accept if degree small relative to data
            if deg <= r+1:
                best = (combo, sp.expand(num))
                return best
    return None

for n in range(0, NMAX+1):
    seq = [A[n][m] for m in range(1,5)]
    if all(v==0 for v in seq):
        print(f"  a_{n}: all zero"); continue
    r = fit_rational(seq)
    if r:
        combo, num = r
        denom = "*".join(f"(1-{j}y)" for j in combo) if combo else "1"
        print(f"  a_{n}(y) = ({num}) / [{denom}]   [verify: matches y^1..y^4]")
        # PREDICT a_n coeff at m=5 and m=6 from this rational fit
        full = sp.expand(num)/ (sp.prod([(1-j*y) for j in combo]) if combo else 1)
        ser = sp.series(full, y, 0, 7).removeO()
        p5 = sp.Poly(ser, y).coeff_monomial(y**5)
        p6 = sp.Poly(ser, y).coeff_monomial(y**6)
        print(f"        PREDICT a_{n}[m=5]={p5}, a_{n}[m=6]={p6}")
    else:
        print(f"  a_{n}: no clean (1-jy)-product denominator found in search range")

print("\nNow connect a_n to P_m TOP-DOWN coefficients.")
print("lead P_m = a_1[m]; we want to know which P_m coeff each a_n controls.")
# Relation: P_m(x), x = s/(1+s).  Near s=-1, x->infty.  Express P_m top-coeffs via a_n.
# Use: T_m = Q_m(s); P_m(x) = Q_m(s)/(s^m (1+s)^{m-1}) with x=s/(1+s).
# coeff of x^{m-1-n} in P_m  <-> a_n  (to be confirmed numerically):
for m in range(2,5):
    Pm = sp.cancel(Q[m]/(s**m*(1+s)**(m-1)))   # = P_m(x) with x=s/(1+s); express in x
    # substitute s = x/(1-x)
    x = sp.Symbol('x')
    Pm_x = sp.cancel(Pm.subs(s, x/(1-x)))
    Pm_poly = sp.Poly(sp.expand(Pm_x), x)
    print(f"  m={m}: P_m(x) = {Pm_poly.as_expr()}  (coeffs low->high {Pm_poly.all_coeffs()[::-1]})")
