#!/usr/bin/env python3
"""
paley_starstar_taylor_pochhammer_monad.py
monad-explorer-2026-06-07 (deep-research, 11th session)

CONJECTURE (new): the s=-1 Taylor coefficients of V(s,y),
   a_n(y) = sum_m [Q_m^{(n)}(-1)/n!] y^m,
are rational with denominator the Pochhammer-like product
   D_n(y) = prod_{j=1}^{n+1} (1 - j y),
so a_n[m] = sum_{j=1}^{n+1} c_{n,j} j^m  (a sum of geometric sequences).
n=0:  -y/(1-y) = -y          (handoff #1)
n=1:  y/((1-y)(1-2y)) -> 2^m-1   (handoff #2)

Test with PROPER rational fits, then cross-check using the PARTIAL Q_5 (which carries
c2,c3,c4) together with the handoff values c4=0, c3=31 -> solve for c2 three ways.
"""
import sympy as sp
import itertools

s, y = sp.symbols('s y')
c2,c3,c4 = sp.symbols('c2 c3 c4')
Q = {
 1: s,
 2: 3*s**2 + 3*s**3,
 3: 13*s**3 + 33*s**4 + 20*s**5,
 4: 69*s**4 + 304*s**5 + 416*s**6 + 181*s**7,
 5: 421*s**5 + 2740*s**6 + (5694+c2)*s**7 + (4852+2*c2+c3)*s**8 + (1477+c2+c3+c4)*s**9,
}

def a_n_m(n, m):
    return sp.expand(sp.diff(Q[m], s, n).subs(s, -1) / sp.factorial(n))

print("a_n[m] table (m=1..4 exact; m=5 carries c2,c3,c4):")
for n in range(0,5):
    row = [a_n_m(n,m) for m in range(1,5)]
    a5 = a_n_m(n,5)
    print(f"  a_{n}: {row}   | a_{n}[5] = {a5}")

print("\n--- PROPER rational fit a_n(y)=N(y)/D(y), D=prod (1-jy), deg N < deg D ---")
def proper_fit(n):
    data = [a_n_m(n,m) for m in range(1,5)]  # m=1..4
    fser = sum(data[m-1]*y**m for m in range(1,5))
    best = None
    for r in range(1,5):
        for combo in itertools.combinations([1,2,3,4,5], r):
            D = sp.prod([(1-j*y) for j in combo])
            prod = sp.series(fser*D, y, 0, 5).removeO()
            P = sp.Poly(sp.expand(prod), y)
            if P.degree() < r:  # PROPER
                # we used 4 data pts; numerator has <= r coeffs (deg<r => <=r coeffs incl const)
                # number of independent checks = 4 - (#nonzero-allowed numerator coeffs)
                best = (combo, sp.expand(prod), r)
                return best
    return None

preds = {}
for n in range(0,5):
    if all(a_n_m(n,m)==0 for m in range(1,5)):
        print(f"  a_{n}: identically 0 on m<=4"); continue
    r = proper_fit(n)
    if r is None:
        print(f"  a_{n}: NO proper (1-jy)-product fit found (search r<=4)")
        continue
    combo, num, deg = r
    Dstr = "".join(f"(1-{j}y)" for j in combo)
    full = sp.expand(num)/sp.prod([(1-j*y) for j in combo])
    ser = sp.series(full, y, 0, 8).removeO()
    p5 = sp.Poly(ser,y).coeff_monomial(y**5)
    p6 = sp.Poly(ser,y).coeff_monomial(y**6)
    p7 = sp.Poly(ser,y).coeff_monomial(y**7)
    preds[n] = p5
    # partial-fraction to show geometric structure
    pf = sp.apart(full, y)
    print(f"  a_{n}(y) = ({num}) / {Dstr}")
    print(f"        PROPER, denom-degree {len(combo)}; predicts a_{n}[5]={p5}, [6]={p6}, [7]={p7}")
    print(f"        partial fractions: {pf}")

print("\n--- CROSS-TEST: pin c2 from a_n[5] predictions (set c4=0, c3=31 from handoffs) ---")
subs0 = {c4:0, c3:31}
for n in [2,3,4]:
    if n not in preds: continue
    lhs = a_n_m(n,5).subs(subs0)   # linear in c2
    eq = sp.Eq(lhs, preds[n])
    solc2 = sp.solve(eq, c2)
    print(f"  a_{n}[5]: {lhs} = {preds[n]}  ->  c2 = {solc2}")
print("\n(If all three give the SAME integer c2, the Pochhammer law + c2 are corroborated;")
print(" the background k=7 run (t(7,5)) provides the independent confirmation.)")
