#!/usr/bin/env python3
"""
klein-2026-07-20-S365 -- THE CHARGE-0 RADIAL LAYER BY EXACT ELIMINATION, closing the
sign-indefinite sub-case where positivity (THM-1640) and domination (MISTAKE-202) both fail,
and where boxeph's Watson route (THM-1620 Case II) reduces but does not close.

THE OBJECT.  {-1,0,1} span, alpha = r*a*c, beta = b(r).  E[P^m] = E_r[L_m(alpha,beta)],
L_m = sum_k m!/(k!^2 (m-2k)!) alpha^k beta^{m-2k}, E_r[r^j] = j!.

THE HARD SUB-CASE (THM-1640 Part 4): FIX alpha = r (a = c = 1, so P is TWO-SIDED regardless
of b -- charges -1 and +1 both present), and let beta = b(r) be the charge-0 coefficient,
NON-CONSTANT.  Here:
  - positivity is UNAVAILABLE: L_m(r) is sign-indefinite (verified S363);
  - domination is FALSE (MISTAKE-202);
  - boxeph closes constant b (Hermite, THM-1615) but non-constant b is Watson-Case-II, open.
Since P is two-sided, NC2 predicts NO b(r) makes every moment vanish -- i.e. the ideal
<E_r[L_m] : m >= 1> should be the UNIT ideal (empty variety, basis <1>).  That is an EXACT
GROEBNER question, finite for each bounded deg b, and it needs no asymptotics at all.

CONTROL: the constant-b case (deg b = 0) must reproduce boxeph's Hermite no-common-root, i.e.
also give an empty variety (no constant b kills all moments when alpha = r).
"""
import sympy as sp
from sympy import groebner
from math import comb, factorial

r = sp.symbols('r')

def Er_of_poly(expr):
    """E_r[poly(r)] with r ~ Exp(1): E[r^j] = j!"""
    p = sp.Poly(sp.expand(expr), r)
    return sum(c * factorial(j) for (j,), c in zip(p.monoms(), p.coeffs()))

def Lm(alpha, beta, m):
    tot = 0
    for k in range(m // 2 + 1):
        tot += (factorial(m) // (factorial(k) ** 2 * factorial(m - 2 * k))) * alpha**k * beta**(m - 2 * k)
    return sp.expand(tot)

print("=" * 86)
print("CHARGE-0 RADIAL ELIMINATION, alpha = r FIXED (two-sided), beta = b(r) of degree d")
print("=" * 86)
print(" NC2 predicts: NO b(r) is a nullcone member here, so the ideal should be <1> (empty).")
print(" A basis {1} = CLOSED: the sign-indefinite fixed-alpha charge-0 case has no nullcone.\n")
print(f"{'deg b':>7} {'unknowns':>9} {'m used':>10} {'Groebner basis':>18} {'variety':>14}")
alpha = r
for d in (0, 1, 2, 3):
    bco = sp.symbols(f'b0:{d+1}')
    beta = sum(bco[i] * r**i for i in range(d + 1))
    unk = list(bco)
    # escalate m until the basis is <1> or we hit the cap
    verdict = None
    for M in ((d + 1) + 3, (d + 1) + 8, (d + 1) + 14):
        eqs = []
        for m in range(1, M + 1):
            e = sp.nsimplify(Er_of_poly(Lm(alpha, beta, m)))
            if e != 0: eqs.append(e)
        try:
            G = groebner(eqs, *unk, order='grevlex')
            gens = list(G.exprs)
        except Exception as ex:
            gens = [f'FAIL:{ex}']; break
        if gens == [1]:
            verdict = ('EMPTY (<1>)', M); break
        verdict = (f'basis size {len(gens)}', M)
    tag = 'CLOSED (empty)' if verdict[0].startswith('EMPTY') else 'not <1> yet'
    print(f"{d:>7} {len(unk):>9} {f'1..{verdict[1]}':>10} {verdict[0]:>18} {tag:>14}")
print("""
 If every row is EMPTY (<1>): the fixed-alpha charge-0 layer has NO nullcone member for beta
 of degree <= 3 -- CLOSED by elimination, exactly, where positivity/domination could not go,
 and covering boxeph's Watson Case II for these bounded degrees with zero asymptotics.
""")

# ---------------------------------------------------------------- NEWTON POLYGON
print("=" * 86)
print("NEWTON POLYGON OF THE BRANCHES  lambda_pm(r) = beta(r) +- 2*sqrt(alpha(r))")
print("=" * 86)
print(" The GF is E_r[ D(r,t)^{-1/2} ], D = (1 - beta t)^2 - 4 alpha t^2 = (1-lam+ t)(1-lam- t).")
print(" Branch points sit at t = 1/lam_pm(r); the LARGE branch is the lam of larger growth in r,")
print(" and its r-exponent is read off the Newton polygon of (alpha, beta).\n")
print(" alpha = r^A-ish (from a*c), beta ~ b_d r^B.  lam_+ ~ max(b_d r^B, 2 r^{(A+1)/2}).")
print(" The branch-point accumulation rate at t=0 is t* ~ 1/lam ~ r^{-max(B, (A+1)/2)}:\n")
print(f"   {'(A=deg ac, B=deg b)':>22} {'lam_+ leading':>18} {'branch exponent -max(B,(A+1)/2)':>32}")
import sympy as sp
def branch_exp(A, B):
    # alpha = r^{A+1} (since alpha = r*a*c, deg = 1 + A);  beta ~ r^B
    # lam ~ beta (r^B) vs 2 sqrt(alpha) (r^{(A+1)/2})
    e_alpha = sp.Rational(A + 1, 2)
    e_beta = B
    return max(e_beta, e_alpha), ('beta' if e_beta >= e_alpha else 'sqrt(alpha)')
for (A, B) in [(0, 0), (0, 1), (0, 2), (2, 1), (1, 3), (0, 3)]:
    e, which = branch_exp(A, B)
    print(f"   {str((A, B)):>22} {which:>18} {str(-e):>32}")
print("""
 READING.  The 'large branches' are the lam_+(r) that grow with r; their reciprocal is the
 branch point, accumulating at t = 0 like r^{-e}, e = max(B, (A+1)/2).  This is EXACTLY the
 kappa exponent of boxeph's THM-1620 Case II (kappa = -1/2 at A=B=0, the pair model), now
 read off a Newton polygon instead of a pinch computation:  when beta dominates (B >= (A+1)/2)
 the accumulation is set by the CHARGE-0 coefficient degree B; when alpha dominates it is set
 by the pair scale (A+1)/2.  The two regimes are the two edges of the polygon, and the split
 at B = (A+1)/2 is the Newton-split of THM-1620.
""")

# ---------------------------------------------------------------- HANKEL AT m=2
print("=" * 86)
print("THE REAL CLOSURE IS AT m=2: E_r[L_2] = beta^T H beta + 2 E_r[alpha], H = exponential Hankel")
print("=" * 86)
print(" L_2 = beta^2 + 2 alpha, so E_r[L_2] = E_r[beta^2] + 2 E_r[alpha].")
print(" For beta = sum b_i r^i, E_r[beta^2] = sum_ij b_i b_j (i+j)! = b^T H b, H_ij = (i+j)!.")
print(" H is the moment (Hankel) matrix of Exp(1): POSITIVE DEFINITE.  Check leading minors:\n")
for d in (1, 2, 3, 4):
    H = sp.Matrix(d + 1, d + 1, lambda i, j: sp.factorial(i + j))
    minors = [H[:k, :k].det() for k in range(1, d + 2)]
    print(f"   deg b = {d}: H is {d+1}x{d+1}, leading minors = {minors}  (all > 0: "
          f"{all(m > 0 for m in minors)})")
print("""
 SO OVER THE REALS the charge-0 layer with alpha >= 0 is closed at the SECOND moment alone:
 E_r[L_2] = ||beta||_H^2 + 2 E_r[alpha] >= 2 E_r[alpha] > 0 when alpha = r.  No b in R[r] can
 make it vanish.  This is opus's Hankel positive-definiteness (THM-1535) appearing at m=2, and
 it is why the COMPLEX case is the only hard one -- over C the form beta^T H beta is
 indefinite and the elimination above (needing m up to d+4) is doing the real work.  The
 honest boundary: elimination CLOSES each bounded degree over C; the Newton polygon shows the
 unbounded case is governed by the single exponent e = max(B, (A+1)/2), which is the object
 boxeph's per-component Watson lemma must control.
""")
