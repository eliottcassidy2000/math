"""JC(2): is W = 0 FORCED?  (the first open question named in THM-3016 section 5)

THM-3016 (R):  J * Jac(W, H) = 0  with  W = P_{n-1} - kappa H^{a-b} Q_{m-1},
deg W = n-1,  H homogeneous of degree g,  n = ga,  m = gb,  gcd(a,b) = 1.
J = Jac(P,Q) is a NONZERO constant, so  Jac(W,H) = 0.

CLAIM.  For binary forms, Jac(W,H) = 0 with W != 0 forces W^g = c H^{n-1},
and then, because gcd(g, n-1) = gcd(g, ga-1) = 1, H must be a PURE POWER OF A
SINGLE LINEAR FORM, H = l^g.  That is K = 1 -- the one-place-at-infinity
configuration.  Hence for any pair with K >= 2 directions at infinity,
W = 0 is FORCED.
"""
import sympy as sp
from math import gcd

x, y = sp.symbols('x y')

def jac(f, g_):
    return sp.expand(sp.diff(f, x) * sp.diff(g_, y) - sp.diff(f, y) * sp.diff(g_, x))

print("STEP 1.  For binary forms, Jac(W,H)=0  <=>  W^deg(H) = c H^deg(W).")
print("   (check both directions on random forms)")
ok = True
for (H, W, tag) in [
    ((x + y) ** 3, (x + y) ** 5, "common linear form"),
    ((x + 2 * y) ** 2 * (x - y) ** 2, ((x + 2 * y) * (x - y)) ** 6, "common quadratic"),
    ((x + y) * (x - y), (x + y) ** 2, "NOT proportional powers"),
]:
    J = sp.expand(jac(W, H))
    dH, dW = sp.total_degree(sp.expand(H)), sp.total_degree(sp.expand(W))
    prop = sp.simplify(sp.expand(W ** dH) / sp.expand(H ** dW))
    isprop = prop.free_symbols == set()
    print(f"   {tag:26s}: Jac=0? {J == 0};  W^dH / H^dW constant? {isprop}")
    ok &= (J == 0) == isprop

print(f"   equivalence held on all samples: {ok}")

print("\nSTEP 2.  The arithmetic: gcd(g, n-1) = gcd(g, ga-1) = 1 always.")
bad = [(g_, a_) for g_ in range(1, 40) for a_ in range(1, 40) if gcd(g_, g_ * a_ - 1) != 1]
print(f"   pairs (g,a) with 1<=g,a<=39 and gcd(g, ga-1) != 1:  {bad if bad else 'NONE'}")
print("   => in W^g = c H^{n-1}, writing H = prod l_i^{e_i} forces g | e_i(n-1),")
print("      and gcd(g,n-1)=1 gives g | e_i; with sum e_i = g and e_i >= 1 there is")
print("      EXACTLY ONE i, with e_i = g.  So H = l^g:  a single direction, K = 1.")

print("\nSTEP 3.  Direct check: for K >= 2 (H with >= 2 distinct linear factors),")
print("   is W = 0 the ONLY homogeneous solution of Jac(W,H)=0 in degree n-1 = ga-1?")
for (Hexpr, g_, a_, tag) in [
    ((x + y) * (x - y), 2, 2, "H=(x+y)(x-y), g=2, a=2 -> deg W=3"),
    ((x + y) * (x - y), 2, 3, "H=(x+y)(x-y), g=2, a=3 -> deg W=5"),
    ((x + y) ** 2 * (x - y), 3, 2, "H=(x+y)^2(x-y), g=3, a=2 -> deg W=5"),
    (x * y * (x + y), 3, 3, "H=xy(x+y), g=3, a=3 -> deg W=8"),
    ((x + y) * (x - y) * (x + 2 * y), 3, 4, "H=(x+y)(x-y)(x+2y), g=3,a=4 -> deg W=11"),
]:
    dW = g_ * a_ - 1
    cs = sp.symbols(f'c0:{dW+1}')
    W = sum(cs[i] * x ** i * y ** (dW - i) for i in range(dW + 1))
    eqs = sp.Poly(jac(W, Hexpr), x, y).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    allzero = all(all(sp.simplify(s.get(c, 0)) == 0 for c in cs) for s in sol) if sol else False
    K = len(sp.factor_list(sp.expand(Hexpr))[1])
    print(f"   {tag:44s} K={K}: only W=0? {allzero}")

print("\nSTEP 4.  Control -- when K = 1 the solution space is genuinely NONZERO:")
for (Hexpr, g_, a_) in [((x + y) ** 2, 2, 2), ((x + y) ** 3, 3, 2)]:
    dW = g_ * a_ - 1
    cs = sp.symbols(f'c0:{dW+1}')
    W = sum(cs[i] * x ** i * y ** (dW - i) for i in range(dW + 1))
    sol = sp.solve(sp.Poly(jac(W, Hexpr), x, y).coeffs(), cs, dict=True)
    Wsol = sp.simplify(W.subs(sol[0])) if sol else 0
    print(f"   H=(x+y)^{g_}, deg W={dW}:  W = {sp.factor(Wsol)}   (nonzero <=> K=1)")
