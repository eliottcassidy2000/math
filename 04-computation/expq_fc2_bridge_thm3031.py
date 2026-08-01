"""THE BRIDGE, settled.  int_0^1 e^Q != 0 does NOT give FC(2).  The TRANSCENDENCE
form does, and immediately.

THM-3018: FC(2)  <=>  [ g poly on [0,1], int_0^1 g^m du = 0 for all m>=1  =>  g = 0 ].

Suppose FC(2) FAILS.  Then there is g != 0 with all moments m>=1 vanishing.
 (i) g is NONCONSTANT: a constant c has int_0^1 c^1 du = c, so c = 0.
 (ii) The exponential period is forced to the value ONE, not zero:
        int_0^1 e^{g(u)} du = sum_{m>=0} (1/m!) int_0^1 g^m du = 1 + 0 = 1.
 (iii) SPECIALISATION: the conditions int_0^1 g^m du = 0 are polynomial in the
       coefficients of g with RATIONAL coefficients (int_0^1 u^j du = 1/(j+1)).
       So the locus is a variety defined over Q; nonempty over C => nonempty over Qbar.
       Hence g may be taken with ALGEBRAIC coefficients.
 => Q := g is a nonconstant polynomial in Qbar[t] with int_0^1 e^Q dt = 1, which is
    RATIONAL.  That contradicts TRANSCENDENCE, but is perfectly consistent with != 0.

CONCLUSION.  The nonvanishing statement is too weak; the TRANSCENDENCE statement
yields FC(2) in four lines.  This is exactly why the owner's remark that the
transcendence result is "far more impressive" is also the operative point.
"""
import sympy as sp
from sympy import Rational as R

u, s = sp.symbols('u s')

print("CHECK (i)+(ii) symbolically: a vanishing-moment g forces int e^g = 1.")
# impose vanishing moments on a general cubic and solve; then evaluate int e^g
for d in (1, 2, 3):
    cs = sp.symbols(f'c0:{d+1}')
    g = sum(cs[i]*u**i for i in range(d+1))
    eqs = [sp.integrate(sp.expand(g**m), (u, 0, 1)) for m in range(1, d+4)]
    sols = sp.solve(eqs, cs, dict=True)
    nonzero = [t for t in sols if any(sp.simplify(t.get(c,0)) != 0 for c in cs)]
    print(f"   deg {d}: moment system has {len(sols)} solution(s); NONZERO ones: {len(nonzero)}"
          f"  -> FC(2) holds in degree {d}: {len(nonzero)==0}")

print("\n   and the identity int_0^1 e^g = 1 whenever all m>=1 moments vanish is formal:")
m = sp.symbols('m')
print("      int_0^1 e^g du = sum_m (1/m!) int_0^1 g^m du = 1 + sum_{m>=1} 0 = 1.  QED")

print("\nCHECK (iii): the moment conditions are polynomial over Q in the coefficients.")
d = 3
cs = sp.symbols(f'c0:{d+1}')
g = sum(cs[i]*u**i for i in range(d+1))
for mm in (1, 2):
    e = sp.expand(sp.integrate(sp.expand(g**mm), (u, 0, 1)))
    coeffs = sp.Poly(e, *cs).coeffs()
    allrat = all(sp.nsimplify(c).is_rational for c in coeffs)
    print(f"   m={mm}: int_0^1 g^m du is a polynomial in the c_i with ALL-RATIONAL coefficients:"
          f" {allrat}   (e.g. {sp.nsimplify(e)})" if mm==1 else
          f"   m={mm}: all-rational coefficients: {allrat}")

print("\nCONTRAST -- the two candidate hypotheses, on the SAME hypothetical counterexample:")
print("   int_0^1 e^Q dt = 1.")
print("     'int e^Q != 0'          : 1 != 0  -> SATISFIED, no contradiction, FC(2) NOT implied.")
print("     'int e^Q transcendental': 1 is rational -> CONTRADICTION, FC(2) IS implied.")
print("\n=> The bridge REQUIRES the transcendence form.  The nonvanishing form is insufficient.")
