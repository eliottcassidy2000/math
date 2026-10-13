"""HYP-9257: the cusp curve of THM-1300's Keller map carries no constituent.

Checks (sympy, exact): the Jelonek surface c_3 = 0 is singular exactly along the cusp curve
C = (4/(27 t^2), 4/(3 t), t); a transversal slice at a point of C is an ordinary cusp (multiplicity
two, double-line tangent cone, Milnor number two), so the local Euler obstruction of the surface
along C is two.  Then the integer bookkeeping of the Kashiwara-Dubson index formula on the fibre
stratification (U, J \\ C, C) with stalk Euler characteristics 3, 1, 0 forces the characteristic
cycle multiplicities (3, 2, 1) and leaves no room for a constituent supported on C:
length_(A_3)(V_F) = 3 on that stratification, not 4.  See
05-knowledge/results/keller_module_length_cusp_20261012.md for the argument and its cited inputs.
Run from the repository root: python3 -B <script>
"""
import sympy as sp

CHECKS = 0


def check(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


a, b, c, t, u, v = sp.symbols('a b c t u v')
c3 = 27*a**2*c**2 - 18*a*b*c + 16*a + b**3*c - b**2          # Jelonek surface J (THM-4621, Gao's c_3)
C = {a: sp.Rational(4, 27)/t**2, b: sp.Rational(4, 3)/t, c: t}  # the cusp curve

check(sp.simplify(c3.subs(C)) == 0, 'C lies on J')
grad = [sp.diff(c3, w) for w in (a, b, c)]
check(all(sp.simplify(g.subs(C)) == 0 for g in grad), 'J is singular along C')
sol = sp.solve([c3]+grad, [a, b, c], dict=True)
check(len(sol) == 1 and sp.simplify(sol[0][a]-sp.Rational(4, 27)/c**2) == 0 and sp.simplify(sol[0][b]-sp.Rational(4, 3)/c) == 0,
      'the singular locus of J is exactly C')

# transversal slice at t = 1 (tangent of C there has a nonzero c-component, so the plane c = const is transversal)
p0 = {a: sp.Rational(4, 27), b: sp.Rational(4, 3), c: 1}
tangent = [sp.diff(C[w], t).subs(t, 1) for w in (a, b, c)]
check(tangent[2] != 0, 'the plane c = 1 is transversal to C at t = 1')
f = sp.expand(c3.subs({a: p0[a]+u, b: p0[b]+v, c: 1}))
P = sp.Poly(f, u, v)
mult = min(sum(m) for m in P.monoms())
check(mult == 2, 'slice curve has multiplicity two at the point of C')
low = sum(coef*u**m[0]*v**m[1] for m, coef in zip(P.monoms(), P.coeffs()) if sum(m) == 2)
check(sp.factor(low) == 3*(3*u-v)**2, 'tangent cone is a double line')
fu, fv = sp.diff(f, u), sp.diff(f, v)
r = sp.Poly(sp.resultant(fu, fv, v), u)
coeffs = r.all_coeffs()[::-1]
check(coeffs[0] == 0 and coeffs[1] == 0 and coeffs[2] != 0, 'Milnor number two: an ordinary cusp (A_2)')

# the index-formula bookkeeping (section 3 of the note)
phi = {'U': 3, 'J-C': 1, 'C': 0}            # stalk Euler characteristics of Rf_* C
Eu_J_on_C = 2                                # local Euler obstruction of a cuspidal edge
m3 = phi['U']
mJ = phi['J-C']-(-m3)*(-1)                   # -m3 + mJ*1 = -phi(J-C)  ->  mJ = m3 - phi(J-C)
mJ = m3-phi['J-C']
mC = mJ*Eu_J_on_C-m3-(-phi['C'])             # -m3 + mJ*Eu - mC = -phi(C)
mC = mJ*Eu_J_on_C-m3+phi['C']
check((m3, mJ, mC) == (3, 2, 1), 'characteristic cycle multiplicities (3, 2, 1)')
# constituents: O (0 on T*C), IC(std) (a_C >= 0), IC_J (b_C >= 1), constituents on C (each >= 1)
b_C_min = 1
check(mC-b_C_min == 0, 'no multiplicity left for a constituent supported on C')
print({'multiplicities': (m3, mJ, mC), 'length_on_fibre_stratification': 3, 'checks': CHECKS})
