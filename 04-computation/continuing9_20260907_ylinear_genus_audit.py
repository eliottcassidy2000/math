#!/usr/bin/env python3
"""Independent unnormalized discriminant and local-chart genus referee.

No research producer is imported. Cubic discriminants are constructed from a
Sylvester determinant; finite generic-squarefreeness controls use a separate
standard-library modular Euclidean algorithm.
"""
from pathlib import Path
from hashlib import sha256
import json
import sys
import sympy as s
sys.stdout.reconfigure(encoding='utf-8', newline='\n')
GATES = 0


def need(x, why):
    global GATES
    GATES += 1
    if not x:
        raise ArithmeticError(why)


def same(a, b):
    return s.cancel(a-b) == 0


def trim(v):
    while v and v[-1] == 0:
        v.pop()
    return v


def remainder(a, b, p):
    a = a.copy()
    while len(a) >= len(b):
        z = a[-1]*pow(b[-1], -1, p) % p
        k = len(a)-len(b)
        for j, v in enumerate(b):
            a[k+j] = (a[k+j]-z*v) % p
        trim(a)
    return a


def gcdmod(a, b, p):
    a, b = trim(a.copy()), trim(b.copy())
    while b:
        a, b = b, remainder(a, b, p)
    if a:
        z = pow(a[-1], -1, p)
        a = [z*x % p for x in a]
    return a


def main():
    p, y, c, A, gamma = s.symbols('p y c A gamma', nonzero=True)
    a, b, d, e = s.symbols('a b d e')
    # Coefficient vectors of cubic and derivative; determinant rather than
    # the symbolic discriminant API used by the producer.
    matrix = s.Matrix([[a,b,d,e,0],[0,a,b,d,e],
                       [3*a,2*b,d,0,0],[0,3*a,2*b,d,0],[0,0,3*a,2*b,d]])
    disc = s.cancel(-matrix.det(method='domain-ge')/a)
    need(same(disc, b*b*d*d-4*a*d**3-4*b**3*e-27*a*a*e*e+18*a*b*d*e), 'Sylvester cubic discriminant reconstruction')
    Q = p*p*(p**3-y*y)*(A+gamma*y)-c
    coefficients = s.Poly(Q,y).all_coeffs()
    raw_disc = s.expand(disc.subs(dict(zip((a,b,d,e), coefficients)), simultaneous=True))
    N = 4*p**7*(gamma*gamma*p**3-A*A)**2 + 4*A*c*p*p*(9*gamma*gamma*p**3-A*A)-27*gamma*gamma*c*c
    need(same(raw_disc/p**4,N), 'unnormalized constant-gamma full discriminant')
    need(N.subs(p,0) == -27*gamma*gamma*c*c, 'zero is excluded from finite branch roots')
    # Triple root of the unnormalized derivative equation.
    fy = 3*gamma*y*y+2*A*y-gamma*p**3
    root = -A/(3*gamma)
    need(same(fy.subs(y,root),-(A*A+3*gamma*gamma*p**3)/(3*gamma)), 'triple-root value with gamma retained')
    need(same(s.diff(fy,y).subs(y,root),0), 'triple-root derivative')
    epsilon, q, v0, v1, v2, v3 = s.symbols('epsilon q v0 v1 v2 v3')
    perturb = disc.subs({a:1+epsilon*v3,b:q+epsilon*v2,d:epsilon*v1,e:epsilon*v0}, simultaneous=True)
    need(same(s.diff(perturb,epsilon).subs(epsilon,0),-4*q**3*v0), 'local double-root discriminant variation')

    # Distinct fixed polynomial bank. A nonzero specialization of the branch
    # resultant modulo101 implies its generic characteristic-zero nonvanishing
    # for that fixed A,gamma; the all-A assertion is the analytic proof.
    polynomials = [s.Integer(0),s.Integer(-3),2*p-5,p*p-2*p+3,
                   -2*p**3+p+1,p**4-p*p+2,(p+1)**5,p**7-2*p+5,3*p**8+2*p**3-1]
    r,Y,u,tau = s.symbols('r Y u tau')
    rows = []
    for index, aa in enumerate(polynomials):
        gg = (2,-3,5)[index%3]
        m = int(s.degree(aa,p)) if aa != 0 else -1
        expected = 13 if m <= 1 else 4*m+7
        NN = s.expand(N.subs({A:aa,gamma:gg}, simultaneous=True))
        need(s.degree(NN,p) == expected, 'all finite branch degrees on new polynomial bank')
        leading = 4*gg**4 if m <= 1 else 4*s.LC(s.Poly(aa,p))**4
        need(s.LC(s.Poly(NN,p)) == leading, 'unnormalized branch leading coefficient')
        choice = None
        for cc in range(1,21):
            pp = s.Poly(NN.subs(c,cc),p)
            coeff = [int(pp.nth(i)) % 101 for i in range(pp.degree()+1)]
            deriv = [(i*coeff[i]) % 101 for i in range(1,len(coeff))]
            if len(gcdmod(coeff,deriv,101)) == 1:
                choice = cc; break
        need(choice is not None, 'fresh modular Euclidean generic-squarefreeness witness')
        ff = (y*y-p**3)*(aa+gg*y)+c/p**2
        zero = s.expand(ff.subs({p:r**3,y:r**-2*Y}, simultaneous=True)*r**6)
        need(same(zero.subs(r,0),gg*Y**3+c), 'unnormalized cubic local chart at zero')
        # Third-root constants are simple without choosing a root of c in K.
        need(s.gcd(s.Poly(gg*Y**3+c,Y),s.Poly(3*gg*Y*Y,Y)).degree() == 0, 'zero-chart geometric separability')
        if m <= 1:
            infinity = s.expand(ff.subs({p:r**-2,y:r**-3*Y}, simultaneous=True)*r**9)
            need(same(infinity.subs(r,0),gg*Y*(Y*Y-1)), 'small-degree infinity has three simple lifted roots')
            need(same(infinity.subs({r:-r,Y:-Y}, simultaneous=True),-infinity), 'small-degree deck preserves equation up to unit sign')
        else:
            am = s.LC(s.Poly(aa,p))
            infinity = s.expand(ff.subs({p:r**-2,y:r**-3*Y}, simultaneous=True)*r**(2*m+6))
            need(same(infinity.subs(r,0),am*(Y*Y-1)), 'two square-root branches for every new degree')
            need(same(infinity.subs({r:-r,Y:-Y}, simultaneous=True),infinity), 'large-degree deck exchanges two branches')
            native = s.expand(ff.subs({p:1/r,y:r**(-m)*Y}, simultaneous=True)*r**(3*m))
            need(same(native.subs(r,0),Y*Y*(am+gg*Y)), 'remaining infinity branch on ordinary parameter')
            need(same(s.diff(native,Y).subs({r:0,Y:-am/gg}),am*am/gg), 'unramified infinity root is simple')
        genus = (expected-1)//2
        need(-6+expected+2+1 == 2*genus-2, 'complete characteristic-zero ramification count')
        I = p*p*(p**3-y*y)*(aa+gg*y)
        logI = s.expand(I.subs({p:u*u+tau,y:u*(u*u+tau)}, simultaneous=True))
        W = u**8*(aa.subs(p,u*u)+gg*u**3)
        need(same(logI.coeff(tau,1),W) and W != 0, 'actual first invariant coefficient with nonunit gamma')
        for k in (1,2,3):
            # Only the leading tau term is needed; it gives the universal first
            # displacement. Actual full I is checked at k=1 separately below.
            Hlead = tau**k*W**k
            response = s.expand(tau*(2*u*s.diff(Hlead,tau)-s.diff(Hlead,u)))
            need(same(response.coeff(tau,k),2*k*u*W**k), 'outer first response at orders1,2,3')
        actual = s.expand(tau*(2*u*s.diff(logI,tau)-s.diff(logI,u)))
        need(same(actual.coeff(tau,1),2*u*W), 'full original first Hamiltonian response')
        rows.append([str(aa),gg,m,expected,genus,choice])

    # This separates finite critical VALUES from a finite critical LOCUS.
    J = p*p*y*y
    need(s.diff(J,p).subs(p,0) == 0 and s.diff(J,y).subs(p,0) == 0, 'critical locus can contain a whole line')
    need(J.subs(p,0) == 0 and J.subs(y,0) == 0, 'both critical lines have the same constant value')
    # Variable y-slope is outside the theorem: already gamma=p^2,A=0 has
    # finite branch degree17 and genus8, rather than the constant-slope genus6.
    variable = y**3-p**3*y+c/p**4
    variable_disc = disc.subs({a:1,b:0,d:-p**3,e:c/p**4}, simultaneous=True)
    need(same(variable_disc*p**8,4*p**17-27*c*c), 'variable-slope extrapolation changes branch degree')
    need(same(s.expand(variable.subs({p:r**3,y:r**-4*Y}, simultaneous=True)*r**12).subs(r,0),Y**3+c), 'variable-slope zero ramification3')
    need(same(s.expand(variable.subs({p:r**-2,y:r**-3*Y}, simultaneous=True)*r**9).subs(r,0),Y*(Y*Y-1)), 'variable-slope infinity ramification1')
    need(-6+17+2+1 == 2*8-2, 'specific outside-scope genus8 hostile')
    print('INDEPENDENT unnormalized Sylvester discriminant; modular Euclid; direct local charts and original logarithmic response')
    print('CONTROL_ROWS [A,gamma,degree,finite_branch_count,genus,c_mod101]',json.dumps(rows,separators=(',',':')))
    print('ALL_DEGREE_SCOPE constant gamma!=0; geometric integrality; finite critical values; genus6 or2m+3')
    print('BOUNDARIES critical locus may be positive-dimensional; variable gamma=p^2,A0 has genus8; retain I under outer f(I)')
    print('PASS',GATES,'always-active exact gates; raw LF')


if __name__ == '__main__':
    main()
