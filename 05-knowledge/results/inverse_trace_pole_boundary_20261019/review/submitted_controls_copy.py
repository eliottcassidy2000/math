#!/usr/bin/env python3
"""Finite exact controls; not a search for or proof of a Keller realization."""
from pathlib import Path
import json
import sympy as S

CHECKS = 0

def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(label)


def eq(a, b, label):
    check(S.cancel(a-b) == 0, label)


z, u, v, q, p, t, w, a, b = S.symbols('z u v q p t w a b')


def newton_sums(poly, n):
    cs = S.Poly(poly, z).all_coeffs()
    d = len(cs)-1
    check(cs[0] == 1, 'monic polynomial')
    sums = [S.Integer(d)]
    for k in range(1, n+1):
        total = S.Integer(0)
        if k <= d:
            for j in range(1, k):
                total += cs[j]*sums[k-j]
            total += k*cs[k]
        else:
            for j in range(1, d+1):
                total += cs[j]*sums[k-j]
        sums.append(S.cancel(-total))
    return sums


def companion(poly):
    cs = S.Poly(poly, z).all_coeffs()
    d = len(cs)-1
    M = S.zeros(d)
    for j in range(d-1):
        M[j+1,j] = 1
    for j in range(d):
        M[j,d-1] = -cs[d-j]
    return M


def val_t(expr):
    expr = S.cancel(expr)
    if expr == 0:
        return S.oo
    num, den = expr.as_numer_denom()
    pn, pd = S.Poly(num,t), S.Poly(den,t)
    return min(k[0] for k,c in pn.terms() if c != 0) - min(k[0] for k,c in pd.terms() if c != 0)


families = []
for d in range(2, 11):
    # Polynomial, non-Keller hostile with one surviving sheet.
    chi = z**d-z/u-v/u
    sums = newton_sums(chi,d+2)
    M = companion(chi)
    power = S.eye(d)
    for k in range(1,d+3):
        power = power*M
        eq(S.trace(power),sums[k],f'companion/Newton d={d}, k={k}')
    for k in range(1,d-1):
        eq(sums[k],0,f'first blind polynomial moment d={d}, k={k}')
    eq(sums[d-1],(d-1)/u,f'first polar polynomial moment d={d}')
    eq(S.det(S.Matrix([p,p*q**d-q]).jacobian([q,p])),1-d*p*q**(d-1),f'non-Keller determinant d={d}')
    eq((p*q**d-q).subs({p:0,q:-v}),v,f'finite source branch d={d}')
    eq((1-d*p*q**(d-1)).subs(p,0),1,f'finite branch unramified d={d}')
    # Exact rational symplectic hostile, with one finite sheet.
    f = q/(1+q**d)
    fp = S.diff(f,q)
    eq(fp,(1-(d-1)*q**d)/(1+q**d)**2,f'cotangent derivative d={d}')
    eq(S.det(S.Matrix([f,p/fp]).jacobian([q,p])),1,f'cotangent determinant d={d}')
    eq((u*z**d-z+u).subs({u:f,z:q}),0,f'cotangent inverse equation d={d}')
    sums_rat = newton_sums(z**d-z/u+1,d)
    for k in range(1,d-1):
        eq(sums_rat[k],0,f'first blind rational moment d={d}, k={k}')
    eq(sums_rat[d-1],(d-1)/u,f'first rational polar moment d={d}')
    # Removing finite-sheet survival raises sharp detection depth from d-1 to d.
    all_escape = newton_sums(z**d-1/t,d)
    for k in range(1,d):
        eq(all_escape[k],0,f'all-escape blindness d={d}, k={k}')
    eq(all_escape[d],d/t,f'all-escape detector d={d}')
    families.append({'d':d,'first_polar_moment_with_survivor':d-1,'value':str(sums[d-1]),'without_survivor':d})

# Exact root-of-unity packets. Use cyclotomic remainder, not floating roots.
for m in range(1,13):
    cyc = S.cyclotomic_poly(m,z)
    for k in range(1,m+1):
        residue = S.rem(sum(z**(j*k) for j in range(m)),cyc,z)
        eq(residue,m if k == m else 0,f'root packet m={m},k={k}')

# Discriminants from symbolic resultants, with finite coefficients also present.
# The coefficient of z has pole a; the constant coefficient has pole b<=a.
discriminant_cases = []
for d in range(2,10):
    # Once per d: derive exact generic trinomial discriminant by resultant.
    disc_generic = S.discriminant(z**d+a*z+b,z)
    expected = (-1)**(d*(d-1)//2)*(d**d*b**(d-1)+(-1)**(d-1)*(d-1)**(d-1)*a**d)
    eq(disc_generic,expected,f'trinomial discriminant d={d}')
    for pole in range(1,5):
        for cpole in range(pole+1):
            disc = disc_generic.subs({a:-t**(-pole),b:-v*t**(-cpole)})
            check(val_t(disc) == -d*pole,f'discriminant pole d={d},a={pole},b={cpole}')
            g = S.gcd(pole,d-1)
            sign_exp = (d-1)-int(g)
            check(sign_exp % 2 == (d*pole) % 2,f'inertia/discriminant parity d={d},a={pole}')
            if d % 2 == 0:
                check(sign_exp % 2 == 0,f'even-degree even inertia d={d},a={pole}')
            discriminant_cases.append({'d':d,'moment_pole':pole,'constant_pole':cpole,'discriminant_order':-d*pole,'inertia_cycle_length':(d-1)//int(g),'inertia_cycle_count':int(g)})

# Generic quartic regular first two coefficients do not change the leading pole.
A1,A2,A3,A4=S.symbols('A1 A2 A3 A4')
D4=S.discriminant(z**4+A1*z**3+A2*z**2+A3*z+A4,z)
for pole in range(1,4):
    D=D4.subs({A3:a/t**pole,A4:b/t**pole})
    check(val_t(D)==-4*pole,f'quartic full coefficient discriminant pole={pole}')
    eq(S.cancel(t**(4*pole)*D).subs(t,0),-27*a**4,f'quartic leading discriminant pole={pole}')

# At d=4, maximal delay permits only identity or a 3-cycle, never transposition.
for pole in range(1,13):
    g=int(S.gcd(pole,3))
    check((3//g,g) in [(1,3),(3,1)],f'quartic packet inertia pole={pole}')

# Distinct leading vectors with first moment zero but second moment nonzero.
lin=[a,b,-a-b]
eq(sum(lin),0,'vector first moment')
eq(sum(x*x for x in lin),2*(a*a+a*b+b*b),'vector second moment')
check(S.expand(sum(x*x for x in lin)) != 0,'noncollinear packet detects degree two')

# Pole escalation and vanishing residues are separate.
for pole in range(1,8):
    for n in range(0,8):
        derivative = S.diff(t**(-pole),t,n)
        expected = (-1)**n*S.rf(pole,n)*t**(-pole-n)
        eq(derivative,expected,f'pole derivative a={pole},n={n}')
        check(val_t(derivative)==-pole-n,f'escalating order a={pole},n={n}')
    eq(S.residue(S.diff(t**(-pole),t),t,0),0,f'exact polar differential residue a={pole}')

# Local 3-cycle discrepancy hostile: boundary blowup chart and target model.
xloc=1/t; yloc=t**4*w
Tloc=t**3; Sloc=-w/3
eq(S.det(S.Matrix([xloc,yloc]).jacobian([t,w])),-t**2,'boundary volume order two')
eq(S.det(S.Matrix([Tloc,Sloc]).jacobian([t,w])),-t**2,'triple ramification volume agreement')

# Genuine polynomial Keller positive control, generic field degree one.
P=q+p**2; Q=p+P**2
ix=u-(v-u**2)**2; iy=v-u**2
eq(S.det(S.Matrix([P,Q]).jacobian([q,p])),1,'Henon-type positive control determinant')
eq(P.subs({q:ix,p:iy},simultaneous=True),u,'positive control first inverse')
eq(Q.subs({q:ix,p:iy},simultaneous=True),v,'positive control second inverse')

report={'status':'PASS','checks':CHECKS,'scope':'Finite exact algebra controls, no Keller realization census or general Jacobian closure','families':families,'discriminant_cases':discriminant_cases,'sympy_version':S.__version__}
Path(__file__).with_name('verification-results.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS: {CHECKS} exact checks')
print('Universe: d=2..10 trace controls; d=2..9, pole=1..4 discriminant controls; packet sizes 1..12.')
print('No global Keller realization or general trace regularity is inferred from these controls.')
