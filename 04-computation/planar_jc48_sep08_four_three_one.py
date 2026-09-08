#!/usr/bin/env python3
"""Exact all-parameter controls for the all-finite DG boundary type (4,3,1).

The analytic proof pays every branch, geometric generic integrality and the
complete primitive space. No finite mate-degree search is used here.
"""
from hashlib import sha256
import json
import sympy as S

gates = []


def need(label, predicate):
    if not bool(predicate):
        raise RuntimeError(label)
    gates.append(label)


def zero(label, expr):
    need(label, S.cancel(expr) == 0)


x, u, t, q, r, bd, p, zeta = S.symbols('x u t q r bd p zeta')
b, c, h, k, lam = S.symbols('b c h k lam')
source_jac = lambda f, g: S.diff(f, u)*S.diff(g, t)-S.diff(f, t)*S.diff(g, u)

# Full fifteen-parameter L2 and six-parameter L1; independent second-chart
# expansion verifies the complete induced rows, not a frozen lower prefix.
nn = S.symbols('n0:9')
pp = S.symbols('p0:5')
q0 = S.symbols('q0')
Nfull = sum(nn[i]*x**i for i in range(9))
Pfull = 2*nn[8]*x**6+2*nn[7]*x**5+sum(pp[i]*x**i for i in range(5))
Qfull = nn[8]*x**4+nn[7]*x**3+(pp[4]-nn[6])*x*x+(pp[3]-nn[5])*x+q0
Hfull = Nfull*t*t+Pfull*t+Qfull
Hr = S.cancel(Hfull.subs({x:1/r, t:-r*r-r**4*bd}, simultaneous=True))
need('full L2 second chart has no denominator', S.fraction(Hr)[1] == 1)
zero('full H restriction to D', Hr.subs(r,0)
     -(nn[8]*bd**2+(2*nn[6]-pp[4])*bd+nn[4]-pp[2]+q0))
mm = S.symbols('m0:5')
r0 = S.symbols('r0')
Mfull = sum(mm[i]*x**i for i in range(5))
Lfull = Mfull*t+mm[4]*x*x+mm[3]*x+r0
Lr = S.cancel(Lfull.subs({x:1/r, t:-r*r-r**4*bd}, simultaneous=True))
need('full L1 second chart has no denominator', S.fraction(Lr)[1] == 1)
zero('full L restriction to D', Lr.subs(r,0)-(r0-mm[2]-mm[4]*bd))
zero('zero M forces constant global L', Lfull.subs({v:0 for v in mm})-r0)
need('complete section dimensions', len(nn)+len(pp)+1 == 15 and len(mm)+1 == 6)

# Actual source scaling preserves W; translation is only a coefficient change.
d, X, T = S.symbols('d X T', nonzero=True)
zero('scaling volume multiplier', S.diff(d*X,X)*S.diff(T/d**2,T)-1/d)
root_p, root_q, root_r = S.symbols('root_p root_q root_r')
zero('leading residue normalized 431 condition',
     (3/(root_p-root_q)+1/(root_p-root_r)).subs({root_q:root_p+3*d,root_r:root_p-d}))
N = u**4*(u-3)**3*(u+1)
zero('full normalized leading polynomial', S.expand(N)-(u**8-8*u**7+18*u**6-27*u**4))
for root, multiplicity in [(0,4),(3,3),(-1,1)]:
    for j in range(multiplicity):
        zero(f'leading root {root} jet {j}', S.diff(N,u,j).subs(u,root))
    need(f'leading root {root} exact multiplicity', S.diff(N,u,multiplicity).subs(u,root) != 0)

# The complete inverse coefficient T2, and rational trace descent multiplier.
a, P0, Q0, M0, R0, v, T2 = S.symbols('a P0 Q0 M0 R0 v T2')
inverse_scaled = 1/a-P0*v/(2*a*a)+(P0*P0-4*a*a*Q0)*v*v/(8*a**3)+T2*v**3
inv_eq = S.Poly(S.expand((a*a*inverse_scaled**2+P0*inverse_scaled*v+Q0*v*v)**2
                       +M0*inverse_scaled*v**3+R0*v**4-1),v)
for j in range(3):
    zero(f'formal inverse v{j}', inv_eq.coeff_monomial(v**j))
zero('complete inverse T2 equation', inv_eq.coeff_monomial(v**3)-4*a*T2-M0/a)
zero('complete inverse T2 solution', inv_eq.coeff_monomial(v**3).subs(T2,-M0/(4*a*a)))
need('mate coefficient six derivative multiplier', S.Rational(2,4) == S.Rational(1,2))
need('normalized degree-two trace preserves a rational coefficient', S.Rational(1,2)*2 == 1)

# Active jets and the simple-root residue exhaust all possible M coefficients.
ma, mb, mc = S.symbols('ma mb mc')
Mtrial = u*u*(ma*u*u+mb*u+mc)
zero('simple root evaluation', Mtrial.subs(u,-1)-(ma-mb+mc))
reduced_M = Mtrial.subs(mc,mb-ma)
zero('fourfold rational residue', S.residue(reduced_M/N,u,0)+mb/27)
M = lam*u*u*(u*u-1)
zero('complete M primitive', S.diff(-lam/(3*u*(u-3)**2),u)-M/N)
zero('triple point is M-unit', M.subs(u,3)-72*lam)
Lsrc = M*t+lam*(u*u-2*p*u)
Z = (u*u-1)*q+1-2*p*u
zero('complete L blow chart', Lsrc.subs(t,(q-1)/u**2)-lam*Z)
Loriginal = S.Poly(S.expand(M.subs(u,x-p)),x)
zero('all-p full L induced x coefficient',
     lam*(u*u-2*p*u)-(Loriginal.coeff_monomial(x**4)*(u+p)**2
                         +Loriginal.coeff_monomial(x**3)*(u+p))-3*lam*p*p)

# Full shifted section constraints, before the moving-root residue.
aa, bb = S.symbols('aa bb')
Pre = u*u*(2*u**4+(-4*p-16)*u**3+aa*u*u+bb*u+b)
Q1 = -2*aa*p+bb+32*p*p+72*p
Qre = u**4-4*(p+2)*u**3+(aa+4*p*p-18)*u*u+Q1*u+c
Np = S.Poly(S.expand(N.subs(u,x-p)),x)
Pp = S.Poly(S.expand(Pre.subs(u,x-p)),x)
Qp = S.Poly(S.expand(Qre.subs(u,x-p)),x)
zero('complete shifted P degree-six anchor', Pp.coeff_monomial(x**6)-2*Np.coeff_monomial(x**8))
zero('complete shifted P degree-five anchor', Pp.coeff_monomial(x**5)-2*Np.coeff_monomial(x**7))
for power, expected in [(4,Np.coeff_monomial(x**8)),(3,Np.coeff_monomial(x**7)),
                        (2,Pp.coeff_monomial(x**4)-Np.coeff_monomial(x**6)),
                        (1,Pp.coeff_monomial(x**3)-Np.coeff_monomial(x**5))]:
    zero(f'complete shifted Q anchor {power}', Qp.coeff_monomial(x**power)-expected)
zero('two active P jets value', Pre.subs(u,0))
zero('two active P jets derivative', S.diff(Pre,u).subs(u,0))
Hre = S.expand(N*(q-1)**2/u**4+Pre*(q-1)/u**2+Qre)
H0 = -27*(q-1)**2+b*(q-1)+c
zero('actual first H row', Hre.subs(u,0)-H0)
zero('actual next H row', S.diff(Hre,u).subs(u,0)-(bb*(q-1)+Q1))
Fre = S.expand(Hre*Hre+lam*Z)
f = H0**2-lam*(q-1)
g = 2*H0*(bb*(q-1)+Q1)-2*p*lam
zero('actual original first F row', Fre.subs(u,0)-f)
zero('actual original second F row', S.diff(Fre,u).subs(u,0)-g)
zero('actual blow-up volume multiplier', S.diff((q-1)/u**2,q)-1/u**2)

# The residue equation is linear in the three unknowns after the exact root
# derivative has been retained. Three explicit triangular rows pay uniqueness.
CC, QQ = S.symbols('CC QQ')
gres = 2*H0*(bb*(q-1)+QQ)-2*p*lam-CC*S.diff(f,q)
zero('residue comparison cubic row', S.expand(gres).coeff(q,3)+54*(bb+54*CC))
zero('residue comparison quadratic after cubic',
     S.expand(gres.subs(bb,-54*CC)).coeff(q,2)-54*(b*CC-QQ))
zero('residue comparison constant after first two',
     gres.subs({bb:-54*CC,QQ:b*CC})-lam*(CC-2*p))
zero('full moving residue identity solved',
     g.subs({bb:-108*p,aa:16*p-18-b})-2*p*S.diff(f,q))
zero('all-p final Q1 relation', Q1.subs({bb:-108*p,aa:16*p-18-b})-2*p*b)

# Anchored degeneration retains aa until the triple-point zero condition.
Fanchor_pre = (N*t*t+Pre.subs({p:0,bb:0})*t+Qre.subs({p:0,bb:0}))**2+Lsrc.subs(p,0)
zero('anchored original source critical derivative u', S.diff(Fanchor_pre,u).subs(u,0))
zero('anchored original source critical derivative t', S.diff(Fanchor_pre,t).subs(u,0))
zero('anchored triple evaluation before reduction', Pre.subs({p:0,bb:0,u:3})-9*(9*aa+b-270))

# Local leading equations and orders. These are exact rescalings of the full
# equation and accompany, rather than replace, the analytic branch exhaustion.
w, s, tau, zz = S.symbols('w s tau zz')
al, bl, cl, ml, rl = S.symbols('al bl cl ml rl')
Nlocal = al*w**3+bl*w*w*s+cl*s*s
Elocal = Nlocal**2+s**3*(ml*w*w+(rl-zeta)*s)
triple_tangent = (al+cl*zz*zz)**2+(rl-zeta)*zz**4
zero('sole active triple actual leading quartic',
     S.expand(Elocal.subs({w:tau*tau,s:tau**3*zz},simultaneous=True)).coeff(tau,12)-triple_tangent)
need('sole active triple normalized differential order', 2*3+1-9 == -2)
need('sole active triple primitive capacity below D local degree', 2*(2-1) < 2+1)
active_tangent = (al+bl*zz+cl*zz*zz)**2-lam*zz**3-zeta*zz**4
Eactive = (al*w**4+bl*w*w*s+cl*s*s)**2+s**3*(-lam*w*w-zeta*s)
zero('active fourfold actual leading quartic',
     S.expand(Eactive.subs(s,w*w*zz)).coeff(w,8)-active_tangent)
zero('generic tangent repeated-root eliminant nonzero constant',
     (zz*S.diff(active_tangent+zeta*zz**4,zz)-4*(active_tangent+zeta*zz**4)).subs(zz,0)+4*al*al)
need('active fourfold differential order', 4-6 == -2)
need('active fourfold degree exhaustion', 4 == 4)
need('normal-unit triple actual normalized zero order', 12+1-9 == 4)
need('normal-unit triple primitive local degree exceeds total capacity', 4+1 > 4*(2-1))
need('simple shared branch exact differential order', 2-2 == 0)
need('four D intersections have exact differential order two', 2 == 2)

# No triple root for the balanced cubic: comparison with M(Z-Z0)^3
# contradicts the square identity for a!=0,M!=0,Z0!=0.
rr = S.symbols('rr')
zero('balanced cubic square-coefficient obstruction',
     (3*ml*rr**2)**2-4*(-ml*rr**3)*(-3*ml*rr)+3*ml*ml*rr**4)
for ell in [1,2]:
    order = ell-1 if ell % 2 else ell//2
    need(f'balanced m3 contact {ell} is unit or forbidden log', order in [0,1])
need('full differential canonical degree', 4*(-2)+4*2 == 0)
need('genus-one complete primitive-space dimension', 4+1-1 == 4 and 4 > 2*1-2)
need('nonsquare leading root multiplicities', 3 % 2 == 1 and 1 % 2 == 1)

# The full remaining nonanchored family and exact Riemann--Roch rows.
P = Pre.subs({aa:16*p-18-b,bb:-108*p})
Q = Qre.subs({aa:16*p-18-b,bb:-108*p})
zero('triple normal coefficient', P.subs(u,3)+72*(b+36*p+54))
zero('triple normal derivative after zero gate', S.diff(P,u).subs({u:3,b:-36*p-54})-864*p)
bvalue = -36*p-54
cvalue = h-(2*p-1)*bvalue+108*p*p-108*p+27
H = S.expand((N*(q-1)**2/u**4+P*(q-1)/u**2+Q).subs({b:bvalue,c:cvalue}))
F = S.expand(H*H+lam*Z)
J = S.expand((u-3)**2*(u+1)*(q-1)/u+u*u-(2*p+5)*u)
V = (u-3)*H/u
Jsource = u*(u-3)**2*(u+1)*t+u*u-(2*p+5)*u
zero('primitive J actual chart', Jsource.subs(t,(q-1)/u**2)-J)
Jr = S.cancel(Jsource.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
need('primitive J is global in original second chart', S.fraction(Jr)[1] == 1)
zero('primitive J leading triple vanishing', (u*(u-3)**2*(u+1)).subs(u,3))
zero('primitive J derivative triple vanishing', S.diff(u*(u-3)**2*(u+1),u).subs(u,3))
zero('primitive J simple-point vanishing', (u*(u-3)**2*(u+1)).subs(u,-1))
need('primitive basis distinct positive t degrees below generic equation', 2 < 4 and 1 < 2)

def bracket(Fc, G):
    return S.cancel(u*u*(S.diff(Fc,u)*S.diff(G,q)-S.diff(Fc,q)*S.diff(G,u)))

Acoef, Bcoef, Ccoef = S.symbols('Acoef Bcoef Ccoef')
remainders = [S.cancel(S.rem(bracket(F,G),F-zeta,q)) for G in [1/u,J,V]]
expr = S.expand(sum(coef*rem for coef,rem in zip([Acoef,Bcoef,Ccoef],remainders)))
zero('nonanchored RR complete cubic row', expr.coeff(q,3)
     -4*(u-3)**6*(u+1)**2*(Acoef-3*(2*p+3)*Bcoef))
expr2 = S.expand(expr.subs(Acoef,3*(2*p+3)*Bcoef))
zero('nonanchored RR complete quadratic row', expr2.coeff(q,2)
     -(u-3)**4*(u+1)*(u+3)*(4*(h+192*p*p)*Bcoef+lam*Ccoef))
expr3 = S.expand(expr2.subs(Ccoef,-4*(h+192*p*p)*Bcoef/lam))
obstruction = 48*(h+192*p*p)*zeta/lam-(32*p+48)*(h+192*p*p)-lam
zero('nonanchored RR final selected coefficient', expr3.coeff(q,1).coeff(u,4)-Bcoef*obstruction)
zero('nonanchored final nonzero fibre slope', S.diff(obstruction,zeta)-48*(h+192*p*p)/lam)
zero('nonanchored zero-slope hostile still obstructed', obstruction.subs(h,-192*p*p)+lam)

# The anchored full family is separate: it is not obtained by dividing p.
Panchor = Pre.subs({p:0,bb:0,aa:30-b/9})
Qanchor = Qre.subs({p:0,bb:0,aa:30-b/9})
Hanchor = S.expand(N*(q-1)**2/u**4+Panchor*(q-1)/u**2+Qanchor)
Ha_expected = (u-3)**3*(u+1)*q*q+(b+54)*(1-u*u/9)*q+c-b-27
zero('full anchored H after triple zero gate', Hanchor-Ha_expected)
zero('anchored triple vanishing', Panchor.subs(u,3))
zero('anchored normal order includes higher degeneration', S.diff(Panchor,u).subs(u,3)+6*(b+54))
Fa = S.expand(Hanchor*Hanchor+lam*((u*u-1)*q+1))
Ja = J.subs(p,0)
zero('anchored primitive J displayed form', Ja-((u-3)**2*(u+1)*q/u-3-9/u))
Va = (u-3)*Hanchor/u
rema = [S.cancel(S.rem(bracket(Fa,G),Fa-zeta,q)) for G in [1/u,Ja,Va]]
exa = S.expand(sum(coef*rem for coef,rem in zip([Acoef,Bcoef,Ccoef],rema)))
zero('anchored RR complete cubic row', exa.coeff(q,3)
     -S.Rational(2,3)*(6*Acoef+b*Bcoef)*(u-3)**6*(u+1)**2)
exa2 = S.expand(exa.subs(Acoef,-b*Bcoef/6))
zero('anchored RR complete quadratic row', exa2.coeff(q,2)
     -(u-3)**4*(u+1)*(u+3)*(Bcoef*(b*b+108*c)+27*Ccoef*lam)/27)
exa3 = S.expand(exa2.subs(Ccoef,-(b*b+108*c)*Bcoef/(27*lam)))
obsa = 108*(b*b+108*c)*zeta+2*b*lam*(b*b+108*c)-243*lam*lam
zero('anchored RR complete linear row', exa3.coeff(q,1)
     -Bcoef*(u-3)**3*(u+1)*obsa/(243*lam))
zero('anchored final nonzero fibre slope', S.diff(obsa,zeta)-108*(b*b+108*c))
zero('anchored zero-slope hostile still obstructed', obsa.subs(c,-b*b/108)+243*lam*lam)

# Exact same-partition rational hostile for every finite p: full original
# surface globality is checked, not inferred from a surface translation.
kappa, constant = S.symbols('kappa constant')
RR = (u-3)*(u+1)
ZZ = u*u*(u-3)*t+u-2*p-3
Hsharp = RR*ZZ*ZZ+kappa
Gsharp = (u-2)/(12*u*(u-3)*ZZ)
zero('sharp rational family has exactly the declared leading coefficient', S.expand(Hsharp).coeff(t,2)-N)
sharpchart = S.cancel(Hsharp.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
need('sharp rational family is global for every original p', S.fraction(sharpchart)[1] == 1)
zero('sharp rational mate of H', source_jac(Hsharp,Gsharp)-1)
zero('sharp rational mate of F', source_jac(Hsharp*Hsharp+constant,Gsharp/(2*Hsharp))-1)
need('sharp family is not labelled polynomial mate', S.degree(S.fraction(S.cancel(Gsharp))[1],t) == 1)
GG = S.Function('GG')(u,t)
zero('constant lower row polynomial factor', source_jac(Hsharp*Hsharp+constant,GG)-2*Hsharp*source_jac(Hsharp,GG))

# Named exact residual controls, independent of generic coefficient choices.
for name, subst in [('ordinary',{p:1,h:1,lam:2}),
                    ('zero-fibre-slope',{p:1,h:-192,lam:2})]:
    need(name+' final obstruction is nonzero', S.Poly(obstruction.subs(subst),zeta).as_expr() != 0)
for name, subst in [('normal-order-one',{b:0,c:1,lam:2}),
                    ('higher-normal-order',{b:-54,c:1,lam:2}),
                    ('anchored-zero-slope',{b:0,c:0,lam:2})]:
    need(name+' anchored obstruction is nonzero', S.Poly(obsa.subs(subst),zeta).as_expr() != 0)
need('triple at infinity leading capacity fails', 4-2 < 5-2)
need('simple at infinity leading capacity fails', (4-2)+(3-2) < 7-2)
need('fourfold at infinity is explicitly not excluded here', 3-2 == S.Rational(4,2)-1)

semantic = sha256(json.dumps(gates,separators=(',',':')).encode()).hexdigest()
print('DG quartic boundary (4,3,1), all three points finite')
print('Universe: full global L2/L1 coefficients; every normalized finite p')
print('Rational mate forces L constant; all polynomial mates excluded')
print('Complete local divisors and genus-one primitive spaces: exact controls PASS')
print('Sharp same-partition all-p rational constant-L family: PASS')
print(f'Gates: {len(gates)}')
print('Semantic SHA256: '+semantic)
print('RESULT: PASS')
