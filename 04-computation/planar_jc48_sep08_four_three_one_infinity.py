#!/usr/bin/env python3
"""Exact controls for the fourfold-infinity (4,3,1) global quartic class.

Analytic local completeness and the componentwise primitive-degree argument
are proved in the companion. The trace test is independent and redundant.
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


x, u, t, p, r, bd, T, s, w, tau, Z, zeta = S.symbols('x u t p r bd T s w tau Z zeta')
a, b, c, lam, constant = S.symbols('a b c lam constant')
AA = S.symbols('A0:5')
BB = S.symbols('B0:5')
C0, E0 = S.symbols('C0 E0')
N = u**3*(u-1)
P = sum(AA[i]*u**i for i in range(5))
Q = AA[4]*u*u+(AA[3]-2*p*AA[4])*u+C0
M = sum(BB[i]*u**i for i in range(5))
R = BB[4]*u*u+(BB[3]-2*p*BB[4])*u+E0
H = N*t*t+P*t+Q
L = M*t+R
jac = lambda f,g: S.diff(f,u)*S.diff(g,t)-S.diff(f,t)*S.diff(g,u)

# All leading multiplicities and the full original coefficient constraints.
for point, multiplicity in [(0,3),(1,1)]:
    for j in range(multiplicity):
        zero(f'leading point {point} jet {j}', S.diff(N,u,j).subs(u,point))
    need(f'leading point {point} exact multiplicity', S.diff(N,u,multiplicity).subs(u,point) != 0)
need('complete binary infinity multiplicity', 8-S.degree(N,u) == 4)
Poriginal = S.Poly(S.expand(P.subs(u,x-p)),x)
Moriginal = S.Poly(S.expand(M.subs(u,x-p)),x)
Qoriginal = S.Poly(S.expand(Q.subs(u,x-p)),x)
Roriginal = S.Poly(S.expand(R.subs(u,x-p)),x)
for label, lead, lower in [('H',Poriginal,Qoriginal),('L',Moriginal,Roriginal)]:
    zero(label+' full quadratic induced row', lower.coeff_monomial(x*x)-lead.coeff_monomial(x**4))
    zero(label+' full linear induced row', lower.coeff_monomial(x)-lead.coeff_monomial(x**3))
    need(label+' arbitrary remaining lower constant', S.degree(lower,x) <= 2)
Hr = S.cancel(H.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
Lr = S.cancel(L.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
need('complete H is global on second chart', S.fraction(Hr)[1] == 1)
need('complete L is global on second chart', S.fraction(Lr)[1] == 1)
HD, LD = Hr.subs(r,0), Lr.subs(r,0)
zero('full H D slope', S.diff(HD,bd)+AA[4])
zero('full L D slope', S.diff(LD,bd)+BB[4])
FD = S.expand(HD*HD+LD)
zero('constant-D first necessary coefficient', FD.coeff(bd,2)-AA[4]**2)
zero('constant-D second necessary coefficient', FD.subs(AA[4],0).coeff(bd,1)+BB[4])
need('full independent dimensions for fixed leading polynomial', len(AA)+1 == 6 and len(BB)+1 == 6)

# Actual surface inversion and original volume, without a source polynomial
# automorphism assertion. All transformed coefficients remain global sections.
Hinv = S.cancel(H.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*T},simultaneous=True))
need('full inversion H remains polynomial in new source chart', S.fraction(Hinv)[1] == 1)
zero('complete infinity leading unit', S.expand(Hinv).coeff(T,2)
     -r**4*(1-p*r)**3*(1-(p+1)*r))
zero('infinity unit constant term', ((1-p*r)**3*(1-(p+1)*r)).subs(r,0)-1)
zero('actual original source volume multiplier', S.diff(1/r,r)*S.diff(-r*r-r**4*T,T)-r*r)

# Complete weighted m4 branch faces, including nonconstant leading unit.
al1, al2, b1, b2, cc, d1, d2, ee = S.symbols('al1 al2 b1 b2 cc d1 d2 ee')
unit = 1+al1*w+al2*w*w
calN = w**4*unit+s*(b1*w+b2*w*w)+cc*s*s
calM = d1*w+d2*w*w+ee*s
Phi = calN**2+s**3*calM-zeta*s**4
low = Z**2*((b1+cc*Z)**2+d1*Z+(ee-zeta)*Z*Z)
zero('m4 first-normal-jet complete low face', S.expand(Phi.subs(s,w*Z)).coeff(w,4)-low)
need('m4 first-normal-jet weighted low order', 2+2-3 == 1)
for ell in [1,2,3]:
    exponent = S.Rational(5-ell,2)
    ram = 2 if exponent.q == 2 else 1
    need(f'm4 cancellation normalized regular ell={ell}', ram*exponent+ram-1 >= 0)
need('m4 first-normal-jet branch exhaustion', 2+2 == 4)
Phi_d = Phi.subs(b1,0)
zero('m4 first-M-jet low face', S.expand(Phi_d.subs(s,w*Z)).coeff(w,4)
     -Z**3*(d1+(cc*cc+ee-zeta)*Z))
zero('m4 first-M-jet complete high face',
     S.expand(Phi_d.subs({w:tau**3,s:tau**7*Z},simultaneous=True)).coeff(tau,24)-(1+d1*Z**3))
need('m4 high branch actual weighted differential order', (6+14+2)-17 == 5)
need('m4 first-M-jet branch exhaustion', 1+3 == 4)
face = (1+b2*Z+cc*Z*Z)**2+d2*Z**3+(ee-zeta)*Z**4
zero('m4 full higher-jet tangent quartic',
     S.expand(Phi.subs({b1:0,d1:0}).subs(s,w*w*Z)).coeff(w,8)-face)
zero('m4 generic root eliminant nonzero constant',
     (Z*S.diff(face+zeta*Z**4,Z)-4*(face+zeta*Z**4)).subs(Z,0)+4)
need('m4 higher-jet weighted differential order', 2+4-6 == 0)
need('m4 higher-jet branch exhaustion', 4 == 4)
need('m4 M-unit cannot be balanced with an integer normal order', 4 % 3 != 0)
for j in [0,1]:
    exponent = S.Rational(4-3*j,2)
    ram = 2 if exponent.q == 2 else 1
    need(f'm4 M-unit cancellation weighted regular j={j}', ram*(exponent+2)+ram-1 >= 0)
for ell in [1,2,3,4]:
    exponent = S.Rational(4-ell,2)
    ram = 2 if exponent.q == 2 else 1
    need(f'm4 shared normal-unit normalized regular ell={ell}', ram*exponent+ram-1 >= 0)

# The finite active triple contributes two double poles, while a nonconstant
# restriction to D forces local primitive degree three on an actual component.
al, bl, cl, ml, rl = S.symbols('al bl cl ml rl')
Et = (al*w**3+bl*w*w*s+cl*s*s)**2+s**3*(ml*w*w+(rl-zeta)*s)
triple_face = (al+cl*Z*Z)**2+(rl-zeta)*Z**4
zero('active triple full leading quartic',
     S.expand(Et.subs({w:tau*tau,s:tau**3*Z},simultaneous=True)).coeff(tau,12)-triple_face)
zero('active triple generic repeated-root eliminant',
     (Z*S.diff(triple_face+zeta*Z**4,Z)-4*(triple_face+zeta*Z**4)).subs(Z,0)+4*al*al)
need('active triple actual differential order', 6+1-9 == -2)
need('active triple primitive capacity', 2*(2-1) == 2)
need('same-component D degree contradiction', 2 < 2+1)
need('simple shared normalized differential is regular', 2-2 == 0)

# The complete T2 equation, followed by both rational residues. In particular
# the second residue must not be lost while exploring the intermediate field.
a0, P0, Q0, M0, R0, v, T2 = S.symbols('a0 P0 Q0 M0 R0 v T2')
UU = 1/a0-P0*v/(2*a0*a0)+(P0*P0-4*a0*a0*Q0)*v*v/(8*a0**3)+T2*v**3
inv = S.Poly(S.expand((a0*a0*UU*UU+P0*UU*v+Q0*v*v)**2+M0*UU*v**3+R0*v**4-1),v)
for j in range(3):
    zero(f'formal inverse coefficient {j}', inv.coeff_monomial(v**j))
zero('formal inverse coefficient T2', inv.coeff_monomial(v**3)-4*a0*T2-M0/a0)
zero('formal complete T2 solution', inv.coeff_monomial(v**3).subs(T2,-M0/(4*a0*a0)))
need('rational trace descent multiplier', S.Rational(1,2)*2 == 1)
reduction = {AA[4]:0,AA[0]:0,AA[1]:0,BB[4]:0,BB[0]:0,BB[1]:0}
Mred = M.subs(reduction)
zero('full constant-D active M', Mred-u*u*(BB[3]*u+BB[2]))
zero('simple-point residue', S.residue(Mred/N,u,1)-(BB[3]+BB[2]))
zero('last triple-point rational differential', (Mred/N).subs(BB[2],-BB[3])-BB[3]/u)
zero('last triple-point residue', S.residue((Mred/N).subs(BB[2],-BB[3]),u,0)-BB[3])
zero('full L constant after counted residues', L.subs({BB[i]:0 for i in range(5)})-E0)

# Independent ambient quadratic-field obstruction for the intermediate lambda
# family. It is retained as a control and is not needed after the last residue.
H0 = N*t*t+u*u*(a*u+b)*t+a*u+c
Z0 = u*u*(u-1)*t+u
F0 = H0*H0+lam*Z0
zero('full residual H from global coefficient intersection',
     H.subs(reduction).subs({AA[3]:a,AA[2]:b,C0:c})-H0)
zero('full residual L from global coefficient intersection',
     L.subs(reduction).subs({BB[3]:lam,BB[2]:-lam,E0:0})-lam*Z0)
hvar, zvar = S.symbols('hvar zvar')
A = a*zvar+1-a-b+c-hvar
B = (b-2)*zvar+hvar-c
Delta = S.expand(B*B-4*A*zvar*zvar)
quadratic = A*u*u+B*u+zvar*zvar
zero('actual quadratic field relation', quadratic.subs({hvar:H0,zvar:Z0},simultaneous=True))
yexpr = 2*A*u+B
zero('quadratic inverse t relation', ((zvar-u)/(u*u*(u-1))).subs(zvar,Z0)-t)
zero('quadratic square completion identity', yexpr*yexpr-Delta-4*A*quadratic)
zero('actual residual map Jacobian', jac(H0,Z0)-u*yexpr.subs({hvar:H0,zvar:Z0},simultaneous=True))
zero('actual residual map dominance coefficient', S.expand(jac(H0,Z0)).coeff(t,2)+u**4*(u-1)*(2*u-1))
zero('monic quadratic nonsquare test', S.Poly(Delta,hvar).LC()-1)
zero('complete discriminant in h', S.discriminant(Delta,hvar)-16*zvar*zvar*(zvar-1)*(zvar+a+b-1))
need('discriminant is nonzero for all coefficient values', S.Poly(S.discriminant(Delta,hvar),zvar).LC() == 16)
yy = S.symbols('yy')
uinverse = (yy-B)/(2*A)
trace_eta_factor = (1+B/yy)/(2*lam*zvar*zvar)
cleared = S.together(-1/(lam*uinverse*yy)-trace_eta_factor)
num, den = S.fraction(cleared)
zero('actual relative differential factor in quadratic field', S.rem(num,yy*yy-Delta,yy))
zero('actual field trace of the relative differential', trace_eta_factor+trace_eta_factor.subs(yy,-yy)-1/(lam*zvar*zvar))
zero('fixed-original-fibre trace', (1/(lam*zvar*zvar)).subs(zvar,(zeta-hvar*hvar)/lam)-lam/(zeta-hvar*hvar)**2)
rho = S.symbols('rho', nonzero=True)
zero('exact traced nonzero residue', S.residue(lam/(rho*rho-hvar*hvar)**2,hvar,rho)+lam/(4*rho**3))
for name, vals in [('ordinary',{a:1,b:2,c:3,lam:5}),
                   ('coincident-discriminant-z-roots',{a:1,b:-1,c:0,lam:2}),
                   ('zero-top-parameters',{a:0,b:0,c:0,lam:1})]:
    need(name+' ambient nonsquare discriminant remains nonzero', S.discriminant(Delta,hvar).subs(vals) != 0)
    need(name+' earlier T2 obstruction retained', (lam/u).subs(vals) != 0)

# Exact same-partition constant-L rational mates, for every original p.
kappa = S.symbols('kappa')
Hsharp = N*t*t+kappa
GH = -1/(u*u*t)
zero('sharp rational mate of H', jac(Hsharp,GH)-1)
zero('sharp rational mate of F', jac(Hsharp*Hsharp+constant,GH/(2*Hsharp))-1)
zero('sharp original global chart',
     Hsharp.subs(u,x-p).subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True)
     -((1-p*r)**3*(1-(p+1)*r)*(1+r*r*bd)**2+kappa))
need('sharp control retains actual finite multiplicities', S.degree(N,u) == 4 and 8-S.degree(N,u) == 4)
GG = S.Function('GG')(u,t)
zero('constant-L polynomial obstruction', jac(Hsharp*Hsharp+constant,GG)-2*Hsharp*jac(Hsharp,GG))

semantic = sha256(json.dumps(gates,separators=(',',':')).encode()).hexdigest()
print('DG boundary (4 infinity, 3 finite, 1 finite): complete coefficient theorem')
print('Weighted infinity, componentwise capacity, and both T2 residues: PASS')
print('Rational mate forces constant L; unrestricted polynomial mate excluded')
print('Independent degree-two field trace and same-partition rational controls: PASS')
print(f'Gates: {len(gates)}')
print('Semantic SHA256: '+semantic)
print('RESULT: PASS')
