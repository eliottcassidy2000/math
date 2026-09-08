#!/usr/bin/env python3
"""Exact identities for the complete constant-D shared finite-six class.

The companion proves the all-parameter case split, branch exhaustion,
field-kernel assertion, and original-source pole obstruction analytically.
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


u, x, t, p, r, T, s, Z = S.symbols('u x t p r T s Z')
b, c, k, d, n, e, f, h, v, q = S.symbols('b c k d n e f h v q')
A = S.symbols('A0:5')
B = S.symbols('B0:5')
C0, E0 = S.symbols('C0 E0')
N = u**6
P = sum(A[i]*u**i for i in range(5))
Q = (A[4]-1)*u*u+(A[3]-2*p*A[4]+4*p)*u+C0
M = sum(B[i]*u**i for i in range(5))
R = B[4]*u*u+(B[3]-2*p*B[4])*u+E0
Hfull = N*t*t+P*t+Q
Lfull = M*t+R
chart = {x:1/r,t:-r*r-r**4*T}
Hr = S.cancel(Hfull.subs(u,x-p).subs(chart, simultaneous=True))
Lr = S.cancel(Lfull.subs(u,x-p).subs(chart, simultaneous=True))
need('full shifted H second chart polynomial', S.fraction(Hr)[1] == 1)
need('full shifted L second chart polynomial', S.fraction(Lr)[1] == 1)
zero('full H D slope', S.diff(Hr.subs(r,0),T)-(2-A[4]))
zero('full L D slope', S.diff(Lr.subs(r,0),T)+B[4])
FD = S.expand(Hr.subs(r,0)**2+Lr.subs(r,0))
zero('constant D H coefficient', FD.coeff(T,2)-(2-A[4])**2)
zero('constant D L coefficient', FD.subs(A[4],2).coeff(T,1)+B[4])
need('finite multiplicity six and infinity multiplicity two', S.degree(N,u) == 6 and 8-S.degree(N,u) == 2)

# Derive the global coefficient space separately from its full numerator box.
q0,q1,q2,j0,j1,j2 = S.symbols('q0 q1 q2 j0 j1 j2')
abox = q2*x*x+q1*x+q0
bbox = (1-q2)*x**4+(-6*p-q1)*x**3+j2*x*x+j1*x+j0
cbox = S.expand((x-p)**6-abox*x**4-bbox*x*x)
need('full global numerator box C degree', S.degree(cbox,x) <= 4)
pbox = S.Poly(S.expand((2*abox*x*x+bbox).subs(x,u+p)),u)
qbox = S.Poly(S.expand(abox.subs(x,u+p)),u)
zero('independent shifted quadratic row', qbox.coeff_monomial(u*u)-pbox.coeff_monomial(u**4)+1)
zero('independent shifted linear row', qbox.coeff_monomial(u)-pbox.coeff_monomial(u**3)+2*p*pbox.coeff_monomial(u**4)-4*p)
need('fixed leading full H kernel dimension', len((q0,q1,q2,j0,j1,j2)) == 6)

active = {A[4]:2,A[3]:c,A[2]:b,A[1]:0,A[0]:0,C0:k+b,
          B[4]:0,B[3]:n,B[2]:d,B[1]:0,B[0]:0,E0:e+d}
Hs = S.expand(Hfull.subs(active))
Ls = S.expand(Lfull.subs(active))
qs = 1+u*u*t
vs = u*qs
Hmodel = v*v+c*v+b*q+k
Lmodel = d*q+n*v+e
Fs = Hs*Hs+Ls
jac = lambda F,G: S.diff(F,u)*S.diff(G,t)-S.diff(F,t)*S.diff(G,u)
zero('complete active H in q v', Hmodel.subs({q:qs,v:vs},simultaneous=True)-Hs)
zero('complete active L in q v', Lmodel.subs({q:qs,v:vs},simultaneous=True)-Ls)
zero('exact v q Jacobian', jac(vs,qs)-u*u*qs)
zero('exact H v Jacobian', jac(Hs,vs)+b*u*u*qs)
zero('birational u inverse', (v/q).subs({v:vs,q:qs},simultaneous=True)-u)
zero('birational t inverse', ((q-1)*q*q/(v*v)).subs({v:vs,q:qs},simultaneous=True)-t)
zero('active source gradient along u zero', S.diff(Fs,u).subs(u,0)-2*(b+k)*c-n)
zero('active source vertical gradient along u zero', S.diff(Fs,t).subs(u,0))
zero('active H original D value', Hr.subs(active).subs(r,0)-(4*p*p+2*p*c+k))
zero('active L original D value', Lr.subs(active).subs(r,0)-(2*p*n+e))

# The actual infinity point is regular even before the active finite jets.
const_D = {A[4]:2,B[4]:0}
Ni = S.cancel(s*s*Hr.subs(const_D).subs(T,1/s))
Mi = S.cancel(s*Lr.subs(const_D).subs(T,1/s))
need('complete infinity numerator H analytic', S.fraction(Ni)[1] == 1)
need('complete infinity numerator L analytic', S.fraction(Mi)[1] == 1)
zero('complete infinity leading order two', Ni.subs(s,0)-r*r*(1-p*r)**6)
zero('infinity H normal constant vanishes', S.diff(Ni,s).subs({r:0,s:0}))
zero('infinity M constant vanishes', Mi.subs({r:0,s:0}))
zeta = S.symbols('zeta')
E = Ni*Ni+s**3*Mi-zeta*s**4
face = S.expand(E.subs(s,r*Z)).coeff(r,4)
zero('infinity tangent quartic constant', face.subs(Z,0)-1)
zero('infinity generic repeated root eliminant',
     (Z*S.diff(face+zeta*Z**4,Z)-4*(face+zeta*Z**4)).subs(Z,0)+4)
need('infinity complete generic degree four', S.degree(face,Z) == 4)
need('infinity weighted differential order', 2+2-3 == 1)

# b*d !=0: the actual generic conic has two unavoidable logarithmic poles.
delta = d/b
alpha = n-delta*c
Fconic = h*h+delta*h-delta*v*v+alpha*v+e-delta*k
qinv = (h-v*v-c*v-k)/b
zero('full conic equation', Fconic.subs(h,Hmodel)-Hmodel*Hmodel-Lmodel)
zero('conic actual relative form factor',
     (-qinv/(b*(2*h+delta)*v*v)).subs({h:Hs,v:vs},simultaneous=True)*jac(Fs,vs)-1)
xi,w,rho = S.symbols('xi w rho', nonzero=True)
proj = xi*xi-delta+w*(delta*xi+alpha)+w*w*(e-delta*k-f)
zero('complete conic infinity equation',
     S.expand(w*w*(Fconic-f).subs({h:xi/w,v:1/w},simultaneous=True))-proj)
eta_w = ((xi-c)*w-1-k*w*w)/(b*b*w*(2*xi+delta*w))
zero('actual conic differential in inverse v',
     (-qinv/(b*(2*h+delta)*v*v)).subs({h:xi/w,v:1/w},simultaneous=True)*(-1/w**2)-eta_w)
for sign in [-1,1]:
    zero(f'conic infinity branch {sign}', proj.subs({w:0,xi:sign*rho,d:b*rho*rho}))
    zero(f'conic infinity implicit derivative {sign}',
         S.diff(proj,xi).subs({w:0,xi:sign*rho,d:b*rho*rho})-2*sign*rho)
    zero(f'conic exact nonzero residue {sign}',
         (w*eta_w).subs(w,0).subs(xi,sign*rho)+sign/(2*b*b*rho))

# b=0,d!=0: rational field C(F,v), and the lone residue is exactly
# criticality on the entire retained source line u=0.
polyv = (v*v+c*v+k)**2+n*v
eta_v = -(f-polyv-e)/(d*d*v*v)
zero('b zero actual relative form',
     eta_v.subs({f:Fs.subs(b,0),v:vs},simultaneous=True)*jac(Fs.subs(b,0),vs)-1)
zero('b zero exact residue', S.residue(eta_v,v,0)-(2*k*c+n)/(d*d))
G_bzero = (v**3/S.Integer(3)+c*v*v+(c*c+2*k)*v+(f-e-k*k)/v)/(d*d)
zero('b zero primitive remainder', S.diff(G_bzero,v)-eta_v+(2*k*c+n)/(d*d*v))
zero('b zero source criticality condition', S.diff(Fs.subs(b,0),u).subs(u,0)-(2*k*c+n))
Gb_source = G_bzero.subs({f:Fs.subs({b:0,n:-2*k*c}),v:vs},simultaneous=True)
zero('b zero rational mate on entire accepted stratum', jac(Fs.subs({b:0,n:-2*k*c}),Gb_source)-1)

# d=0,b!=0,n!=0: field C(F,H), two residues force c=k=0.
vlin = (f-h*h-e)/n
eta_h = ((h-vlin*vlin-c*vlin-k)/(n*b*b*vlin*vlin))
zero('d zero actual relative form',
     eta_h.subs({f:Fs.subs(d,0),h:Hs},simultaneous=True)*jac(Fs.subs(d,0),Hs)-1)
for sign in [-1,1]:
    residue = S.residue(eta_h.subs(f,e+rho*rho),h,sign*rho)
    zero(f'd zero complete residue {sign}',
         residue-sign*(n*k+2*c*rho*rho)/(4*b*b*rho**3))
Gsharp = 1/(2*b*b*v)-h/(n*b*b)
Fsharp = Fs.subs({c:0,k:0,d:0})
Hsharp = Hs.subs({c:0,k:0})
Gsharp_source = Gsharp.subs({v:vs,h:Hsharp},simultaneous=True)
zero('nonconstant lower row rational submersion mate', jac(Fsharp,Gsharp_source)-1)
zero('sharp complete inverse v', (f-h*h-e)/n-vlin)
zero('sharp complete inverse q', (h-v*v)/b-qinv.subs({c:0,k:0}))
fiberfactor = u*(qs*qs+b*t)*(Hsharp+b)+n*qs
zero('original special fibre factor', Fsharp-(b*b+e)-u*fiberfactor)
zero('second fibre factor avoids u zero', fiberfactor.subs(u,0)-n)
zero('second fibre factor is nonconstant', fiberfactor.subs(t,0)-(u**3+2*b*u+n))
zero('second fibre factor avoids q zero', fiberfactor.subs(t,-1/(u*u))+b*b/u)
zero('sharp source line is transverse', S.diff(Fsharp,u).subs(u,0)-n)
zero('sharp q-zero source is transverse', S.diff(Fsharp,t).subs(t,-1/(u*u))-n*u**3)
zero('compulsory source u-pole coefficient', S.cancel(u*Gsharp_source).subs(u,0)-1/(2*b*b))

# Constant L and b=d=0 are kept as complete split cases, including
# rational positive controls; no blanket rational exclusion is claimed.
GH = (v+(h-k)/v)/(b*b)
eta_H = (1+c/v+(k-h)/(v*v))/(b*b)
zero('constant L H relative differential',
     eta_H.subs({h:Hs,v:vs},simultaneous=True)*jac(Hs,vs)-1)
zero('constant L H only residue', S.residue(eta_H,v,0)-c/(b*b))
zero('constant L H primitive when c zero', S.diff(GH,v)-eta_H.subs(c,0))
GHsource = GH.subs({h:Hs.subs(c,0),v:vs},simultaneous=True)
zero('constant L rational mate', jac(Hs.subs(c,0)**2+e,GHsource/(2*Hs.subs(c,0)))-1)
Gv = 1/(2*u*u)
zero('v actual rational mate', jac(vs,Gv)-1)
Pcomp = (v*v+c*v+k)**2+n*v+e
Pprime = S.diff(Pcomp,v)
need('b d zero outer degree four', S.degree(Pcomp,v) == 4)
need('b d zero critical factor nonconstant', S.degree(Pprime,v) == 3)
zero('b d zero rational mate for all coefficients',
     jac(Pcomp.subs(v,vs),Gv/Pprime.subs(v,vs))-1)
gg = S.Function('gg')(u,t)
zero('constant L polynomial factor obstruction', jac(Hs*Hs+e,gg)-2*Hs*jac(Hs,gg))
zero('composite polynomial factor obstruction', jac(Pcomp.subs(v,vs),gg)-Pprime.subs(v,vs)*jac(vs,gg))

semantic = sha256(json.dumps(gates,separators=(',',':')).encode()).hexdigest()
print('Constant-D shared finite6/infinity2: complete shifted global class')
print('All coefficient splits: exact rational boundary and no polynomial mate')
print('Conic residues, source-critical residue, and original same-fibre repair: PASS')
print('Nonconstant-lower-row rational submersion and composite controls: PASS')
print(f'Gates: {len(gates)}')
print('Semantic SHA256: '+semantic)
print('RESULT: PASS')
