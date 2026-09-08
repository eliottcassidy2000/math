#!/usr/bin/env python3
"""Exact family intersection for translated finite-six constant-D quartics.

No surface automorphism or prior mathematical implementation is imported.
The rational-mate exclusion is an analytic consumer of the fixed-zero proof.
"""
import hashlib
import json
import sympy as S

w, x, z, s, t, p, r, bd = S.symbols('w x z s t p r bd')
a, b, c0, c1, c2, beta0, beta1 = S.symbols('a b c0 c1 c2 beta0 beta1')
m0, m1, d, n0, n1, n2 = S.symbols('m0 m1 d n0 n1 n2')
gates = 0
records = {}

def check(ok, label):
    global gates
    gates += 1
    if not ok:
        raise RuntimeError(label)

def eq(left, right, label):
    check(S.cancel(left-right) == 0, label)

def coefficient(poly, var, degree):
    return S.expand(poly).coeff(var, degree)

def jac(left, right, u, v):
    return S.diff(left,u)*S.diff(right,v)-S.diff(left,v)*S.diff(right,u)

graph = (w+p)**2
C = c0+c1*w+c2*w*w
B = beta0+beta1*w+b*w*w+(c1+2*p*(c2-a))*w**3+(a+c2)*w**4
Ngeneral = a*w**6+s*B+s*s*C
Qg = S.expand(Ngeneral.subs(s,z-graph))
check(S.Poly(Qg,w,z).degree(w)<=4, 'full preconstant N x box')
check(S.Poly(Qg,w,z).degree(z)<=2, 'full preconstant N z box')
eq(Qg.subs(z,graph),a*w**6,'prescribed translated octic')
Hg = S.expand(Ngeneral.subs(s,1/t)*t*t)
Hginf = S.cancel(Hg.subs({w:1/r-p,t:-r*r-r**4*bd},simultaneous=True))
check(Hginf.is_polynomial(r,bd),'preconstant actual infinity regularity')
eq(coefficient(Hginf.subs(r,0),bd,1),a-c2,'constant D forces c2=a')

N = S.expand(Ngeneral.subs(c2,a))
N0 = a*w**6+s*(beta0+beta1*w+b*w*w+c1*w**3+2*a*w**4)+s*s*(c0+c1*w+a*w*w)
eq(N,N0,'whole final numerator independent of point p')
Q = S.expand(N.subs(s,z-graph))
Mgeneral = m0+m1*w+d*w*w+(n1+2*p*n2)*w**3+n2*w**4+s*(n0+n1*w+n2*w*w)
Mg = S.expand(Mgeneral.subs(s,z-graph))
check(S.Poly(Mg,w,z).degree(w)<=2,'full preconstant M x box')
check(S.Poly(Mg,w,z).degree(z)<=1,'full preconstant M z box')
Lg = S.expand(t*Mgeneral.subs(s,1/t))
Lginf = S.cancel(Lg.subs({w:1/r-p,t:-r*r-r**4*bd},simultaneous=True))
check(Lginf.is_polynomial(r,bd),'preconstant actual L infinity regularity')
eq(coefficient(Lginf.subs(r,0),bd,1),-n2,'constant L on D forces n2=0')
M = S.expand(Mgeneral.subs(n2,0))
eq(M,m0+m1*w+d*w*w+n1*w**3+s*(n0+n1*w),'whole final M independent of p')
eq(M.subs({w:0,s:0}),m0,'translated M-unit value')
Qm = S.expand(M.subs(s,z-graph))

# Exact ranks for every complex p, paid by constant nonzero minors.
mon = [w**i*z**j for i in range(5) for j in range(3)]
rows = [[coefficient(mm.subs(z,graph),w,k) for mm in mon]
        for k in range(9) if k!=6]
rows.append([coefficient(mm.subs(z,graph),w,6)-
             coefficient(coefficient(mm,z,2),w,2) for mm in mon])
constraints = S.Matrix(rows)
cols = constraints.subs(p,0).rref()[1]
check(len(cols)==9,'nine independent constraints at zero')
minor = S.factor(constraints[:,list(cols)].det())
check(minor in (-1,1),'same constraint minor nonzero at every p')
params = [a,c1,c0,beta0,beta1,b]
columns = S.Matrix([[coefficient(coefficient(S.diff(Q,par),w,i),z,j)
                     for par in params] for i in range(5) for j in range(3)])
check(constraints*columns==S.zeros(9,6),'complete translated basis satisfies constraints')
selected = [mon.index(mm) for mm in (w*w*z*z,w*z*z,z*z,z,w*z,w*w*z)]
eq(columns[selected,:].det(),1,'complete six-parameter basis for every p')
check(constraints.cols-constraints.rows==6,'no additional N parameter')
mparams = [n0,n1,m0,m1,d]
mmon = [w**i*z**j for i in range(3) for j in range(2)]
mcols = S.Matrix([[coefficient(coefficient(S.diff(Qm,par),w,i),z,j)
                   for par in mparams] for i in range(3) for j in range(2)])
eq(coefficient(coefficient(Qm,w,2),z,1),0,'complete M constant-D constraint')
mselected = [mmon.index(mm) for mm in (z,w*z,1,w,w*w)]
eq(mcols[mselected,:].det(),1,'complete five-parameter M basis for every p')

H0 = S.expand(t*t*N.subs(s,1/t))
L0 = S.expand(t*M.subs(s,1/t))
F0 = S.expand(H0*H0+L0)
Hp = H0.subs(w,x-p)
Lp = L0.subs(w,x-p)
eq(Hp.subs(x,w+p),H0,'exact source H translation')
eq(Lp.subs(x,w+p),L0,'exact source L translation')
eq(jac(w+p,t,w,t),1,'source translation preserves volume form')
Hinf = S.cancel(Hp.subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
Linf = S.cancel(Lp.subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
check(Hinf.is_polynomial(r,bd) and Linf.is_polynomial(r,bd),'all translated functions global')
eq(Hinf.subs(r,0),c0-b+2*p*c1+4*a*p*p,'actual translated H boundary constant')
eq(Linf.subs(r,0),n0-d+2*p*n1,'actual translated L boundary constant')
eq(Hinf.subs({r:0,p:0}),c0-b,'auxiliary fixed-zero H boundary')
eq(Linf.subs({r:0,p:0}),n0-d,'auxiliary fixed-zero L boundary')
eq(H0.coeff(t,2),a*w**6,'exact source degree two preserved')
eq(coefficient(F0,t,4),a*a*w**12,'exact quartic source layer preserved')

# A direct rational chain-rule control retains an unrestricted rational mate.
testG=(x*x+t)/(x+t+1)
testF=S.expand(Hp*Hp+Lp)
eq(jac(testF,testG,x,t).subs(x,w+p),
   jac(F0,testG.subs(x,w+p),w,t),'literal rational Jacobian transport')

# This plane translation does NOT extend to W for p nonzero.
boundary_coordinate = -(x-p)**2-(x-p)**4*t
translated_bd = S.cancel(boundary_coordinate.subs({x:1/r,t:-r*r-r**4*bd},simultaneous=True))
eq(translated_bd,-2*p/r+5*p*p-4*p**3*r+p**4*r*r+(1-p*r)**4*bd,
   'exact nonextendable boundary coordinate')
eq(S.cancel(r*translated_bd).subs(r,0),-2*p,'unremovable simple boundary pole when p nonzero')
eq(translated_bd.subs(p,0),bd,'identity point control')

h=x*x+x**4*t
eq(jac(h**4+h,1/(3*x**3*(4*h**3+1)),x,t),1,
   'genuine rational mate outside the six-plus-two boundary pattern')

# Symbolic identities retain lower jets; named controls also exercise them.
for lower in (beta0,beta1):
    eq(Hinf.subs({p:1,lower:1}).subs(r,0),c0-b+2*c1+4*a,
       'named nonzero lower normal jet '+str(lower))

# The second-kind subfamily is transported at real and nonreal points.
for pp in (S.Integer(0),S.Integer(1),S.Rational(-3,2),S.I):
    data={a:S.Rational(-4,27),b:1,c0:0,c1:0,beta0:0,beta1:0,
          m0:-1,m1:0,d:S.Rational(8,9),n0:0,n1:0,p:pp}
    eq(Hinf.subs(data).subs(r,0),-1-S.Rational(16,27)*pp**2,
       'named translated second-kind H boundary '+str(pp))
    eq(Linf.subs(data).subs(r,0),-S.Rational(8,9),
       'named translated second-kind L boundary '+str(pp))
    check(S.Poly(Q.subs(data),w,z).degree(w)<=4,'named complete global numerator '+str(pp))
    if pp!=0:
        check(S.cancel(r*translated_bd.subs(p,pp)).subs(r,0)!=0,
              'named nonextendability hostile '+str(pp))

records['constant_constraint_minor']=str(minor)
records['complete_dimensions']=[6,5]
records['boundary_values']=['c0-b+2*p*c1+4*a*p^2','n0-d+2*p*n1']
records['actual_map']='x=w+p; t fixed; determinant 1; no W automorphism'
records['scope']='finite p arbitrary; boundary octic a*(x-p)^6, a*m0!=0; F_D constant'
semantic=hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest()
print('PASS finite-point constant-D family translation:',gates,'always-active exact gates')
print('Complete section spaces for every complex point: dimensions 6 and 5')
print('Exact plane symplectic transport into fixed-zero excluded family')
print('Nonextendable surface-translation hostile retained; no new surface automorphism')
print('Semantic SHA256:',semantic)
