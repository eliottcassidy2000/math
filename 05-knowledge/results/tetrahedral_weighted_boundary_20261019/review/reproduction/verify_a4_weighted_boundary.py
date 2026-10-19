#!/usr/bin/env python3
"""Exact finite checks accompanying A4_WEIGHTED_BOUNDARY.md.

Uses explicit exceptions, not Python assertions. Written geometric proofs are
not replaced by the finite algebra and group checks in this companion.
"""
from fractions import Fraction as F
from collections import Counter
from pathlib import Path
import json
import sympy as S
from sympy.matrices.normalforms import smith_normal_form

checks = []
def check(name, condition):
    if not condition:
        raise RuntimeError('FAIL: ' + name)
    checks.append(name)
def zero(name, expression):
    check(name, S.cancel(S.expand(expression)) == 0)

# Q(delta), delta^2=-3, represented without floating point or radicals.
def k(a=0,b=0): return (F(a),F(b))
def add(x,y): return (x[0]+y[0],x[1]+y[1])
def neg(x): return (-x[0],-x[1])
def sub(x,y): return add(x,neg(y))
def mul(x,y): return (x[0]*y[0]-3*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def inv(x):
    d=x[0]*x[0]+3*x[1]*x[1]
    if d == 0: raise ZeroDivisionError
    return (x[0]/d,-x[1]/d)
def mm(a,b):
    return tuple(add(mul(a[2*i],b[j]),mul(a[2*i+1],b[2+j])) for i in range(2) for j in range(2))
def det(a): return sub(mul(a[0],a[3]),mul(a[1],a[2]))
def mi(a):
    d=inv(det(a)); return tuple(mul(d,t) for t in (a[3],neg(a[1]),neg(a[2]),a[0]))
def mv(a,v): return (add(mul(a[0],v[0]),mul(a[1],v[1])),add(mul(a[2],v[0]),mul(a[3],v[1])))
I=(k(1),k(),k(),k(1)); minusI=tuple(neg(t) for t in I)
omega=k(F(-1,2),F(1,2)); omega2=mul(omega,omega)
r=(k(1),k(),k(),omega)
dinv=k(0,F(-1,3))
q=tuple(mul(dinv,k(t)) for t in (1,2,1,-1))
def closure(gens):
    found={I}; todo=[I]
    while todo:
        a=todo.pop()
        for b in gens:
            c=mm(a,b)
            if c not in found:
                found.add(c);todo.append(c)
                if len(found)>200: raise RuntimeError('unexpected nonfinite closure')
    return found
G=closure([r,q]); check('G4 order 24',len(G)==24)
check('q squared is minus identity',mm(q,q)==minusI)
check('r cubed is identity',mm(mm(r,r),r)==I)
check('reflection determinant order three',det(r)==omega and det(q)==k(1))
for i,g in enumerate(sorted(G)):
    check(f'inverse {i}',mm(g,mi(g))==I)
    check(f'determinant cube {i}',mul(mul(det(g),det(g)),det(g))==k(1))
vertices=[(k(1),k()),(k(1),k(1)),(omega,k(1)),(omega2,k(1))]
def perm(g):
    out=[]
    for v in vertices:
        w=mv(g,v)
        hits=[i for i,z in enumerate(vertices) if sub(mul(w[0],z[1]),mul(w[1],z[0]))==k()]
        if len(hits)!=1: raise RuntimeError('projective vertex failure')
        out.append(hits[0])
    return tuple(out)
P={perm(g) for g in G}
check('projective group order 12',len(P)==12)
for p in sorted(P):
    inversions=sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))
    check('even permutation '+str(p),inversions%2==0)
check('r has unique fixed sheet',sum(i==j for i,j in enumerate(perm(r)))==1)
check('q double transposition',sum(i==j for i,j in enumerate(perm(q)))==0 and perm(mm(q,q))==(0,1,2,3))
Q={g for g in G if det(g)==k(1)}
check('Q8 has order eight',len(Q)==8)
def order(g):
    z=I
    for n in range(1,25):
        z=mm(z,g)
        if z==I:return n
    raise RuntimeError('order exceeded')
check('Q8 element orders',Counter(order(g) for g in Q)==Counter({1:1,2:1,4:6}))
comm=closure([mm(mm(mm(a,b),mi(a)),mi(b)) for a in Q for b in Q])
check('Q8 commutator is central C2',comm=={I,minusI})
for g in Q-{I}:
    # No eigenvalue one: det(g-I) != 0.
    check('free Q8 element '+str(g),det(tuple(sub(t,s) for t,s in zip(g,I)))!=k())
H=closure([r,minusI]); check('H order six',len(H)==6)
core=set(H)
for g in G: core &= {mm(mm(g,h),mi(g)) for h in H}
check('H core is central C2',core=={I,minusI})
normal_r=closure([mm(mm(g,r),mi(g)) for g in G])
check('one cubic inertia normally generates projective A4',len({perm(g) for g in normal_r})==12)
# Conjugation cyclically permutes the three nonzero Q8/center classes.
cosets=[]
for g in sorted(Q):
    c=frozenset((g,mm(minusI,g)))
    if c not in cosets:cosets.append(c)
nontrivial=[c for c in cosets if I not in c]
conj=[frozenset(mm(mm(r,g),mi(r)) for g in c) for c in nontrivial]
check('triality moves all three nonzero classes',all(a!=b for a,b in zip(nontrivial,conj)) and set(conj)==set(nontrivial))

# Polynomial identities and invariant rings.
x,y,a,b,c,u,v,h,e,z,p,qq=S.symbols('x y a b c u v h e z p qq')
d=S.sqrt(-3); om=(-1+d)/2
AA=x**3*y-y**4; BB=x**4+8*x*y**3; CC=x**6-20*x**3*y**3-8*y**6
zero('tetrahedral cubic relation',CC**2-BB**3+64*AA**3)
for label,f in [('A',AA),('B',BB),('C',CC)]:
    zero('q invariant '+label,f.subs({x:(x+2*y)/d,y:(x-y)/d},simultaneous=True)-f)
    zero('r character '+label,f.subs(y,om*y)-(om*f if label=='A' else f))
U=a*a+8*b; V=a**3-20*a*b-8*c
f=27*a**4-18*u*a*a-8*v*a-u*u
zero('quartic elimination modulo cone',f.subs({u:U,v:V},simultaneous=True)-64*(a*c-b*b))
bb=(u-a*a)/8; cc=(7*a**3-5*a*u-2*v)/16
zero('inverse b presentation',U.subs(b,bb)-u)
zero('inverse c presentation',V.subs({b:bb,c:cc},simultaneous=True)-v)
zero('cone relation recovers quartic',64*(bb**2-a*cc)+f)
check('quartic discriminant',S.factor(S.discriminant(f,a))==-80621568*(u**3-v**2)**2)
zero('generic cusp factorization',f.subs({u:z*z,v:z**3})-(a-z)*(3*a+z)**3)
E=a**3-2*a*b+c
zero('square identity for deleted divisor',(a*a-b)**2-a*E-(b*b-a*c))
# R=D(e) coordinate isomorphism, checked after embedding in C[x,y].
subcone={a:x*x,b:x*y**3,c:y**6}
eh=S.factor(E.subs(subcone)); hh=S.cancel((a*a-b).subs(subcone)/eh)
zero('e is square of cubic',eh-(x**3-y**3)**2)
for label,old,new in [('a',x*x,hh**2*eh),('b',x*y**3,hh**4*eh**2-hh*eh),('c',y**6,eh*(hh**3*eh-1)**2)]:
    zero('inverse localized chart '+label,old-new)
ap=h*h*e; bp=h**4*e**2-h*e; cp=e*(h**3*e-1)**2
zero('localized chart cone equation',bp*bp-ap*cp)
zero('localized chart e recovery',E.subs({a:ap,b:bp,c:cp},simultaneous=True)-e)
zero('localized chart h recovery',ap*ap-bp-h*e)
UU=9*h**4*e**2-8*h*e; VV=-27*h**6*e**3+36*h**3*e**2-8*e
zero('target formula u',U.subs({a:ap,b:bp,c:cp},simultaneous=True)-UU)
zero('target formula v',V.subs({a:ap,b:bp,c:cp},simultaneous=True)-VV)
J=S.diff(UU,h)*S.diff(VV,e)-S.diff(UU,e)*S.diff(VV,h)
zero('Jacobian is 64e',J-64*e)
zero('target cusp pullback',UU**3-VV**2-64*e*e*(h**3*e-1))
zero('inverse h equation',27*(UU**3-VV**2)*h**4+64*VV*h-64*UU)
zero('primitive h recovery of e',16*(9*h*h*UU-8)*e-27*h**3*(UU**3-VV**2)-16*VV)
res=S.resultant(f,h*(27*a**3-9*u*a-2*v)-2*(9*a*a-u),a)
zero('exact inverse resultant',res+314928*(u**3-v**2)*(27*(u**3-v**2)*h**4+64*v*h-64*u))
D=u**3-v**2
coeff=[1,0,0,64*v/(27*D),-64*u/(27*D)]
tr=[S.Integer(4)]
for n in range(1,5):
    value=-sum(coeff[j]*tr[n-j] for j in range(1,n))-n*coeff[n]
    tr.append(S.cancel(value))
for n,want in enumerate([0,0,-64*v/(9*D),256*u/(27*D)],1):zero('trace moment '+str(n),tr[n]-want)
# Dual derivations and all frame identities.
uh=(-81*h**6*e**2+72*h**3*e-8)/(64*e)
ue=(81*h**5*e**2-54*h*h*e)/32
vh=-9*h**4/S.Integer(32)+h/(8*e)
ve=(9*h**3*e-2)/16
def du(t):return S.cancel(uh*S.diff(t,h)+ue*S.diff(t,e))
def dv(t):return S.cancel(vh*S.diff(t,h)+ve*S.diff(t,e))
zero('du u = 1',du(UU)-1);zero('du v = 0',du(VV))
zero('dv u = 0',dv(UU));zero('dv v = 1',dv(VV)-1)
zero('commutator on h',du(dv(h))-dv(du(h)))
zero('commutator on e',du(dv(e))-dv(du(e)))
zero('du h polar coefficient',S.limit(e*uh,e,0)+S.Rational(1,8))
zero('dv h polar coefficient',S.limit(e*vh,e,0)-h/8)
Ur=9*qq**4/p**2-8*qq;Vr=-27*qq**6/p**3+36*qq**3/p-8*p
zero('rational constant Jacobian',S.diff(Ur,qq)*S.diff(Vr,p)-S.diff(Ur,p)*S.diff(Vr,qq)-64)
zero('rational map agrees u',Ur.subs({qq:h*e,p:e},simultaneous=True)-UU)
zero('rational map agrees v',Vr.subs({qq:h*e,p:e},simultaneous=True)-VV)
# Weighted graph checks, not a substitute for the written class-group proof.
M=S.Matrix([[-2,1,1,1],[1,-2,0,0],[1,0,-2,0],[1,0,0,-2]])
check('D4 determinant four',M.det()==4)
SN=smith_normal_form(M,domain=S.ZZ)
check('D4 Smith invariants 1 1 2 2',[abs(SN[i,i]) for i in range(4)]==[1,1,2,2])
check('A1 determinant two',abs(S.Matrix([[-2]]).det())==2)
T=S.Matrix([[1,0,0,0],[0,0,0,1],[0,1,0,0],[0,0,1,0]])
check('triality preserves D4 matrix',T.T*M*T==M)
ker=[S.Matrix(w) for w in __import__('itertools').product(range(2),repeat=4) if all(int(t)%2==0 for t in M*S.Matrix(w))]
check('D4 mod two kernel dimension two',len(ker)==4)
nonzero=[w for w in ker if any(w)]
check('triality acts freely on nonzero kernel',all(any(int(t)%2 for t in T*w-w) for w in nonzero))

out={'status':'PASS','checks':len(checks),'scope':'exact finite algebra/group checks; geometric arguments are in the paper','check_names':checks,
     'group_orders':{'G4':len(G),'A4':len(P),'Q8':len(Q),'H':len(H)},
     'quartic_discriminant':str(S.factor(S.discriminant(f,a))),
     'trace_moments':[str(t) for t in tr[1:]],'D4_smith':[1,1,2,2],'A1_smith':[2]}
Path(__file__).with_name('verification-results.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
print('PASS: A4 weighted-boundary exact companion')
print('checks:',len(checks))
print('groups: G4=24, A4=12, Q8=8, H=6; core(H)=C2')
print('quartic normalization, global hostile, trace packet and dual derivations: PASS')
print('weighted graphs: A1 determinant 2; D4 Smith 1,1,2,2')
print('scope: finite identities checked; geometric proofs remain separately stated')
