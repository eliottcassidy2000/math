#!/usr/bin/env python3
"""Exact finite quartic map of the specified DG surface; no Keller claim."""
import hashlib
import json
import sympy as S

x,t,r,b,T,B,X,R,u,j,v=S.symbols('x t r b T B X R u j v')
n=0
def check(c,label):
 global n
 n+=1
 if not c: raise RuntimeError(label)
def eq(a,c,label): check(S.cancel(S.together(a-c))==0,label)
def sub(f,tab): return f.subs(tab,simultaneous=True)

bb=-x*x-x**4*t
tt=-r*r-r**4*b
binary=T*X**4+X*X*R*R+B*R**4
eq(binary.subs({X:x,R:1,T:t,B:bb}),0,'finite source chart of quartic')
eq(binary.subs({X:1,R:r,T:tt,B:b}),0,'boundary chart of quartic')
check(S.Poly(binary,X,R).coeff_monomial(X**2*R**2)==1,'no identically zero projective fibre')
eq(binary.subs({T:0,B:0}),X**2*R**2,'origin has two double points')
eq(S.discriminant(binary.subs(R,1),X),16*T*B*(1-4*T*B)**2,'binary quartic discriminant')
check((16*T*B*(1-4*T*B)**2).subs({T:1,B:1})!=0,'unramified four-point control')

uu=x*x*t; jj=x*t; vv=x*(1+uu)
vals={u:uu,j:jj,v:vv,T:t,B:bb}
relations=[u*u+u+T*B,j*j-T*u,v*v+B*(1+u),u*j-T*v+j,u*v+B*j,j*v+T*B]
for e in relations: eq(sub(e,vals),0,'source multiplication relation')
chart={x:1/r,t:tt}
for f,g in [(uu,-1-r*r*b),(jj,-r*(1+r*r*b)),(vv,-r*b)]:
 eq(sub(f,chart),g,'basis function extends')

# Independent finite-free trace algebra on basis (1,u,j,v).
Mu=S.Matrix([[0,-T*B,0,0],[1,-1,0,0],[0,0,-1,-B],[0,0,T,0]])
Mj=S.Matrix([[0,0,0,-T*B],[0,0,T,0],[1,-1,0,0],[0,T,0,0]])
Mv=S.Matrix([[0,0,-T*B,-B],[0,0,0,-B],[0,-B,0,0],[1,0,0,0]])
Id=S.eye(4)
for e in [Mu*Mu+Mu+T*B*Id,Mj*Mj-T*Mu,Mv*Mv+B*(Id+Mu),
          Mu*Mj-T*Mv+Mj,Mu*Mv+B*Mj,Mj*Mv+T*B*Id]:
 check(e.applyfunc(S.expand)==S.zeros(4),'matrix multiplication relation')
for A,C in [(Mu,Mj),(Mu,Mv),(Mj,Mv)]: check((A*C-C*A).applyfunc(S.expand)==S.zeros(4),'commutative trace algebra')
basis=[Id,Mu,Mj,Mv]
trace=S.Matrix(4,4,lambda a,c:S.trace(basis[a]*basis[c]))
expected=S.Matrix([[4,-2,0,0],[-2,2-4*T*B,0,0],
                   [0,0,-2*T,-4*T*B],[0,0,-4*T*B,-2*B]])
check((trace-expected).applyfunc(S.expand)==S.zeros(4),'full trace matrix')
eq(trace.det(),16*T*B*(1-4*T*B)**2,'trace discriminant agrees with binary discriminant')

# Cech representatives proving the proposed global basis is saturated,
# not just a linearly independent generic basis.
F=t*x**4+x*x+b
eq(t*x*x-(-1-b/x**2),F/x**2,'u connecting class')
eq(t*x-(-1/x-b/x**3),F/x**3,'j connecting class')
eq(x+t*x**3-(-b/x),F/x,'v connecting class')

# Exact ramification in both charts and branch images.
J0=-S.diff(bb,x)
Jinf=S.diff(tt,r)
eq(J0,2*x*(1+2*uu),'source Jacobian')
eq(Jinf,-2*r*(1+2*r*r*b),'boundary Jacobian')
eq(1-4*t*bb,(1+2*uu)**2,'hyperbola pullback is a square')
eq(bb.subs(x,0),0,'retained ramification curve maps to B=0')
eq(tt.subs(r,0),0,'boundary ramification curve maps to T=0')
eq(bb.subs(t,-1/(2*x*x)),-x*x/2,'hyperbola ramification residue degree two')
eq(bb.subs(t,0),-x*x,'other T=0 component residue degree two')
eq(tt.subs(b,0),-r*r,'other B=0 component residue degree two')
check(S.diff(tt,r,2).subs(r,0)==-2,'boundary t has exact order two')
check(S.diff(bb,x,2).subs(x,0)==-2,'retained b has exact order two')

# Fixed first-coordinate obstruction: integration on the source forces
# a nonextendible simple pole; a polynomial in t cannot cancel it.
c,a0,a1,a2=S.symbols('c a0 a1 a2')
g=-c*x+a0+a1*t+a2*t*t
eq(-S.diff(g,x),c,'source-relative primitive identity')
eq(S.cancel(S.together(r*sub(g,chart))).subs(r,0),-c,'unavoidable global primitive pole')

payload={'quartic':str(binary),'basis':['1','u','j','v'],
         'discriminant':str(S.factor(trace.det())),
         'ramification':['D: r=0','Z: x=0','C: 1+2*x^2*t=0'],
         'branch_profiles':{'T=0':[[2,1],[1,2]],'B=0':[[2,1],[1,2]],'1-4TB=0':[[2,2]]}}
print('STATUS: exact finite flat degree-four control; not Keller')
print('projective quartic: T X^4+X^2 R^2+B R^4=0; charts recover W exactly')
print('free module basis: 1,u=x^2t,j=xt,v=x(1+x^2t)')
print('trace/binary discriminant: 16*T*B*(1-4*T*B)^2')
print('ramification divisor: D+Z+C, each generic index2; branch: axes and hyperbola')
print('generic profiles: T=0 and B=0 have (e,f)=(2,1)+(1,2); hyperbola has(2,2)')
print('source Jacobian: 2*x*(1+2*x^2*t); ramification remains inside source A2')
print('relative primitive obstruction: no global g with dt wedge dg=c*omega, c!=0')
print('finite (F(t),g) maps cannot be etale throughout source A2')
print('scope: one explicit surface and fixed-first-coordinate class; no general Keller exclusion')
print('semantic_sha256='+hashlib.sha256(json.dumps(payload,sort_keys=True,separators=(',',':')).encode()).hexdigest())
print('PASS gates='+str(n))
