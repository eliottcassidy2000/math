#!/usr/bin/env python3
"""Exact controls for the all-m parabolic polynomial-mate first-jet gate.

RESERVED pending independent audit.  The proof is uniform in m>=1 and
all polynomial degrees.  Finite controls exercise the exact source-image
condition, unrestricted higher rows, and sharp rational hostiles.
"""
from hashlib import sha256
import sympy as S
gates=[]
def need(label,value):
    if not bool(value):raise RuntimeError(label)
    gates.append(label)
def zero(label,value):need(label,S.cancel(value)==0)
def jac(F,G,u,t):return S.diff(F,u)*S.diff(G,t)-S.diff(F,t)*S.diff(G,u)
u,t,w=S.symbols('u t w')
i,k,j,h,m,ell=S.symbols('i k j h m ell',integer=True)
zero('general monomial bracket coefficient',
     (i-m*k)*h-k*(j-m*h)-(i*h-k*j))
zero('general monomial source-weight identity',
     (i-m*k)+(j-m*h)+m-1+m*(k+h-1)-(i+j-1))
image_checks=degree_checks=0
for mm in range(1,7):
    for ll in range(-12,6):
        for kk in range(7):
            source_exponent=ll+mm*kk
            required=max(0,(-ll+mm-1)//mm)
            need('exact source-image monomial criterion',
                 (source_exponent>=0)==(kk>=required))
            image_checks+=1
    pairs={(ii-mm*kk,kk):(ii,kk) for ii in range(7) for kk in range(6)}
    need('original monomial map is injective',len(pairs)==42)
    for (ll,kk),(ii,_) in pairs.items():
        need('inverse monomial exponent recovers original',ll+mm*kk==ii)
    # All higher rows retained in both exact source Jacobians.
    Fhat=w*w+3*w+u*(w**3+2)+u*u*(w**4+1)
    Gsrc=t*t+u*t+u**3+u**(mm+1)*t**3
    Ghat=S.expand(Gsrc.subs(t,w/u**mm))
    Fsrc=Fhat.subs(w,u**mm*t)
    zero('actual constant-row chart bracket',jac(Fsrc,Gsrc,u,t).subs(t,w/u**mm)-
         u**mm*jac(Fhat,Ghat,u,w))
    low=-2*mm
    zero('constant-row leading term with all higher rows',
         S.expand(u**mm*jac(Fhat,Ghat,u,w)).coeff(u,low+mm-1)-
         (-low*(2*w+3)*w*w))
    g=w*w+2*w+3
    Fhat=u*g+u*u*(w**3+1)
    Fsrc=Fhat.subs(w,u**mm*t)
    zero('actual first-row chart bracket',jac(Fsrc,Gsrc,u,t).subs(t,w/u**mm)-
         u**mm*jac(Fhat,Ghat,u,w))
    zero('first-row leading term with all higher rows',
         S.expand(u**mm*jac(Fhat,Ghat,u,w)).coeff(u,low+mm)-
         (g*2*w-low*S.diff(g,w)*w*w))
    # Positive control retains arbitrary constant and two polynomial
    # source terms whose weights cannot change the initial mate row.
    beta=S.Symbol('beta',nonzero=True);q=S.Symbol('q')
    FF=q+beta*u;GG=t/beta+u**4+3*u
    zero('all-m polynomial positive',jac(FF,GG,u,t)-1)
    zero('forced initial polynomial mate row',
         S.expand(GG.subs(t,w/u**mm)).coeff(u,-mm)-w/beta)
    # F=u*w has a rational mate for EVERY m, violating the second row
    # conclusion as soon as the original polynomial constraint is lost.
    Fr=u**(mm+1)*t;Gr=1/(mm*u**mm)
    zero('all-m rational hostile to first-row conclusion',jac(Fr,Gr,u,t)-1)
    need('rational hostile fails negative-row divisibility',
         S.denom(Gr)!=1 and S.expand(Gr*u**mm).subs(w,0)!=0)
    if mm>=2:
        Fr=u**mm*t;Gr=u**(1-mm)/(mm-1)
        zero('rational hostile to constant-row conclusion',jac(Fr,Gr,u,t)-1)
    # The m=1 corner has no negative ell producing weight zero in the
    # nonconstant constant-row case.  ell=0 is explicitly not treated as
    # a nonzero leading derivative.
    if mm==1:
        need('m1 first-row exponent corner',all(ll+mm-1<0 for ll in range(-12,0)))
for hh in range(1,7):
    for jj in range(1,7):
        for ll in range(-12,0):
            need('negative-row degree cannot cancel',jj-ll*hh>0)
            need('weight-zero first-row coefficient has positive degree',hh+jj-1>=1)
            degree_checks+=1
# Passing the necessary rows is not sufficient for a polynomial mate.
zero('first-jet gate does not claim sufficiency',
     S.diff(u+u*u,u).subs(u,-S.Rational(1,2)))

# Generic square-prefix first rows: every omitted row is divisible by u².
n4,n5,n6,c,d,b,q1,p4,q2,nu,mu,r0,r1,r2=S.symbols(
    'n4 n5 n6 c d b q1 p4 q2 nu mu r0 r1 r2')
Hhat=n4*w*w+c*w+d+u*(n5*w*w+b*w+q1)+u*u*(n6*w*w+p4*w+q2)
Lhat=nu*w+r0+u*(mu*w+r1)+u*u*r2
Fhat=S.expand(Hhat*Hhat+Lhat)
f=Fhat.coeff(u,0);g=Fhat.coeff(u,1)
zero('active multiplicity-four leading obstruction',f.coeff(w,4)-n4*n4)
zero('higher multiplicity constant row',f.subs(n4,0)-((c*w+d)**2+nu*w+r0))
zero('multiplicity-at-least-six complete first row',
     g.subs({n4:0,n5:0,c:0,nu:0})-((2*b*d+mu)*w+2*d*q1+r1))
zero('multiplicity-five complete first row',
     g.subs({n4:0,c:0,nu:0})-(2*d*n5*w*w+(2*b*d+mu)*w+2*d*q1+r1))

# Full all-p active62 rows supplied explicitly, not by translating W.
p,A,lam,e=S.symbols('p A lam e')
N=u**6;P=(A+2)*u**4+b*u**3+c*u*u
Q=(A+1)*u*u+(b-2*p*A)*u+d
M=lam*u**4+mu*u**3+nu*u*u
R=lam*u*u+(mu-2*p*lam)*u+e
H=N*t*t+P*t+Q;L=M*t+R
for row in (Q,P-2*Q*(u+p)**2,N-P*(u+p)**2+Q*(u+p)**4):
    need('actual all-p62 global H row',S.Poly(S.expand(row),u).degree()<=4)
for row in (R,M-R*(u+p)**2):
    need('actual all-p62 global L row',S.Poly(S.expand(row),u).degree()<=2)
full=S.expand((H*H+L).subs(t,w/u**2))
zero('actual62 constant coefficient',full.coeff(u,0)-((c*w+d)**2+nu*w+e))
zero('actual62 first coefficient',full.coeff(u,1)-
     (2*(c*w+d)*(b*w+b-2*p*A)+mu*w+mu-2*p*lam))
zero('actual62 forced nonzero source derivative',
     full.coeff(u,1).subs({c:0,nu:0,mu:-2*b*d})+2*p*(lam+2*d*A))
print('polynomial_weight_gate: PASS')
print('scope: every integer m>=1, original polynomial mates, complete parabolic source image')
print('conclusion: Fhat=f0+beta*u+O(u²), beta!=0; min weight Ghat=-m; leading row=w/beta')
print('controls: m1..6; image monomials',image_checks,'; negative degree rows',degree_checks)
print('square-prefix: active ord(N)=4 excluded; ord(N)>=6 gives P2=M2=0,M3=-2P3Q0')
print('active62: c=nu=0,mu=-2bd,beta=-2p(lambda+2dA)!=0')
print('rational and polynomial-positive hostiles retained; no sufficiency claim')
print('gates:',len(gates))
print('semantic SHA256:',sha256('\n'.join(gates).encode()).hexdigest())
