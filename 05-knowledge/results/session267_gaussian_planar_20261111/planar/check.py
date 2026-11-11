#!/usr/bin/env python3
"""Exact controls for the proved Gaussian/rotation/flag bridge; no source lift."""
import sympy as s
x,G,H=s.symbols('x G H')
r6=s.sqrt(6)
U,V=H+s.I*G,H-s.I*G

def T(m,f):
    return s.Poly(f,x).terms() and s.expand(sum(c*r6**j*G**j*H**(2*m-j) for (j,),c in s.Poly(f,x).terms()))
def R(p): return s.expand(H*s.diff(p,G)-G*s.diff(p,H))
def moment(k): return 0 if k%2 else s.factorial2(k-1) if k else s.Integer(1)
def E(p): return s.expand(sum(c*moment(a)*moment(b) for (a,b),c in s.Poly(s.expand(p),G,H).terms()))
def L(m,f): return -x*f+(x*x+6)*s.diff(f,x)/(2*m)
def require(condition):
    if not condition:
        raise AssertionError("Exact control failed")

checks=0
for m in range(1,9):
    for j in range(2*m+1):
        require(s.expand(R(T(m,x**j))-2*m/r6*T(m,s.expand(L(m,x**j))))==0)
        mu=0 if j%2 else 6**(j//2)*moment(j)*moment(2*m-j)/moment(2*m)
        require(s.simplify(E(T(m,x**j))/moment(2*m)-mu)==0)
        checks+=2
    for k in range(-m,m+1):
        p=s.expand(U**(m+k)*V**(m-k))
        require(s.expand(R(p)-2*s.I*k*p)==0)
        require(E(p)==(2**m*s.factorial(m) if k==0 else 0))
        checks+=2
for m in range(1,4):
    for n in range(1,4):
        a={k:s.Integer(k*k+2*k+3) for k in range(-m,m+1)}
        b={k:s.Integer(2*k*k-k+1) for k in range(-n,n+1)}
        p=s.expand(sum(v*U**(m+k)*V**(m-k) for k,v in a.items()))
        q=s.expand(sum(v*U**(n+k)*V**(n-k) for k,v in b.items()))
        expected=2**(m+n)*s.factorial(m+n)*sum(v*b.get(-k,0) for k,v in a.items())
        require(E(p*q)==expected)
        checks+=1
require(E(U**2)==E(V**2)==0 and E(U**2*V**2)==8)
require(E((U**2+V**2)**2)==16)
# Gaussian vanishing does not pay the row16 high coefficient x^11.
require(E(T(16,x**11))==0)
m=16; d=(m-1)//2+2
require(d==9 and 11>=d+2)
# For this flag all L-images have degree <=10, proving x^11 unrepairable.
require(all(s.degree(L(m,x**j),x)<=10 for j in range(d+1)))
checks+=5
print(f'PASS: {checks} exact Gaussian/rotation/convolution/flag controls.')
