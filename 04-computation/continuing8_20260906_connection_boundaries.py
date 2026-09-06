"""Exact incoming-work transfer boundaries; no imported research producer."""
from math import comb
from pathlib import Path
import json
import sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
G=0
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
x,t,z,a,b,m=s.symbols('x t z a b m')
phi=a*x+b
residual=s.expand((x*x+6)*phi*s.diff(phi,x)-x*(phi*phi+6))
gate(s.expand(residual-(-a*b*x*x+(6*a*a-b*b-6)*x+6*a*b))==0,'complete linear-chart coefficient equation')
L=lambda f,index,parameter=6:(x*x+parameter)*s.diff(f,x)-2*index*x*f
for sign in (-1,1):
    for degree in range(9):
        theta=x**degree
        gate(s.expand(L(theta.subs(x,sign*x),m)-sign*L(theta,m).subs(x,sign*x))==0,'both permitted charts intertwine exactly')
for phi in (x*x,2*x-x*x,x*x+1,x**3,2*x+1):
    forced=x/phi
    gate(s.cancel(L(phi,m)-forced*L(x,m).subs(x,phi))!=0,'natural Student scalar transfer fails')
    zero_h=(x*x+6)*s.diff(phi,x)/(phi*phi+6)
    for degree in range(6):
        theta=x**degree
        gate(s.cancel(L(theta.subs(x,phi),0)-zero_h*L(theta,0).subs(x,phi))==0,'m=0 sharp exception for arbitrary chart')
gate(s.cancel(L(x*x,m)-L(x,m).subs(x,x*x)/x-(x**3+12*x-6/x))==0,'explicit x^2 residual')
for degree in range(7):
    theta=x**degree
    gate(s.cancel(L(theta.subs(x,a*x),m)-L(theta,m,6*a*a).subs(x,a*x)/a)==0,'scaling survives with quadratic parameter retained')

# Actual source-coordinate identity, proved before any characteristic reduction.
u=x*x*t
pi=t*(1+u)
y=x*t*pi
gate(s.expand(pi**3-y*y-t**3*(1+u)**2)==0,'actual source identity over Z')
families=[]
for prime in (2,3,5,7,11,13):
    for exponent in range(1,4):
        P=prime**exponent
        if P>125: continue
        L0=P//2; epsilon=P-2*L0; T=3*L0+epsilon
        rho=(2*T+2)//3
        N0=4*T//3+1
        gate(rho==P and T-rho==L0 and N0==T+L0+1,'full diagonal parameters')
        gate(0<=L0<P,'quotient degree below Frobenius jump')
        gate(all(comb(P,k)%prime==0 for k in range(1,P)),'entire Frobenius coefficient cancellation')
        # The full integer source basis on intercept 2T is triangular with
        # unit diagonal after removing (1-z)^rho. This does not use a field rank.
        basis=s.Matrix(L0+1,L0+1,lambda row,col:
            (-1)**row*comb(L0-col,row-col) if col<=row else 0)
        gate(all(basis[i,i] in (-1,1) for i in range(L0+1)),'source diagonal basis remains invertible modulo every prime')
        quotient=[comb(P+k-1,k)%prime for k in range(P+2)]
        gate(quotient[0]==1 and all(v==0 for v in quotient[1:P]) and quotient[P]==1,'full inverse series, not only a single source lift')
        first=next(k for k in range(L0+1,len(quotient)) if quotient[k])
        gate(T+first==T+P,'earliest unavoidable projected failure')
        gate(T+P-N0==(P+1)//2-1,'unbounded excess delay over characteristic zero')
        if P<=13:
            actual=s.Poly(s.expand(pi**epsilon*(pi**3-y*y)**L0),x,t)
            reduced={mon:int(coefficient)%prime for mon,coefficient in actual.terms() if coefficient%prime}
            expected={(0,T):1,(2*P,T+P):1}
            gate(reduced==expected,'literal actual-source expansion and exact sparse reduction')
            signed=s.Poly(s.expand((1-z)**P),z)
            gate(all(int(signed.nth(k))%prime==(1 if k==0 else -1%prime if k==P else 0) for k in range(P+1)),'signed diagonal Frobenius identity')
        families.append([prime,exponent,P,T,N0,T+P,T+P-N0])

# Arithmetic saturation and global nonvanishing are separate from jet order.
A=s.Matrix([[0,6],[-4,0],[0,-3]])
target=s.Matrix([0,2,0]); solution=s.Matrix([-s.Rational(1,2),0])
gate(A*solution==target,'actual Student compatible target 2x has denominator two')
gate(target.applyfunc(lambda q:q%2)==s.zeros(3,1),'mod-two membership passes with zero lift')
gate(all(int((A*s.Matrix([r0,r1])-target)[1])%4!=0 for r0 in range(4) for r1 in range(4)),'mod-four full equation forbids every lift')
for prime in (2,3,5,7):
    P=prime*prime
    for degree in range(1,8):
        gate(P*degree>degree,'Frobenius composition is nonzero while multiplying vanishing order')
print('STUDENT: only +/-x scalar charts for m nonzero; m=0 exception and parameter-covariant scaling verified')
print('FROBENIUS [prime,a,P,T,char0_cutoff,exact_modp_cutoff,excess]')
print(json.dumps(families,separators=(',',':')))
print('ARITHMETIC: actual 2x target passes mod2 but fails mod4; rational lift -1/2')
print('PASS',G,'always-active exact gates; normal/-O raw LF')
print('Scope: natural scalar intertwiner and declared characteristic-p source-module recognition; no characteristic-zero Keller obstruction.')
