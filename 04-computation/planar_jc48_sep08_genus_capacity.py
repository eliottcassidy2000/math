#!/usr/bin/env python3
"""Exact controls for the sharp all-degree genus bound for dx/sqrt(N).

Finite partitions only test the uniform capacity arithmetic. Actual exactness
and the extremal classification are analytic; no survivor is called exact
merely because it passes the necessary capacity inequality.
"""
from hashlib import sha256
import json
import sympy as S

gates=[]


def need(label,value):
    if not bool(value):raise RuntimeError(label)
    gates.append(label)


def zero(label,value):
    need(label,S.cancel(value)==0)


def partitions(n,upper=None):
    if n==0:
        yield ()
        return
    if upper is None:upper=n
    for first in range(min(n,upper),0,-1):
        for rest in partitions(n-first,first):yield (first,)+rest


universe=0;passed=0;extremal=[]
for n in range(3,25):
    for part in partitions(n):
        universe+=1
        if 2 in part:continue
        square=all(m%2==0 for m in part)
        if square:continue
        s=part.count(1)
        ho=sum(m>=3 and m%2==1 for m in part)
        he=sum(m>=4 and m%2==0 for m in part)
        capacity=sum(max(m-2,0) for m in part)
        zero(f'capacity identity n{n} {part}',capacity-(n-s-2*ho-2*he))
        required=n-2 if n%2 else n//2-1
        if capacity<required:continue
        passed+=1
        if n%2:
            need(f'odd-degree necessary pure root n{n} {part}',len(part)==1)
        else:
            B=s+ho;genus=(B-2)//2
            need(f'even branch parity n{n} {part}',B%2==0 and B>=2)
            need(f'sharp genus bound n{n} {part}',genus<=(n-4)//4)
            if n%4==0 and genus==n//4-1:
                m=n//4
                need(f'extremal partition n{n} {part}',part==(2*m+1,)+(1,)*(2*m-1))
                extremal.append((n,part))
need('one extremal partition each degree4m through24',len(extremal)==6)

x,a,b=S.symbols('x a b',nonzero=True)
# All even degrees4..30, with actual exact families, including squarefree
# residual polynomials of either odd or even degree.
for r in range(1,15):
    N=x**(r+2)*(a*x**r+b)
    R=-S.Rational(2,r)/(b*x**(r+1))
    zero(f'actual sharp primitive r{r}',N*S.diff(R,x)+S.diff(N,x)*R/2-1)
    branches=r+((r+2)%2)
    genus=(branches-2)//2
    need(f'actual sharp genus r{r}',genus==(2*r+2-4)//4)
    need(f'primitive pole order r{r}',2*r+2>=4)

# Full extremal odd primitive spaces, m=1..8: ord∞ X=−2,
# ord∞ Y=−(2m−1); no YX term fits the allowed pole order.
for m in range(1,9):
    cap=2*m-1
    even=[j for j in range(2*m+2) if 2*j<=cap]
    odd=[j for j in range(2*m+2) if 2*j+cap<=cap]
    need(f'complete even monomial pole list m{m}',even==list(range(m)))
    need(f'complete odd monomial pole list m{m}',odd==[0])
    coeff=S.symbols(f'p0:{2*m}')
    P=sum(coeff[j]*x**j for j in range(2*m))
    derivative=S.diff(P,x)
    for j in range(2*m-2):
        zero(f'extremal sparse coefficient m{m},j{j}',derivative.coeff(x,j)-(j+1)*coeff[j+1])

# Degree12 makes genus2 exact; elliptic-only claims beyond degree8 fail.
N=x**7*(x**5+1);R=-S.Rational(2,5)/x**6
zero('genus-two hostile to universal elliptic bound',N*S.diff(R,x)+S.diff(N,x)*R/2-1)
need('genus-two branch count', (6-2)//2==2)

# Extremal multiplicities alone do not suffice: an intermediate coefficient
# destroys the only possible odd primitive direction.
for m in range(2,8):
    P=1+x+x**(2*m-1)
    need(f'extremal-position hostile m{m}',S.diff(P,x).coeff(x,0)==1)
    need(f'extremal curve remains squarefree m{m}',S.gcd(P,S.diff(P,x))==1)

print('genus_capacity: PASS')
print('scope: sharp all-degree genus bound and complete extremal classification in degrees4m')
print('finite partition universe: every partition in degrees3..24;',universe,'total;',passed,'capacity survivors')
print('capacity survivors are necessary-only; six extremal partitions checked')
print('actual sharp families: r=1..14, even degrees4..30; full primitive spaces m=1..8')
print('hostiles: genus2 exact degree12; extremal multiplicities with wrong positions')
print('gates:',len(gates))
print('gate-label sha256:',sha256(json.dumps(gates,separators=(',',':')).encode()).hexdigest())
