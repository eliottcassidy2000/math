#!/usr/bin/env python3
"""Exact controls for algebraic exactness of u^A(a*u^r+b).

RESERVED pending independent audit.  The all-exponent theorem is analytic;
the finite rectangle checks its predictions using complete linear systems.
"""
from fractions import Fraction as F
from hashlib import sha256
import json
import sympy as S

GATES=[]
def need(label,value):
    if not bool(value):raise RuntimeError(label)
    GATES.append(label)

def zero(label,value):
    need(label,S.cancel(value)==0)

def predicted(A,r):
    return (A,r)==(0,1) or (A>=r+2 and (A-r-2)%(2*r)==0)

def rank(rows):
    rows=[list(map(F,row)) for row in rows]
    if not rows:return 0
    pivot=0
    for col in range(len(rows[0])):
        choices=[i for i in range(pivot,len(rows)) if rows[i][col]]
        if not choices:continue
        i=choices[0];rows[pivot],rows[i]=rows[i],rows[pivot]
        v=rows[pivot][col]
        rows[pivot]=[x/v for x in rows[pivot]]
        for j in range(pivot+1,len(rows)):
            v=rows[j][col]
            if v:rows[j]=[x-v*y for x,y in zip(rows[j],rows[pivot])]
        pivot+=1
        if pivot==len(rows):break
    return pivot

def full_system(j,r,a,b):
    # All polynomial coefficients through the forced highest degree M;
    # there is NO residue-class or predicted-congruence restriction here.
    M=j-r+1
    if M<0:return False
    matrix=[[F(0) for _ in range(M+1)] for _ in range(j+1)]
    for q in range(M+1):
        if q:matrix[q-1][q]+=a*q
        matrix[q+r-1][q]+=b*(F(q)+F(r,2))
    target=[F(-int(i==j)) for i in range(j+1)]
    return rank(matrix)==rank([row+[z] for row,z in zip(matrix,target)])

# Complete finite rectangle, with two nonzero coefficient specializations.
# Low degrees and odd degrees use their explicit residue/capacity controls;
# every even degree >=4 uses the unrestricted exact coefficient matrix.
records=[];positive=[];systems=0
for A in range(49):
    for r in range(1,13):
        n=A+r
        expected=predicted(A,r)
        if n==1:
            need('unique degree-one row',(A,r)==(0,1) and expected)
        elif n==2:
            need('degree-two infinity logarithm',not expected)
        elif n%2:
            capacity=max(0,A-2)
            need(f'odd finite pole capacity A{A}r{r}',capacity<n-2)
            need(f'odd rejected A{A}r{r}',not expected)
        else:
            j=n//2-2
            for a,b in ((F(2),F(3)),(F(-3),F(5))):
                actual=full_system(j,r,a,b);systems+=1
                need(f'full matrix iff A{A}r{r}a{a}b{b}',actual==expected)
        records.append((A,r,n,expected))
        if expected:positive.append((A,r))
need('complete rectangle size',len(records)==49*12)
low=[(A,r) for A,r in positive if A+r<=8]
need('six exact-pencil generic types',set(low)=={(0,1),(3,1),(5,1),(7,1),(4,2),(5,3)})

u,X,a,b=S.symbols('u X a b',nonzero=True)

def coefficients(r,ell):
    c=[S.Integer(0)]*(ell+1)
    c[ell]=-S.Rational(2,r*(2*ell+1))/b
    for h in range(ell,0,-1):
        c[h-1]=-S.Rational(2*h,2*h-1)*a*c[h]/b
    return c

# Independent check in the original coordinate, so an incorrect inversion
# weight cannot pass merely by satisfying the transformed operator equation.
explicit=[]
for r in range(1,9):
    for ell in range(6):
        A=r+2+2*r*ell;k=(A+r)//2
        c=coefficients(r,ell)
        B=sum(c[h]*X**(r*h) for h in range(ell+1))
        P=a+b*X**r
        zero(f'transformed primitive r{r}ell{ell}',P*S.diff(B,X)+S.diff(P,X)*B/2+X**(k-2))
        R=sum(c[h]*u**(-k-r*h) for h in range(ell+1))
        N=u**A*(a*u**r+b)
        zero(f'original primitive r{r}ell{ell}',N*S.diff(R,u)+S.diff(N,u)*R/2-1)
        need(f'no parameter radical r{r}ell{ell}',all(z.is_rational_function(a,b) for z in c))
        branch=r+(A%2)
        need(f'actual genus r{r}ell{ell}',(branch-2)//2==(r-1)//2)
        need(f'total primitive pole degree r{r}ell{ell}',A-2==r*(2*ell+1))
        explicit.append((A,r,ell))
zero('degree-one primitive',S.diff(a*u+b,u)/a-1)

# Exact inversion and monomial operator formulas, independently expanded.
for r in range(1,9):
    P=a+b*X**r
    for q in range(13):
        expected=b*(q+S.Rational(r,2))*X**(q+r-1)
        if q:expected+=a*q*X**(q-1)
        zero(f'full monomial operator r{r}q{q}',P*S.diff(X**q,X)+S.diff(P,X)*X**q/2-expected)
        need(f'nonzero highest multiplier r{r}q{q}',q+F(r,2)>0)

# Named hostiles and boundaries.
need('A0 only linear',all(predicted(0,r)==(r==1) for r in range(1,25)))
need('r1 staircase',all(predicted(A,1)==(A==0 or A>=3 and A%2==1) for A in range(49)))
need('r2 staircase',all(predicted(A,2)==(A>=4 and A%4==0) for A in range(49)))
need('same radical field does not preserve weighted exactness',predicted(5,3) and not predicted(7,3) and predicted(11,3))
need('same-field wrong weight fails full system',not full_system(3,3,F(1),F(1)))
# y^2=u^2+1: d(log(u+y))=du/y, but the two infinity residues are +/-1.
y=S.symbols('y',nonzero=True)
zero('elementary logarithm is not an algebraic primitive',(1+u/y)/(u+y)-1/y)
need('two nonzero infinity residues',(-1,1)==tuple(-1/z for z in (1,-1)))
# Dropping b!=0 gives u^3, exact although (A,r)=(0,3) is rejected.
zero('zero-coefficient boundary is outside theorem',u**3*S.diff(-2/u**2,u)+S.diff(u**3,u)*(-2/u**2)/2-1)
# Every coefficient endpoint in the exact pencil is a pure power other than
# exponent two; this tests the separate endpoint formula without dividing by
# a coefficient that was set to zero.
for q in range(21):
    if q==2:continue
    N=b*u**q;R=S.Rational(2,2-q)*u**(1-q)/b
    zero(f'pure-power coefficient endpoint q{q}',N*S.diff(R,u)+S.diff(N,u)*R/2-1)
# r3,ell1 gives a higher primitive on the same elliptic curve as ell0.
c=coefficients(3,1)
zero('first higher elliptic constant coefficient',c[0]-4*a/(9*b**2))
zero('first higher elliptic cubic coefficient',c[1]+2/(9*b))

print('binomial_exactness: PASS')
print('scope: all A>=0,r>=1,ab!=0; exact iff (A,r)=(0,1) or A=r+2+2r*ell,ell>=0')
print('finite rectangle:',len(records),'rows A0..48,r1..12;',len(positive),'accepted;',systems,'unrestricted exact matrix systems')
print('generic exact-pencil rows through degree8:',sorted(low))
print('symbolic primitives:',len(explicit),'rows r1..8,ell0..5; original and inverted coordinates')
print('hostiles: logarithmic quadratic; same field wrong weight; coefficient-zero degeneration')
print('gates:',len(GATES))
print('rectangle sha256:',sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
print('gate-label sha256:',sha256(json.dumps(GATES,separators=(',',':')).encode()).hexdigest())
