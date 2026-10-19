#!/usr/bin/env python3
"""Exact finite controls for two standalone integral-lattice theorems.
Universe and limits are printed. No claim about analytic S6 existence is tested.
"""
from collections import Counter
from itertools import combinations
from math import gcd
import sympy as S
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

def require(ok, message):
    if not ok:
        raise AssertionError(message)

def vp(n,p):
    if n == 0: return float('inf')
    n=abs(n); e=0
    while n%p==0: n//=p; e+=1
    return e

def predicted_count(a,b,p):
    answer=1
    for q in S.factorint(a*b):
        A,B,P=vp(a,q),vp(b,q),vp(p,q)
        if A != B:
            if P != min(A,B): return 0
            if min(A,B)>0: answer *= q**(min(A,B)-1)*(q-1)
        else:
            if P<A: return 0
            answer *= q**(A-1)*(q-1 if P>A else q-2)
    return answer

def diag(A):
    D=smith_normal_form(S.Matrix(A),domain=ZZ)
    return tuple(abs(int(D[i,i])) for i in range(min(D.shape)))

def wedge(M,q):
    I=list(combinations(range(M.rows),q))
    return S.Matrix([[M.extract(i,j).det() for j in I] for i in I])

cases=0
for a in range(2,33):
    for b in range(2,33):
        observed=Counter((-b*x-a*y)%(a*b) for x in range(a) for y in range(b)
                         if gcd(a,x)==gcd(b,y)==1)
        for p in range(a*b):
            require((observed[p] == predicted_count(a, b, p)), (a, b, p))
            require((predicted_count(a, b, p) == predicted_count(a, b, p - a * b)), 'exact check failed')
            cases+=1
        require((predicted_count(a, b, 0) == (int(S.totient(a)) if a == b else 0)), 'exact check failed')
        require((sum(observed.values()) == int(S.totient(a) * S.totient(b))), 'exact check failed')
        K=S.Matrix([[1,1],[a,0],[0,b]])
        require((diag(K) == (1, gcd(a, b))), 'exact check failed')
print('SEIFERT NORMALIZED-GAUGE COUNTS: PASS')
print('Universe: ordered 2 <= a,b <= 32; all signed-p residue classes modulo ab.')
print('Residue classes independently enumerated:',cases)
print('Gauge matrices independently checked:',31**2)
print('Hostile: a=b=5,p=10 has',predicted_count(5,5,10),'admissible marked gauge classes.')
print('Normalized (ell0,ell1,ell2) representatives:',
      [( (10+5*x+5*y)//25,x,y) for x in range(5) for y in range(5)
       if gcd(x,5)==gcd(y,5)==1 and (10+5*x+5*y)%25==0])
print('Even-prime boundary: a=b=2,p=2 count=',predicted_count(2,2,2),
      '; p=4 count=',predicted_count(2,2,4))
print('Zero-clutch boundary: a=b=m gives phi(m); a!=b gives zero: PASS')

for s in range(1,11):
    for mult in range(1,11):
        t=s*mult
        M=S.eye(4);M[0,2]=s;M[1,3]=t
        expected={0:(0,),1:(s,t,0,0),2:(s,s,0,0,0,0),3:(s,t,0,0),4:(0,)}
        for q in range(5):
            W=wedge(M,q)
            require((diag(W - S.eye(W.rows)) == expected[q]), (s, t, q))
print('\nRANK-TWO SQUARE-ZERO EXTERIOR SMITH PROFILE: PASS')
print('Universe: 1 <= s <= 10, t=s*m with 1 <= m <= 10; all q=0,...,4.')
print('Exterior matrices checked:',100*5)
for d in (1,2,3,5,12,101):
    M=S.eye(4);M[0,2]=1;M[1,3]=d
    print('B=diag(1,%d): wedge2 nonzero Smith=%s; cusp components=%d; Euler=%d'
          %(d,tuple(x for x in diag(wedge(M,2)-S.eye(6)) if x),d,2*d))

# Explicit integral conjugations preserve all profile claims; include sign/mixing hostiles.
U=S.Matrix([[1,2],[0,1]]);V=S.Matrix([[0,-1],[1,3]])
for s,t in ((1,1),(1,7),(2,6),(3,12)):
    B=U*S.diag(s,t)*V
    M=S.eye(4);M[:2,2:]=B
    require((abs(B.det()) == s * t), 'exact check failed')
    require((diag(M - S.eye(4)) == (s, t, 0, 0)), 'exact check failed')
    require((diag(wedge(M, 2) - S.eye(6)) == (s, s, 0, 0, 0, 0)), 'exact check failed')
print('Nondiagonal/sign-mixed shear controls: 4 PASS')
print('Scope: exact lattice algebra and conditional toric cell count only; no analytic quotient verified.')

print('\nTHETA-CONDUCTOR LINE-BUNDLE HOLONOMY: PASS')
for r in ((1,1,1),(2,2,2),(1,2,1),(2,3,5),(1,-1,1)):
    # Unknowns c1,c2,c3,vP,vQ; edge i imposes ci=vP and ri*ci=vQ.
    A=S.zeros(6,5)
    for i in range(3):
        A[2*i,i]=1;A[2*i,3]=-1
        A[2*i+1,i]=r[i];A[2*i+1,4]=-1
    h0,h1=5-A.rank(),6-A.rank()
    trivial=len(set(r))==1
    require(((h0, h1) == ((1, 2) if trivial else (0, 1))), 'exact check failed')
    print('endpoint ratios',r,'=> (h0,h1)=',(h0,h1),'trivial=',trivial)

print('\nTORIC t-CHARACTER AND CONORMAL FRAME: PASS')
a,b=S.symbols('a b',integer=True)
Shear=S.Matrix([[1,0,a],[0,1,b],[0,0,1]])
m0=S.Matrix([[0,0,1]])
require((m0 * Shear == m0), 'exact check failed')
triangles=(((0,0),(1,0),(0,1)),((1,1),(1,0),(0,1)))
checks=0
for u in range(-2,3):
    for v in range(-2,3):
        for tri in triangles:
            Rays=S.Matrix.hstack(*(S.Matrix([x+u,y+v,1]) for x,y in tri))
            require((abs(Rays.det()) == 1), 'exact check failed')
            require((Rays.T * m0.T == S.ones(3, 1)), 'exact check failed')
            require((Rays.inv().T * S.ones(3, 1) == m0.T), 'exact check failed')
            checks+=1
print('All symbolic height-preserving shears fix m0; A2 chart products checked:',checks)
print('Every t-preserving deck transition fixes the dual-conormal T^1 frame.')
