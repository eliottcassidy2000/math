#!/usr/bin/env python3
"""Independent exact prime-density and unbiased path-Hessian checks."""
from itertools import permutations,combinations
from fractions import Fraction
from math import factorial,isqrt
import json
import sympy as s
out={}
N=10**7
prime=bytearray(b'\1')*N;prime[:2]=b'\0\0'
for p in range(2,isqrt(N-1)+1):
 if prime[p]:prime[p*p:N:p]=b'\0'*(((N-1-p*p)//p)+1)
counts=[sum(prime[10**(k-1):10**k]) for k in range(1,8)]
assert counts==[4,21,143,1061,8363,68906,586081]
a=[Fraction(c,9*10**i) for i,c in enumerate(counts)]
out['decimal_prime_shells']=[{'digits':i+1,'count':counts[i],'density':str(a[i])} for i in range(7)]
# These are only necessary finite moment tests, never a realization proof.
out['density_hankel_tests']=[]
for shift in range(2):
 for size in range(1,(len(a)-shift+1)//2+1):
  M=s.Matrix(size,size,lambda i,j:s.Rational(a[i+j+shift].numerator,a[i+j+shift].denominator))
  det=M.det();out['density_hankel_tests'].append({'shift':shift,'size':size,'determinant':str(det),'positive':bool(det>0)})
out['hessian_path_audits']=[]
for n in range(2,8):
 edges=list(combinations(range(n),2)); idx={e:i for i,e in enumerate(edges)};m=len(edges)
 M=[[0]*m for _ in range(m)]
 for path in permutations(range(n)):
  es=[(idx[tuple(sorted((u,v)))],1 if u<v else -1) for u,v in zip(path,path[1:])]
  for (i,a1),(j,a2) in combinations(es,2):M[i][j]+=a1*a2;M[j][i]+=a1*a2
 B=s.zeros(n,m)
 for j,(u,v) in enumerate(edges):B[u,j]=1;B[v,j]=-1
 H=s.Matrix(M)/2**(n-1);c=s.Rational(factorial(n-2),2**(n-2)); expected=c*(2*s.eye(m)-B.T*B)
 assert H==expected
 eig=H.eigenvals();expected_eig={} if n==2 else {-(n-2)*c:n-1,2*c:(n-1)*(n-2)//2}
 if n==2:expected_eig={s.Integer(0):1}
 assert eig==expected_eig
 out['hessian_path_audits'].append({'n':n,'paths':factorial(n),'arcs':m,'eigenvalues':{str(k):v for k,v in eig.items()}})
# A symmetric function need not have only uniform critical points.
t=s.symbols('t');f=t**4-s.Rational(1,8)*t**2
assert s.diff(f,t).subs(t,s.Rational(1,4))==0
assert s.diff(f,t).subs(t,-s.Rational(1,4))==0
out['symmetric_critical_hostile']={'simplex':'x1+x2=1; t=x1-1/2','function':'t^4-t^2/8','interior_critical_t':['-1/4','0','1/4']}
print(json.dumps(out,indent=2))
