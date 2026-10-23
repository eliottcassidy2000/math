#!/usr/bin/env python3
"""Independent audit: divisor recurrence and exact Q(zeta_24) arithmetic.
No producer source imports; producer certificate is compared only as data.
"""
from fractions import Fraction as F
from pathlib import Path
import json,mpmath as mp
N=30
sigma=[0]+[sum(d for d in range(1,n+1) if n%d==0) for n in range(1,N+1)]
p=[1]
for n in range(1,N+1):
 v=sum(sigma[k]*p[n-k] for k in range(1,n+1));assert v%n==0;p.append(v//n)
coeff=[]
for a in (0,1):
 h=F(a,4);v=[0]*(N+1)
 for n in range(-10,11):
  norm=F((2*n+a)**2,4);level=norm-h;assert level.denominator==1
  if level<=N:
   for j in range(int(level),N+1):v[j]+=p[j-int(level)]
 coeff.append(v)
expected=json.loads((Path(__file__).resolve().parents[1]/'vertex/LATTICE_VOA_CERTIFICATE.json').read_text())
assert coeff==[expected['vacuum_coefficients'],expected['nonvacuum_coefficients']]
# Q[z]/(Phi_24(z)) with Phi_24=z^8-z^4+1. No radical simplifier.
def elt(a):return tuple(map(F,a))+tuple(F(0) for _ in range(8-len(a)))
zero=elt([]);one=elt([1]);z=elt([0,1])
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,s):return tuple(s*x for x in a)
def mul(a,b):
 c=[F(0)]*15
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 for k in range(14,7,-1):c[k-4]+=c[k];c[k-8]-=c[k]
 return tuple(c[:8])
def power(a,n):
 v=one
 for _ in range(n):v=mul(v,a)
 return v
assert power(z,24)==one
root2=add(power(z,3),power(z,21));assert mul(root2,root2)==scale(one,2)
s=scale(root2,F(1,2));i=power(z,6);t=power(z,23)
S=[[s,s],[s,scale(s,-1)]];theta=[[one,zero],[zero,i]];T=[[t,zero],[zero,mul(t,i)]];identity=[[one,zero],[zero,one]]
def mm(A,B):return [[add(mul(A[a][0],B[0][b]),mul(A[a][1],B[1][b])) for b in range(2)] for a in range(2)]
def cube(A):return mm(mm(A,A),A)
assert mm(S,S)==identity and cube(mm(S,T))==identity
anomaly=power(z,3);assert cube(mm(S,theta))==[[anomaly,zero],[zero,anomaly]]
# Since inverse S0j is sqrt2, no field division needed.
fusion=[]
for a in range(2):
 for b in range(2):
  for c in range(2):
   v=zero
   for j in range(2):v=add(v,mul(root2,mul(S[a][j],mul(S[b][j],S[c][j]))))
   assert v==(one if c==(a^b) else zero);fusion.append([a,b,c,int(c==(a^b))])
# Distinct analytic control uses library theta/q-Pochhammer, not producer finite sums.
mp.mp.dps=90
E=lambda x:mp.exp(2j*mp.pi*x)
def chars(tau):
 q=E(tau);eta=E(tau/24)*mp.qp(q)
 return [mp.jtheta(3,0,q)/eta,mp.jtheta(2,0,q)/eta]
checks=[]
for tau in [mp.mpc('0.11','0.91'),mp.mpc('-0.07','1.31'),mp.mpc('0','0.58')]:
 v=chars(tau);w=chars(-1/tau);u=chars(tau+1)
 # jtheta(2) fractional nome branch must be tracked across T; evaluate via S only.
 err=max(abs(w[0]-(v[0]+v[1])/mp.sqrt(2)),abs(w[1]-(v[0]-v[1])/mp.sqrt(2)))
 assert err<mp.mpf('1e-75');checks.append({'tau':str(tau),'S_error':mp.nstr(err,8)})
result=dict(status='PASS',coefficient_method='divisor-sum recurrence, independent of producer partition and pentagonal methods',coefficient_level=N,coefficients=coeff,exact_field='Q(zeta_24), modulus x^8-x^4+1',S_square='I',character_ST_cube='I',ribbon_ST_cube='zeta_24^3 I',fusion=fusion,numerical_controls=checks,numerical_status='not interval certified; analytic proof separately audited')
Path(__file__).with_name('VOA_ADDENDUM_AUDIT_CERTIFICATE.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS: independent exact coefficients, cyclotomic matrices, Verlinde, and numerical S controls')
