#!/usr/bin/env python3
"""Independent arithmetic, exact Fourier, and alternate Euler certificate audit.
Does not import producer code; no assert statements; normal and -O equivalent.
"""
from pathlib import Path
from math import comb, gcd, isqrt, factorial
from fractions import Fraction
from collections import Counter
import json
from mpmath import iv, mp
ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/'actual-manifest-copy/results.json').read_text())
checks=0

def require(ok,msg):
 global checks
 checks+=1
 if not ok: raise RuntimeError(msg)

def prime(n):
 return n>1 and all(n%d for d in range(2,isqrt(n)+1))

def hasse_E(p):
 m=(p-1)//2
 return sum(comb(m,i)*comb(m-i,i)*pow(4,m-2*i,p)*pow(2,i,p) for i in range(m//2+1))%p

def hasse_short(p,a,b):
 m=(p-1)//2
 return sum(comb(m,i)*comb(m-i,2*m-3*i)*pow(a,2*m-3*i,p)*pow(b,2*i-m,p) for i in range((m+1)//2,2*m//3+1))%p

E=[]
for row in DATA['CM_curve']['rows']:
 p=row['p']; H=hasse_E(p)
 require(prime(p),'prime input')
 require(H==row['trace']%p,'Hasse E equals trace mod p')
 require((H==0)==row['supersingular'],'Hasse E agrees with supersingularity')
 require((H==0)==(p%8 in [5,7]),'Hasse E agrees with CM character')
 m=(p-1)//2
 Hminus=sum(comb(m,i)*comb(m-i,i)*pow(-4,m-2*i,p)*pow(2,i,p) for i in range(m//2+1))%p
 require(Hminus==((-1)**m*H)%p,'Bolza E-minus quadratic twist Hasse relation')
 E.append({'p':p,'hasse':H})
all_j=[]
for row in DATA['Ogg_finite_certificate']:
 if row['p']>107: continue
 p=row['p']; ss=[]
 for j in range(p):
  if j==0: a,b=0,1
  elif j==1728%p: a,b=1,0
  else:
   k=j*pow((1728-j)%p,-1,p)%p; a,b=3*k%p,2*k%p
  require((4*a**3+27*b*b)%p!=0,'nonsingular short model')
  require(1728*4*a**3*pow((4*a**3+27*b*b)%p,-1,p)%p==j,'model j computed explicitly')
  if hasse_short(p,a,b)==0: ss.append(j)
 require(ss==row['rational_supersingular_j'],'Hasse all-j vs point-count census')
 all_j.append({'p':p,'supersingular_j':ss})

# Exact group-algebra Fourier filter via Ramanujan coefficients, no floats.
# c_15(j)=sum_{d|gcd(15,j)} d*mu(15/d).
def ramanujan(j):
 return sum(d*mu for d,mu in [(1,1),(3,-1),(5,-1),(15,1)] if j%d==0)
# Sum_j c15(j)*zeta^(jn) evaluated exactly modulo Phi15(z).
phi=[1,-1,0,1,-1,1,0,-1,1]
def mod_phi(poly):
 a=poly[:]
 while len(a)>=len(phi):
  c=a[-1]; offset=len(a)-len(phi)
  for i,v in enumerate(phi): a[offset+i]-=c*v
  while a and a[-1]==0: a.pop()
 return a+[0]*(8-len(a))
for n in range(1,121):
 a=[0]*15
 for j in range(15): a[(j*n)%15]+=ramanujan(j)
 require(mod_phi(a)==[15*int(gcd(n,15)==1)]+[0]*7,'root-of-unity projection exact')
require([ramanujan(j) for j in range(15)]==[8,1,1,-2,1,-4,-2,1,1,-2,-4,1,-2,1,1],'Ramanujan coefficients')

# Independent finite-difference Euler realization, using a different split and
# truncation. Uses explicit repeated differences, not producer binomial weights.
iv.dps=155
M,N=512,128
bounds={}
for d in (1,3,5,15):
 f=lambda n:iv.exp(iv.log(iv.mpf(d)*n)/(d*n))-1
 prefix=sum(((-1)**n*f(n) for n in range(1,M)),iv.mpf(0))
 diffs=[f(M+k) for k in range(N)]
 tail=iv.mpf(0)
 for k in range(N):
  tail+=((-1)**k)*diffs[0]/iv.mpf(2)**(k+1)
  diffs=[diffs[i+1]-diffs[i] for i in range(len(diffs)-1)]
 L=iv.log(iv.mpf(3)*d*M/2)+iv.pi/6
 B=2*factorial(N)*iv.exp(2*L/(d*M))/d*(L/iv.mpf(M)**(N+1)+(L/N+iv.mpf(1)/(N*N))/iv.mpf(M)**N)
 bounds[d]=prefix+tail+iv.mpf([-1,1])*B
 # Producer intervals are decimal renderings, parse directed endpoints anew.
 text=DATA['periodic_filter']['M_d_intervals'][str(d)]
 lo,hi=text[1:-1].split(', ')
 require(bool(bounds[d].a>iv.mpf(lo)) and bool(bounds[d].b<iv.mpf(hi)),'alternate interval strictly inside producer interval')
C=bounds[1]; U=C-bounds[3]-bounds[5]+bounds[15]
D=2/iv.pi*iv.exp(-iv.log(8)/8)*iv.sin(iv.pi/8)-C
require(bool(U.a>iv.mpf('.53758442393739534102470436896889584658743294509034679972563051925142108122520')),'independent U lower digits')
require(bool((U-8*C/15).a>iv.mpf('.437')),'independent density rejection')
require(bool(D.a>0) and bool((iv.exp(-4*iv.pi*iv.sqrt(2))-D).a>0),'independent defect/theta rejection')

# Direct pairing, with alternating-series remainders, checks convergence without
# using Euler acceleration at all. Decimal intervals only establish consistency.
mp.dps=40
for d in (1,3,5,15):
 K=4000
 f=lambda n:mp.expm1(mp.log(d*n)/(d*n))
 S=sum((-1)**n*f(n) for n in range(1,K+1))
 # K is even; next term is negative, so S-f(K+1) <= M(d) <= S.
 lo=mp.mpf(str(bounds[d].a)[1:].split(',')[0])
 require(S-f(K+1) < lo < S,'ordinary partial-sum alternating enclosure')

report={'status':'PASS','checks':checks,'E_hasse_rows':E,'all_j_hasse_rows':all_j,
        'ramanujan15':[ramanujan(j) for j in range(15)],
        'alternate_euler':{'M':M,'N':N,'interval_dps':iv.dps,'M_d':{str(d):str(x) for d,x in bounds.items()},'unit_MRB':str(U),'D':str(D)}}
(ROOT/'independent_results.json').write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
print('PASS:',checks,'independent active checks')
print('Hasse E:',len(E),'primes; all-j Hasse census:',len(all_j),'primes')
print('Root-of-unity projection: exact modulo Phi15, n=1..120')
print('Alternate M=512, N=128 interval certificates nested in all producer M(d) intervals')
print('Direct ordinary alternating partial sums confirm all four M(d) independently')
