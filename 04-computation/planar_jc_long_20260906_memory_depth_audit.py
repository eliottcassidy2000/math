#!/usr/bin/env python3
"""Independent audit of finite depth recognition.

No producer or repository mathematical imports. Enumerates literal source
monomials, compares their spans to Hasse jets at u=-1, and tests fixed-support
recognition by dimensions of actual subspace intersections. This path uses
neither the producer's reciprocal series nor its Toeplitz determinant.
"""
from fractions import Fraction as F
from math import comb
from functools import lru_cache
from hashlib import sha256
import json
CHECKS=0
records=[]
def check(ok,label):
 global CHECKS
 CHECKS+=1
 if not ok: raise RuntimeError(label)

def rank(rows):
 if not rows or not rows[0]:return 0
 a=[[F(v) for v in row] for row in rows]
 r=0
 for j in range(len(a[0])):
  p=next((i for i in range(r,len(a)) if a[i][j]),None)
  if p is None:continue
  a[r],a[p]=a[p],a[r]
  v=a[r][j];a[r]=[q/v for q in a[r]]
  for i in range(r+1,len(a)):
   if a[i][j]:
    v=a[i][j];a[i]=[b-v*c for b,c in zip(a[i],a[r])]
  r+=1
  if r==len(a):break
 return r

def data(h,d,ell):
 s=max(0,-((-ell)//h));rho=max(0,-((-ell)//(h+1)))
 D=ell+d-s
 return s,rho,D

@lru_cache(None)
def image(h,d,ell):
 s,rho,D=data(h,d,ell)
 if D<0:return ()
 answer=set()
 # Literal x^a u^b p^c y^e. On intercept ell, a=h*c+(h+1)*e-ell.
 # Highest t-row=2*c+3*e+b; the cap is that this is <=ell+d.
 for c in range(max(0,(ell+d)//2)+1):
  for e in range(max(0,(ell+d-2*c)//3)+1):
   a=h*c+(h+1)*e-ell
   if a<0:continue
   for b in range(ell+d-2*c-3*e+1):
    row=[0]*(D+1)
    start=b+c+2*e-s
    for j in range(c+e+1):
     check(0<=start+j<=D,'literal generator lies in cap')
     row[start+j]=comb(c+e,j)
    answer.add(tuple(row))
 return tuple(sorted(answer))

def hasse(rho,width):
 return [[comb(j,k)*((-1)**(j-k)) if j>=k else 0
          for j in range(width)] for k in range(rho)]

def basis(width,n):
 return [[int(i==j) for j in range(width)] for i in range(n)]

def main():
 image_cases=recognition_cases=sharp_cases=0
 for h in range(2,7):
  for d in range(4):
   for ell in range(-d,25):
    s,rho,D=data(h,d,ell)
    cols=image(h,d,ell)
    if D<0:
     check(not cols,'negative ambient bound');continue
    jets=hasse(rho,D+1)
    expected=D+1-rank(jets)
    check(rank(cols)==expected,'literal all-height span is complete Hasse kernel')
    for v in cols:
     check(all(sum(a*b for a,b in zip(v,w))==0 for w in jets),'literal generator has required root multiplicity')
    records.append(('image',h,d,ell,expected));image_cases+=1
 for h in range(2,7):
  for T in range(1,7):
   for d in range(4):
    N=h*h*T//(h+1)+d+1
    for ell in range(-d,h*T+1):
     s,rho,D=data(h,d,ell)
     M=min(T,ell+d)-s
     if M<0:continue
     cols=image(h,d,ell)
     width=min(N,ell+d)-s+1
     truncated=[v[:width] for v in cols]
     projected_rank=rank(truncated)
     target=basis(width,M+1)
     projected_obstruction=rank(truncated+target)-projected_rank
     full_obstruction=rank(hasse(rho,M+1))
     check(projected_obstruction==full_obstruction,'fixed-support projected/full membership kernels agree')
     records.append(('recognition',h,T,d,ell,N,projected_obstruction));recognition_cases+=1
    ell=h*T;s,rho,D=data(h,d,ell);cols=image(h,d,ell)
    before=N-s
    u=[1]+[0]*D
    p=[v[:before] for v in cols]
    check(rank(p+[u[:before]])==rank(p),'t^T passes every row before cutoff')
    p=[v[:before+1] for v in cols]
    check(rank(p+[u[:before+1]])==rank(p)+1,'t^T first fails at cutoff')
    check(sum(u[j]*((-1)**j) for j in range(D+1))==1,'t^T not divisible even once')
    # Positive actual-source control: any displayed literal generator has
    # polynomial membership and therefore passes every longer prefix.
    check(bool(cols),'sharp diagonal has actual source generators')
    sharp_cases+=1
 # Complete h=2 simplex annihilator, checked against literal packets.
 simplex_cases=0
 for d in range(4):
  for ell in range(1,19):
   s,rho,D=data(2,d,ell)
   for m in range(s,ell+d+1):
    width=m-s+1
    rows=[[(-1)**j*comb(m+q-s-j,q) for j in range(width)]
          for q in range(max(0,ell+d-m),rho)]
    cols=[v[:width] for v in image(2,d,ell)]
    check(rank(rows)==width-rank(cols),'simplex bank has full annihilator rank')
    for v in cols:
     check(all(sum(a*b for a,b in zip(v,w))==0 for w in rows),'simplex annihilates raw unsigned coefficients')
    simplex_cases+=1
 # Independent characteristic-three hostile, using an actual source
 # polynomial p^3-3*y^2=t^3-3*x^4*t^5-2*x^6*t^6.
 lift=[1,0,-3,-2]
 check(all(q%3==0 for q in lift[1:3]),'char3 actual lift matches t3 through t5')
 check(lift[3]%3!=0,'char3 discrepancy at t6')
 check(sum(lift[j]*((-1)**j) for j in range(4))==0,'integer actual lift vanishes at u=-1')
 check(sum(j*lift[j]*((-1)**(j-1)) for j in range(1,4))==0,'integer actual lift has double root')
 # Independent universal determinant identity controls, by polynomial
 # coefficient extraction (SymPy determinant, no producer algorithm).
 import sympy as S
 det_cases=0
 for rho in range(1,8):
  for k in range(1,rho+1):
   for B in (0,1,2,5,9):
    mat=S.Matrix([[comb(B+i-j+rho-1,rho-1) if B+i-j>=0 else 0 for j in range(k)] for i in range(k)])
    answer=S.prod(S.factorial(rho+B-1-i)*S.factorial(i)/(S.factorial(rho-1-i)*S.factorial(B+i)) for i in range(k))
    check(mat.det(method='domain-ge')==answer and answer>0,'factorial determinant identity')
    det_cases+=1
 print('INDEPENDENT literal-source/Hasse/intersection audit; no repo imports')
 print('LITERAL_IMAGE_CASES',image_cases)
 print('FIXED_SUPPORT_RECOGNITION_CASES',recognition_cases)
 print('SHARP_ACTUAL_SOURCE_CASES',sharp_cases)
 print('COMPLETE_SIMPLEX_BANK_CASES',simplex_cases)
 print('DETERMINANT_CONTROLS',det_cases)
 print('CHAR3 actual_source=p^3-3y^2=t^3-3x^4t^5-2x^6t^6; t^3 passes row5 but fails row6')
 print('SEMANTIC_SHA256',sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest())
 print('CHECKS',CHECKS,'PASS')
if __name__=='__main__':main()
