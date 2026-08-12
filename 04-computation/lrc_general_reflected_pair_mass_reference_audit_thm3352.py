#!/usr/bin/env python3
"""General exact floor-moment reflected pair overlap, arbitrary level ratio.

This extends THM-3349's low-ratio two-lift triangle evaluator.  When the
periodized tent radius spans several residue moduli, every positive and
negative lift is summed separately; each lift is still one residue prefix or
suffix and is evaluated by Euclidean floor moments.

Proved scope: L >= 168 with 14 | L, integral cell, 1 <= e,f <= 14, and
positive integral levels p,q.  The implementation also works more broadly
whenever its displayed denominators are positive, but that is not claimed.
"""
from fractions import Fraction as F

def floor_moments(n,m,a,b):
 if n==0:return 0,0,0
 s1=n*(n-1)//2;s2=n*(n-1)*(2*n-1)//6
 qa,a0=divmod(a,m);qb,b0=divmod(b,m)
 base0=qa*s1+qb*n;base1=qa*s2+qb*s1;base2=qa*qa*s2+2*qa*qb*s1+qb*qb*n
 if a0==0:return base0,base1,base2
 height=(a0*(n-1)+b0)//m
 if height==0:return base0,base1,base2
 u0,u1,u2=floor_moments(height,a0,m,m-b0+a0-1)
 r0=n*height-u0;r1=height*s1-(u2-u0)//2;r2=n*height*height-2*u1-u0
 return base0+r0,base1+r1,base2+2*qa*r1+2*qb*r0+r2

def residue_prefix(n,m,a,b,threshold,base):
 """Count and sum residues (at+b mod m) strictly below threshold."""
 if threshold<=0:return 0,0
 if threshold>=m:
  total=a*n*(n-1)//2+b*n-m*base[0]
  return n,total
 shifted=floor_moments(n,m,a,b+m-threshold)
 d0,d1=shifted[0]-base[0],shifted[1]-base[1]
 y0d=(shifted[2]-base[2]-d0)//2
 high_sum=a*d1+b*d0-m*y0d
 total=a*n*(n-1)//2+b*n-m*base[0]
 return n-d0,total-high_sum

def triangle_sum_general(n,m,a,b,peak,L,base,total):
 if peak<=0:return 0
 radius=(peak-1)//L;answer=0
 # Nonnegative lifts: |r+k m|=r+k m.
 for k in range(radius//m+1):
  threshold=min(m,radius-k*m+1)
  count,summation=residue_prefix(n,m,a,b,threshold,base)
  answer+=(peak-L*k*m)*count-L*summation
 # Negative lifts -k: |r-km|=km-r, k>=1.
 for k in range(1,(radius+m-1)//m+1):
  threshold=max(0,k*m-radius)
  before_count,before_sum=residue_prefix(n,m,a,b,threshold,base)
  count,summation=n-before_count,total-before_sum
  answer+=(peak-L*k*m)*count+L*summation
 return answer

def mass(L,cell,e,p,f,q):
 if L*p-e>L*q-f:return mass(L,cell,f,q,e,p)
 z,w=L*p-e,L*q-f
 r,s=e*cell%L,f*cell%L
 determinant=r*w-s*z
 if determinant%L:raise RuntimeError(('nonintegral phase',L,cell,e,p,f,q,determinant))
 b,a=(determinant//L)%z,w%z
 base=floor_moments(p,z,a,b)
 total=a*p*(p-1)//2+b*p-z*base[0]
 unit=L//14
 outer,inner=unit*(z+w),unit*(w-z)
 numerator=triangle_sum_general(p,z,a,b,outer,L,base,total)
 numerator-=triangle_sum_general(p,z,a,b,inner,L,base,total)
 # The floor moments sum nominal indices k=0,...,p-1.  The actual clipped
 # teeth can also include negative/upper boundary indices, and nominal
 # endpoint indices can be empty. Replace the nominal index set by the exact
 # support test for the first clause.  Only O(1) indices differ.
 def one_tooth(k,peak):
  residue=(a*k+b)%z
  return sum(max(0,peak-L*abs(residue+lift*z))
             for lift in range(-(peak//(L*z)+2),peak//(L*z)+3))
 # Boundary teeth require clipping against [0,1]; the unbounded trapezoid
 # expression is invalid there.  Remove every nominal boundary tooth and add
 # the exact literal clipped intersection for every actual boundary tooth.
 radius=F(L,14*z);radius2=F(L,14*w)
 actual=[];lo=(-r)//L-1;hi=(z-r)//L+1
 for k in range(lo,hi+1):
  centre=F(r+k*L,z)
  if centre+radius>0 and centre-radius<1:actual.append((k,centre))
 nominal=set(range(p))
 boundary_nominal={k for k,c in actual if k in nominal and (c-radius<0 or c+radius>1)}
 boundary_actual=[(k,c) for k,c in actual if k not in nominal or k in boundary_nominal]
 for k in (nominal-{k for k,c in actual})|boundary_nominal:
  numerator-=one_tooth(k,outer)-one_tooth(k,inner)
 for k,centre in boundary_actual:
  left,right=max(F(0),centre-radius),min(F(1),centre+radius)
  literal=F(0)
  # Conservative covering range for every second tooth meeting [left,right].
  # The exact minimal integer range is floor(A)+1,...,ceil(B)-1 for
  # A=(w*(left-radius2)-s)/L and B=(w*(right+radius2)-s)/L.
  nlo=((w*(left-radius2)-s)//L)-1
  nhi=((w*(right+radius2)-s)//L)+2
  for n in range(nlo,nhi+1):
   c2=F(s+n*L,w)
   literal+=max(F(0),min(right,c2+radius2)-max(left,c2-radius2))
  numerator+=literal*z*w
 return F(numerator,z*w)

if __name__=='__main__':
 print('general reflected pair floor-moment engine; import mass(L,cell,e,p,f,q)')
