#!/usr/bin/env python3
"""Complete fixed-place obstruction to simultaneous marked retractions.

Exact rational intervals certify geometric parameters with 0<N<4.
All finite-field masks are complete over F29² and checked by literal pairs.
"""
import json
from fractions import Fraction as Q
from pathlib import Path
from itertools import combinations
import math


def need(ok,msg):
 if not ok:raise RuntimeError(msg)
p=29

def add(a,b):return ((a[0]+b[0])%p,(a[1]+b[1])%p)
def neg(a):return (-a[0]%p,-a[1]%p)
def mul(a,b):return ((a[0]*b[0]+2*a[1]*b[1])%p,(a[0]*b[1]+a[1]*b[0])%p)
def norm(a):return (a[0]*a[0]-2*a[1]*a[1])%p
def power(a,n):
 r=(1,0)
 for _ in range(n):r=mul(r,a)
 return r

def red(v):
 # Basis1,j,h,hj; h=10t,j=3t,t²=2, so hj=2 mod29.
 roots=[(1,0),(0,3),(0,10),(2,0)];r=(0,0)
 for c,x in zip(v,roots):
  q=c.numerator*pow(c.denominator,-1,p)%p;r=add(r,(q*x[0]%p,q*x[1]%p))
 return r

# Rational interval arithmetic, used only for the selected real embedding.
def iadd(a,b):return (a[0]+b[0],a[1]+b[1])
def ineg(a):return (-a[1],-a[0])
def imul(a,b):
 v=[x*y for x in a for y in b];return min(v),max(v)
def idiv(a,b):
 need(b[0]>0 or b[1]<0,'nonzero interval denominator');return imul(a,(1/b[1],1/b[0]))
def iconst(q):return (Q(q),Q(q))
def sqrt_interval(n,bits=40):
 den=1<<bits;k=math.isqrt(n*den*den);return Q(k,den),Q(k+1,den)
def w_interval():
 lo,hi=Q(1),Q(2)
 poly=lambda w:w**3+2*w*w-4*w-4
 for _ in range(50):
  mid=(lo+hi)/2
  if poly(mid)<0:lo=mid
  else:hi=mid
 need(poly(lo)<0<poly(hi),'cubic root bracket; derivative positive on[1,2]')
 return lo,hi

def complex_interval(v):
 # a+bj+ch+dhj, with h=i sqrt3,j=i sqrt11.
 a,b,c,d=v;s3=sqrt_interval(3);s11=sqrt_interval(11);s33=sqrt_interval(33)
 return iadd(iconst(a),imul(iconst(-d),s33)),iadd(imul(iconst(b),s11),imul(iconst(c),s3))
def cmul(a,b):return iadd(imul(a[0],b[0]),ineg(imul(a[1],b[1]))),iadd(imul(a[0],b[1]),imul(a[1],b[0]))
def cadd(a,b):return iadd(a[0],b[0]),iadd(a[1],b[1])
def unit_interval():
 h=iadd(w_interval(),iconst(3));h2=imul(h,h);den=iadd(h2,iconst(11))
 return idiv(iadd(h2,iconst(-11)),den),idiv(imul(iconst(-2),imul(h,sqrt_interval(11))),den)

def main():
 source=Path('04-computation/experiments/marked_contact_current_20261017e/catalog.json');data=json.loads(source.read_text())
 V=[tuple(Q(c) for c in a) for a in data['differences']];U=sorted((a,b) for a in range(p) for b in range(p) if norm((a,b))==1)
 need(len(U)==30,'complete norm-one group')
 diff={add(a,neg(b)) for a in U for b in U};maskrows=[];candidates=[];rot=unit_interval()
 for i,row in enumerate(data['rows']):
  d,e=red(V[row['d']]),red(V[row['e']]);mask=0
  for ri,r in enumerate(U):
   z=add(d,mul(r,e));literal=any(norm(add(z,t))==1 for t in U)
   need(literal==(z in diff),'literal two-circle enumeration equals difference-set criterion')
   if literal:mask|=1<<ri
  need(mask!=0,'single coupler has reduction of a known specialization')
  maskrows.append(mask)
  z=cadd(complex_interval(V[row['d']]),cmul(rot,complex_interval(V[row['e']])))
  N=iadd(imul(z[0],z[0]),imul(z[1],z[1]))
  # Interval arithmetic may overestimate the square across zero; the selected
  # lower endpoint here is positive, so it is nevertheless a valid certificate.
  if N[0]>0 and N[1]<4:candidates.append((i,mask,N))
 # Exact candidates already have rational real-embedding certificates.
 unique={}
 for i,mask,N in candidates:unique.setdefault(mask,(i,N))
 need(all(a&b for a,b in combinations(unique,2)),'all locally admissible pairs have a common residue rotation')
 # Prefer pairwise characteristic-zero compatibility, not just residue
 # compatibility. This makes the obstruction genuinely ternary for the
 # displayed contact equations. It does not classify other joint edges.
 chosen=None
 for ix,(i,im,_) in enumerate(candidates):
  ib=data['rows'][i]['allowed_rotation_mask']
  for jx in range(ix+1,len(candidates)):
   j,jm,_=candidates[jx];jb=data['rows'][j]['allowed_rotation_mask']
   if not ib&jb:continue
   for k,km,_ in candidates[jx+1:]:
    kb=data['rows'][k]['allowed_rotation_mask']
    if ib&kb and jb&kb and not im&jm&km:chosen=[i,j,k];break
   if chosen:break
  if chosen:break
 need(chosen is not None,'find a simultaneous local obstruction')
 selected=[];common=(1<<30)-1
 for i in chosen:
  row=data['rows'][i];N=next(N for j,mask,N in candidates if i==j);common&=maskrows[i]
  selected.append(dict(row=i,d=row['d'],e=row['e'],d_coefficients=data['differences'][row['d']],e_coefficients=data['differences'][row['e']],
    residue_d=red(V[row['d']]),residue_e=red(V[row['e']]),allowed_mask=maskrows[i],allowed_rotations=[r for j,r in enumerate(U) if maskrows[i]>>j&1],
    norm_interval=[str(q) for q in N],norm_interval_decimal=[float(q) for q in N]))
 need(common==0,'no common rotation anywhere in the full finite norm-one group')
 pair_witnesses=[]
 for i,j in combinations(chosen,2):
  shared=data['rows'][i]['allowed_rotation_mask']&data['rows'][j]['allowed_rotation_mask']
  ri=(shared&-shared).bit_length()-1
  witnesses=[next(w for w in data['rows'][k]['rotation_witnesses'] if w[0]==ri) for k in (i,j)]
  pair_witnesses.append(dict(rows=[i,j],rotation=ri,witnesses=witnesses))
 # This also authenticates the extension of the fixed place from the
 # biquadratic coefficient field to K=Q(zeta15,sqrt(-11)).
 zeta3=(14,5) # (-1+h)/2, h=10t.
 zeta15=next(z for z in U if power(z,15)==(1,0) and power(z,5)==zeta3 and power(z,3)!=(1,0))
 need(power(zeta15,29)==(zeta15[0],-zeta15[1]%p),'Frobenius equals complex conjugation')
 # The local obstruction is global for a K-linear port-fixing graph retraction:
 # norm-one r,t are units at every conjugation-fixed place.
 return dict(status='PROVED finite obstruction + exact rational geometric intervals; joint forcing OPEN',prime=p,
  residue_field='F29[t]/(t²-2)',base_place=dict(sqrt_minus3=[0,10],sqrt_minus11=[0,3],zeta15=zeta15),
  complete_unit_group=U,parameters=len(maskrows),real_certified_parameters=len(candidates),distinct_masks=len(unique),
  allowed_masks=maskrows,selected=selected,common_mask=common,pairwise_characteristic_zero_witnesses=pair_witnesses,
  pairwise_scope='displayed contact equations only; no claim of pairwise full unit-graph retractions',
  restriction='rules out simultaneous K-linear port-fixing graph retractions; no claim about nonlinear colourings or missing six-colour boundaries')

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
