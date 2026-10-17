#!/usr/bin/env python3
"""Exact finite atlas of spindle difference contact parameters.

A missing witness means only missing from the declared fourteen-direction
specialization bank. It does not assert no rational point on a linkage.
"""
import json
from fractions import Fraction as Q
from escaped_ports_20261017d import Field,need


def sign33(a,b):
 # sign(a+b sqrt33), exact over Q.
 if not a:return (b>0)-(b<0)
 if not b:return (a>0)-(a<0)
 if (a>0)==(b>0):return (a>0)-(a<0)
 q=a*a-33*b*b
 return ((a>0)-(a<0))*((q>0)-(q<0))

def main():
 f=Field([[3,0,1],[11,0,1]]);h,j=f.gen(0),f.gen(1)
 cj=f.substitution([f.scale(h,-1),f.scale(j,-1)]);bar=lambda x:f.apply(x,cj);norm=lambda x:f.mul(x,bar(x))
 eta=f.scale(f.add(f.one,h),Q(1,2));rho=f.scale(f.add(f.scalar(5),j),Q(1,6));d=f.add(f.one,eta)
 M=[f.zero,f.one,eta,d,rho,f.mul(rho,eta),f.mul(rho,d)]
 diffs=sorted({f.sub(a,b) for a in M for b in M if a!=b});ns=[norm(a) for a in diffs]
 U=[a for a,n in zip(diffs,ns) if n==f.one]
 need(len(U)==14,'fourteen oriented spindle unit directions')
 def sign(v):
  need(v[1]==v[2]==0,'real quadratic norm coordinate');return sign33(v[0],-v[3])
 unit_diff={f.sub(w,s):(i,k) for i,s in enumerate(U) for k,w in enumerate(U)}
 rows=[]
 for di,a in enumerate(diffs):
  for ei,e in enumerate(diffs):
   nd,ne=ns[di],ns[ei]
   if nd==ne:continue
   # Some complex unit rotation can make |d+u e| <=2 iff
   # (sqrt(nd)-sqrt(ne))²<=4. Keep exact endpoint inequalities.
   q=f.sub(f.add(nd,ne),f.scalar(4))
   if sign(q)>0 and sign(f.sub(f.mul(q,q),f.scale(f.mul(nd,ne),4)))>0:continue
   witnesses=[]
   for ri,r in enumerate(U):
    t=f.add(a,f.mul(r,e));witness=None
    if t in unit_diff:
     si,wi=unit_diff[t];witness=[ri,si,wi]
    else:
     for si,s in enumerate(U):
      if norm(f.add(t,s))==f.one:witness=[ri,si,None];break
    if witness is not None:
     rr,ss=U[witness[0]],U[witness[1]]
     need(norm(rr)==norm(ss)==norm(f.add(a,f.add(f.mul(rr,e),ss)))==f.one,'port-fixing specialization certificate')
     witnesses.append(witness)
   rows.append(dict(d=di,e=ei,witness=witnesses[0] if witnesses else None,
                    allowed_rotation_mask=sum(1<<w[0] for w in witnesses),rotation_witnesses=witnesses))
 # Compatibility concerns the SAME first rotation for multiple couplers.
 masks={}
 for i,row in enumerate(rows):masks.setdefault(row['allowed_rotation_mask'],i)
 disjoint=[]
 for a,i in sorted(masks.items()):
  for b,j in sorted(masks.items()):
   if i<j and a&b==0:disjoint.append([i,j])
 # Prefer a pair with the same first displacement d, to keep its source port.
 same_first=[p for p in disjoint if rows[p[0]]['d']==rows[p[1]]['d']]
 chosen=(same_first or disjoint)[:1]
 compatibility=dict(distinct_rotation_masks=len(masks),disjoint_mask_pairs=len(disjoint),
                    same_first_displacement_pairs=len(same_first),witness_row_pairs=chosen,
                    warning='No common specialization in the fourteen-direction bank, not in the full base field')

 return dict(status='FINITE-EXACT declared parameter/specialization atlas',seed_vertices=7,
  differences=[[str(c) for c in a] for a in diffs],norms=[[str(c) for c in a] for a in ns],
  units=[[str(c) for c in a] for a in U],eligible=len(rows),specialized=sum(r['witness'] is not None for r in rows),
  residual=[r for r in rows if r['witness'] is None],rows=rows,
  compatibility=compatibility,scope='nonzero spindle differences, unequal norms, geometric contact feasible; r,s from fourteen spindle unit directions; output any exact base-field unit')

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
