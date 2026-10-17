#!/usr/bin/env python3
"""Independent four-coordinate reconstruction of every specialization mask."""
import json,argparse
from fractions import Fraction as F
from itertools import product
from pathlib import Path


def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def add(a,b):return tuple(x+y for x,y in zip(a,b))
def neg(a):return tuple(-x for x in a)
def sub(a,b):return add(a,neg(b))
def scale(a,c):return tuple(x*c for x in a)
def mul(a,b):
 out=[F(0)]*4
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   common=i&j;k=i^j
   out[k]+=x*y*((-11) if common&1 else 1)*((-3) if common&2 else 1)
 return tuple(out)
def norm(a):return mul(a,(a[0],-a[1],-a[2],a[3]))
def sign(a):
 need(a[1]==a[2]==0,'real element')
 x,y=a[0],-a[3]
 if x==0:return (y>0)-(y<0)
 if y==0 or x*y>0:return (x>0)-(x<0)
 return ((x>0)-(x<0))*((x*x>33*y*y)-(x*x<33*y*y))


def main(path):
 d=json.loads(Path(path).read_text());one=(F(1),F(0),F(0),F(0));zero=scale(one,0)
 eta=(F(1,2),F(0),F(1,2),F(0));rho=(F(5,6),F(1,6),F(0),F(0));a=add(one,eta)
 M=[zero,one,eta,a,rho,mul(rho,eta),mul(rho,a)]
 D=sorted({sub(x,y) for x in M for y in M if x!=y});U=[x for x in D if norm(x)==one]
 unpack=lambda xs:[tuple(F(c) for c in x) for x in xs]
 need(D==unpack(d['differences']) and U==unpack(d['units']),'independently reconstructed ordered banks')
 universe=[];four=scale(one,4)
 for i,x in enumerate(D):
  for j,y in enumerate(D):
   A,B=norm(x),norm(y)
   if A==B:continue
   q=sub(add(A,B),four)
   if sign(q)>0 and sign(sub(mul(q,q),scale(mul(A,B),4)))>0:continue
   universe.append((i,j))
 need(universe==[(r['d'],r['e']) for r in d['rows']],'complete filtered parameter universe')
 norm_checks=0
 for row in d['rows']:
  x,y=D[row['d']],D[row['e']];mask=0
  for ri,r in enumerate(U):
   candidates=[norm(add(x,add(mul(r,y),s)))==one for s in U];norm_checks+=len(U)
   if any(candidates):mask|=1<<ri
  need(mask==row['allowed_rotation_mask'] and mask!=0,'every positive and negative mask bit exact')
  for ri,si,wi in row['rotation_witnesses']:
   out=add(x,add(mul(U[ri],y),U[si]));need(norm(out)==one,'stored source-preserving witness')
   if wi is not None:need(out==U[wi],'stored output direction')
 masks={}
 for i,row in enumerate(d['rows']):masks.setdefault(row['allowed_rotation_mask'],i)
 pairs=[(i,j) for a,i in masks.items() for b,j in masks.items() if i<j and a&b==0]
 c=d['compatibility']
 need(len(masks)==c['distinct_rotation_masks'] and len(pairs)==c['disjoint_mask_pairs'],'exact compatibility quotient counts')
 need(sum(d['rows'][i]['d']==d['rows'][j]['d'] for i,j in pairs)==c['same_first_displacement_pairs'],'same marked first displacement')
 for i,j in c['witness_row_pairs']:need((i,j) in pairs,'declared incompatible pair')
 return dict(status='FINITE-EXACT independent full mask verification',parameters=len(universe),norm_checks=norm_checks,
             distinct_rotation_masks=len(masks),incompatible_representative_pairs=len(pairs))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('input');a=p.parse_args();print(json.dumps(main(a.input),sort_keys=True))
