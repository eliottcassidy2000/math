#!/usr/bin/env python3
"""Literal F8 / order21 semilinear action and failed evaluation quotient."""
from pathlib import Path
from hashlib import sha256
import json
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def mul(x,y):
 out=0
 while y:
  if y&1:out^=x
  y>>=1;x<<=1
  if x&8:x^=11
 return out
def power(x,n):
 out=1
 for _ in range(n):out=mul(out,x)
 return out
def main():
 V=list(range(8));G=[(a,b) for b in range(3) for a in range(7)]
 S={(0,1),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2),(3,0),(4,1),(4,2)}
 def product(g,h):return ((g[0]+pow(2,g[1],7)*h[0])%7,(g[1]+h[1])%3)
 def action(g,x):return mul(power(2,g[0]),power(x,2**g[1]))
 need({power(2,a) for a in range(7)}==set(range(1,8)),'theta has order seven')
 maps={g:tuple(action(g,x) for x in V) for g in G}
 need(len(set(maps.values()))==21,'faithful action')
 for g in G:
  need(set(maps[g])==set(V),'invertible maps')
  for x in V:
   for y in V:need(action(g,x^y)==action(g,x)^action(g,y),'F2-linearity')
  for h in G:
   for x in V:need(action(g,action(h,x))==action(product(g,h),x),'composition law')
 out={u:{product(u,s) for s in S} for u in G}
 for u in G:
  for v in G:
   if u!=v:need((v in out[u])!=(u in out[v]),'actual tournament pair')
 u=(0,0);up=(0,1);v=(1,1)
 need(action(u,1)==action(up,1)==1,'same evaluation fibre')
 need(v in out[u] and up in out[v],'opposite orientations to same target')
 need(product((0,2),v)==(4,0),'relative element of phase-one source')
 fibre0=[g for g in G if g[0]==0];fibre1=[g for g in G if g[0]==1]
 cross=sum(vv in out[uu] for uu in fibre0 for vv in fibre1)
 need(cross==8,'nonuniform eight versus one seam')
 need([pow(2,b,7) for b in range(3)]==[1,2,4],'Frobenius exponent orbit')
 return {'status':'PASS; complete finite action, not an arithmetic conjecture',
  'checks':checks,'field_modulus_bits':11,'primitive_element_bits':2,
  'group_elements':21,'field_elements':8,'nonzero_points':7,
  'frobenius_exponents':[1,2,4],
  'evaluation_quotient_hostile':{'same_fibre':[u,up],'target':v,'cross_orientation_counts':[cross,9-cross]},
  'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
