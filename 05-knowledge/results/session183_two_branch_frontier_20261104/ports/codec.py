"""Recompute inherited field codecs, never infer them from source phase.
Codec copied mathematically from pinned return6325_completion/transport.py.
Input counts independently supplied by ports/joint.cpp.
"""
import json,hashlib
from pathlib import Path
HERE=Path(__file__).parent

def need(x):
 if not x:raise RuntimeError('codec check failed')
def mul(a,b):
 c=[0]*5
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 for i in (4,3):c[i-1]-=c[i];c[i-3]-=c[i]
 return tuple(x%7 for x in c[:3])
def power(a,k):
 out=(1,0,0)
 while k:
  if k&1:out=mul(out,a)
  a=mul(a,a);k//=2
 return out
def chart(g):
 h,x,y=g;z=tuple(v%7 for v in g);a,b,c=(x+4*(h-y))%7,y%7,(h-y)%7
 need(power(z,343)==z);need(((b+c)%7,(a-4*c)%7,b)==z)
 return {'observer_mod7':z,'F343_vector_codec':[z,power(z,7),power(z,49)],'F343_codec_degree':next(i for i in (1,3) if power(z,7**i)==z),'F49_charge':[a,b,c],'F49_frobenius':[(a+4*b)%7,(-b)%7,c]}
