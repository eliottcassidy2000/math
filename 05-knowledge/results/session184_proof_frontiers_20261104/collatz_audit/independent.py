"""Independent native lift-tree audit, without importing the author's code."""
from pathlib import Path
from itertools import product
import json

def need(x,msg):
 if not x: raise RuntimeError(msg)
def orbit(n,w):
 x=n;first=0
 for i,a in enumerate(w,1):
  z=3*x+1;b=0
  while z%2==0:z//=2;b+=1
  if b!=a:return None
  x=z
  if not first and x<n:first=i
 return x,first
# Each child residue is discovered by literal native-map testing of lifts,
# rather than by the affine inverse used in the proposed theorem.
levels=[((),1,2,0,1,1)]
words=sources=0
for length in range(1,7):
 nxt=[]
 for w,r,mod,B,P,Q in levels:
  for a in range(1,5):
   ww=w+(a,);mm=mod*2**a
   good=[r+k*mod for k in range(2**a) if orbit(r+k*mod,ww) is not None]
   need(len(good)==1,'unique native residue')
   rr=good[0];BB=3*B+Q;PP=P*3;QQ=Q*2**a
   need((PP*rr+BB)%(2*QQ)==QQ,'claimed residue matches literal lift')
   if PP<QQ:
    # Only retain words whose first multiplicative crossing is their end.
    threshold=BB//(QQ-PP)+1
    n=rr+max(0,(threshold-rr+mm-1)//mm)*mm
    for t in range(3):
     result=orbit(n+t*mm,ww)
     need(result is not None and result[1]==length,'exact first descent')
     sources+=1
    words+=1
   else:nxt.append((ww,rr,mm,BB,PP,QQ))
 levels=nxt
# Confirm displayed examples using literal odd trajectories.
for q,r,n,y in [(1,1,43,37),(5,3,720891,548047)]:
 need(orbit(n,(1,2)*q+(2,)*r)==(y,2*q+r),'displayed example')
# Rotation existence independently checks every expanding word through length8.
rotations=0
for L in range(1,9):
 for w in product(range(1,4),repeat=L):
  if 3**L <= 2**sum(w):continue
  candidates=[w[k:]+w[:k] for k in range(L)]
  need(any(all(3**j>=2**sum(v[:j]) for j in range(1,L+1)) for v in candidates),'rotation')
  rotations+=1
print(json.dumps({'status':'PASS','slope_first_crossing_words_through_length6_letters1to4':words,'native_sources':sources,'expanding_rotation_words_through_length8_letters1to3':rotations},sort_keys=True))
