"""Independent literal graph checks for fixed-C5 tail activation."""
import contextlib,io
with contextlib.redirect_stdout(io.StringIO()):
 from independent import subset_count
from itertools import combinations
import json

def padded(l,r,coremask,flipleft=(),flipright=()):
 n=l+5+r;G=[[False]*(n+1) for _ in range(n+1)]
 for i in range(n):
  out=(i>=l+5) or (l<=i<l+5 and bool(coremask>>(i-l)&1))
  if i in flipleft or i in flipright:out=not out
  G[0][i+1]=out;G[i+1][0]=not out
  for j in range(i+1,n):
   f=(j-i)%5<=2 if l<=i<j<l+5 else True
   G[i+1][j+1]=f;G[j+1][i+1]=not f
 return G
checks=0
for mask in range(32):
 k=mask.bit_count();h=subset_count(padded(0,0,mask))
 for d in range(4):
  L=subset_count(padded(d+1,0,mask,flipleft=(0,)))
  R=subset_count(padded(0,d+1,mask,flipright=(5+d,)))
  assert L==h+(120-24*k if d==0 else 2**d*(120-21*k)),(mask,d,'L',h,L)
  assert R==h+(24*k if d==0 else 2**d*(21*k+15)),(mask,d,'R',h,R)
  checks+=2
# Complete nearest 1+1 padding boundary cube; no producer cut forms.
vals=[subset_count(padded(1,1,(m>>1)&31,flipleft=(0,) if m&1 else (),flipright=(6,) if not(m>>6&1) else ())) for m in range(128)]
island=[m for m in range(128) if not(m&1) and m>>6&1 and ((m>>1)&31).bit_count() in (2,3)]
allmin=999;tailmin=999
for x in island:
 for y in range(128):
  if y in island or not(1<=(x^y).bit_count()<=2):continue
  d=abs(vals[y]-vals[x]);allmin=min(allmin,d)
  if (x^y)&65:tailmin=min(tailmin,d)
assert (allmin,tailmin)==(8,14),(allmin,tailmin)
print(json.dumps({'status':'PASS','single_anomaly_checks':checks,'whole_128_cut_boundary':{'minimum_exit':allmin,'minimum_tail_exit':tailmin}},indent=2))
