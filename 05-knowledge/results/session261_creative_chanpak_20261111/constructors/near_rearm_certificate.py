#!/usr/bin/env python3
from module_tools import *
import json
HERE=Path(__file__).resolve().parent
T=(6,4,8,19,103,79,15);O=20
if m.data(T,O)!=(25,24,146):raise RuntimeError('balanced seed')
R=tuple(sum(1<<i for i in range(7) if T[i]>>j&1) for j in range(7))
if m.data(R,O)[2]!=146:raise RuntimeError('balanced endpoint response')
f=m.make_family('q146_binary_debt',T,O,5,'none',lambda k:(10*2**k+5,9-4*2**k,74*2**k-2))
reverse_family=m.make_family('q146_opposite_endpoint',R,O,5,'none',lambda k:(10*2**k+5,9-4*2**k,24*3**k+38*2**k-2))
def tail(B,O,t):
 n=len(B);out=list(B)+[O]+[0]*t
 for i in range(n):
  out[i]|=((1<<t)-1)<<(n+1)
  if not O>>i&1:out[i]|=1<<n
 for j in range(t):out[n+1+j]=(1<<n)|sum(1<<(n+1+k) for k in range(j+1,t))
 return tuple(out)
controls=[]
for mm in range(3):
 B,P=m.inflate(T,O,5,mm)
 for tt in range(1,4):
  Z=tail(B,P,tt);hh=sum(m.endpoints(Z)[-1]);expected=(74*2**tt-4)*2**mm+9-2**(tt+1)
  if hh!=expected:raise RuntimeError('literal tail mismatch')
  controls.append({'m':mm,'t':tt,'order':len(Z),'H':hh})
out={'status':'PASS','seed_H_W_U':[25,24,146],'family':f,'opposite_endpoint_family':reverse_family,'literal_tail_controls':controls,'scope':'Physical binary +2 response with nonunit defect after first clone, converted by a literal tail to the stated two-parameter Hamilton family; not a legal unit rearm.'}
HERE.joinpath('near_rearm_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
