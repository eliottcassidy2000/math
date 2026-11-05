#!/usr/bin/env python3
"""Exact all-order tensor-annihilator certificate; no fitted-recurrence inference."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import hashlib,json
from verify import TRANS,INIT,ACCEPT,require,tensor_path_moments

def multiply(a,b):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]+=x*y
 return out

P=[1]
for root,mult in [(0,3),(36,1),(24,2),(16,2),(12,2),(8,4),(4,3)]:
 for _ in range(mult):P=multiply(P,[-root,1])
kernel=[]
for i in range(36):
 for j in range(36):
  counts=defaultdict(int)
  for c in range(4):
   for a in TRANS[c][i]:
    for b in TRANS[c][j]:counts[a*36+b]+=1
  kernel.append(list(counts.items()))
v=[a*b for a in INIT for b in INIT]
res=[P[0]*x for x in v]
for k in range(1,len(P)):
 nxt=[0]*1296
 for i,x in enumerate(v):
  if x:
   for j,c in kernel[i]:nxt[j]+=x*c
 v=nxt
 for j,x in enumerate(v):res[j]+=P[k]*x
require(not any(res),'full 1296-coordinate Krylov annihilator')

def formula(m):
 return (F(841,25600)*9**m-F(5*(10*m-109),3072)*6**m
  +F(14255*m-71716,345600)*4**m-F(3*m-89,144)*3**m
  -F(4*m**3-24*m**2+68*m-687,3072)*2**m
  +F(1536*m**2-5656*m+9365,27648))
values=tensor_path_moments(100)
for m in range(3,101):require(formula(m)==F(values[m],4**m),'spectral closed formula')
# Fourteen consecutive values identify the positive-root solution after the 3-head.
# The vector annihilator above, not this numerical fit, proves continuation.
receipt={'status':'exact full-vector annihilator certificate',
 'operator':'K=sum_c M_c tensor M_c; vectors have 1296 integer coordinates',
 'polynomial_ascending':P,'degree':len(P)-1,'residual_nonzero_coordinates':sum(bool(x) for x in res),
 'initial_EH_squared':[str(F(values[m],4**m)) for m in range(3)],
 'closed_formula_checked_orders':[3,100],'tail_solution_dimension':14,
 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'transfer_source_sha256':hashlib.sha256(Path(__file__).with_name('verify.py').read_bytes()).hexdigest()}
s=json.dumps(receipt,indent=2,sort_keys=True)+'\n'
Path(__file__).with_name('spectral_receipt.json').write_text(s)
print(s,end='')
