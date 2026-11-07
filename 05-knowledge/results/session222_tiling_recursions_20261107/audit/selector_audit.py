#!/usr/bin/env python3
from itertools import product
from functools import lru_cache
from pathlib import Path
import json

def need(x):
 if not x:raise RuntimeError('audit failure')
@lru_cache(None)
def H(w):
 n=len(w)+1
 E={(i,j) for i in range(n-1) for j in range(i+1,n-1)}|{(i,n-1) if b=='1' else (n-1,i) for i,b in enumerate(w)}
 @lru_cache(None)
 def f(rem,last):
  if not rem:return 1
  return sum(f(rem^(1<<v),v) for v in range(n) if rem>>v&1 and (last,v) in E)
 return sum(f(((1<<n)-1)^(1<<v),v) for v in range(n))
count=0;strong=0
for m in range(3,9):
 for b in product('01',repeat=m):
  w=''.join(b)
  for i in range(1,m-1):
   if w[i:i+2]!='10' or '1' in w[i+3:]:continue
   z=w[:i]+'01'+w[i+2:]
   u=w[:i-1]+'0'+w[i-1:];v=z[:i-1]+'0'+z[i-1:]
   N=1+sum((w[j]=='0')*2**(i-j-2) for j in range(i-1))
   need(H(z)-H(w)==2*N);need(H(v)-H(u)==4*N)
   need(u[:i-1]+u[i:]==w and v[:i-1]+v[i:]==z)
   count+=1
   strong+=int(w[0]==z[0]=='0' and w[-1]==z[-1]=='1')
out={'status':'PASS','universe':'all apex words lengths3..8 and all marked10 positions i>=1 with R=0','cases':count,'strong_cases':strong,'N1_special_counts':[H(s) for s in ['0101','0011','00101','00011']]}
Path(__file__).with_name('SELECTOR_RESULT.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
