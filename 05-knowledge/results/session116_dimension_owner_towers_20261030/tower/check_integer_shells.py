"""Exact finite shell audit: integer zero divisors in W have count r3(m)r4(m)."""
from itertools import product
from collections import Counter
from math import isqrt
from pathlib import Path
import json
LIMIT=128;gates=0
def need(v,s):
 global gates
 gates+=1
 if not v:raise RuntimeError(s)
def reps(k):
 C=Counter()
 for v in product(range(-isqrt(LIMIT),isqrt(LIMIT)+1),repeat=k):
  n=sum(t*t for t in v)
  if n<=LIMIT:C[n]+=1
 return C
r3,r4=reps(3),reps(4)
rows=[]
for m in range(1,LIMIT+1):
 q=m
 while q%4==0:q//=4
 need((r3[m]>0)==(q%8!=7),'Legendre shell existence')
 jacobi=8*sum(d for d in range(1,m+1) if m%d==0 and d%4)
 need(r4[m]==jacobi,'Jacobi shell count')
 rows.append({'m':m,'norm_squared':2*m,'r3':r3[m],'r4':r4[m],'zero_divisor_elements_in_W':r3[m]*r4[m]})
need(rows[0]['zero_divisor_elements_in_W']==48,'first shell')
need(rows[6]['zero_divisor_elements_in_W']==0,'norm14 hostile')
out={'status':'PASS','gates':gates,'shell_limit':LIMIT,'counts_include_signs_and_coordinate_order':True,'rows':rows}
Path(__file__).with_name('integer_shell_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='rows'},sort_keys=True))
