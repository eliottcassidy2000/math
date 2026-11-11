#!/usr/bin/env python3
from module_tools import *
from itertools import permutations
import json
HERE=Path(__file__).resolve().parent

def literal_tensor(T,O):
 n=len(T);out=[[0]*4 for _ in range(n+1)]
 for p in permutations(range(n)):
  bad=sum(1<<i for i in range(n-1) if not T[p[i]]>>p[i+1]&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad!=bad:continue
   bounds=[0]+[i+1 for i in range(n-1) if cuts>>i&1]+[n];r=len(bounds)-1
   st=bool(O>>p[0]&1);en=bool(O>>p[bounds[1]-1]&1)
   two=r>=2 and en and bool(O>>p[bounds[1]]&1)
   for j,z in enumerate((1,st,en,two)):out[r][j]+=z
 return out

def main():
 T=(2,4,9,19,7);left=9;right=18
 a=literal_tensor(T,left);b=literal_tensor(T,right)
 rho=tuple(reversed(range(5)))
 if not all((T[i]>>j&1)==(T[rho[j]]>>rho[i]&1) for i in range(5) for j in range(5) if i!=j):raise RuntimeError('anti-automorphism')
 for p in permutations(range(5)):
  if sum(1<<p[i] for i in range(5) if left>>i&1)==right and all((T[i]>>j&1)==(T[p[i]]>>p[j]&1) for i in range(5) for j in range(5)):
   raise RuntimeError('unexpected marked automorphism')
 if a!=b:raise RuntimeError('marked profile mismatch')
 if m.data(T,left)!=(9,8,32) or m.data(T,right)!=(9,8,32):raise RuntimeError('base seed')
 # Exact finite skeleton proof of the inherited mixed-port formula in this gauge.
 fam=m.make_family('reselected_q32_mixed_port',T,right,0,'none',lambda k:(3*2**k+3,1,4*3**k+9*2**k+2) if k else (5,-1,15))
 old_binary=m.make_family('balanced_old_binary_response',T,left,4,'none',lambda k:(3*2**k+3,1,16*2**k) if k else (5,-1,15))
 hostile=[]
 for k in (1,2,3):
  Z,P=m.inflate(T,left,0,k,'first');h,w,u=m.data(Z,P);hostile.append((k,h,h-w,u))
 if hostile[-1][2]!=3:raise RuntimeError('selected-port boundary')
 out={'status':'PASS','rows':T,'old_selected_mask':left,'new_selected_mask':right,'tensor_rows_C_S_E_K':a,'base_H_W_U':[9,8,32],'activated_port':fam,'old_selected_first_controls':hostile,'old_binary_response':old_binary,'reselection_anti_automorphism':rho,'scope':'Full marked-tensor equality gives all-host replacement; the mixed ladder identity is inherited and certified in this vertex gauge.'}
 HERE.joinpath('marked_reselection_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
