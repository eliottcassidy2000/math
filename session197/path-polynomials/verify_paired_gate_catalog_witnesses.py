#!/usr/bin/env python3
from verify import *
from math import isqrt
from pathlib import Path

def least_factor(x):
 if x<2:return 0
 return next((d for d in range(2,isqrt(x)+1) if x%d==0),1)
def opposite(a):return tuple(v^((1<<len(a))-1)^(1<<i) for i,v in enumerate(a))
def strong(a):
 for v in range(len(a)):
  mask=1<<v
  while True:
   nxt=mask
   for i,row in enumerate(a):
    if mask>>i&1:nxt|=row
   if nxt==mask:break
   mask=nxt
  if mask!=(1<<len(a))-1:return False
 return True

def run():
 data=json.loads(Path(__file__).with_name('paired_gate_catalog_normal.json').read_text());need=0;direct=[]
 for item in data['values']:
  a=tuple(item['input_rows']);validate(a);v=(cover_polynomial(a)+[0,0,0])[1:4]
  require(v==item['covers'],'literal witness covers')
  h,u,w=v;H=5*h*h+12*h*u+28*u*u+48*u*w+72*w*w
  require(H==item['H'] and H<data['exclusive_limit'],'exact gate count')
  require(least_factor(H)==item['least_divisor_or_one'],'integer output primality field')
  require((0 if H%7 else least_factor(H//7))==item['seven_quotient_least_divisor_or_one'],'integer quotient primality field');need+=1
  if item['least_divisor_or_one']==1 or item['seven_quotient_least_divisor_or_one']==1:
   g=cyclic(a,opposite(a),transitive(2));validate(g);n=len(a);tau=list(range(n,2*n))+list(range(n))+[2*n+1,2*n]
   require(all(((g[i]>>j)&1)==((g[tau[j]]>>tau[i])&1)for i,j in combinations(range(len(g)),2)),'actual anti-involution')
   require(strong(g) and h_subsets(g)[-1]==H,'direct strong physical prime/7p count')
   direct.append({'H':H,'order':len(g),'rows':g,'seven_prime':item['seven_quotient_least_divisor_or_one']==1})
 require(need==119 and len(direct)==11,'complete catalog physical controls')
 print(json.dumps({'status':'PASS','input_cover_controls':need,'direct_output_controls':direct},sort_keys=True,indent=2))
if __name__=='__main__':run()
