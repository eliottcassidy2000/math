"""Square-sum edge inventory: independent enumeration, standard library."""
from math import isqrt
import json
def need(v,msg):
 if not v:raise RuntimeError(msg)
def inv(n):
 oo=ee=cross=0
 for i in range(1,n+1):
  for j in range(i+1,n+1):
   s=i+j
   if isqrt(s)**2==s:
    if i%2 and j%2:oo+=1
    elif i%2==0 and j%2==0:ee+=1
    else:cross+=1
 return oo,ee,cross
def predicted(n):return isqrt((n-1)//2) if n%2 else isqrt((n+1)//4)
rows=[]
for n in range(1,501):
 oo,ee,x=inv(n);d=oo-ee
 need(d==predicted(n),(n,oo,ee,d,predicted(n)))
 if n in [1,2,15,16,17,24,25,26,99,100,101,500]:rows.append((n,oo,ee,x,d))
print(json.dumps({'status':'PASS','n_through':500,'samples':rows},indent=2))
