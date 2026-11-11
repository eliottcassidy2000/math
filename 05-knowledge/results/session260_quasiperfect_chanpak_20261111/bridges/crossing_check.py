from functools import lru_cache
from itertools import combinations
from pathlib import Path
import json

def graph(a,b):
 n=a+b+1; out=[0]*n
 for i,j in combinations(range(n),2):
  if j==n-1 and i<a:out[j]|=1<<i
  else:out[i]|=1<<j
 return out

def h(out):
 n=len(out);full=(1<<n)-1
 @lru_cache(None)
 def f(mask,v):
  if mask==full:return 1
  avail=out[v]&~mask;ans=0
  while avail:
   bit=avail&-avail;avail-=bit;ans+=f(mask|bit,bit.bit_length()-1)
  return ans
 return sum(f(1<<i,i) for i in range(n))

def response(a,b,i,j):
 A=2**(a-1);B=2**(b-1);d=a-i;e=j-1;M=(A-1)*(B-1)
 L=2*(A>>d)*(B>>e)-(d==0)*(B>>e)-(e==0)*(A>>d)+2*(d==e==0)
 X=2**(i-1)*3**d-2**d;Y=2**(b-j)*3**e-2**e
 G=2*X*Y+2**(d+e+1)-(j==b)*2**d-(i==1)*2**e
 return L,G

def check(q,s):
 if not q:raise RuntimeError(s)
checks=0;rows=[]
for a in range(2,6):
 for b in range(2,6):
  out=graph(a,b);H=h(out);negative=[]
  for i in range(1,a+1):
   for j in range(1,b+1):
    alt=out.copy();u=i-1;v=a+j-1;alt[u]^=1<<v;alt[v]^=1<<u
    L,G=response(a,b,i,j);delta=h(alt)-H
    check(delta==G-L,(a,b,i,j,delta,G-L));checks+=1
    check((delta<0)==(i==a and j==1),('sign',a,b,i,j,delta))
    if delta<0:negative.append({'i':i,'j':j,'decrease':-delta})
  rows.append({'a':a,'b':b,'H':H,'decreasing_crossings':negative})
result={'status':'FINITE-EXACT','literal_subset_DP_crossing_cases':checks,'range':'2<=a,b<=5; every crossing edge','rows':rows}
Path(__file__).with_name('crossing_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
