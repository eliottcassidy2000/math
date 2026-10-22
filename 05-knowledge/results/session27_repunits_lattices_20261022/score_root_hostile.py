import json
import itertools,json
from collections import defaultdict

def hamilton(out):
 n=len(out); dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,1<<n):
  for v,z in enumerate(dp[s]):
   if z:
    nxt=out[v]&~s
    while nxt:
     bit=nxt&-nxt;u=bit.bit_length()-1;nxt-=bit;dp[s|bit][u]+=z
 return sum(dp[-1])
def Qout(k,mask):
 out=[0]*k
 for i,(u,v) in enumerate(itertools.combinations(range(k),2)):
  a,b=(u,v) if mask>>i&1 else (v,u);out[a]|=1<<b
 return out
def response(out,word):
 k=len(out);o=out[:]+[0]*len(word)
 for i,x in enumerate(word):
  for j in range(k):
   if x>>j&1:o[j]|=1<<(k+i)
   else:o[k+i]|=1<<j
  for j in range(i+1,len(word)):o[k+i]|=1<<(k+j)
 return hamilton(o)

q=7;R={1,2,4};a=[sum(1<<v for v in range(q) if (v-u)%q in R) for u in range(q)];b=a[:]
for u,v in [(0,1),(1,3),(3,0)]:b[u]^=1<<v;b[v]^=1<<u
z=[]
for out in [a,b]:
 full=out[:]+[0]
 for u in range(7):
  if 7>>u&1:full[u]|=1<<7
  else:full[7]|=1<<u
 z.append({'out_masks':full,'H':hamilton(full),'scores':[x.bit_count() for x in full],'core_H':hamilton(out)})
print(json.dumps(z,indent=2))
