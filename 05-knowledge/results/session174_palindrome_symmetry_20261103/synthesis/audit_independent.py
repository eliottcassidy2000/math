from itertools import combinations
import json

def paths(g):
 n=len(g);d={}
 for j in range(n):d[1<<j,j]=1
 totals=[0]*(1<<n);totals[0]=1
 for s in range(1,1<<n):
  for j in range(n):
   v=d.get((s,j),0);totals[s]+=v
   if v:
    for k in range(n):
     if not(s>>k&1) and g[j]>>k&1:d[s|1<<k,k]=d.get((s|1<<k,k),0)+v
 return totals

def wrapper(g):
 n=len(g);a=[0]*(2*n+1)
 for i in range(n):
  for j in range(n):
   if g[i]>>j&1:a[i]|=1<<j;a[n+j]|=1<<(n+i)
  for j in range(n,2*n):a[i]|=1<<j
  a[n+i]|=1<<(2*n);a[2*n]|=1<<i
 return a
count=0;values={}
for n in range(1,5):
 E=list(combinations(range(n),2));seen=set()
 for mask in range(1<<len(E)):
  g=[0]*n
  for k,(i,j) in enumerate(E):
   if mask>>k&1:g[i]|=1<<j
   else:g[j]|=1<<i
  hs=paths(g);allmask=(1<<n)-1;a=hs[-1];b=sum(hs[s]*hs[allmask^s] for s in range(1<<n))
  h=paths(wrapper(g))[-1]
  if h!=(b-a)**2+2*a*a or h%8!=3 or h in (7,21):raise RuntimeError((g,a,b,h))
  seen.add(h);count+=1
 values[n]=sorted(seen)
print(json.dumps({'result':'PASS','all_labelled_tournaments_through_order4':count,'wrapper_H_by_seed_order':values},indent=2))
