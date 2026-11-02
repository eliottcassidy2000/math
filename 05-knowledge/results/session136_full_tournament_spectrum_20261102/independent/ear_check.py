"""Independent even-path network DP and vertex-transitive band controls."""
from kernel import *
import json,hashlib,pathlib

def even_path_network(a):
 n=len(a);full=(1<<n)-1
 @lru_cache(None)
 def ending(s,v):
  rem=s^(1<<v)
  if not rem:return 1
  return sum(ending(rem,u)for u in range(n)if rem>>u&1 and a[u]>>v&1)
 hs=[1]+[sum(ending(s,v)for v in range(n)if s>>v&1)for s in range(1,full+1)]
 r=[[0]*n for _ in a]
 for start in range(n):
  d=[{}for _ in range(1<<n)];d[1<<start][start]=1
  for mask,row in enumerate(d):
   for end,count in row.items():
    if not mask.bit_count()%2:r[start][end]+=count*hs[full^mask]
    targets=a[end]&(full^mask)
    while targets:
     bit=targets&-targets;targets-=bit;new=mask|bit;v=bit.bit_length()-1
     d[new][v]=d[new].get(v,0)+count
 return r,hs

# Literal validation on every parent through4; strong n5 singleton boundary.
count=0;singletonfails=[]
for n in range(2,6):
 for code in range(1<<(n*(n-1)//2)):
  a=tournament(n,code);r,hs=even_path_network(a);full=(1<<n)-1
  fs={s:hs[-1]+2*sum(r[i][j]for i in range(n)if s>>i&1 for j in range(n)if not s>>j&1)for s in range(1,full)}
  if n<=4:
   for s,f in fs.items():require(h_literal(extend(a,s))==f,('network literal',n,code,s));count+=1
  if len(scc(a))==1:
   if min(f for s,f in fs.items()if s.bit_count() in(1,n-1))!=min(fs.values()):singletonfails.append((n,code))
require(not singletonfails,'small singleton hostile boundary')

controls=[]
for n in (5,7,9,11):
 a=tuple(sum(1<<((i+j)%n)for j in range(1,(n+1)//2))for i in range(n));r,hs=even_path_network(a);full=(1<<n)-1
 rho=min(hs[full^(1<<i)^(1<<j)]for i in range(n)for j in range(i+1,n))
 fs={s:hs[-1]+2*sum(r[i][j]for i in range(n)if s>>i&1 for j in range(n)if not s>>j&1)for s in range(1,full)}
 low=min(fs.values());cuts=[s for s,f in fs.items()if f<=low+rho-1];lowvals=sorted({fs[s]for s in cuts})
 bound=0;p=1
 while p*3<=n:p*=3;bound+=1
 require(len(lowvals)<=bound,'VT low band')
 if n in(5,7,11):
  require(len(cuts)==2*n and len(lowvals)==1,'primitive low band')
 for s in (1,2,3,full//3,full-1):require(h_recursive(extend(a,s))==fs[s],('direct cyclic ear',n,s))
 controls.append({'n':n,'H':hs[-1],'rho':rho,'minimum_ear':low,'maximum_ear':max(fs.values()),'low_cut_count':len(cuts),'low_values':lowvals,'exact_missing_low_band':(rho+1)//2-len(lowvals)})
print(json.dumps({'status':'PASS','literal_child_checks':count,'no_singleton_minimum_failure_all_strong_through5':True,'cyclic_controls':controls,'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()},sort_keys=True,indent=2))
