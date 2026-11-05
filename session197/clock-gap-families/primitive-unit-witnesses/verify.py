import itertools,json,math
from pathlib import Path
P=Path(__file__).resolve().parent
def need(b,s):
 if not b:raise RuntimeError(s)
def graph(j,m):
 w=[1,0,1]+[0]*j+[2,0,2];n=3+len(w)+m;r=[0]*n
 for a,b in ((0,1),(0,2),(2,1)):r[a]|=1<<b
 for i,t in enumerate(w):
  v=3+i
  for u in range(3):
   if u<2 and (t>>u)&1:r[u]|=1<<v
   else:r[v]|=1<<u
  for k in range(i+1,len(w)):r[v]|=1<<(3+k)
  for k in range(m):r[v]|=1<<(3+len(w)+k)
 for k in range(m):
  v=3+len(w)+k;r[v]|=7
  for z in range(k+1,m):r[v]|=1<<(3+len(w)+z)
 return r
def validate(r):
 n=len(r);need(all(0<=x<1<<n and not x>>i&1 for i,x in enumerate(r)),'rows')
 need(all((r[i]>>j&1)+(r[j]>>i&1)==1 for i in range(n) for j in range(i+1,n)),'tournament')
def H(r):
 n=len(r);N=1<<n;d=[[0]*n for _ in range(N)]
 for v in range(n):d[1<<v][v]=1
 for mask in range(1,N):
  for v in range(n):
   x=d[mask][v]
   if not x:continue
   w=r[v]&~mask
   while w:
    b=w&-w;w-=b;d[mask|b][b.bit_length()-1]+=x
 return sum(d[-1]) if n else 1
def runs(r,root):
 v=[i for i in range(len(r)) if i!=root];n=len(v);N=1<<n;d=[[[0]*5 for _ in range(n)] for _ in range(N)]
 for i,u in enumerate(v):d[1<<i][i][0]=1;d[1<<i][i][1]=r[root]>>u&1
 for mask in range(1,N):
  for i,u in enumerate(v):
   if not mask>>i&1:continue
   for j,z in enumerate(v):
    if mask>>j&1:continue
    a=d[mask][i];b=d[mask|1<<j][j]
    if r[u]>>z&1:
     for k in range(5):b[k]+=a[k]
    if (r[u]>>root&1) and (r[root]>>z&1):
     for k in range(4):b[k+1]+=a[k]
 a=[0]*5
 for i,u in enumerate(v):
  for k in range(5):
   a[k]+=d[-1][i][k]
   if k<4 and r[u]>>root&1:a[k+1]+=d[-1][i][k]
 return a
def module_masks(r):
 n=len(r);out=[]
 for mask in range(1,1<<n):
  if not 1<mask.bit_count()<n:continue
  vs=[v for v in range(n) if mask>>v&1]
  if all(all((r[v]>>u&1)==(r[vs[0]]>>u&1) for v in vs) for u in range(n) if not mask>>u&1):out.append(mask)
 return out
def automorphisms(r):
 n=len(r);degrees=[x.bit_count() for x in r];groups=[[i for i in range(n) if degrees[i]==d] for d in sorted(set(degrees))];out=[];attempts=0
 for pieces in itertools.product(*(itertools.permutations(g) for g in groups)):
  f=[0]*n
  for g,p in zip(groups,pieces):
   for u,v in zip(g,p):f[u]=v
  attempts+=1
  if all((r[i]>>j&1)==(r[f[i]]>>f[j]&1) for i in range(n) for j in range(n)):out.append(f)
 return out,attempts
def strong(r):
 n=len(r);F=(1<<n)-1
 for v in range(n):
  reach=1<<v
  while True:
   q=reach
   for u in range(n):
    if reach>>u&1:reach|=r[u]
   if q==reach:break
  if reach!=F:return False
 return True
def prime_certificate(p):
 need(p>=2,'prime candidate');divisors=list(range(2,math.isqrt(p)+1));need(all(p%d for d in divisors),'primality')
 return {'value':p,'trial_divisor_max':math.isqrt(p),'tested_divisors':len(divisors)}
def main():
 cases=[]
 for j,m in ((0,1),(1,1),(0,2),(0,3),(0,6)):
  r=graph(j,m);validate(r);h=H(r);A=260*2**j;beta=-1 if j==0 else 1
  need(h==A*3**m+beta,'literal Hamilton count');need(strong(r),'strong witness')
  z={'j':j,'m':m,'rows':r,'H':h,'A':A,'beta':beta,'strong':True}
  if m==1:
   a=runs(r,len(r)-1);need(a==[A+beta,3*A+beta,3*A,A,0],'full root run profile')
   mods=module_masks(r);need(not mods,'module indecomposable source');aut,attempts=automorphisms(r);need(aut==[list(range(len(r)))],'rigid source')
   z.update(run=a,nontrivial_modules=mods,automorphisms=aut,degree_preserving_permutations=attempts)
  if (j,m) in ((0,2),(0,3)):z['prime']=prime_certificate(h)
  if (j,m) in ((1,1),(0,6)):
   need(h%7==0,'seven times prime');z['seven_times_prime']=prime_certificate(h//7)
  cases.append(z)
 need(cases[0]['rows']==[46,320,2,1014,999,966,901,775,517,7],'explicit order-ten rows')
 result={'status':'PASS','cases':cases,'scope':'Two explicit primitive rigid unit-clock sources and finite prime/7p witnesses; no unbounded primality claim.'}
 text=json.dumps(result,sort_keys=True,indent=2)+'\n';print(text,end='')
if __name__=='__main__':main()
