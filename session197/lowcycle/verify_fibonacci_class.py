"""Exact canonical rooted class controls. No removable assertions."""
import itertools as it,json,pathlib
from functools import lru_cache
from verify import canonical,extend_classes,strong,count_cycles,require
P=pathlib.Path(__file__).resolve().parent

def rooted_canonical(rows,root):
 n=len(rows);colors=[(int(u==root),r.bit_count()) for u,r in enumerate(rows)]
 pal={v:i for i,v in enumerate(sorted(set(colors)))};colors=[pal[x] for x in colors]
 while True:
  sig=[(colors[u],tuple(sorted(colors[v] for v in range(n) if rows[u]>>v&1))) for u in range(n)]
  pal={v:i for i,v in enumerate(sorted(set(sig)))};nc=[pal[s] for s in sig]
  if len(set(nc))==len(set(colors)):break
  colors=nc
 cells=[[root]]+[[u for u in range(n) if nc[u]==c and u!=root] for c in sorted(set(nc))]
 cells=[x for x in cells if x];best=None
 for ps in it.product(*(it.permutations(c) for c in cells)):
  p=sum(ps,());pos={u:i for i,u in enumerate(p)}
  r=tuple(sum(1<<pos[v] for v in range(n) if rows[u]>>v&1) for u in p)
  if best is None or r<best:best=r
 return best

def starts(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for mask in range(1,1<<n):
  for v in range(n):
   if mask>>v&1:
    rest=mask^(1<<v)
    d[mask][v]=sum(d[rest][u] for u in range(n) if r[v]>>u&1) if rest else 1
 return d[-1]

def build(U,w,bits):
 m=len(U);n=m+w+1;root=n-1;r=list(U)+[0]*(w+1)
 for u in range(m):
  for x in range(m,m+w):r[u]|=1<<x
  if u:r[u]|=1<<root
  else:r[root]|=1<<u
 for i in range(w):
  x=m+i
  for j in range(i+1,w):r[x]|=1<<(m+j)
  if bits>>i&1:r[root]|=1<<x
  else:r[x]|=1<<root
 return rooted_canonical(r,root)

def induced(r,keep):
 return [sum(1<<j for j,v in enumerate(keep) if r[u]>>v&1) for u in keep]

def main():
 F=[0,1]
 for _ in range(202):F.append(F[-1]+F[-2])
 b=[0,1,0]
 for n in range(3,101):
  b.append(b[-1]+sum(2**(n-m-2)*b[m] for m in range(1,n-1)))
  require(b[-1]==F[2*n-4],'all-order recurrence control')
 bank={1:{(0,)},2:set()};stats=[];physical=0;twobijections=0
 for n in range(1,10):
  if n>=3:
   bank[n]=set()
   for m in range(1,n):
    w=n-m-1
    if w==0 and m==1:continue
    for U in bank[m]:
     for bits in range(1<<max(0,w-1)):
      r=build(U,w,bits)
      require(r not in bank[n],'canonical decomposition collision')
      bank[n].add(r)
   require(len(bank[n])==b[n],'rooted class count')
  unrooted={canonical(r) for r in bank[n]}
  c=sum(r[0].bit_count()==1 for r in bank[n]);d=0
  if n>=3:
   require(c==sum(b[1:n]),'outdegree rooted count')
   require(len(unrooted)==(1 if n==3 else b[n-1]),'unrooted count')
  for r in bank[n]:
   require(strong(r) and starts(r)[0]==1,'constructed start-one graph')
   require(n==1 or count_cycles(r)==1,'constructed cycle count');physical+=1
  if n>=4:
   images=set()
   for S in bank[n-1]:
    if S[0].bit_count()!=1:continue
    r=list(S)+[1]
    for v in range(1,n-1):r[v]|=1<<(n-1)
    images.add(canonical(r))
   for r in unrooted:
    low=[v for v,row in enumerate(r) if row.bit_count()==1]
    require(1<=len(low)<=2,'lowdegree alternatives')
    if len(low)==2:
     d+=1;u,v=low
     if not(r[u]>>v&1):u,v=v,u
     keep=[i for i in range(n) if i!=u];S=induced(r,keep)
     key=rooted_canonical(S,keep.index(v))
     require(key in bank[n-1] and key[0].bit_count()==1,'inverse two-vertex map')
     twobijections+=1
   target={r for r in unrooted if sum(x.bit_count()==1 for x in r)==2}
   require(images==target,'two-vertex image equality')
   require(d==sum(b[1:n-1]),'two-vertex count')
  stats.append({'n':n,'rooted':len(bank[n]),'outdegree_one_rooted':c,'unrooted':len(unrooted),'two_outdegree_one':d})
 ts={(0,)};complete=[]
 for n in range(1,8):
  if n>1:ts=extend_classes(ts)
  actual=set()
  for r in ts:
   if not strong(r):continue
   for v,h in enumerate(starts(r)):
    if h==1:actual.add(rooted_canonical(r,v))
  require(actual==bank[n],'complete tournament universe mismatch')
  complete.append({'n':n,'all_tournament_classes':len(ts),'rooted_start_one':len(actual)})
 result={'status':'PASS','physical_rooted_graphs':physical,'two_vertex_inverse_checks':twobijections,'recurrence_max_order':100,'class_inventory':stats,'independent_complete_inventory':complete}
 text=json.dumps(result,indent=2,sort_keys=True)+'\n';(P/'fibonacci_class_certificate.json').write_text(text);print(text,end='')
if __name__=='__main__':main()
