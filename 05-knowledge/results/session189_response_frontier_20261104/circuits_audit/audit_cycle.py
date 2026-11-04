"""Independent labeled DAG incidence census; constant-region topology only."""
import itertools,json

def cc(v,edges):
 a={x:set() for x in v}
 for x,y in edges:
  if x in a and y in a:a[x].add(y);a[y].add(x)
 rem=set(v);k=0
 while rem:
  k+=1;stack=[rem.pop()]
  while stack:
   for y in a[stack.pop()]&rem:rem.remove(y);stack.append(y)
 return k

def beta(v,e):return len(e)-len(v)+cc(v,e) if v else 0
R=0;by={}
for n in range(2,7):
 count=0;opts=[list(itertools.chain.from_iterable(itertools.combinations(range(j),k) for k in range(min(2,j)+1))) for j in range(n)]
 for parents in itertools.product(*opts):
  edges=[(u,v) for v in range(n) for u in parents[v]]
  live={n-1}
  for v in range(n-1,-1,-1):
   if v in live:live.update(parents[v])
  if len(live)!=n:continue
  out=[sum(u==i for u,v in edges) for i in range(n)]
  for mask in range(1,1<<(n-1)):
   C={i for i in range(n) if mask>>i&1};V=set(range(n));rem=V-C
   if any(parents[v] and len(parents[v])==2 and set(parents[v])<=C for v in rem):continue
   kap=cc(rem,edges);a=sum(u in rem and v in C for u,v in edges);E=sum(out[v]-1 for v in C)
   diff=beta(V,edges)-beta(rem,[(u,v) for u,v in edges if u in rem and v in rem]);Q=E+a-kap+1
   if diff!=Q or Q<E or E<0:raise RuntimeError((n,parents,C,diff,Q,E))
   count+=1
 by[n]=count;R+=count
print(json.dumps({'status':'PASS','regions':R,'by_order':by},sort_keys=True))
