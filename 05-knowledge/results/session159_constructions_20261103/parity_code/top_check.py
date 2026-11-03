from itertools import combinations,permutations,product

def direct(n,edges):
 adj=[[i<j for j in range(n)] for i in range(n)]
 for a,b in edges:adj[a][b]=False;adj[b][a]=True
 cyc=[]
 for m in range(3,n+1,2):
  for s in combinations(range(n),m):
   for p in permutations(s[1:]):
    p=(s[0],)+p
    if all(adj[p[i]][p[(i+1)%m]] for i in range(m)):
     cyc.append(sum(1<<v for v in s))
 def pack(i,mask,k):
  if k==0:return 1
  return sum(pack(j+1,mask|cyc[j],k-1) for j in range(i,len(cyc)) if not mask&cyc[j])
 return pack(0,0,len(edges))
def fourier(n,edges):
 endpoints={v for e in edges for v in e}
 total=0
 for signs in product((-1,1),repeat=len(edges)):
  p=1
  for e in signs:p*=e
  for v in range(n):
   if v not in endpoints:p*=1+sum(signs[i] for i,(a,b) in enumerate(edges) if a<v<b)
  total+=p
 if total % (1<<len(edges)):raise ValueError('not divisible')
 return total//(1<<len(edges))
cases=[(5,[(0,3),(1,4)]),(7,[(0,5),(1,6)]),(8,[(0,4),(1,7),(2,6)]),(9,[(0,6),(1,7),(2,8)]),(9,[(0,4),(2,6),(5,8)]),(9,[(0,8),(1,5),(3,7)]),(8,[(0,2),(3,5)]),(8,[(0,7)]),(9,[(0,2),(3,5),(6,8)])]
for n,edges in cases:
 a,b=direct(n,edges),fourier(n,edges)
 if a!=b:raise ValueError((n,edges,a,b))
 print(n,edges,'top coefficient',a)
print('PASS general overlap cases:',len(cases))
