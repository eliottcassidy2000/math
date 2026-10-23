from itertools import combinations, permutations
from collections import Counter,defaultdict
import json

def matrix(n,m):
 M=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(combinations(range(n),2)):
  M[i][j]=1 if m>>k&1 else -1; M[j][i]=-M[i][j]
 return M

def Kkey(M):
 n=len(M); z=n-1
 return tuple(M[i][j]-M[i][z]+M[j][z] for i in range(z) for j in range(i+1,z))
def hamilton(M):
 n=len(M); dp=[{} for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for v,c in dp[mask].items():
   for u in range(n):
    if not mask>>u&1 and M[v][u]==1:
     q=dp[mask|1<<u];q[u]=q.get(u,0)+c
 return sum(dp[-1].values())
def scc(M):
 n=len(M); reach=[sum(1<<j for j in range(n) if M[i][j]>0)|(1<<i) for i in range(n)]
 for k in range(n):
  for i in range(n):
   if reach[i]>>k&1:reach[i]|=reach[k]
 comps=[];todo=set(range(n))
 while todo:
  i=min(todo); c=tuple(j for j in todo if reach[i]>>j&1 and reach[j]>>i&1);comps.append(c);todo.difference_update(c)
 all_comps=tuple(comps)
 comps.sort(key=lambda c: -sum(M[c[0]][d[0]]>0 for d in all_comps if c!=d))
 return comps

def rootpair(M):
 n=len(M);roots=list(permutations(range(n),2));out=Counter()
 for i,j in roots:
  for k,l in roots:
   g=(i==k)-(i==l)-(j==k)+(j==l)
   a=M[i][k]-M[i][l]-M[j][k]+M[j][l]
   out[g,a]+=1
 return tuple(sorted(out.items()))
def paley8():
 C=[[0]*8 for _ in range(8)]
 for i in range(8):
  for j in range(i+1,8):
   C[i][j]=1 if i==0 or (j-i)%7 in {1,2,4} else -1; C[j][i]=-C[i][j]
 return C
if __name__=='__main__':
 results={}
 for n in range(1,7):
  fibres=defaultdict(list); strong=0
  for m in range(1<<(n*(n-1)//2)):
   M=matrix(n,m);c=scc(M);h=hamilton(M)
   fibres[Kkey(M)].append((m,len(c),h))
  for key,f in fibres.items():
   if len({x[2] for x in f})!=1:raise Exception(('H collision',n,key,f))
   if any(x[1]!=len(f) for x in f):raise Exception(('fibre size',n,key,f))
  results[n]={'tournaments':1<<(n*(n-1)//2),'forms':len(fibres),'fibre_sizes':dict(sorted(Counter(map(len,fibres.values())).items()))}
  print(n,results[n],flush=True)
 C=paley8();fp=defaultdict(list)
 for m in range(128):
  d=[1]+[-1 if m>>(i-1)&1 else 1 for i in range(1,8)]
  M=[[d[i]*C[i][j]*d[j] for j in range(8)] for i in range(8)]
  fp[rootpair(M)].append((m,hamilton(M)))
 print('P8 root-pair classes',len(fp),flush=True)
 results['conference_rootpair']=[{'size':len(v),'Hs':dict(Counter(h for _,h in v)),'rep':v[0][0],'data':k} for k,v in fp.items()]
 for r in results['conference_rootpair']:print(r,flush=True)
 with open('session36-lattice/explore.json','w') as f:json.dump(results,f,indent=2)
