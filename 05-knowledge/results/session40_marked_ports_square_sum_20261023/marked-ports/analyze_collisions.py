import itertools,json,collections
from pathlib import Path
P=Path(__file__).resolve().parent
def matrix(code,n=8):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)): A[i][j]=code>>(n*(n-1)//2-1-k)&1;A[j][i]=1-A[i][j]
 return A
def refine(A,s,t):
 n=len(A);c=[0 if v==s else 1 if v==t else 2 for v in range(n)];history=[]
 while True:
  k=[(c[v],tuple(sum(A[v][w] for w in range(n) if c[w]==i) for i in range(max(c)+1))) for v in range(n)]
  history.append(tuple(sorted(k)));values={z:i for i,z in enumerate(sorted(set(k)))};d=[values[z] for z in k]
  if len(set(c))==len(set(d)):break
  c=d
 cells=[[v for v in range(n) if d[v]==i] for i in range(max(d)+1)]
 return history,cells

def sign(A,cells):
 triples=[b for b in cells if len(b)==3]
 if len(triples)!=2:return None
 X,Y=triples;q=sum(A[X[0]][v] for v in Y)
 if q in (0,3):return None
 pi={u:next(v for v in Y if A[u][v]==(q==1)) for u in X}
 u=X[0];w=next(v for v in X if A[u][v]);return 1 if A[pi[u]][pi[w]] else -1

def path(A,s,t):
 n=len(A)
 def rec(p,mask):
  if len(p)==n:return p if p[-1]==t else None
  for v in range(n):
   if not(mask>>v&1) and A[p[-1]][v] and (v!=t or len(p)==n-1):
    r=rec(p+[v],mask|1<<v)
    if r:return r
 return rec([s],1<<s)
if __name__=="__main__":
 rows=[]
 for line in (P/'endpoint_scan.txt').read_text().splitlines():
  if not line.startswith('collision '):continue
  G=[]
  for spec in line.split()[1:]:
   code,s,t=map(int,spec.split(':'));A=matrix(code);h,c=refine(A,s,t);p=path(A,s,t)
   G.append({'code':code,'start':s,'end':t,'cells':c,'sign':sign(A,c),'path':p})
  rows.append(G)
 print('groups',len(rows),'signs',[sorted(set(x['sign'] for x in g)) for g in rows]);print(rows[1])
 (P/'collision_analysis.json').write_text(json.dumps(rows,indent=2)+'\n')
