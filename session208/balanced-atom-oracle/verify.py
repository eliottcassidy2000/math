import itertools,json
from pathlib import Path
R=Path(__file__).resolve().parent
n=13;code=241005792788926351623915
incoming=[0]*n;outgoing=[0]*n;reach=[1<<i for i in range(n)]
for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
 u,v=(i,j) if code>>k&1 else(j,i)
 incoming[v]|=1<<u;outgoing[u]|=1<<v;reach[u]|=1<<v
for k in range(n):
 for i in range(n):
  if reach[i]>>k&1:reach[i]|=reach[k]
assert all(x==(1<<n)-1 for x in reach)
N=1<<n;E=[[0]*n for _ in range(N)];B=[[0]*n for _ in range(N)]
for i in range(n):E[1<<i][i]=B[1<<i][i]=1
for S in range(1,N):
 if S&(S-1)==0:continue
 for v in range(n):
  if S>>v&1:
   P=S^(1<<v)
   E[S][v]=sum(E[P][u] for u in range(n) if P>>u&1 and incoming[v]>>u&1)
   B[S][v]=sum(B[P][u] for u in range(n) if P>>u&1 and outgoing[v]>>u&1)
assert sum(E[-1])==sum(B[-1])==759375==15**5
normalization=0
for u,v in itertools.combinations(range(n),2):
 uv=vu=0;rest=(N-1)^(1<<u)^(1<<v);s=rest
 while True:
  A=s|(1<<u);uv+=E[A][u]*B[(N-1)^A][v]
  A=s|(1<<v);vu+=E[A][v]*B[(N-1)^A][u]
  if s==0:break
  s=(s-1)&rest
 normalization+=uv if outgoing[u]>>v&1 else vu
assert normalization==(n-1)*759375
result={'n':n,'code':str(code),'strong':True,'count':759375,'power_of_15':5,'global_upper_rate':'13/5','end_counts':E[-1],'start_counts':B[-1],'arc_normalization':normalization,'status':'PASS','scope':'Explicit witness and global upper bound only.'}
(R/'RESULT.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
