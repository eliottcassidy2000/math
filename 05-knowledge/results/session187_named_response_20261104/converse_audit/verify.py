import json,itertools,sys
from pathlib import Path
BASE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'collision'))
from verify import polynomial

def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def count(rows):
 n=len(rows);pred=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)];dp=[[0]*n for _ in range(1<<n)]
 for s in range(1,1<<n):
  for v in range(n):
   if not s>>v&1:continue
   rem=s^(1<<v)
   if not rem:dp[s][v]=1;continue
   p=pred[v]&rem;value=0
   while p:
    b=p&-p;p-=b;value+=dp[rem][b.bit_length()-1]
   dp[s][v]=value
 return sum(dp[-1])
def extension(rows,a,c):
 n=len(rows);mat=[[False]*(n+2) for _ in range(n+2)]
 for u in range(n):
  for v in range(n):mat[u][v]=bool(rows[u]>>v&1)
  for v,m in [(n,a),(n+1,c)]:mat[u][v]=bool(m>>u&1);mat[v][u]=not mat[u][v]
 mat[n][n+1]=True
 return [sum(1<<v for v in range(n+2) if mat[u][v]) for u in range(n+2)]
def reachable(rows,u):
 seen={u};todo=[u]
 while todo:
  v=todo.pop()
  for w in range(len(rows)):
   if rows[v]>>w&1 and w not in seen:seen.add(w);todo.append(w)
 return seen
r=json.loads((BASE/'converse/receipt.json').read_text());checks=0
for rec in r['records']:
 q=rec['quotient_rows'];m=rec['core_rows'];a=rec['ambient_rows'];b=rec['active_rows'];n=len(m)
 reach=[reachable(q,u) for u in range(len(q))];K=[u for u in range(len(q)) if u in reach[0] and 0 in reach[u]]
 keep=list(range(n))+[n+u-1 for u in K if u];incoming=sum(1<<(n+u-1) for u in range(1,len(q)) if u not in K and q[u]&1)
 need(len(keep)==len(b),'active size');factor=rec['factor']
 for x,y in [(0,0),((1<<len(b))-1,0),(0,(1<<len(b))-1),(1,2),(3,3)]:
  lift=lambda mask:incoming|sum(1<<keep[v] for v in range(len(b)) if mask>>v&1)
  need(count(extension(a,lift(x),lift(y)))==factor*count(extension(b,x,y)),'independent SCC query');checks+=1
h=json.loads((BASE/'converse/hostile.json').read_text())
for rec in h['sixteen_literal_queries']:
 v=count(extension(h['ambient_rows'],rec['ambient_a'],rec['ambient_c']))-count(extension(h['ambient_reversed_module_rows'],rec['ambient_a'],rec['ambient_c']))
 need(v==rec['difference'],'hostile independent query');checks+=1
b=h['active_rows'];br=h['active_reversed_rows'];n=len(b);bo=[((1<<n)-1)^(1<<u)^b[u] for u in range(n)]
p=polynomial(b)
for rows,key in [(br,'complete_active_module_reversal_polynomial'),(bo,'complete_active_global_converse_polynomial')]:
 p2=polynomial(rows);want={x|(y<<n):v for x,y,v in h[key]};got={s:p.get(s,0)-p2.get(s,0) for s in p.keys()|p2.keys() if p.get(s,0)!=p2.get(s,0)}
 need(want==got,'independent complete word polynomial')
need(all(len(reachable(b,u))==n for u in range(n)),'active strong')
print(json.dumps({'status':'PASS','independent_query_checks':checks,'complete_word_polynomials_checked':2,'active_strong':True},indent=2,sort_keys=True))
