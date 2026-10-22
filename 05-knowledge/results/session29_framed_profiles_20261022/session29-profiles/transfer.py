import itertools as it,json
from pathlib import Path
import seamed as s
STATES=['S','E','I','DE','DB']
def classify(paths,a):
 if len(paths)==1:
  p=paths[0]
  return 'S' if p[0]==a else 'E' if p[-1]==a else 'I'
 p,q=paths
 if p[-1]!=a:return None # canonical order: active component first
 return 'DB' if len(p)==1 else 'DE'
def countcovers(n,adj,a):
 ans=[0]*5
 for p in it.permutations(range(n)):
  bad=[i+1 for i in range(n-1) if not adj[p[i]][p[i+1]]]
  if not bad:ans[STATES.index(classify([p],a))]+=1
  if len(bad)>1:continue
  cuts=bad if bad else range(1,n)
  for cut in cuts:
   key=classify([p[:cut],p[cut:]],a)
   if key:ans[STATES.index(key)]+=1
 return ans

def transition(tile,state):
 q=1 if state in ('S','E','I') else 2;n=q+5
 start_a=state in('S','DB');end_a=state in('E','DE','DB')
 adj=[[False]*n for _ in range(n)]
 for i in range(n):
  for j in range(n):
   if i==j:continue
   if i>=q and j>=q:adj[i][j]=bool(s.arc(tile,5,i-q,j-q))
   elif i<q and j>=q:adj[i][j]=(i==0 and end_a and j==q)
   elif i>=q and j<q:adj[i][j]=not(i==q and j==0 and start_a)
 return countcovers(n,adj,n-1)
def matvec(v,M):return [sum(v[i]*M[i][j] for i in range(5)) for j in range(5)]
T={t:s.complement(t,5) for m in(40,41) for t in s.tilelist(m)}
initial={};mat={}
for t,x in T.items():
 adj=[[False if i==j else bool(s.arc(x,5,i,j)) for j in range(5)] for i in range(5)]
 initial[t]=countcovers(5,adj,4)
 mat[t]=[transition(x,state) for state in STATES]
for row in s.rows:
 t,u=row['tiles'];v=matvec(initial[t],mat[u]);s.require(sum(v[:3])==row['H'],'transfer vs324DP')
# independent three-block tests, including two-word types and all mirror statuses.
checks=[]
letters=list(T)
for inds in ((0,0,0),(0,5,0),(5,0,5),(17,17,17),(1,8,16),(4,12,7)):
 ts=[letters[i] for i in inds];v=initial[ts[0]]
 for t in ts[1:]:v=matvec(v,mat[t])
 actual=s.hdata(s.complement(s.join(ts),15),15)[0][-1]
 s.require(sum(v[:3])==actual,'three-block transfer vs32768subsetDP')
 checks.append({'tiles':ts,'H':actual})
V=[sum(v[j] for v in initial.values()) for j in range(5)]
M=[[sum(m[i][j] for m in mat.values()) for j in range(5)] for i in range(5)]
v=V;seq=[]
for r in range(1,13):seq.append(sum(v[:3]));v=matvec(v,M)
import sympy as sp
z=sp.symbols('z');F=sp.cancel(z*(sp.Matrix([V])*(sp.eye(5)-z*sp.Matrix(M)).inv()*sp.Matrix([1,1,1,0,0]))[0])
result={'states':STATES,'initial':initial,'matrices':mat,'sum_initial':V,'sum_matrix':M,'H_mass_coefficients_r1_to12':seq,'H_mass_GF':str(F),'three_block_checks':checks}
Path(__file__).with_name('transfer-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in('initial','matrices')},indent=2))
# Independent subset-DP validation of the path-cover transition enumerator.
def independent(adj,a):
 n=len(adj);full=(1<<n)-1
 def dpcount(adj):
  dp=[[0]*n for _ in range(1<<n)]
  for i in range(n):dp[1<<i][i]=1
  for sub in range(1,1<<n):
   for i in range(n):
    if not dp[sub][i]:continue
    for j in range(n):
     if not(sub>>j)&1 and adj[i][j]:dp[sub|1<<j][j]+=dp[sub][i]
  return dp
 e=dpcount(adj);start=dpcount(list(map(list,zip(*adj))));h=[sum(row) for row in e]
 S=start[full][a];E=e[full][a];Q=sum(e[sub][a]*h[full^sub] for sub in range(1,full) if(sub>>a)&1);DB=h[full^(1<<a)]
 return [S,E,h[full]-S-E,Q-DB,DB]
for t,x in T.items():
 for state in STATES:
  q=1 if state in('S','E','I') else 2;n=q+5
  start_a=state in('S','DB');end_a=state in('E','DE','DB')
  adj=[[False]*n for _ in range(n)]
  for i in range(n):
   for j in range(n):
    if i==j:continue
    if i>=q and j>=q:adj[i][j]=bool(s.arc(x,5,i-q,j-q))
    elif i<q and j>=q:adj[i][j]=(i==0 and end_a and j==q)
    elif i>=q and j<q:adj[i][j]=not(i==q and j==0 and start_a)
  s.require(independent(adj,n-1)==mat[t][STATES.index(state)],'independent subsetDP transition')
 s.require(mat[t][3]==mat[t][4],'DE DB exact lumping')
print('INDEPENDENT TRANSITION DP PASS all90rows; DE/DB lumping PASS')
