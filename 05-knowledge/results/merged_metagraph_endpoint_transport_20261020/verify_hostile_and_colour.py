import json,collections as co
from fractions import Fraction as F
from pathlib import Path
from probe import tiles

def matrix_sum(A,B):return [[a+b for a,b in zip(x,y)] for x,y in zip(A,B)]
def leak(A):
 d=list(map(sum,A));q=len(A)
 return [[F(d[i]*(i==j))-sum(F(A[i][k]*A[k][j],d[k]) for k in range(q) if d[k]) for j in range(q)] for i in range(q)]
def scorekey(s):return min(tuple(sorted(s)),tuple(sorted(len(s)-1-a for a in s)))
reports=[]
for n in range(3,8):
 d=json.load(open(f'n{n}.json'));q=d['q'];B=d['Ablue'];K=d['Ablack'];A=matrix_sum(B,K)
 # Independent HP-orbit classifier cross-check from a separate implementation, under explicit bit gauge.
 p=Path(f'crosscheck_census_n{n}.json')
 cross=None
 if p.exists():
  other=json.loads(p.read_text());theirtiles=[(x-1,y-1) for x,y in other['tiles']];ourtiles=tiles(n);newbits=[ourtiles.index(t) for t in theirtiles]
  labelmap={}
  for t,u in enumerate(d['pi']):
   tt=sum((1-((t>>k)&1))<<j for j,k in enumerate(newbits));v=other['pi'][tt]
   assert labelmap.setdefault(u,v)==v
  assert len(set(labelmap.values()))==q
  cross={'independent_classifiers_agree':True,'node_count':q}
 if n<=6:
  L=leak(A);LB=leak(B);LK=leak(K);b=list(map(sum,B));k=list(map(sum,K));dd=list(map(sum,A));mixed=[u for u in range(q) if b[u]*k[u]]
  for i in range(q):
   for j in range(q):
    rhs=LB[i][j]+LK[i][j]+sum(F(b[u]*k[u],dd[u])*(F(B[u][i],b[u])-F(K[u][i],k[u]))*(F(B[u][j],b[u])-F(K[u][j],k[u])) for u in mixed)
    assert L[i][j]==rhs
 reports.append({'n':n,'colour_matrix_identity':n<=6,'cross_classifier':cross})

x=json.load(open('n5.json'));pi=x['pi'];N=64;q=x['q'];h=x['c3'];C=[t^63 for t in range(N)];C[0x29]=0x1a;C[0x1a]=0x29;C[0x16]=0x25;C[0x25]=0x16
assert all(C[C[t]]==t and C[t]!=t for t in range(N))
assert all(x['reflection'][C[t]]==C[x['reflection'][t]] for t in range(N))
Bnew=[[0]*q for _ in range(q)];Knew=[[0]*q for _ in range(q)]
for t in range(N):(Bnew if x['reflection'][t]==t else Knew)[pi[t]][pi[C[t]]]+=1
assert Knew==x['Ablack']
assert list(map(sum,Bnew))==list(map(sum,x['Ablue']))
assert [Bnew[i][i] for i in range(q)]==[x['Ablue'][i][i] for i in range(q)]
def joint(B,key):
 c=co.Counter()
 for u in range(q):
  for v in range(q):c[key(u),key(v)]+=B[u][v]
 return +c
assert joint(Bnew,lambda i:h[i])==joint(x['Ablue'],lambda i:h[i])
assert joint(Bnew,lambda i:scorekey(x['scores'][i]))==joint(x['Ablue'],lambda i:scorekey(x['scores'][i]))
assert all(not Bnew[u][v] or x['Ablue'][u][v] for u in range(q) for v in range(q))
changes=[]
for u in range(q):
 before=co.Counter();after=co.Counter()
 for v in range(q):before[h[v]-h[u]]+=x['Ablue'][u][v];after[h[v]-h[u]]+=Bnew[u][v]
 if +before!=+after:changes.append({'node':u,'before':dict(+before),'after':dict(+after)})
assert [c['node'] for c in changes]==[2,3]
out={'controls':reports,'hostile':{'involution':True,'commutes_with_reflection':True,'colour_masses_preserved':True,'diagonals_preserved':True,'global_joint_triangle_counts_preserved':True,'global_joint_score_counts_preserved':True,'only_previously_supported_pairs_used':True,'support_itself_preserved':False,'violating_rows':changes},'failures':0}
Path('hostile_and_colour_audit.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
