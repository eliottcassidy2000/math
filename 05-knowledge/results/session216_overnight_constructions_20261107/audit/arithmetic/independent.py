from pathlib import Path
import json,itertools
P=Path(__file__).parent;W=json.loads((P/'frozen/interval_witnesses.json').read_text())
def add(a,b):
 c=[0]*max(len(a),len(b))
 for i,x in enumerate(a):c[i]+=x
 for i,x in enumerate(b):c[i]+=x
 return c
def mul(a,b):
 c=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c

def coeff(a,v):
 # Independently contract a visit into each external transition, with optional
 # prefix/suffix visits. No enumeration of external permutations or gap subsets.
 ex=[u for u in range(len(a)) if u!=v];n=len(ex);N=1<<n
 dp=[[[] for j in ex] for S in range(N)]
 for j,u in enumerate(ex):dp[1<<j][j]=[1,int(bool(a[v]>>u&1))]
 for S in range(1,N):
  for j,u in enumerate(ex):
   if not dp[S][j]:continue
   for l,w in enumerate(ex):
    if not S>>l&1:
     t=[int(bool(a[u]>>w&1)),int(bool(a[u]>>v&1) and bool(a[v]>>w&1))]
     dp[S|1<<l][l]=add(dp[S|1<<l][l],mul(dp[S][j],t))
 out=[]
 for j,u in enumerate(ex):out=add(out,mul(dp[-1][j],[1,int(bool(a[u]>>v&1))]))
 return (out+[0]*(len(a)+1))[:len(a)+1]

def counts(a):
 n=len(a);out=[]
 for S in range(1<<n):
  vertices=[j for j in range(n) if S>>j&1]
  out.append(sum(all(a[u]>>v&1 for u,v in zip(p,p[1:])) for p in itertools.permutations(vertices)))
 return out

def tournaments(n):
 edges=list(itertools.combinations(range(n),2))
 for bits in itertools.product((0,1),repeat=len(edges)):
  a=[0]*n
  for (i,j),bit in zip(edges,bits):
   if bit:a[i]|=1<<j
   else:a[j]|=1<<i
  yield a

def substitute(a,v,U):
 # Place module FIRST to avoid copying producer's vertex convention.
 m=len(U);ex=[j for j in range(len(a)) if j!=v];out=list(U)+[0]*len(ex)
 for i,u in enumerate(ex):
  for j,w in enumerate(ex):
   if a[u]>>w&1:out[m+i]|=1<<(m+j)
  for j in range(m):
   if a[u]>>v&1:out[m+i]|=1<<j
   else:out[j]|=1<<(m+i)
 return out
seeds=[a for n in range(1,5) for a in tournaments(n)]+[[sum(1<<j for j in range(i+1,n)) for i in range(n)] for n in range(1,7)]
records=[];graphs=[];empty=[]
for w in W:
 a=w['rows'];v=w['root'];n=len(a);cy=w['hamilton_cycle']
 if not (n==w['n']<=8 and sorted(cy)==list(range(n))):raise RuntimeError('dimensions')
 if any(x<0 or x>>n or x>>i&1 for i,x in enumerate(a)):raise RuntimeError('row domain')
 if any((a[i]>>j&1)+(a[j]>>i&1)!=1 for i in range(n) for j in range(i)):raise RuntimeError('tournament')
 if not all(a[u]>>z&1 for u,z in zip(cy,cy[1:]+cy[:1])):raise RuntimeError('cycle')
 c=coeff(a,v)
 if not(c[1]==w['h'] and c[2]==36 and not any(c[3:])):raise RuntimeError(c)
 records.append({'h':w['h'],'coefficients':c})
 empty.append({'d':w['h']-72,'actual_empty_output':c[0],'invalid_extension':36+w['h']-72,'coincides':c[0]==w['h']-36})
 for U in seeds:
  hs=counts(U);q=sum(x*hs[(len(hs)-1)^s] for s,x in enumerate(hs));target=36*q+(w['h']-72)*hs[-1]
  graph=substitute(a,v,U);graphs.append((target,graph))
if [w['h']-72 for w in W]!=list(range(-7,14,2)):raise RuntimeError('offset set')
if records!=json.loads((P/'frozen/VERIFY.json').read_text())['coefficient_certificates']:raise RuntimeError('receipt')
with (P/'GRAPHS.txt').open('w') as f:
 for h,a in graphs:print(h,len(a),*a,file=f)
print(json.dumps({'coefficient_certificates':records,'literal_graphs_prepared':len(graphs),'max_graph_order':max(len(a) for h,a in graphs),'seed_count':75,'empty_boundary':empty},indent=2))
