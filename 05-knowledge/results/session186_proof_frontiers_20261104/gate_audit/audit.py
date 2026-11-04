import itertools,json
from pathlib import Path
SEEDS=((16,5,33,39,46,3),(16,5,49,39,42,3),(18,4,49,39,42,3))
def adjacency(seed,r):
 n=r+5
 return [[False if i==j else (i<j if i>=5 and j>=5 else bool(seed[min(i,5)]>>min(j,5)&1)) for j in range(n)] for i in range(n)]
def reverse(a,internal=False):
 return [[not a[i][j] if i!=j and (not internal or i>=5 and j>=5) else a[i][j] for j in range(len(a))] for i in range(len(a))]
def perm_c(a):
 n=len(a);C=[[0]*n for _ in range(n)];H=0
 for p in itertools.permutations(range(n)):
  bad=[k for k in range(n-1) if not a[p[k]][p[k+1]]]
  H+=not bad
  if len(bad)>1:continue
  for k in (bad if bad else range(n-1)):C[p[k]][p[k+1]]+=1
 return H,C
records=[]
for j,seed in enumerate(SEEDS):
 a=adjacency(seed,3);want=2*(4*3+9*2-3+2*j)
 for label,b,sign in [('original',a,1),('whole',reverse(a),-1),('module',reverse(a,True),-1)]:
  h,c=perm_c(b);L=c[6][7]-c[5][7]+c[7][5]-c[7][6]
  if L!=sign*want or h!=4*27+9*8-3+2*j:raise RuntimeError((j,label,h,L))
  records.append(dict(state=j,case=label,H=h,L=L))
# Independent Hamilton recurrence uses subset-deletion recursion, adjacency booleans.
from functools import lru_cache
def hamilton(a):
 @lru_cache(None)
 def end(mask,v):
  rem=mask^(1<<v)
  return 1 if not rem else sum(end(rem,u) for u in range(len(a)) if rem>>u&1 and a[u][v])
 return sum(end((1<<len(a))-1,v) for v in range(len(a)))
q=json.loads((Path(__file__).resolve().parent.parent/'gate'/'queries.json').read_text())
for row in q['literal_two_port_witnesses']:
 a=adjacency(SEEDS[row['state']],row['r']);n=len(a)
 for key,base in [('F',a),('F_op',reverse(a))]:
  out=[r+[False,False] for r in base]+[[False]*(n+2),[False]*(n+2)];out[n][n+1]=True
  for u in range(n):
   for p,m in [(n,row['a']),(n+1,row['c'])]:out[u][p]=bool(m>>u&1);out[p][u]=not out[u][p]
  if hamilton(out)!=row[key]:raise RuntimeError(('query',row,key))
# Ordered path-cover profile, computed by permutations and all cut sets.
profiles=[]
for rows in [(4,1,2,7,15),(10,0,3,6,15)]:
 counts=[0]*5
 for p in itertools.permutations(range(5)):
  bad=sum(1<<i for i in range(4) if not(rows[p[i]]>>p[i+1]&1))
  for cuts in range(16):
   if cuts&bad==bad:counts[cuts.bit_count()]+=1
 profiles.append(counts)
if profiles!=[[3,42,162,240,120]]*2:raise RuntimeError(('hostile profile',profiles))
pairs=[];values=[]
for rows in [(4,1,2,7,15),(10,0,3,6,15)]:
 pairs.append([(i,j) for i,j in itertools.combinations(range(5),2) if all(bool(rows[i]>>k&1)==bool(rows[j]>>k&1) for k in range(5) if k not in (i,j))])
 a=[[bool(rows[i]>>j&1) for j in range(5)] for i in range(5)]
 out=[r+[False,False] for r in a]+[[False]*7,[False]*7];out[5][6]=True
 for i in range(5):out[5][i]=True;out[i][6]=(i==0);out[6][i]=(i!=0)
 values.append(hamilton(out))
if pairs!=[[(3,4)],[]] or values!=[15,9]:raise RuntimeError(('hostile',pairs,values))
ambient=[]
for j,seed in enumerate(SEEDS):
 basecounts=[];portcounts=[]
 for rows in [(4,1,2,7,15),(10,0,3,6,15)]:
  module=[[bool(rows[i]>>k&1) for k in range(5)] for i in range(5)]
  extended=[r+[False,False] for r in module]+[[False]*7,[False]*7];extended[5][6]=True
  for i in range(5):extended[5][i]=True;extended[i][6]=(i==0);extended[6][i]=(i!=0)
  for m,dest in [(module,basecounts),(extended,portcounts)]:
   n=len(m);a=[[False]*(n+5) for _ in range(n+5)]
   for u in range(n+5):
    for v in range(n+5):
     if u<n and v<n:a[u][v]=m[u][v]
     elif u>=n and v>=n:a[u][v]=bool(seed[u-n]>>(v-n)&1)
     elif u<n:a[u][v]=bool(seed[5]>>(v-n)&1)
     else:a[u][v]=bool(seed[u-n]>>5&1)
   # Independent transitive closure checks the strong lift.
   reach=[r[:] for r in a]
   for k in range(n+5):
    for u in range(n+5):
     for v in range(n+5):reach[u][v]|=reach[u][k] and reach[k][v]
   if not all(all(r) for r in reach):raise RuntimeError('strong')
   dest.append(hamilton(a))
 if basecounts!=[1611+6*j]*2 or portcounts!=[21015+30*j,17793+18*j]:raise RuntimeError(('ambient',j,basecounts,portcounts))
 ambient.append(dict(state=j,H=basecounts,port_values=portcounts))
print(json.dumps(dict(status='PASS',permutation_controls=records,two_port_queries=24,hostile_profiles=profiles,hostile_pairs=pairs,hostile_queries=values,strong_contexts=ambient),indent=2))
