"""Independent literal two-apex DP and all-order GF coefficient controls."""
import json,itertools,math,collections,contextlib,io
from pathlib import Path
with contextlib.redirect_stdout(io.StringIO()):
 from independent import subset_count
P=Path(__file__).resolve().parents[1]/'wildcard'
D=json.loads((P/'COUNT_DFA21.json').read_text());F=json.loads((P/'FIBRES63.json').read_text())['fibres']
def graph(w):
 n=len(w)+2;u=n-2;v=n-1;G=[[False]*n for _ in range(n)];G[u][v]=True
 for i,x in enumerate(w):
  for j in range(i+1,n-2):G[i][j]=True
  for z,b in ((u,x&1),(v,x>>1)):
   G[z][i]=bool(b);G[i][z]=not b
 return G

def expected(h,m):
 R=F[str(h)]
 return sum(int(c)*(m-int(j)+1) for j,c in R['two_apex_strong'].items() if m>=int(j))+sum(int(c)*math.comb(m-int(j)+2,2) for j,c in R['separate_apex_numerator'].items() if m>=int(j))
trials=0;gfs=0
for m in range(6):
 counts=collections.Counter()
 for w in itertools.product(range(4),repeat=m):
  h=subset_count(graph(w));counts[h]+=1;s=D['initial']
  for c in w:s=D['transitions'][s][c]
  assert D['outputs'][s]==(h if h<=21 else -1),(w,h,D['outputs'][s]);trials+=1
 for h in range(1,64,2):
  assert expected(h,m)==counts[h],(h,m,expected(h,m),counts[h]);gfs+=1
print(json.dumps({'status':'PASS','independent_literal_words':trials,'full_GF_coefficient_checks':gfs,'lengths':'0..5','targets':'all odd1..63'},indent=2))
