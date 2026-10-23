from explore import hamilton,scc,rootpair
from itertools import product,combinations,permutations
from math import gcd
import json
from pathlib import Path

def require(q,msg):
 if not q:raise RuntimeError(msg)
def near_transitive(a):
 n=a+2;b=a;c=a+1;M=[[0]*n for _ in range(n)]
 for i in range(a):
  for j in range(i+1,a):M[i][j]=1;M[j][i]=-1
  M[i][b]=1;M[b][i]=-1;M[c][i]=1;M[i][c]=-1
 M[b][c]=1;M[c][b]=-1
 return M
results={'near_transitive':[],'gcd_checks':0}
for a in range(1,13):
 h=hamilton(near_transitive(a));require(h==2**a+1,('near_transitive',a,h));results['near_transitive'].append([a,h])
ps=[5,7,11,13,17,19,23,29,31,37]
for p in ps:
 require((2**p+1)%3==0 and (2**p+1)%9!=0,('v3',p))
for p,q in combinations(ps,2):
 require(gcd(2**p+1,2**q+1)==3,('gcd',p,q));require(gcd((2**p+1)//3,(2**q+1)//3)==1,('coprime',p,q));results['gcd_checks']+=1
seen=set()
for exps in product(range(4),repeat=5):
 h=1
 for p,e in zip(ps,exps):h*=(2**p+1)**e
 require(h not in seen,('injectivity',exps));seen.add(h)
results['distinct_H_products']=len(seen)
def compositions(N):
 if N==0:yield [];return
 for a in range(2,N+1):
  for tail in compositions(N-a):yield [a]+tail
def interval(parts):
 N=sum(parts);M=[[0 if i==j else (1 if i<j else -1) for j in range(N+1)] for i in range(N+1)]
 l=0
 for a in parts:
  r=l+a;M[l][r]=-1;M[r][l]=1;l=r
 return M
count=0
for N in range(2,13):
 for parts in compositions(N):
  M=interval(parts);want=1
  for a in parts:want*=2**(a-1)+1
  require(hamilton(M)==want,('interval H',parts))
  require(len(scc(M))==1,('interval strong',parts))
  scores=[sum(x==1 for x in row) for row in M]
  require(scores==[N-1]+[N-t for t in range(1,N)]+[1],('interval scores',parts,scores))
  count+=1
results['interval_compositions_checked']=count
hostiles=[];pair_reference=None
for parts in [[2,2,4],[2,3,3]]:
 M=interval(parts);n=len(M);pairs=rootpair(M)
 if pair_reference is None:pair_reference=pairs
 require(pairs==pair_reference,('explicit genus2 hostile',parts))
 c3=sum(all(M[p[i]][p[(i+1)%3]]==1 for i in range(3)) for p in permutations(range(n),3))//3
 c4=sum(all(M[p[i]][p[(i+1)%4]]==1 for i in range(4)) for p in permutations(range(n),4))//4
 hostiles.append({'parts':parts,'n':n,'H':hamilton(M),'c3':c3,'c4':c4,'scores':[sum(x==1 for x in row) for row in M],'outmasks':[sum(1<<j for j in range(n) if M[i][j]==1) for i in range(n)]})
require([x['H'] for x in hostiles]==[81,75],'explicit H values')
require(all(x['c3']==7 and x['c4']==10 for x in hostiles),'explicit cycle values')
results['strong_fixed_score_genus2_hostile']=hostiles
results['status']='PASS' 
print(json.dumps(results,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(results,indent=2)+'\n')
