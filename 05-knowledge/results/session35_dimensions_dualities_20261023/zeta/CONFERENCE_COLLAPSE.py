from itertools import permutations,combinations
from collections import Counter,defaultdict
from itertools import product
import sympy as s
N=8
C=[[0]*N for _ in range(N)]
for j in range(1,N): C[0][j]=1;C[j][0]=-1
for i in range(1,N):
 for j in range(1,N):
  if i!=j:C[i][j]=1 if (j-i)%7 in {1,2,4} else -1
assert s.Matrix(C)**2==-7*s.eye(N)
def hp(A):
 dp=[[0]*N for _ in range(1<<N)]
 for j in range(N):dp[1<<j][j]=1
 for mask in range(1,1<<N):
  for j in range(N):
   if mask>>j&1:
    dp[mask][j]+=sum(dp[mask^(1<<j)][i] for i in range(N) if mask>>i&1 and A[i][j])
 return sum(dp[-1])
def odd_cycles(A):
 cyc={}
 for r in [3,5,7]:
  for vs in combinations(range(N),r):
   cyc[sum(1<<v for v in vs)]=sum(all(A[p[i]][p[(i+1)%r]] for i in range(r)) for tail in permutations(vs[1:]) for p in [(vs[0],)+tail])
 c={r:sum(v for k,v in cyc.items() if k.bit_count()==r) for r in [3,5,7]}
 c33=sum(v*w for k,v in cyc.items() for l,w in cyc.items() if k<l and k.bit_count()==l.bit_count()==3 and not k&l)
 c35=sum(v*w for k,v in cyc.items() for l,w in cyc.items() if k.bit_count()==3 and l.bit_count()==5 and not k&l)
 return (c[3],c[5],c[7],c33,c35)
def maps(A,B):
 da=list(map(sum,A));db=list(map(sum,B))
 if sorted(da)!=sorted(db):return
 groups=[([i for i in range(N) if da[i]==d],[j for j in range(N) if db[j]==d]) for d in sorted(set(da))]
 for choices in product(*(permutations(v) for _,v in groups)):
  p=[None]*N
  for (us,_),vs in zip(groups,choices):
   for u,v in zip(us,vs):p[u]=v
  if all(A[i][j]==B[p[i]][p[j]] for i in range(N) for j in range(N)):yield p
def iso(A,B):return next(maps(A,B),None) is not None
classes=[]; buckets=defaultdict(list)
for mask in range(128):
 d=[1]+[(-1 if mask>>(i-1)&1 else 1) for i in range(1,N)]
 A=[[int(d[i]*C[i][j]*d[j]>0) for j in range(N)] for i in range(N)]
 G=A
 scores=tuple(sorted(map(sum,A)))
 cl=next((i for i,(g,sc,a,m) in enumerate(classes) if scores==sc and iso(g,G)),None)
 if cl is None:cl=len(classes);classes.append((G,scores,A,mask))
 buckets[cl].append(mask)
for i,(G,sc,A,mask) in enumerate(classes):
 H=hp(A);cs=odd_cycles(A);assert H==1+2*sum(cs[:3])+4*sum(cs[3:])
 aut=sum(1 for _ in maps(G,G))
 print('class',i,'representative mask',mask,'orbit size',len(buckets[i]),'scores',sc,'H',H,'odd(c3,c5,c7,c33,c35)',cs,'Aut',aut,flush=True)
 print('masks',buckets[i],flush=True)
 print('reverse class',next(j for j,(g,*_) in enumerate(classes) if iso([list(x) for x in zip(*G)],g)),flush=True)
assert len(classes)==4
print('H histogram',dict(Counter({h:sum(len(buckets[i]) for i,(_,_,A,_) in enumerate(classes) if hp(A)==h) for h in [189,333]})))
# Universal-vertex switching uniqueness: for each v and desired sign, solve d_v*d_j*C_vj=sign with d_0=1.
uni=[]
for v in range(N):
 for sign in [-1,1]:
  if v==0:d=[1]+[sign*C[0][j] for j in range(1,N)]
  else:
   dv=sign*C[v][0];d=[1]+[(dv if j==v else sign*dv*C[v][j]) for j in range(1,N)]
  mask=sum((d[j]<0)<<(j-1) for j in range(1,N));uni.append(mask)
assert len(set(uni))==16
assert set(uni)=={m for i,(_,_,A,_) in enumerate(classes) if hp(A)==189 for m in buckets[i]}
print('Exactly16 universal-vertex switches verified; H189 iff universal vertex.')
source=set(buckets[0]);sink=set(buckets[3]);sink_complement={127^m for m in sink}
for code in [source,sink_complement]:
 assert len(code)==8 and 0 in code and all(x^y in code for x in code for y in code)
 assert Counter(x.bit_count() for x in code)=={0:1,4:7}
 lines=[127^x for x in code if x]
 assert all(sum((line>>i&1) and (line>>j&1) for line in lines)==1 for i,j in combinations(range(7),2))
 print('Simplex [7,3,4] code',sorted(code),'Fano line masks',sorted(lines))
assert source!=sink_complement
