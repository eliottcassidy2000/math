def require(x):
 if not x: raise AssertionError("recursion check failed")
from itertools import permutations, combinations
n=5; mask=8
pairs=list(combinations(range(n),2)); bits={p:mask>>k&1 for k,p in enumerate(pairs)}
def beat(a,b):return bool(bits[a,b]) if a<b else not bits[b,a]
def ispath(p):return all(beat(a,b) for a,b in zip(p,p[1:]))
H=[]; M=[[0]*n for _ in range(n)]
for p in permutations(range(n)):
 if ispath(p):H.append(p)
 for k in range(1,n):
  if ispath(p[:k]) and ispath(p[k:]):M[p[k-1]][p[k]]+=1
start=[sum(p[0]==i for p in H) for i in range(n)]; end=[sum(p[-1]==i for p in H) for i in range(n)]
print('H',len(H),'start',start,'end',end,'M',M)
print('W',[[M[i][j]+M[j][i] for j in range(n)] for i in range(n)])
vals={}
for m in range(1<<n):
 S={i for i in range(n) if m>>i&1}
 formula=len(H)+sum(end[i]-start[i] for i in S)+sum(M[i][j] for i in S for j in range(n) if j not in S)
 direct=sum(ispath(p[:p.index(n)]) and ispath(p[p.index(n)+1:]) and (p.index(n)==0 or p[p.index(n)-1] in S) and (p.index(n)==n or p[p.index(n)+1] not in S) for p in permutations(range(n+1)))
 if formula!=direct: raise AssertionError((S,formula,direct))
 vals[m]=formula
print('ear values by mask',list(vals.values()))
print('proper nonempty image',sorted(set(vals[m] for m in range(1,31))))
print('pair complement sums',[(m,vals[m]+vals[31^m]) for m in range(16)])
w=[[M[i][j]+M[j][i] for j in range(n)] for i in range(n)]
a=[(end[i]-start[i]+sum(M[i]))//2 for i in range(n)]
b=[[w[i][j]//2 for j in range(n)] for i in range(n)]
d=[2*a[i]-sum(b[i]) for i in range(n)]
print('alpha',a,'beta',b,'drift',d)
require(sum(d)==0)
for m in range(32):
 S={i for i in range(n) if m>>i&1}
 require(vals[m]==len(H)+sum(d[i] for i in S)+sum(b[i][j] for i in S for j in range(n) if j not in S))
autos=[p for p in permutations(range(n)) if all(beat(p[i],p[j])!=beat(i,j) for i,j in combinations(range(n),2))]
print('anti-autos',autos)
for p in autos:
 require(all(vals[m]==vals[sum(1<<p[i] for i in range(n) if not m>>i&1)] for m in range(32)))
from collections import Counter
print('value multiplicities',sorted(Counter(vals.values()).items()))
print('by cardinality',[(k,sorted(Counter(vals[m] for m in range(32) if m.bit_count()==k).items())) for k in range(6)])
