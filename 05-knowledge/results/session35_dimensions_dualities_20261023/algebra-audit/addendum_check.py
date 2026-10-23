from itertools import permutations, combinations
from collections import defaultdict
import json
import sympy as s

def require(x,label):
 if not x: raise RuntimeError(label)
def hp(out,subset=None):
 vertices=range(len(out)) if subset is None else subset
 return sum(all(out[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(vertices))
def matrix(out):
 n=len(out)
 return [[0 if i==j else 1 if out[i]>>j&1 else -1 for j in range(n)] for i in range(n)]
def signature(out):
 M=matrix(out);n=len(M);P=[[int(i==j) for j in range(n)] for i in range(n)];key=[]
 for k in range(n):
  P=[[sum(P[i][h]*M[h][j] for h in range(n)) for j in range(n)] for i in range(n)]
  key += [sum(P[i][i] for i in range(n)),sum(map(sum,P))]
 return tuple(key)
def make(n,mask,conditioned):
 edges=[(i,j) for i in range(n) for j in range(i+1,n) if not(conditioned and j==i+1)]
 out=[0]*n
 if conditioned:
  for i in range(n-1):out[i]|=1<<(i+1)
 for k,(i,j) in enumerate(edges):
  a,b=(i,j) if mask>>k&1 else (j,i);out[a]|=1<<b
 return out
summaries=[]
for n in range(1,7):
 conditioned=n==6
 total=1<<((n-1)*(n-2)//2 if conditioned else n*(n-1)//2)
 groups=defaultdict(set)
 for mask in range(total):
  out=make(n,mask,conditioned);groups[signature(out)].add(hp(out))
 row={'n':n,'range':'Hamiltonian-order presentations' if conditioned else 'all labelled tournaments','count':total,'classes':len(groups),'H_ambiguous_classes':sum(len(v)>1 for v in groups.values())}
 summaries.append(row)
 require(row['H_ambiguous_classes']==(3 if n==6 else 0),('census ambiguity',row))
 require(row['classes']==[1,1,2,3,9,28][n-1],('census classes',row))
rows=[[58,4,9,18,38,14],[22,12,8,17,38,15]]
z=s.symbols('z'); expected_p=z**6+15*z**4+47*z**2+1;expected_q=2*z*(3*z**4+34*z**2+59)
witness=[]
for out,H in zip(rows,[29,25]):
 require(hp(out)==H,'direct permutation HP')
 M=s.Matrix(matrix(out)); require(M==-M.T,'valid skew tournament')
 # Numeric determinants at seven points determine a degree<=6 polynomial.
 samples=[(k,(k*s.eye(6)-M).det()) for k in range(7)]
 require(s.expand(s.interpolate(samples,z)-expected_p)==0,'interpolated p')
 # Numerator deg<=5: six positive z values, no poles for real skew M.
 qs=[]
 for k in range(1,7):
  B=k*s.eye(6)-M
  response=sum(B.LUsolve(s.ones(6,1)))
  qs.append((k,s.cancel(response*expected_p.subs(z,k))))
 require(s.expand(s.interpolate(qs,z)-expected_q)==0,'interpolated q')
 cycles={}
 for k in [3,5]:
  cycles[k]=[]
  for sub in combinations(range(6),k):
   for rest in permutations(sub[1:]):
    c=(sub[0],)+rest
    if all(out[c[i]]>>c[(i+1)%k]&1 for i in range(k)):cycles[k].append(frozenset(c))
 pairs=sum(not(a&b) for a,b in combinations(cycles[3],2))
 require(1+2*(len(cycles[3])+len(cycles[5]))+4*pairs==H,'odd-cycle expansion control')
 witness.append({'H':H,'c3':len(cycles[3]),'c5':len(cycles[5]),'disjoint_triangle_pairs':pairs})
diffs={k:0 for k in range(1,7)}
for k in range(1,7):
 for sub in combinations(range(6),k):
  diffs[k]+=hp(rows[0],sub)!=hp(rows[1],sub)
require(list(diffs.values())==[0,0,8,9,3,1],'all subset count differences')
print(json.dumps({'census':summaries,'witness':witness,'subset_coefficient_differences':diffs,'p':str(expected_p),'q':str(expected_q)},indent=2))
