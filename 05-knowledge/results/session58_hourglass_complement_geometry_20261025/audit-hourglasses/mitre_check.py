from itertools import combinations,permutations,product
from collections import Counter
import json

def require(x,msg):
 if not x:raise RuntimeError(msg)
P=range(9)
# Independent coordinate arithmetic with integer-coded points, no producer imports.
def third(a,b):return (-(a//3+b//3))%3*3+(-(a%3+b%3))%3
L={frozenset((a,b,third(a,b))) for a,b in combinations(P,2)}
require(len(L)==12 and all(len(t)==3 for t in L),'lines')
mitres={}
for blocks in combinations(L,5):
 support=set().union(*blocks)
 if len(support)!=7:continue
 degree=Counter(v for line in blocks for v in line)
 if sorted(degree.values())!=[2]*6+[3]:continue
 c=next(v for v in degree if degree[v]==3)
 outer=[line for line in blocks if c not in line]
 if len(outer)!=2 or outer[0]&outer[1]:continue
 mitres[frozenset(blocks)]=c
require(len(mitres)==36,'all five-line subsets mitres')
patterns=Counter(frozenset(a|b) for a,b in combinations(mitres,2) if mitres[a]==mitres[b])
require(len(patterns)==54 and set(patterns.values())=={1},'same-center pairs')
# Exhaustive 9! labelled arrays; test all eight original line sums directly.
idx=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
layouts=Counter()
for perm in permutations(P):
 if all(third(perm[i],perm[j])==perm[k] for i,j,k in idx):
  layouts[frozenset(frozenset(perm[i] for i in t) for t in idx)]+=1
require(sum(layouts.values())==432 and set(layouts.values())=={8},'permutation layout count')
require(set(layouts)==set(patterns),'pattern equality')
# Symbolic coefficient triples for E,x,z; verify ordered omitted line defects.
vectors=[(1,1,0),(1,-1,-1),(1,0,1),(1,-1,1),(1,0,0),(1,1,-1),(1,0,-1),(1,1,1),(1,-1,0)]
omitted=['AFH','BDI','BFG','CDH'];defects=[]
for line in omitted:
 c=tuple(sum(vectors[ord(t)-65][j] for t in line)-(3 if j==0 else 0) for j in range(3))
 defects.append(c)
require(defects==[(0,3,0),(0,-3,0),(0,0,-3),(0,0,3)],'ordered exact defects')
# Normalize E=1 for nonzero square centers: scaling by inverse center preserves residues/distinctness.
first=[]
for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
 Q={a*a%p for a in range(p)};bad=[]
 for A,C in product(Q,repeat=2):
  h=[1,A,C,(3-A-C)%p,(2-A)%p,(2-C)%p,(A+C-1)%p]
  D=(1+C-A)%p;F=(1+A-C)%p
  if len(set(h))==7 and set(h)<=Q and D not in Q and F not in Q:bad.append((A,C))
 first.append((p,len(bad)))
require(all(n==0 for p,n in first if p<37) and first[-1][1]>0,'normalized prime minimality')
try:require(False,'intentional')
except RuntimeError:pass
else:raise RuntimeError('disabled guards')
print(json.dumps({'status':'PASS','mitres':len(mitres),'patterns':len(patterns),'layouts':sum(layouts.values()),'ordered_defects':defects,'normalized_hostile_counts':first},sort_keys=True))
