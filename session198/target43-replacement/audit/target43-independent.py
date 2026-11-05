"""No imports from the packet: raw-coordinate and literal-word audit."""
from itertools import product,permutations
from functools import lru_cache
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
R=json.loads((HERE/'target43-replay/VALIDATION_NORMAL.json').read_text())
def move(s,g):
 D,E,U,d=s
 return {'O':(D,0,2*U,d),'Z':(D-E+U,d,2*U-d,d),'J':(D-E+U,D-E,D+U,D)}[g]
def word(s,w):
 for c in w:s=move(s,c)
 return s
@lru_cache(None)
def constants(L):
 if L==1:return 1,1,1
 if L%2==0:
  g,q,b=constants(L//2);return g,q,2*b
 a=(L+1)//2;g,q,b=constants(a);return g+a-a%2,q+b,2*b

def inventory(cap):
 bank={(1,1,1,1)};queue=list(bank)
 for s in queue:
  D,E,U,d=s
  candidates=[(D,0,2*U,d)]
  for L in range(1,cap//U+1):
   g,q,b=constants(L);candidates.append((g*D-(L%2)*E+L*U,D-(L==1)*E,q*D+b*U,D))
  for t in candidates:
   if t[2]<=cap and t not in bank:bank.add(t);queue.append(t)
 return bank
JM=[((0,0,0,3),(0,0,0,-3),(0,0,12,-3),(3,0,-12,3)),((0,0,0,3),(0,0,0,-2),(0,0,12,-6),(3,-1,-9,5)),((0,0,0,3),(0,0,0,-1),(0,0,12,-9),(3,-2,-6,5)),((0,0,0,3),(0,0,0,0),(0,0,12,-12),(3,-3,-3,3))]
def value(f,v,w):
 D,E,U,d=v;g,e,L,h=w
 if f.startswith('J'):return sum(v[i]*JM[int(f[1:])][i][j]*w[j]for i in range(4)for j in range(4))
 rows={'N':(D-E+U,D-E+d,D+3*U-d,U+E),'M':(D+2*U,U+d,2*D-E+6*U-d,U+E),'A':(D-E+3*U-2*d,D-E+2*U-2*d,3*D-2*E+7*U-2*d,2*D-2*E+2*U),'B':(D-E+3*U-2*d,U,2*D-E+8*U-4*d,U+E),'C':(D+4*U,U+d,4*D-E+14*U-d,U+E+2*d),'T':(D+2*U,D+U+d,3*D-E+8*U-d,2*D+4*U)}
 a,b,c,d=rows[f];return a*g-b*e+c*L-d*h
one=(1,1,1,1);t=(1,0,2,1);w=(3,0,6,1)
b7=inventory(7)-{one};b14=inventory(14)
assert b7=={tuple(s)for s in R['C2_sources']} and len(b14)==20
hits=[]
for f in ['N','M']:
 for v in b7:
  for z in b7:
   if value(f,v,z)==43:hits.append((f,v,z))
assert set(hits)=={('N',t,w),('N',w,t)}
support=set();admissible=0
for f in ['A','B','C','T','J0','J1','J2','J3']:
 for v in b14:
  for z in b14:
   yes=(z!=one and (f not in ['A','B']or v!=one))if not f.startswith('J')else((v!=one and z!=one)or(f=='J0'and v!=one and z==one)or(f=='J3'and v==one and z!=one))
   if yes:
    admissible+=1;n=value(f,v,z)
    if n<=43:support.add(n)
assert admissible==2964 and support=={15,27,29,31,33}
words=['']+[''.join(p)for k in [1,2]for p in product('OZJ',repeat=k)]
reach=set();trials=0
for a in words:
 for b in words:
  if len(a)+len(b)>2:continue
  for f in ['N','M','A','B','C','T','J0','J1','J2','J3']:
   trials+=1;n=value(f,word(t,a),word(w,b))
   if n<=87:reach.add(n)
assert trials==340 and reach=={43,69,77,79,81,83}
def count(A):
 n=len(A)
 H=sum(all((A[a]>>b)&1 for a,b in zip(p,p[1:]))for p in permutations(range(n)))
 C=sum(all((A[a]>>b)&1 for a,b in zip((0,)+p,p+(0,)))for p in permutations(range(1,n)))
 return H,C
for row in R['literal_graphs']:assert count(row['rows'])==(row['H'],row['C'])
parent=R['literal_graphs'][0]['rows'];n=len(parent);births=[0]*(1<<n)
# One independent exhaustive set of new-vertex Hamilton words, classified by
# the constraints each word places on the incoming/outgoing birth mask.
for p in permutations(range(n+1)):
 force_one=force_zero=0
 for a,b in zip(p,p[1:]):
  if a==n:force_one|=1<<b
  elif b==n:force_zero|=1<<a
  elif not(parent[a]>>b&1):break
 else:
  for mask in range(1<<n):
   if mask&force_one==force_one and not mask&force_zero:births[mask]+=1
assert births==R['birth_counts'] and not{85,87}&set(births)
# Exhaust all labelled rooted tournaments of orders <=4 to check that the
# source profiles used at43 occur only at the stated orders and graph types.
rooted=[]
for order in range(1,5):
 edges=[(i,j)for i in range(order)for j in range(i+1,order)]
 for mask in range(1<<len(edges)):
  A=[0]*order
  for k,(i,j)in enumerate(edges):
   if mask>>k&1:A[i]|=1<<j
   else:A[j]|=1<<i
  paths={subset:[]for subset in range(1<<order)}
  for k in range(order+1):
   for p in permutations(range(order),k):
    if all(A[a]>>b&1 for a,b in zip(p,p[1:])):paths[sum(1<<v for v in p)].append(p)
  full=(1<<order)-1
  for c in range(order):
   if sum(bool(p)and p[0]==c for p in paths[full])!=1:continue
   D=len(paths[full]);E=sum(p[-1]==c for p in paths[full]);d=len(paths[full^(1<<c)])
   U=sum(len(paths[full^S])*sum(bool(p)and p[0]==c for p in paths[S])for S in paths)
   if (D,E,U,d)in(t,w):rooted.append((order,(D,E,U,d)))
assert {o for o,s in rooted if s==t}=={2} and {o for o,s in rooted if s==w}=={4}
out={'status':'PASS','alternative_source_grammar':'indexed M_L and O','C2_profiles_cap7':len(b7),'C3_profiles_cap14':len(b14),'C3_admissible_pairs':admissible,'C2_H43_presentations':hits,'C3_support_through43':sorted(support),'descendant_tests':trials,'descendant_support_through87':sorted(reach),'literal_graphs_permutation_counted':len(R['literal_graphs']),'all_birth_masks_independently_classified':len(births),'small_source_rooted_occurrences':len(rooted)}
(HERE/'TARGET43_INDEPENDENT.json').write_text(json.dumps(out,indent=2,sort_keys=True)+'\n');print(json.dumps(out,sort_keys=True))
