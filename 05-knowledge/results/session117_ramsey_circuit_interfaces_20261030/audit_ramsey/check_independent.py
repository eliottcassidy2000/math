import itertools,json,math
from collections import Counter

def need(test,msg):
 if not test: raise RuntimeError(msg)
def tournament(n,bits):
 a=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  a[i][j]=bool(bits>>k&1);a[j][i]=not a[i][j]
 return a
def double(a):
 n=len(a);b=[[False]*(2*n) for _ in range(2*n)]
 for i,j in itertools.permutations(range(2*n),2):
  u,c=divmod(i,2);v,d=divmod(j,2)
  b[i][j]=(c==0) if u==v else (a[v][u] if c==d==1 else a[u][v])
 return b
def trans(a,s):
 # Independent triangle criterion (not producer's induced-degree test).
 return not any((a[i][j] and a[j][k] and a[k][i]) or (a[j][i] and a[k][j] and a[i][k]) for i,j,k in itertools.combinations(s,3))
def subsets(n):
 for m in range(1<<n):yield tuple(i for i in range(n) if m>>i&1)
def poly(a):return Counter(len(s) for s in subsets(len(a)) if trans(a,s))
def verify(a):
 n=len(a);b=double(a);p=poly(a);q=Counter();z=Counter();lifts=Counter()
 for s in subsets(2*n):
  if not trans(b,s):continue
  q[len(s)]+=1
  support=tuple(sorted(set(v//2 for v in s)))
  twins=len(s)-len(support)
  if twins:
   need(twins==1,'two twins');need(trans(a,support),'projection not transitive');lifts[support]+=1
  else:z[len(s)]+=1
 for s in subsets(n):
  expected=len(s)*(len(s)+1)//2 if trans(a,s) else 0
  need(lifts[s]==expected,'lift multiplicity')
 for k,c in p.items():z[k+1]+=k*(k+1)//2*c
 need(+z==+q,'derivative polynomial')
 need(max(q)<=min(2*max(p),max(n,max(p)+1)),'max bound')
 for bits in itertools.product((0,1),repeat=n):
  pred=True
  for tri in itertools.combinations(range(n),3):
   if not trans(a,tri):pred &= sum(bits[i] for i in tri)==2
   else:
    order=sorted(tri,key=lambda i:-sum(a[i][j] for j in tri))
    pred &= tuple(bits[i] for i in order)!=(1,0,1)
  need(pred==trans(b,tuple(2*i+bits[i] for i in range(n))),'3CNF predicate')
 return q

def alt(l):
 # u indices 0..l-1, w indices l..2l.
 n=2*l+1;a=[[False]*n for _ in range(n)]
 for i,j in itertools.combinations(range(n),2):
  a[i][j]=(i<j) if j<l else ((i>j) if i>=l else not (j-l<=i))
  a[j][i]=not a[i][j]
 return a
count=0
for n in range(1,6):
 for bits in range(1<<(n*(n-1)//2)):verify(tournament(n,bits));count+=1
rows=[]
for l in range(1,5):
 a=alt(l);p=poly(a);expected=Counter()
 for i in range(l):
  for k,c in enumerate((1,l+2,math.comb(l+2,2))):expected[i+k]+=math.comb(l-1,i)*c
 need(p==expected,'alternating base polynomial')
 q=verify(a) if l<=2 else poly(double(a))
 need(max(q)==2*l+1,'alternating maximum');need(q[max(q)]==(12 if l==1 else 1),'max count')
 rows.append({'l':l,'maximum':max(q),'multiplicity':q[max(q)]})
# Native cyclic-triple equations are necessary; their solution set can be
# independently enumerated without doubling all 34 vertices at l=8.
full=[]
for l in range(1,9):
 solutions=[]
 for bits in itertools.product((0,1),repeat=2*l+1):
  if all(bits[l+a]+bits[j]+bits[l+b]==2 for j in range(l) for a in range(j+1) for b in range(j+1,l+1)):
   need(trans(double(alt(l)),tuple(2*i+bits[i] for i in range(2*l+1))),'necessary equations not sufficient here')
   solutions.append(bits)
 need(len(solutions)==(3 if l==1 else 1),'full support count')
 full.append(len(solutions))
print(json.dumps({'status':'PASS','all_labeled_tournaments_n1_to5':count,'alternating_rows':rows,'full_support_counts_l1_to8':full,'method':'independent triangle scans; no imports from producer; explicit exception gates'},indent=2))
