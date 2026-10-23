"""Independent orbit-volume chain and explicit Walsh eigenvector verification."""
import itertools,math,json
from fractions import Fraction
from collections import Counter
checks=0
def need(ok,why):
 global checks
 checks+=1
 if not ok:raise RuntimeError(why)
reports=[]
for n in range(2,6):
 edges=list(itertools.combinations(range(n),2));lookup={e:k for k,e in enumerate(edges)};M=len(edges);size=1<<M;full=size-1
 permutations=list(itertools.permutations(range(n)));tables=[];sourceflips=[];targetflips=[]
 for p in permutations:
  images=[lookup[tuple(sorted((p[u],p[v])))] for u,v in edges]
  sf=sum(1<<i for i,(u,v) in enumerate(edges) if p[u]>p[v]);tf=sum(1<<images[i] for i in range(M) if sf>>i&1)
  table=[0]*size
  for x in range(1,size):
   low=x&-x;table[x]=table[x-low]|(1<<images[low.bit_length()-1])
  tables.append(table);sourceflips.append(sf);targetflips.append(tf)
 gcan=[min(table[x] for table in tables) for x in range(size)]
 tcan=[min(table[x]^flip for table,flip in zip(tables,targetflips)) for x in range(size)]
 reps=sorted(set(gcan));orbit={h:set(table[h] for table in tables) for h in reps}
 if n==4:
  def cubic(x):
   vals=[]
   for v in range(n):
    sign=1
    for k,(a,b) in enumerate(edges):
     if v in (a,b):sign*=(-1 if (x>>k&1) else 1)*(1 if v==a else -1)
    vals.append(sign)
   return sum(vals)
  need(set(map(cubic,range(size)))=={-4,0,4},'cubic observable values')
  for g in set(gcan):need(sum(cubic(x) for x in range(size) if gcan[x]==g)==0,'cubic erased by graph projection')
  for t in set(tcan):need(len({cubic(x) for x in range(size) if tcan[x]==t})==1,'cubic tournament invariant')
 mu={h:Fraction(sum((-1)**((h&flip).bit_count()) for flip in sourceflips),math.factorial(n)) for h in reps}
 need(max(abs(mu[h]) for h in reps if h)<=Fraction(1,3),'signed moment bound')
 for paired in (False,True):
  glabel=[min(gcan[x],gcan[full^x]) if paired else gcan[x] for x in range(size)]
  tlabel=[min(tcan[x],tcan[full^x]) if paired else tcan[x] for x in range(size)]
  gs=sorted(set(glabel));ts=sorted(set(tlabel));counts=Counter(zip(glabel,tlabel));wg=Counter(glabel);wt=Counter(tlabel)
  K=[[sum((Fraction(counts[g,t]*counts[h,t],wg[g]*wt[t]) for t in ts),Fraction()) for h in gs] for g in gs]
  need(all(sum(row)==1 for row in K),'stochastic')
  need(all(wg[g]*K[i][j]==wg[h]*K[j][i] for i,g in enumerate(gs) for j,h in enumerate(gs)),'volume detailed balance')
  shapes=[h for h in reps if not paired or h.bit_count()%2==0]
  need(len(shapes)==len(gs),'quotient Fourier dimension')
  for h in shapes:
   vector=[sum((-1)**((a&g).bit_count()) for a in orbit[h]) for g in gs]
   need(any(vector),'nonzero Walsh orbit vector')
   for i in range(len(gs)):
    need(sum(K[i][j]*vector[j] for j in range(len(gs)))==mu[h]**2*vector[i],'exact eigenvector equation')
  positive_rank=sum(bool(mu[h]) for h in shapes)
  need(positive_rank=={2:1,3:2,4:3,5:10}[n],'exact rank')
  nontrivial=[mu[h]**2 for h in shapes if h]
  second=max(nontrivial,default=Fraction())
  need(second==(Fraction(1,9) if n>=3 else 0),'sharp second eigenvalue')
  need(all(K[i][i]>0 for i in range(len(gs))),'aperiodicity')
  reachable={0};stack=[0]
  while stack:
   i=stack.pop()
   for j,x in enumerate(K[i]):
    if x and j not in reachable:reachable.add(j);stack.append(j)
  need(len(reachable)==len(gs),'irreducibility')
  reports.append({'n':n,'paired':paired,'graph_states':len(gs),'tournament_states':len(ts),'second_eigenvalue':str(second),'gap':str(1-second),'eigenvectors_verified':len(shapes),'positive_rank':positive_rank,'invisible_tournament_dimensions':len(ts)-positive_rank,'invisible_graph_dimensions':len(gs)-positive_rank,'stationary_measure':'labelled orbit volume / 2^binom(n,2)'})
# Direct moment induction checks on all labelled n<=5 graphs; recurrence is signed.
for n in range(2,6):
 edges=list(itertools.combinations(range(n),2));perms=list(itertools.permutations(range(n)))
 for mask in range(1,1<<len(edges)):
  E=[e for k,e in enumerate(edges) if mask>>k&1]
  direct=Fraction(sum((-1)**sum(p[a]>p[b] for a,b in E) for p in perms),math.factorial(n))
  terms=[]
  for v in range(n):
   remaining=[u for u in range(n) if u!=v];er=[(a,b) for a,b in E if v not in (a,b)]
   vals=[]
   for order in itertools.permutations(remaining):
    pos={a:i for i,a in enumerate(order)};vals.append((-1)**sum(pos[a]>pos[b] for a,b in er))
   sign=(-1)**sum(b==v for a,b in E)
   terms.append(sign*Fraction(sum(vals),math.factorial(n-1)))
  need(direct==sum(terms)/n,'first-vertex signed recurrence')
  need(abs(direct)<=Fraction(1,3),'all labelled moment bound')
  degrees=sorted(sum(v in e for e in E) for v in range(n))
  positive=[d for d in degrees if d]
  shape=positive in ([1,1,2],[2,2,2,2])
  need((abs(direct)==Fraction(1,3))==shape,'independent equality classification')
print(json.dumps({'status':'PASS','checks':checks,'reports':reports,'scope':'exact operator eigenvectors on quotient carriers n2..5 plus labelled moment recurrence; all-size theorem audited separately'},indent=2))
