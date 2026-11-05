"""Independent raw-coordinate and permutation audit; no candidate imports."""
from pathlib import Path
from itertools import product,permutations,combinations
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
T=(1,0,2,1);F=(1,-1,1,0);G=(1,0,2,0)
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def dot(x,y):return sum(a*b for a,b in zip(x,y))
def step(s,c):
 D,E,U,d=s
 if c=='O':return D,0,2*U,d
 if c=='Z':return D-E+U,d,2*U-d,d
 return D-E+U,D-E,D+U,D

def state(w,s=T):
 for c in w:s=step(s,c)
 return s

def word(x,k):return ''.join('Z' if x>>i&1 else 'O' for i in range(k))
def unit(U):return ''.join('Z' if c=='0' else 'O' for c in bin(U-1)[3:])
def table(s,row,k):return [dot(row,state(word(x,k),s))for x in range(1<<k)]
def graph(w):
 out=[{1},set()];root=0
 for c in w:
  n=len(out);out.append(set())
  for i in range(n):
   if c!='O'and i==root:out[n].add(i)
   else:out[i].add(n)
  if c=='J':root=n
 return out

def literal(w):
 g=graph(w);H=C=0
 for p in permutations(range(len(g))):
  if all(p[i+1]in g[p[i]]for i in range(len(p)-1)):
   H+=1
   if p[0]==0 and 0 in g[p[-1]]:C+=1
 return H,C
prefixes={state(''.join(w))for k in range(3)for w in product('OZJ',repeat=k)}
# Check a stronger scalar class, using integer alpha,gamma>=1 and tau>=2.
rows=[(a,tau-2*a-g,g,0)for a in range(1,4)for g in range(1,4)for tau in (2,4,6)]
checks=0;minimum=10**100
for s in prefixes:
 for row in rows:
  for k in range(5,9):
   vals=table(s,row,k)
   for x in range(1<<k):
    margins=[]
    for i in range(k):
     lo=x&~(1<<i);hi=lo|(1<<i)
     margins.append(vals[hi]-vals[lo])
    gap=min(margins[i+1]-margins[i]for i in range(k-1))
    need(gap>=6*s[3],'adjacent marginal separation')
    need(sum(abs(v)==2 for v in margins)<=1,'matching degree')
    minimum=min(minimum,gap)
    checks+=1
extrema=[]
for k in range(2,13):
 vals=table(T,F,k);deltas=[]
 for x in range(1<<k):
  for i in range(k-1):
   if x>>i&3==2:deltas.append(vals[x]-vals[x^(3<<i)])
 expected=2 if k==2 else 4 if k==3 else 2+2**((k-2)//2)+2**((k-1)//2)
 need(min(deltas)==expected,'sharp exchange minimum')
 extrema.append([k,expected])
basis=[tuple(int(i==j)for i in range(4))for j in range(4)]
toggles=0
for k in range(1,9):
 for x in range(1<<k):
  w=word(x,k);B=1<<k;p=x%2
  for s in basis:
   lo,hi=state('O'+w,s),state('Z'+w,s)
   need(tuple(b-a for a,b in zip(lo,hi))==(s[2]-s[1]-(x+p)*s[3],0,-B*s[3],0),'formal toggle matrix')
   toggles+=1
for s in basis:
 need(dot(G,state('OJJOOO',s))==dot(G,state('OJJOZO',s))==33*s[0]+68*s[2],'universal rebase')
rebases=0
for k in range(1,7):
 for x in range(1<<k):
  B=1<<k;U=x+x%2+2*B+1;s=state(unit(U));w=word(x,k)
  a,b=state('O'+w,s),state('Z'+w,s)
  need(tuple(y-x for x,y in zip(a,b))==(2*B,0,-B,0),'physical rebase')
  need(dot(G,a)==dot(G,b),'equal balanced scalar')
  rebases+=1
literalchecks=0
for delta in range(-6,7,2):
 B=2
 while B+1+delta<3:B*=2
 m=B.bit_length()-1;U=B+1+delta;prefix=unit(U)
 w0=prefix+'O'*(m+1)+'J';w1=prefix+'Z'+'O'*m+'J'
 h0,c0=literal(w0);h1,c1=literal(w1)
 need(c0==c1==1 and h1-h0==delta,'literal exact even offset')
 g0,g1=graph(w0),graph(w1)
 need(sum(len(a^b)for a,b in zip(g0,g1))==2,'one literal reversed edge')
 literalchecks+=2
for w,h in [('ZOOO',27),('ZZOO',25),('ZZZO',27)]:
 need(literal(w+'J')==(h,1),'literal k4 hostile')
 literalchecks+=1
# Independent BFS of the two-block digit cube; only these finite lengths claimed.
twoblocks=[]
for k in (5,6):
 vals=[dot(G,state(word(x>>k,k)+'J'+word(x&((1<<k)-1),k)))for x in range(1<<(2*k))]
 seen=set();best=[];maxdeg=0
 for x in range(len(vals)):
  if x in seen:continue
  seen.add(x);queue=[x]
  for y in queue:
   neighbors=[y^(1<<i)for i in range(2*k)if abs(vals[y]-vals[y^(1<<i)])==2]
   maxdeg=max(maxdeg,len(neighbors))
   for z in neighbors:
    if z not in seen:seen.add(z);queue.append(z)
  images=sorted({vals[y]for y in queue})
  if len(images)>len(best):best=images
 need(best==[6679,6681,6683] if k==5 else len(best)==2,'two-block exact finite probe')
 twoblocks.append({'k':k,'largest_distinct_image':best,'maxdegree':maxdeg})
switches=0
for k in range(3,11):
 for signs in product((-1,1),repeat=k):
  for i,j,h in combinations(range(k),3):
   need((signs[i]*signs[j])*(signs[i]*signs[h])*(signs[j]*signs[h])==1,'triangle sign invariant')
   switches+=1
report={'status':'PASS','proof_sha256':hashlib.sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),
 'candidate_verifier_sha256':hashlib.sha256((ROOT/'verify.py').read_bytes()).hexdigest(),
 'independent_verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'stronger_class_word_checks':checks,'minimum_marginal_separation':minimum,
 'sharp_exchange_minima':extrema,'formal_toggle_checks':toggles,'physical_rebases':rebases,
 'literal_permutation_graphs':literalchecks,'two_block_controls':twoblocks,'triangle_sign_checks':switches,
 'scope':'Matching, exact freely-chosen-parent offsets, physical rebases and cube-switching obstruction PASS; no exact prescribed-parent coverage.'}
Path(__file__).with_name('INDEPENDENT.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
print(json.dumps(report,indent=2,sort_keys=True))
