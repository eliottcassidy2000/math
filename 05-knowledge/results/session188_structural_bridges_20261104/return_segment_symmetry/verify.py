from functools import lru_cache
from itertools import combinations
from collections import defaultdict,Counter
from math import factorial
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
def need(c,m):
 if not c:raise RuntimeError(m)
def quotient(m,code):
 T=[0]*m
 for b,(i,j) in enumerate(combinations(range(m),2)):
  u,v=(i,j) if code>>b&1 else(j,i);T[u]|=1<<v
 return T
def augment(Q,a,c):
 m=len(Q);T=Q[:]+[1<<(m+1),0]
 for i in range(m):
  if a>>i&1:T[i]|=1<<m
  else:T[m]|=1<<i
  if c>>i&1:T[i]|=1<<(m+1)
  else:T[m+1]|=1<<i
 return T
def count(T,k):
 @lru_cache(None)
 def f(left,last):
  if not any(left):return 1
  total=0
  for j,n in enumerate(left):
   if n and (last<0 or T[last]>>j&1):
    nxt=list(left);nxt[j]-=1;total+=f(tuple(nxt),j)
  return total
 return f(tuple(k)+(1,1),-1)
def words(T,k):
 out=[]
 def go(left,word):
  if not any(left):out.append(tuple(word));return
  for j,n in enumerate(left):
   if n and(not word or T[word[-1]]>>j&1):
    nxt=left[:];nxt[j]-=1;go(nxt,word+[j])
 go(list(k)+[1,1],[]);return out
def pieces(w,i):
 p=[j for j,x in enumerate(w) if x==i]
 return w[:p[0]],tuple(w[a+1:b] for a,b in zip(p,p[1:])),w[p[-1]+1:]
def returns(Q,k,i,R):
 caps=[x//R for x in k];caps[i]=0;out=[]
 def go(left,w):
  if w and Q[w[-1]]>>i&1:out.append(tuple(w))
  for j,n in enumerate(left):
   if n and Q[i if not w else w[-1]]>>j&1:
    nxt=left[:];nxt[j]-=1;go(nxt,w+[j])
 go(caps,[]);return out
def triangle(Q,k,i,R):
 H=[j for j,x in enumerate(k) if j!=i and x>=R]
 return any(Q[i]>>u&1 and Q[u]>>v&1 and Q[v]>>i&1 for u in H for v in H)
def parity_record(Q,a,c,k,i):
 R=k[i]-1;need(R>=2 and not R&(R-1),'power-of-two return count');T=augment(Q,a,c);v=count(T,k);terms=[]
 for W in returns(Q,k,i,R):
  kk=k[:];kk[i]=1
  for j in W:kk[j]-=R
  need(min(kk)>=0,'nonnegative contracted multiplicities');terms.append([list(W),kk,count(T,kk)])
 rhs=sum(t[2] for t in terms);tri=triangle(Q,k,i,R)
 need((v-rhs)%2==0,'exact parity contraction');need(tri or v%2==0,'heavy-label triangle obstruction')
 return [Q,a,c,k,i,v,terms,tri]
def main():
 records=[]
 # Complete three-child queries, all marked labels and a varied multiplicity domain.
 for code in range(8):
  Q=quotient(3,code)
  for a in range(8):
   for c in range(8):
    for i in range(3):
     others=[j for j in range(3) if j!=i]
     for uv in [(1,1),(2,1),(1,2),(2,2),(3,2),(2,3),(3,3),(4,2)]:
      k=[0]*3;k[i]=3
      for j,z in zip(others,uv):k[j]=z
      records.append(parity_record(Q,a,c,k,i))
 # Four-return controls include both repeated-segment sharp cases and zero residuals.
 for code in range(8):
  Q=quotient(3,code)
  for a,c in [(0,0),(0,7),(7,0),(7,7),(1,3),(2,6),(5,1),(3,4)]:
   for i in range(3):
    others=[j for j in range(3) if j!=i]
    for uv in [(1,1),(4,1),(4,4),(5,4),(8,4)]:
     k=[0]*3;k[i]=5
     for j,z in zip(others,uv):k[j]=z
     records.append(parity_record(Q,a,c,k,i))
 # Exact explicit orbit census, with both free and nonfree actions.
 orbit_records=[];explicit_words=0
 for code in range(8):
  Q=quotient(3,code)
  for a,c in [(0,7),(0,0),(7,0),(1,3),(2,6),(5,1)]:
   for k in ([3,2,2],[3,3,3],[5,4,4]):
    ww=words(augment(Q,a,c),k);need(len(ww)==count(augment(Q,a,c),k),'explicit words versus recursion');groups=defaultdict(list)
    for w in ww:
     P,L,S=pieces(w,0);groups[(P,tuple(sorted(L)),S)].append(w)
    hist=Counter()
    for (P,L,S),vs in groups.items():
     multiplicities=Counter(L);size=factorial(k[0]-1)
     for n in multiplicities.values():size//=factorial(n)
     need(len(vs)==size,'exact return-permutation orbit');need((size%2==1)==(len(multiplicities)==1),'power-of-two odd-orbit characterization');hist[size]+=1
    explicit_words+=len(ww);orbit_records.append([Q,a,c,k,len(ww),sorted(hist.items())])
 # Full free-action samples with one other repeated child.
 factorial_records=[]
 for code in range(64):
  Q=quotient(4,code)
  for a,c in [(0,15),(1,3),(6,9),(15,0)]:
   for r in (3,4,5):
    k=[r,r-1,1,1];v=count(augment(Q,a,c),k);need(v%factorial(r-1)==0,'factorial divisibility');factorial_records.append([Q,a,c,k,v])
 sharp=[]
 for R in range(1,8):
  # i->j->V->i, V transitive; source/sink ports.
  Q=[2,sum(1<<j for j in range(2,R+2))]+[1|sum(1<<j for j in range(i+1,R+2)) for i in range(2,R+2)]
  k=[R+1,R]+[1]*R;v=count(augment(Q,0,(1<<len(Q))-1),k);need(v==factorial(R),'sharp factorial control');sharp.append([Q,k,v])
 parity_sharp=[]
 for R in (2,4,8):
  Q=[2,4,1];k=[R+1,R,R];v=count(augment(Q,0,7),k);need(v==1,'sharp odd triangle control');parity_sharp.append(parity_record(Q,0,7,k,0))
 nonsufficient=next((r for r in records if r[-1] and r[5]%2==0),None);need(nonsufficient is not None,'triangle is only necessary')
 pins={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py']};dep=ROOT.parent/'coarse_mod8_affine/PROOF.md'
 return dict(status='PASS',source_pins=pins,dependency_pin=hashlib.sha256(dep.read_bytes()).hexdigest(),parity_controls=len(records),complete_parity_records=records,orbit_controls=orbit_records,explicit_words=explicit_words,factorial_controls=factorial_records,sharp_factorial=sharp,sharp_odd=parity_sharp,triangle_not_sufficient=nonsufficient,scope='All-order return-word symmetry and necessary parity obstruction, not a response separator or an ordered-cover theorem.')
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
