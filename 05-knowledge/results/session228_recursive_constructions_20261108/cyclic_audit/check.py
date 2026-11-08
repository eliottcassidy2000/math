"""Independent literal permutation/cut audit; no producer imports."""
import itertools as it, math, json, hashlib
from pathlib import Path

def require(v,msg):
 if not v: raise RuntimeError(msg)
def graphs(n):
 edges=list(it.combinations(range(n),2))
 for mask in range(1<<len(edges)):
  yield tuple(tuple(i!=j and ((mask>>edges.index(tuple(sorted((i,j))))&1)==(i<j)) for j in range(n)) for i in range(n))
def permstats(g,root=None):
 n=len(g); F=[0]*(n+1); E=U=hs=0
 for p in it.permutations(range(n)):
  bad=[i for i in range(n-1) if not g[p[i]][p[i+1]]]; b=len(bad)
  for d in range(b+1,n+1): F[d]+=math.comb(n-1-b,d-1-b)
  if root is not None:
   if not bad:
    E+=p[-1]==root
    hs+=p[0]==root
   for cut in range(1,n):
    if p[cut]==root and all(i==cut-1 for i in bad): U+=1
 return tuple(F),E,U+hs,hs

def cyc(A,B):
 na=len(A);nb=len(B);n=na+nb+1
 return tuple(tuple(False if i==j else A[i][j] if i<na and j<na else B[i-na][j-na] if na<=i<n-1 and na<=j<n-1 else (i<na and na<=j<n-1) or (na<=i<n-1 and j==n-1) or (i==n-1 and j<na) for j in range(n)) for i in range(n))
def state(a,r,b,s):
 d=a*b; D=3*d+a*s+b*r+r*s
 return D,d,D+b*r,d

def inverse(X,physical):
 D,E,U,d=X;t=U-D;z=2*D-U-3*d
 if d<=0 or E!=d or t<0 or z<0:return set()
 out=set()
 for b in range(1,d+1):
  if d%b or t%b or b*z%(d+t):continue
  a=d//b;r=t//b;s=b*z//(d+t)
  if (a,r) in physical and (b,s) in physical:out.add((a,r,b,s))
 return out

def Q(f,k):
 return sum(math.comb(k,d)*f[d] for d in range(1,min(k,len(f)-1)+1))

def main():
 banks={n:list(graphs(n)) for n in range(1,6)}; pairs={}; floors=0
 for n,gs in banks.items():
  for g in gs:
   f,*_=permstats(g);pair=(f[1],f[2] if n>=2 else 0);pairs[g]=pair
   require(pair[1]>=2**n-2,'F2 floor');floors+=1
 counts=0;images={}
 for na in range(1,5):
  for nb in range(1,5):
   if na+nb>6:continue
   for A in banks[na]:
    for B in banks[nb]:
     C=cyc(A,B);f,e,u,hs=permstats(C,len(C)-1);a,r=pairs[A];b,s=pairs[B]
     X=(f[1],e,u,a*b)
     fa=permstats(A)[0];fb=permstats(B)[0]
     for k in range(1,len(C)+1):
      pred=k*(Q(fa,k)*Q(fb,k)+(Q(fa,k+1)-Q(fa,k))*(Q(fb,k+1)-Q(fb,k))+(Q(fa,k)-Q(fa,k-1))*(Q(fb,k)-Q(fb,k-1)))
      require(Q(f,k)==pred,'full Q identity')
     require(X==state(a,r,b,s),'state');require(hs==a*b,'root start');counts+=1
     require(f[1]>=2**(na+nb-1)+1,'order cap')
     images.setdefault(X,set()).add((a,r,b,s))
 physical=set(pairs.values()); triples=0
 for p in physical:
  for q in physical:
   X=state(*p,*q);require((*p,*q) in inverse(X,physical),'inverse completeness');triples+=1
   for witness in inverse(X,physical):require(state(*witness)==X,'inverse soundness')
 require(inverse((3,1,3,1),physical)=={(1,0,1,0)},'singleton boundary')
 require(not inverse((7,1,7,1),physical),'formal false positive')
 hostile=tuple(tuple(bool(row>>j&1) for j in range(5)) for row in (20,1,2,7,14))
 hf,he,hu,_=permstats(hostile,1)
 deleted=tuple(tuple(hostile[i][j] for j in (0,2,3,4)) for i in (0,2,3,4))
 hd=permstats(deleted)[0][1]
 hx=(hf[1],he,hu,hd)
 require(hx==(9,3,9,3) and not inverse(hx,physical),'physical hostile')
 for v in range(5):
  reached={v}
  for _ in range(5):reached|={j for i in reached for j in range(5) if hostile[i][j]}
  require(len(reached)==5,'hostile strongly connected')
 scalar=sorted({state(*p,*q)[0] for p in physical for q in physical if state(*p,*q)[0]<=64})
 require(13 not in scalar,'13 obstruction')
 producer=json.loads((Path(__file__).resolve().parent/'producer_CHECK.json').read_text())
 require(scalar==[x for x in producer['complete_binary_cyclic_counts_le_128'] if x<=64],'independent scalar image')
 print(json.dumps({'physical_hostile_state':hx,'physical_hostile_strong':True,'complete_scalar_image_le_64':scalar,'status':'PASS','literal_factor_tournaments':floors,'literal_cyclic_pairs':counts,'distinct_tested_states':len(images),'physical_pair_bank':sorted(physical),'inverse_pair_combinations':triples,'hostile_formal_state':[7,1,7,1],'U_convention':'empty first component allowed; second nonempty and root-start'},sort_keys=True,indent=2))
if __name__=='__main__': main()
