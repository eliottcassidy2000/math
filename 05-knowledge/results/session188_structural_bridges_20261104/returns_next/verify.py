"""Exact simultaneous-module and pair-hull checks; standard library only."""
import itertools,json,hashlib
from functools import lru_cache
from pathlib import Path
from probe import hp,ext

def need(x,msg):
 if not x:raise RuntimeError(msg)
def rows_from_mask(n,mask):
 R=[0]*n
 for k,(u,v) in enumerate(itertools.combinations(range(n),2)):
  if mask>>k&1:R[u]|=1<<v
  else:R[v]|=1<<u
 return tuple(R)
def induced(R,S):
 vs=[v for v in range(len(R)) if S>>v&1]
 return tuple(sum(1<<j for j,v in enumerate(vs) if R[u]>>v&1) for u in vs)
def is_module(R,S):
 return all((R[v]&S) in (0,S) for v in range(len(R)) if not S>>v&1)
def modules(R):
 return [S for S in range(1,1<<len(R)) if S.bit_count()>=2 and is_module(R,S)]
def hull(R,S):
 while True:
  add=sum(1<<v for v in range(len(R)) if not S>>v&1 and (R[v]&S) not in (0,S))
  if not add:return S
  S|=add

def reverse(R,sets):
 ans=list(R)
 for u,v in itertools.combinations(range(len(R)),2):
  if sum(bool(S>>u&1 and S>>v&1) for S in sets)%2:ans[u]^=1<<v;ans[v]^=1<<u
 return tuple(ans)
@lru_cache(None)
def profile(R):
 n=len(R);N=(1<<n)-1;H=[hp(induced(R,S)) for S in range(N+1)];Q=[[0]*(N+1) for _ in range(n+1)];Q[0][0]=1
 for k in range(1,n+1):
  for S in range(1,N+1):
   B=S
   while B:Q[k][S]+=H[B]*Q[k-1][S^B];B=(B-1)&S
 return tuple(Q[k][N] for k in range(n+1))
def quotient(R,S):
 rep=(S&-S).bit_length()-1;return induced(R,((1<<len(R))-1)^S|(1<<rep))
def rank(vectors):
 basis={}
 for v in vectors:
  while v:
   p=v.bit_length()-1
   if p in basis:v^=basis[p]
   else:basis[p]=v;break
 return len(basis)
@lru_cache(None)
def core_query(R):
 O=reverse(R,[(1<<len(R))-1])
 for a,c in itertools.product(range(1<<len(R)),repeat=2):
  if (hp(ext(R,a,c))-hp(ext(O,a,c)))%4==2:return a,c
 raise RuntimeError('incoming E3 control failed')
def lifted(R,M,a,c):
 vs=[u for u in range(len(R)) if M>>u&1];rep=vs[0]
 A=sum(1<<u for k,u in enumerate(vs) if a>>k&1);C=sum(1<<u for k,u in enumerate(vs) if c>>k&1)
 for u in range(len(R)):
  if not M>>u&1 and R[u]>>rep&1:A|=1<<u;C|=1<<u
 return A,C

def substitute(Q,blocks):
 starts=[0]
 for B in blocks:starts.append(starts[-1]+len(B))
 R=[0]*starts[-1];masks=[]
 for i,B in enumerate(blocks):
  masks.append(((1<<len(B))-1)<<starts[i])
  for j,row in enumerate(B):
   R[starts[i]+j]=row<<starts[i]
   for k in range(len(blocks)):
    if Q[i]>>k&1:R[starts[i]+j]|=((1<<len(blocks[k]))-1)<<starts[k]
 return tuple(R),masks

def run():
 counts={'labelled_tournaments':0,'pair_hulls':0,'modules':0,'disjoint_pairs':0,'overlapping_pairs':0,'uniform_nested_pairs':0,'ordered_quotient_controls':0}
 examples={}
 for n in range(1,6):
  for code in range(1<<(n*(n-1)//2)):
   R=rows_from_mask(n,code);Ms=modules(R);counts['labelled_tournaments']+=1;counts['modules']+=len(Ms)
   Hs={}
   for u,v in itertools.combinations(range(n),2):
    S=hull(R,(1<<u)|(1<<v));need(is_module(R,S),'hull is module');Hs.setdefault(S,(u,v));counts['pair_hulls']+=1
   need(set(Hs)==set(Ms),'every module has a generating pair')
   edges=list(itertools.combinations(range(n),2));vectors=[sum(1<<k for k,(u,v) in enumerate(edges) if M>>u&1 and M>>v&1) for M in Ms]
   need(rank(vectors)==len(Ms),'module incidence full rank')
   for M,(u,v) in Hs.items():
    for N in Ms:need(bool(N>>u&1 and N>>v&1)==(M&N==M),'inclusion zeta minor')
   h=hp(R);need(len(Ms)<=n*(n-1)//2-(2 if h>1 else 0),'sharp module-count bound')
   for A,B in itertools.combinations(Ms,2):
    T=reverse(R,[A,B]);I=A&B;P=A&~B;Z=B&~A
    if not I:
     aa,cc=core_query(induced(R,A));a,c=lifted(R,A,aa,cc)
     val=hp(ext(R,a,c));both=hp(ext(T,a,c));other=hp(ext(reverse(R,[B]),a,c));single=hp(ext(reverse(R,[A]),a,c))
     need(val==other and both==single,'other disjoint module cancels exactly')
     need((val-both)%4==2,'disjoint actual query separates')
     counts['disjoint_pairs']+=1
    elif P and Z:
     U=A|B
     need(all(is_module(R,S) for S in (P,I,Z,U)),'overlap part modules')
     before=profile(induced(R,U));after=profile(induced(T,U));need(all(b>=a for a,b in zip(before,after)),'overlap every cover channel nondecreasing')
     need(after[1]>=3*before[1],'overlap triple Hamilton floor')
     bound=2*hp(quotient(R,U))*before[1];need(hp(T)-h>=bound,'ambient overlap bound')
     counts['overlapping_pairs']+=1
     if 'overlap' not in examples:examples['overlap']={'rows':R,'A':A,'B':B,'target':T,'H':[h,hp(T)],'bound':bound}
    else:
     small,big=(A,B) if A&B==A else (B,A);C=big^small;u=(small&-small).bit_length()-1
     cross=R[u]&C
     if cross not in (0,C):continue
     need(all((R[v]&C)==cross for v in range(n) if small>>v&1),'nested uniform cross direction')
     # Place x between the two blocks in the original direction.
     a=small if cross==C else C;rep=(big&-big).bit_length()-1
     for v in range(n):
      if not big>>v&1 and R[v]>>rep&1:a|=1<<v
     c=(1<<n)-1
     diff=hp(ext(T,a,c))-hp(ext(R,a,c));bound=2*hp(quotient(R,big))*hp(induced(R,small))*hp(induced(R,C))
     need(diff>=bound,'uniform nested actual two-port bound');counts['uniform_nested_pairs']+=1
 # Arbitrary quotient controls on nontrivial blocks, independently of module-pair enumeration.
 bank=((0,),(2,0),(2,4,1))
 for s in (2,3):
  L=tuple(sum(1<<j for j in range(i+1,s)) for i in range(s))
  for blocks in itertools.product(bank,repeat=s):
   if sum(map(len,blocks))>6:continue
   original,_=substitute(L,blocks);orig=profile(original)
   for qcode in range(1<<(s*(s-1)//2)):
    Q=rows_from_mask(s,qcode)
    for flags in (0,(1<<s)-1):
     Bs=[reverse(B,[(1<<len(B))-1]) if flags>>i&1 else B for i,B in enumerate(blocks)];target,_=substitute(Q,Bs);p=profile(target)
     need(all(b>=a for a,b in zip(orig,p)),'ordered quotient monotonicity')
     need(p[1]>=hp(Q)*orig[1],'quotient Hamilton floor')
     if hp(Q)==1:need(p==orig,'transitive quotient profile invariance')
     counts['ordered_quotient_controls']+=1
 # Exact coarse-query hostile, retaining every actual query rather than a truncated expansion.
 hostiles=[]
 for m in (2,3):
  block=tuple(sum(1<<j for j in range(i+1,m)) for i in range(m));Q=(2,4,1);Qo=(4,1,2)
  R,blocks=substitute(Q,[block]*3);T,_=substitute(Qo,[block]*3);vals=[]
  for a,c in itertools.product(range(8),repeat=2):
   aa=sum(blocks[i] for i in range(3) if a>>i&1);cc=sum(blocks[i] for i in range(3) if c>>i&1)
   vals.append([a,c,hp(ext(R,aa,cc))-hp(ext(T,aa,cc))])
  need(all(v[2]%4==0 for v in vals) and any(v[2] for v in vals),'uniform quotient parity hostile')
  need(max(abs(v[2]) for v in vals)==(64 if m==2 else 2508),'exact hostile maximum')
  hostiles.append({'module_size':m,'queries':vals})
 result={'status':'PASS','counts':counts,'examples':examples,'coarse_parity_hostiles':hostiles,'source_sha256':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py','probe.py')}}
 Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,sort_keys=True))
if __name__=='__main__':run()
