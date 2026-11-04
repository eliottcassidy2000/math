"""Independent ordered-block audit: literal words/cuts and physical port DP."""
from itertools import combinations,permutations,product
from math import comb,prod
from functools import lru_cache
import json

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def rows(n,code):
 a=[0]*n
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(j,i) if code>>bit&1 else(i,j);a[u]|=1<<v
 return tuple(a)

def induced(a,mask):
 names=[i for i in range(len(a)) if mask>>i&1]
 return tuple(sum(1<<j for j,v in enumerate(names) if a[u]>>v&1) for u in names)

def op(a):return tuple(((1<<len(a))-1)^r^(1<<i) for i,r in enumerate(a))

def module(a,M):return all(r&M in(0,M) for i,r in enumerate(a) if not M>>i&1)

def flip(a,M):return tuple(r^(M^(1<<i)) if M>>i&1 else r for i,r in enumerate(a))

def disjoint_or_nested(A,B):return A&B in(0,A,B)

@lru_cache(None)
def q_words(a):
 n=len(a)
 if not n:return (1,)
 q=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(not(a[u]>>v&1) for u,v in zip(w,w[1:]));good=n-1-bad
  for k in range(bad+1,n+1):q[k]+=comb(good,k-1-bad)
 return tuple(q)

@lru_cache(None)
def hp(a):
 n=len(a);full=(1<<n)-1
 @lru_cache(None)
 def finish(mask,last):
  if mask==full:return 1
  out=0;rem=a[last]&~mask
  while rem:
   b=rem&-rem;rem-=b;out+=finish(mask|b,b.bit_length()-1)
  return out
 return sum(finish(1<<i,i) for i in range(n)) if n else 1

def blow(Q,blocks):
 offsets=[0]
 for b in blocks:offsets.append(offsets[-1]+len(b))
 a=[0]*offsets[-1]
 for i,b in enumerate(blocks):
  for u,r in enumerate(b):
   a[offsets[i]+u]=r<<offsets[i]
   for j,c in enumerate(blocks):
    if Q[i]>>j&1:a[offsets[i]+u]|=((1<<len(c))-1)<<offsets[j]
 return tuple(a)

def quotient(a,blocks):
 reps=[(M&-M).bit_length()-1 for M in blocks]
 return tuple(sum(1<<j for j,v in enumerate(reps) if a[u]>>v&1) for u in reps)

def scc(a):
 reach=[r|(1<<i) for i,r in enumerate(a)]
 for k in range(len(a)):
  for i in range(len(a)):
   if reach[i]>>k&1:reach[i]|=reach[k]
 todo=set(range(len(a)));out=[]
 while todo:
  v=min(todo);S={u for u in todo if reach[u]>>v&1 and reach[v]>>u&1}
  todo-=S;out.append(sum(1<<u for u in S))
 Q=quotient(a,out)
 return [out[i] for i in sorted(range(len(out)),key=lambda i:-Q[i].bit_count())]

def extend(a,x,c):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for v in range(n):
  if x>>v&1:b[v]|=1<<n
  else:b[n]|=1<<v
  if c>>v&1:b[v]|=1<<(n+1)
  else:b[n+1]|=1<<v
 return tuple(b)

def ordered_cases():
 pool=(rows(1,0),rows(2,0),rows(3,0),(2,4,1))
 cases=coordinates=lifts=0;records=[]
 for blocks in product(pool,repeat=3):
  if sum(map(len,blocks))>6:continue
  src=blow(rows(3,0),blocks);mask_blocks=[];off=0
  for B in blocks:mask_blocks.append(((1<<len(B))-1)<<off);off+=len(B)
  for pi in permutations(range(3)):
   if pi==(0,1,2):continue
   targetQ=tuple(sum(1<<j for j in range(3) if pi.index(i)<pi.index(j)) for i in range(3))
   dst=blow(targetQ,tuple(op(b) for b in blocks))
   t=next(t for t in(1,2) if set(range(t))!=set(pi[:t]));c=sum(mask_blocks[:t])
   X,Y=extend(src,0,c),extend(dst,0,c)
   if len(X)<=7:
    qx,qy=q_words(X),q_words(Y)
    need(all(v>=u for u,v in zip(qx,qy)),"ordered-slot profile domination")
    coordinates+=len(qx)
   gap=hp(Y)-hp(X);lower=2*prod(hp(b) for b in blocks)
   need(gap>=lower,"physical cut strict bound")
   T=blow((2,4,1),(src,(0,),(0,)));U=blow((2,4,1),(dst,(0,),(0,)))
   n=len(src);ext=1<<(n+1)
   ambient_gap=hp(extend(U,ext,c|ext))-hp(extend(T,ext,c|ext))
   need(ambient_gap>=3*lower,"strong ambient one-run lower bound")
   lifts+=1;cases+=1
   if len(records)<4:records.append({"blocks":blocks,"pi":pi,"cut":t,"gap":gap,"ambient_gap":ambient_gap,"lower":lower})
 return {"cases":cases,"literal_cover_coordinates":coordinates,"strong_ambient_lifts":lifts,"examples":records}

def family_cases():
 families=hereditary=port_cases=module_steps=0;by_n=[]
 for n in range(1,6):
  nf=nh=np=0
  for code in range(1<<(n*(n-1)//2)):
   T=rows(n,code);mods=[M for M in range(1,1<<n) if M.bit_count()>1 and module(T,M)]
   for choice in range(1,1<<len(mods)):
    F=[M for j,M in enumerate(mods) if choice>>j&1]
    if not all(disjoint_or_nested(A,B) for A,B in combinations(F,2)):continue
    U=T
    for M in reversed(F):
     need(module(U,M),"laminar operation lost current modularity");U=flip(U,M);module_steps+=1
    families+=1;nf+=1
    for S in range(1,1<<n):
     need(q_words(induced(T,S))==q_words(induced(U,S)),"hereditary cover invariance")
     hereditary+=1;nh+=1
    maxima=[B for B in F if not any(B!=D and B&D==B for D in F)]
    for B in maxima:
     inside=[A for A in F if A!=B and A&B==A]
     children=[A for A in inside if not any(A!=D and A&D==A for D in inside)]
     if sum(children)!=B:continue
     K=quotient(T,children);components=scc(K)
     if len(components)<2:continue
     R=[sum(children[i] for i in range(len(children)) if comp>>i&1) for comp in components]
     need(all(module(T,M) and module(U,M) for M in R),"component blocks not modules")
     need(quotient(U,R)==op(quotient(T,R)),"component order not reversed")
     need(all(q_words(induced(T,M))==q_words(induced(U,M)) for M in R),"component profiles changed")
     b=(B&-B).bit_length()-1
     ext=sum(1<<v for v in range(n) if not B>>v&1 and T[v]>>b&1)
     a,c=ext,ext|R[0]
     ambient_blocks=[B]+[1<<v for v in range(n) if not B>>v&1]
     Hroot=hp(quotient(T,ambient_blocks))
     lower=2*Hroot*prod(hp(induced(T,M)) for M in R)
     gap=hp(extend(U,a,c))-hp(extend(T,a,c))
     need(gap>=lower,"covered laminar physical cut or exterior cancellation")
     port_cases+=1;np+=1
  by_n.append({"n":n,"laminar_families":nf,"hereditary_subsets":nh,"covered_nonstrong_maximal_cuts":np})
 return {"families":families,"hereditary_subsets":hereditary,"physical_port_cuts":port_cases,"current_module_steps":module_steps,"by_n":by_n}

if __name__=='__main__':
 print(json.dumps({"status":"PASS","ordered_blocks":ordered_cases(),"laminar_families":family_cases(),"scope":"All-order proof audit plus specified finite literal controls, not a universal response reconstruction theorem."},indent=2,sort_keys=True))
